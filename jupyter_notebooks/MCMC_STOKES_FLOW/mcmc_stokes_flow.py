import os, sys
import re
import subprocess
import numpy as np
import random
from matplotlib import pyplot as plt
import pandas as pd
import pymc as pm
import arviz as az
import time
from shutil import rmtree
from IPython.core.pylabtools import figsize
from scipy.stats.mstats import mquantiles

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(ASPECT_LAB_DIR)
from shilofue.VtkPp import VTKP
from shilofue.ParsePrm import ParseFromDealiiInput, ParseToDealiiInput

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

mcmc_dir = "/mnt/lochy/ASPECT_DATA/MCMC_STOKES_FLOW"
assert(os.path.isdir(mcmc_dir))

# Define some utilities functions
def MakeNewPrm(prm_path_ori, prm_path, sinker_viscosity, mantle_sinker_viscosity_ratio, mantle_sinker_density_difference, **kwargs):
    """
    Create a new parameter (.prm) file based on the original file, 
    modifying specific values related to viscosity and density.

    Args:
        prm_path_ori (str): Path to the original parameter file.
        prm_path (str): Path where the new parameter file will be saved.
        sinker_viscosity (float): Viscosity of the sinker material.
        mantle_sinker_viscosity_ratio (float): Ratio of mantle viscosity to sinker viscosity.
        mantle_sinker_density_difference (float): Density difference (sinker density - mantle density).
        **kwargs: Additional keyword arguments. Accepts:
            - output_directory (str): Directory for output files (default: "output-stokes").

    Raises:
        AssertionError: If the original prm file does not exist or the target directory is invalid.
    """
    assert(os.path.isfile(prm_path_ori))
    assert(os.path.isdir(os.path.dirname(prm_path)))

    output_directory = kwargs.get("output_directory", "output-stokes")

    # Read and modify the original prm file
    with open(prm_path_ori, 'r') as fin:
        inputs = ParseFromDealiiInput(fin)

    inputs["Output directory"] = os.path.join(os.path.dirname(prm_path), output_directory)
    inputs["Material model"]["Simple model"]["Viscosity"] = "%.4e" % sinker_viscosity
    inputs["Material model"]["Simple model"]["Composition viscosity prefactor"] = "%.4e" % mantle_sinker_viscosity_ratio
    inputs["Material model"]["Simple model"]["Density differential for compositional field 1"] = "%.4e" % (-mantle_sinker_density_difference)

    # Write the modified inputs to the new prm file
    with open(prm_path, "w") as fout:
        ParseToDealiiInput(fout, inputs)
    
    assert(os.path.isfile(prm_path))


def RunMainASPECT(aspect_executable, prm_path):
    """
    Execute the ASPECT case using the specified executable and parameter file.

    Args:
        aspect_executable (str): Path to the ASPECT executable file.
        prm_path (str): Path to the parameter (.prm) file.

    Returns:
        CompletedProcess: The result of the subprocess run, including stdout and stderr.
    
    Raises:
        AssertionError: If the expected output format is not matched or if there is any error output.
    """
    completed_process = subprocess.run([aspect_executable, prm_path], capture_output=True, text=True)
    stdout = completed_process.stdout
    stderr = completed_process.stderr

    # Ensure the output contains the expected "Total wallclock" line.
    assert(re.match(".*Total wallclock", stdout.split('\n')[-6]))
    # Verify that there are no errors in stderr.
    assert(stderr == "")

    return completed_process


def AnalyzeVy(case_dir):
    """
    Analyze the sinking velocity (Vy) from the output files in the specified case directory.

    Args:
        case_dir (str): Path to the case directory containing the output files.

    Returns:
        float: The analyzed sinking velocity Vy at the specified location.

    Raises:
        FileExistsError: If the expected output file does not exist.
    """
    fields = ["T", "density", "viscosity", "C_1", "velocity"]
    filein = os.path.join(case_dir, "output-stokes", "solution", "solution-00000.pvtu")

    # Assert that the input file exists
    Utilities.my_assert(os.path.isfile(filein), FileExistsError, "%s: %s doesn't exist" % (Utilities.func_name(), filein))

    # Initialize VTKP and read the output file
    VtkP = VTKP(geometry="box")
    VtkP.ReadFile(filein, quiet=True)

    # Construct polygonal data from the fields
    VtkP.ConstructPolyData(fields, quiet=True)

    # Interpolate to find the velocity at the center point
    target_points_np = np.array([[1.445e6, 1.445e6, 0.0]])
    o_poly_data, points_found, interpolated_data, interpolated_vector = VtkP.InterpolateDomain(target_points_np, fields=fields, quiet=True)

    # Extract and return the positive sinking velocity Vy
    v_y = -interpolated_vector[0][0][1] # Negate to get the positive value
    return v_y


def StokesSinkASPECT2dNp(parent_dir, prm_path_ori, drho, log_eta):
    """
    Wraps the parameter file creation, ASPECT execution, and analysis into a single workflow.
    Wrap the follow function of StokesSinkASPECT2d

    Args:
        parent_dir (str): Directory to store the case results.
        prm_path_ori (str): Path to the original parameter file.
        drho (float): Density difference (mantle - sinker), should be negative.
        log_eta (float): Logarithm of mantle viscosity.

    Returns:
        float: The sinking velocity (Vy) obtained from the analysis.

    Raises:
        OSError: If the case directory cannot be created.
    """
    if type(drho) == np.ndarray and type(log_eta) == np.ndarray:
        assert(drho.shape == log_eta.shape)
        v_y = np.zeros(drho.shape)
        for i in range(drho.size):
            v_y[i] = StokesSinkASPECT2d(parent_dir, prm_path_ori, drho[i], log_eta[i])
    elif type(drho) == np.ndarray and type(log_eta) == float:
        v_y = np.zeros(drho.shape)
        for i in range(drho.size):
            v_y[i] = StokesSinkASPECT2d(parent_dir, prm_path_ori, drho[i], log_eta)
    else:
        v_y = StokesSinkASPECT2d(parent_dir, prm_path_ori, drho, log_eta)
    return v_y


def StokesSinkASPECT2d(parent_dir, prm_path_ori, drho, log_eta):
    """
    Wraps the parameter file creation, ASPECT execution, and analysis into a single workflow.

    Args:
        parent_dir (str): Directory to store the case results.
        prm_path_ori (str): Path to the original parameter file.
        drho (float): Density difference (mantle - sinker), should be negative.
        log_eta (float): Logarithm of mantle viscosity.

    Returns:
        float: The sinking velocity (Vy) obtained from the analysis.

    Raises:
        OSError: If the case directory cannot be created.
    """
    # Create the case directory based on input parameters
    case_dir = os.path.join(parent_dir, "stokes_mcmc_2d_drho_%.1f_eta_%.4f" % (drho, log_eta))
    if not os.path.isdir(case_dir):
        os.mkdir(case_dir)
        
    prm_path = os.path.join(case_dir, "stokes.prm")
    
    # Create the parameter file
    log_sinker_viscosity = 24.0
    mantle_sinker_viscosity_ratio = 10**(log_eta - 24.0) # we assume the body to be 10^24
    MakeNewPrm(prm_path_ori, prm_path, 10**log_sinker_viscosity, mantle_sinker_viscosity_ratio, drho)
    
    # Execute the ASPECT case
    completed_process = RunMainASPECT(aspect_executable, prm_path)
    
    # Analyze the sinking velocity
    v_y = AnalyzeVy(case_dir)
    return v_y


def StokesSinkAnalytic2d(parent_dir, prm_path_ori, drho, log_eta):
    """
    Calculate the analytic sinking velocity (Vy) for a 2D Stokes sinker.

    Args:
        parent_dir (str): Directory for storing case results (not used in calculation).
        prm_path_ori (str): Path to the original parameter file (not used in calculation).
        drho (float): Density difference (mantle - sinker), should be negative.
        log_eta (float): Logarithm of mantle viscosity.

    Returns:
        float: The calculated sinking velocity (Vy) in meters per second.
    """
    # Constants
    year = 365 * 24 * 3600.0
    g = 10.0
    R = 200e3

    # Calculate mantle viscosity
    eta = 10**log_eta

    # Compute sinking velocity using the Stokes flow formula
    v_y = 2.0 / 9.0 * (drho * g * R**2.0) / eta * year
    return v_y


stokes_sinker_approach = 1 # 0 - analytic solution; 1 - aspect

if stokes_sinker_approach == 0:
    mcmc_case_dir = os.path.join(mcmc_dir, "analytic2d")
    StokesSinker = StokesSinkAnalytic2d
elif stokes_sinker_approach == 1:
    mcmc_case_dir = os.path.join(mcmc_dir, "aspect2d")
    StokesSinker = StokesSinkASPECT2dNp
if not os.path.isdir(mcmc_case_dir):
    os.mkdir(mcmc_case_dir)

img_dir = os.path.join(mcmc_case_dir, "img")
if not os.path.isdir(img_dir):
    os.mkdir(img_dir)

aspect_executable = "/home/lochy/Softwares/aspect/build_main/aspect"

# Start timing the execution
start = time.time()
print("Intializing")

# Constants
year = 365 * 24 * 3600.0
g = 10.0
R = 200e3

# Define the path to the original parameter file
prm_path_ori = os.path.join(mcmc_dir, "stokes_0", "stokes.prm")

# Assert that the original parameter file exists
assert(os.path.isfile(prm_path_ori))

# Create a parent directory for the new case results
parent_dir = os.path.join(mcmc_case_dir, "stokes_mcmc_2d")

# If the parent directory already exists, remove it
if os.path.isdir(parent_dir):
    rmtree(parent_dir)

# Create the new parent directory
os.mkdir(parent_dir)

# Initialize the workflow by loading data from a CSV file 
# and visualizing the dataset. This dataset is generated 
# from a series of 2D Stokes sinker cases with a mantle 
# viscosity of 1e22. It includes various density differences 
# between the sinker and the mantle.

# Define the path to the input data file
file_in = os.path.join(RESULT_DIR, "sink_v_data.csv")
assert(os.path.isfile(file_in))
sinkV_data = pd.read_csv(file_in)
print("Readed in dataset %s" % file_in)

# Extract the density difference and sinking velocity columns
drhos = sinkV_data["density_difference"]
sinkV = sinkV_data["sinking_velocity"] # take positive values
print("length of data:", len(drhos))

# Create a plot for the sinking velocity data
fig, ax = plt.subplots()
ax.plot(drhos, sinkV, label="perturbation")
ax.plot(drhos, StokesSinker(parent_dir, prm_path_ori, drhos.to_numpy(), 22.0), label="solution, 1e22")
ax.plot(drhos, StokesSinker(parent_dir, prm_path_ori, drhos.to_numpy(), 26.0), label="solution, 1e26")
ax.set_xlabel("Density Difference (kg/m^3)")
ax.set_ylabel("Sinking Velocity (m/yr)")
fig_path = os.path.join(img_dir, "StokesSinkingData.png")
fig.savefig(fig_path)
print("Saved figure: %s" % fig_path)
end = time.time()
print("Intializing takes %.2f" % (end-start))
start = end
print("\n")

# Define the PyMC model context
print("Build PyMC model")
with pm.Model() as aspect_model:
    # Define prior distribution
    log_eta = pm.Normal("log_eta", mu=21.0, sigma=3.0)  # Intercept prior
    sigmaV = pm.HalfNormal("sigmaV", sigma=0.01) # Standard deviation prior

    # Approach: use deterministic variable
    # This section handles two tasks:
    # 1. Extract each of the density difference inputs (drhos) and flatten the corresponding Tensor log_eta 
    #    to pass only scalar values to the StokesSinkASPECT2d function.
    # 2. Collect the outputs from this function into a list, then stack it into a Tensor for PyMC.

    # Uncomment the following section to use 2d analytical solution approach
    if stokes_sinker_approach == 0:
        print("Using 2d analytical solution approach")
        sinkV_model = pm.Deterministic("sinkV_model",
            pm.math.stack([StokesSinkAnalytic2d(parent_dir, prm_path_ori, drhos[i], log_eta) for i in range(len(drhos))]))

    # Note: The shape of sinkV_model will not be printed as expected because it hasn't been sampled yet.
    # print("shape:", sinkV_model.shape)

    # Uncomment the following section to use ASPECT in a black box approach
    if stokes_sinker_approach == 1:
        print("using ASPECT in a black box approach")
        sinkV_model = pm.Deterministic("sinkV_model",
            pm.math.stack([StokesSinkASPECT2d(parent_dir, prm_path_ori, drhos[i], log_eta.eval()) for i in range(len(drhos))]))

    # Define the likelihood of the observed sinking velocity
    sinkV_obs = pm.Normal("sinkV_obs", mu=sinkV_model, sigma=sigmaV, observed=sinkV)

# Record the end time for building the model
end = time.time()
print("Build PyMC model takes %.2f" % (end-start))
start = end # Reset the start time for the next operation
print("\n")

# Define the context for the PyMC model
print("MCMC")
with aspect_model:
    # Perform MCMC sampling using the NUTS (No-U-Turn Sampler) algorithm
    # `return_inferencedata=True` returns the results as an InferenceData object,
    # which is useful for compatibility with ArviZ for analysis and visualization.
    # Note: If you want to use the trace as a dictionary, set `return_inferencedata=False`,
    # but this may not work seamlessly with summary functions.
    trace = pm.sample(1000, return_inferencedata=True)
end = time.time()
print("MCMC takes %.2f" % (end-start))
start = end
print("\n")

print("Visualize Summary and Trace")
# Summarize the trace and posterior distributions using ArviZ
# Ensure that `return_inferencedata=False` is not set, as we are using the trace directly

# Save the summary figure
az.summary(trace)
fig_path = os.path.join(img_dir, "AzSummary.png")
plt.savefig(fig_path)
print("Saved figure: %s" % fig_path)

# Plot the posterior distributions for all parameters
az.plot_posterior(trace)
fig_path = os.path.join(img_dir, "AzPosterior.png")
plt.savefig(fig_path)
print("Saved figure: %s" % fig_path)

# Plot the trace for the specified variable (log_eta)
az.plot_trace(data=trace, var_names=["log_eta"], figsize=(16, 4))
fig_path = os.path.join(img_dir, "AzTrace.png")
plt.savefig(fig_path)
print("Saved figure: %s" % fig_path)

# Plot the posterior for the specific variable (log_eta)
az.plot_posterior(data=trace, var_names=["log_eta"], figsize=(16, 4))
fig_path = os.path.join(img_dir, "AzPosteriorLogEta.png")
plt.savefig(fig_path)
print("Saved figure: %s" % fig_path)

# Plot the autocorrelation for the specified variable (log_eta)
az.plot_autocorr(data=trace, var_names=["log_eta"], figsize=(16, 4))
fig_path = os.path.join(img_dir, "AzAutoCorr.png")
plt.savefig(fig_path)
print("Saved figure: %s" % fig_path)
end = time.time()
print("Visualize Summary and Trace takes %.2f" % (end-start))
start = end
print("\n")

print("Visualize Model Prediction")
# Extract samples from the trace for log_eta and sigmaV, excluding the first 800 samples for burn-in
log_eta_samples = np.concatenate(trace.posterior.log_eta.data[:,800:])[:, None]
sigmaV_samples = np.concatenate(trace.posterior.sigmaV.data[:,800:])[:, None]

# Define a range of density differences for predictions
drho_t = np.arange(drhos.min()-5.0, drhos.max()+5.0, 1.0)[:, None]

if stokes_sinker_approach == 0:
    # Calculate the sinking velocity predictions using the analytical model
    sinkV_t = StokesSinker(parent_dir, prm_path_ori, drho_t.T, log_eta_samples)

    # Compute the mean sinking velocity and the 95% credible interval (quantiles)
    mean_sinkV_t = sinkV_t.mean(axis=0)
    qs = mquantiles(sinkV_t, [0.025, 0.975], axis=0)

    # Plot the synthetic data and fill the credible interval
    figsize(12.5, 4)
    plt.plot(drhos, sinkV, label="synthetic data")
    plt.fill_between(drho_t[:, 0], *qs, alpha=0.7, color="#7A68A6")
    plt.plot(drho_t, mean_sinkV_t, lw=1, ls="--", color='k', label="average posterior \nsinking \
    velocity")
    fig_path = os.path.join(img_dir, "MCMC_Prediction.png")
    plt.savefig(fig_path)
    print("Saved figure: %s" % fig_path)
    end = time.time()
    print("Visualize Model Prediction takes %.2f" % (end-start))
else:
    # For the aspect black box approach, this needs to be modified
    pass