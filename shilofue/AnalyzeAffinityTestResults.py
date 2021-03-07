r"""Analyze affinity test results

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m shilofue.AnalyzeAffinityTestResults analyze_affinity_test_results
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/rene_affinity_test/results/spherical_shell_expensive_solver/peloton-ii-32tasks-core-openmpi-4.0.1

To this end, this scripts required an organized results of affinity tests to run. Upon this, this scripts pull out the time results
and plot figures.
""" 
import numpy as np
import sys, os, argparse
import pathlib, subprocess
import numpy as np
from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']

# analyze test result
# todo
def analyze_affinity_test_results(test_results_dir):
    '''
    analyze affinity test results
    '''
    total_wall_clock = []
    assemble_stokes_system = []
    solve_stokes_system = []
    cores = []
    resolutions = []
    setups = []
    # go into sub dirs
    temp_file = os.path.join(ASPECT_LAB_DIR, 'temp')  # file to save partial results
    path_obj = pathlib.Path(test_results_dir).rglob("output*")
    for _path in path_obj:
        output_file = str(_path)
        output_path = os.path.join(test_results_dir, output_file)
        patterns = output_file.split('_')
        # append data
        subprocess.run("%s/bash_scripts/parse_block_output.sh  analyze_affinity_test_results %s %s" 
                       % (ASPECT_LAB_DIR, output_file, temp_file), shell=True)
        try:
            data = np.genfromtxt(temp_file)
            total_wall_clock.append(data[0, -1])
            assemble_stokes_system.append(data[1, -1])
            solve_stokes_system.append(data[2, -1])
        except Exception:
            pass
        else:
            setups.append(int(patterns[-1]))
            resolutions.append(int(patterns[-2]))
            cores.append(int(patterns[-3]))

    setups = np.array(setups)
    resolutions = np.array(resolutions)
    cores = np.array(cores)
    total_wall_clock = np.array(total_wall_clock)
    assemble_stokes_system = np.array(assemble_stokes_system)
    solve_stokes_system = np.array(solve_stokes_system)
    # rearrange and sort data
    print("Affinity Test Results:")
    print('setups:', setups)
    print('resolutions', resolutions)
    print("cores:", cores)
    print("total wall clocks:", total_wall_clock)
    print("assemble_stokes_system:", assemble_stokes_system)
    print("solve_stokes_system:", solve_stokes_system)
    
    # plot via matplotlib
    odir = os.path.dirname(output_file)  # save in the save dir as the output file

    resolution_options = []
    for resolution in resolutions:
        if resolution not in resolution_options:
            resolution_options.append(resolution)

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    for i in range(len(resolution_options)):
        resolution = resolution_options[i]
        plot_indexes = (resolutions == resolution) 
        xs = cores[plot_indexes]
        sequential_indexes = np.argsort(xs)
        ys1 = total_wall_clock[plot_indexes]
        axs[0].plot(xs[sequential_indexes], ys1[sequential_indexes], ".-", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label='Total Wall Clock(resolution=%d)' % resolution)
        ys2 = assemble_stokes_system[plot_indexes]
        axs[0].plot(xs[sequential_indexes], ys2[sequential_indexes], ".--", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label='Assemble Stokes System(resolution=%d)' % resolution)
        ys3 = ys2 / ys1
        axs[1].plot(xs[sequential_indexes], ys3[sequential_indexes], ".--", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label='Assemble Stokes System(resolution=%d)' % resolution)
        ys4 = solve_stokes_system[plot_indexes]
        axs[0].plot(xs[sequential_indexes], ys4[sequential_indexes], ".-.", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label='Solve Stokes System(resolution=%d)' % resolution)
        ys5 = ys4 / ys1
        axs[1].plot(xs[sequential_indexes], ys5[sequential_indexes], ".-.", color=cm.gist_rainbow(1.0*i/len(resolution_options)),
                label='Solve Stokes System(resolution=%d)' % resolution)
    axs[0].set_xlabel('Cores')
    axs[0].set_ylabel('Time [s]')
    axs[0].grid()
    axs[0].set_title('Wall Clock')
    axs[0].legend()
    axs[1].set_xlabel('Cores')
    axs[1].set_ylabel('Percentage')
    axs[1].grid()
    axs[1].set_title('Percentage of Each Part')
    fig.tight_layout()
    filepath='%s/job_time.png' % odir
    print("output file generated: ", filepath)
    plt.savefig(filepath)

    pass

def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    _commend = sys.argv[1]
    # parse options
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'analyze_affinity_test_results':
        # example:
        # python -m shilofue.AnalyzeAffinityTestResults analyze_affinity_test_results
        # -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/rene_affinity_test/results/spherical_shell_expensive_solver/peloton-ii-32tasks-core-openmpi-4.0.1
        analyze_affinity_test_results(arg.inputs)

# run script
if __name__ == '__main__':
    main()