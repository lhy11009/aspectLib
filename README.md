[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lhy11009/aspectLib/HEAD)

# Logistics

I store all my cases in a structured file systems and wrote my scripts in a hierachical package. 

## Structured file systems

## Hierachical layout of packages

On the first order, there are scripts below 'python_scripts(shilofue)' that should be applicable for all the models I run.

On the second order, there are subfolders under the 'python_scripts' folder (e.g. TwoDSubduction0). Those contains scripts for that specific program.

In the scripts, I define "parental" classes and class functions in the first order scripts and make "daughter" classes that inherit there parents in the second order scripts.

A good example could be seen in the script of (todo).py

(todo) use route like '.'


# Topics

## Run on server

The basic idea is to only use bash to nevigate around the package managing issue of python.

### Command alias

### Case managements

The entire workflow here is:
* read configuration from a json file
* periodically check the running case to update the run time file
* perform operations like restart / terminate/ download a running case when needed.
* put up a report since last check

#### Reading configuration from a json file

script: 
	bash_scripts/parse_case.sh
	utilities/bash_scripts/utilities.sh
	utilities/bash_scripts/JSON.sh


#### Save case runtime info in a file

script: 
	bash_scripts/parse_case.sh
	bash_scripts/awk_states/parse_block_output
The second script is used to parse block output with awk. The functions in the first script simply takes care of the output.

Here, I put case run time info into a separte file "run_time.log",

The information I include:
* case name
* state: state of the case (running, stopped, terminated) (todo)
* step : step in model
* last_restart_step : last step it restarted
* time : time in model
* last_restart_time : last time in model it restarted
* wallclock : wallclock after last restart
* total_wallclock : total wallclock
* last_update: time of last update (todo)

Format: 

* header: 
	# 1 step
	# 2 last_restart_step
	...
* data:
	1000	600 ...

Operations:
* a. Parse runtime output for a single case
* b. Parse runtime output for multiple cases within the same directory.
* c. Assemble this information projectwise by combine these outputs (todo)

run with
	Lib_parse_case -h for more details

##### get the state of a case

determine wheter a case is running (R), stopped (S), terminated (T) or having issue (I)
* R: This case runtime information is found in slurm
* S: This case is stoped and the end time is not reached
* T: This case is terminated manually
* E: This case has run to an end.
* I: This case may run into an issue

#### Restart a previous case (to be tested)

script: 
	bash_scripts/parse_case.sh

Operations:
* a. restart a case if it doesn't reach the end time I set here (e.g. due to a maintainance of the server).
* b. loop in a directory and restart cases where "End time" is not reached (todo)


## Prepare cases

### Create new case

I put the configurations in a json file and wrote a script in each project to create new cases

These options are read into a CASE_OPT class, While the operations are defined in a CASE class.
In all the operations, I took the approach of first importing a prm file and then substitute the entries.

Note: I defined a clase CASE here for preparing and creating cases, while I have another class named CASE_OPTIONS
to deal with post-process.

An example usage (TwoDSubduction)
	Lib_TwoDSubduction0_Cases create_with_json -j \
        /home/lochy/ASPECT_PROJECT/TwoDSubduction/wb_create_test/configure_1.json \n\

(todo) fix the change_plate_ages test, there is still a temperature discontinuity in the prescribing-temperature corner region.


### Create new group

Here, the task is to create mutiple cases at a time and hold them in a single folder.

I made use of the previous script I used for create a single case (Cases.py), and write a new script Group.py.
In this new scripts, I read in a json file with the configurations and write a bundle of json case-configuring files for interfaces defined in Cases.py.

## Cooperate with other models / analytic results / laws

### Cooperate with Hefesto

I wrote a python script "PostHefesto.py" 

The job it does is: 

Check file format (for the first and second entries, I.e. P and T) 

Convert Hefesto file to perplex table 

Resample in terms of P and T 

Fix other fields if they are missing (but not P and T) 

This is a test called test_post_hefesto.py in the integration folder. Look at that to see how this work.

### Work with flow laws



## Post-process

Clearify the usage of prm:

I use case.prm to create cases as well as analytical analysis, while using output/original.prm for post-process.

### Interfaces to substitute words in a script (including a usage for visit)

A CaseOptions class is defined. There, I stored information in a dictionary, then I use the values to substitute
the keys in a script.

#### Automize visit plotting

Software like visit and paraview could be automized by using python scripts.
While there could be the issue of the software only well support python2, this problem of using different 
version of python could be bridged if I use raw python2 scripts with key word substitution and only worki
on the values of these keys in python3.

#### Automize vtk post-processing

This follows similar logic as before with the exception of this being a cpp package.


### Linear plots from fields in Aspect's outputs

For this, I defined a LINEARPLOT class to read in the default file format in aspect.
The format of the file looks like:
	# 1: Time step number
	# 2: Time (years)
	# 3: Time step size (years)
	0.0 0.0 0.0
	...
These also work as a standard interface to parse data from these files


#### Plot results in "statistic" file

(todo) Change the length of number label on axis, currently they are too long so parts are blocked.

There are interfaces for "GetStep", "GetTime" as well

#### Plot Newton solver results

* Parse the outputs from the stdout file (bash_scripts/awk_states/parse_block_newton)
* Plot the outputs

There is option to plot a single step and plot a combined results by combine the steps along the x-axis.
See the example in section 2 in Menno's thesis.

#### Plot time run

* Parse the outputs from the stdout file (bash_scripts/awk_states/parse_block_output)
* Plot the outputs

### Combine figures

Combine plots in different cases
Some notes:
* font size in case names are adjusted to fill the individual boxes

### Generate a journal of the pictures plot within the dya
(todo)

### Using the Image module to prepare results

The task to handle here is combining results from visit, matplotlib and use the right frames that fit.

I implement an IMAGE_OPT module in the "Utilities.py" script which takes a json file and combines results accordingly
For detail, see the test "test_img_operation.py"

For this to work, I need to generate a tranperant frame in AI, and merge that on top of a plot from visit

I also defined a PREPARE_RESULTS class in the file PlotCombine.py for this. This class would take a template (in this case a json file)
and substitute the key words with values (e.g. The directory to saved images.)

For an example usage, see:
	/home/lochy/ASPECT_PROJECT/aspectLib/shilofue/TwoDSubduction0/PlotCase.py


#### work flow of generating the frame

* open the raw figure from visit in Adobe Illustrator.
* pin and paint the frame, open a new file with it.
* Lay a bigger while rectangle and crop out what's inside the frame.
* Change background to transparent (view->show transperant grid)
* Save using the "Export->Save for web" option, choose png24 as format.

### Genearate animation

* loop for steps and use the PIL module to prepare results
* Use the imageio module to generate animation
* There is option to generate animation for one case or for cases in a directory

For example, see:
	/home/lochy/ASPECT_PROJECT/aspectLib/shilofue/TwoDSubduction0/PlotCase.py


## work with prm

(todo) put utilitie functions in the ParsePrm.py, and put interfaces into the Case.py

## Case management

### pinpoint snapshots(restarting file)

(todo) add solution.00000.*.vtu to the file list

# Projects - TwoDSubduction

## Cases

Related script:

	TwoDSubduction/Cases.py

An interfaces for creating new cases.

In this script, I inherit from the parental class (Case.py) and create one for operation and one for handling and passing around variables.

	class CASE(CasesP.CASE):

	class CASE_OPT(CasesP.CASE_OPT):

In this first one, I basically coded sections in the prm file to substitute a existing prm file.
These operations are done in two different ways:
* Substitute by key word and value.
* Rewrite a big chunk of settings directly.
The first one is the standard, while the second one is needed for changing some major setting (sph vs cart).

The idea of working like this are:
* Be able to reproduce all the major settings
* Note down all the major steps that I don't rework them (unless some specific tests are needed.).

### Code logistic~/ASPECT_PROJECT/TwoDSubduction/EBA_CDPTs

With the "class CASE", there are two major operation, namely:

	configure_prm()

and

	configure_wb()

"configure_prm" contains these few sections:
* geometry (i.e. chunk or box)
* material model (e.g. shear zone, yielding, upper to lower mantle)
(todo, write this cleanly like above, i.e, put everything in a "configure geometry function")

## Post process

### trench motion

I use the vtk package in c++ to anaylyze trench motion here:
* a, convert results into inputs for the vtk script
* b, run the vtk script and loop for steps. Combine results into a new file
* c, plot the results in this file.