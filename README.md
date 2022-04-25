[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lhy11009/aspectLib/HEAD)

# Logistics

I store all my cases in a structured file systems and wrote my scripts in a hierachical package. 

## Structured file systems

### Hierachical layout of packages

On the first order, there are scripts below 'python_scripts(shilofue)' that should be applicable for all the models I run.

On the second order, there are subfolders under the 'python_scripts' folder (e.g. TwoDSubduction0). Those contains scripts for that specific program.

In the scripts, I define "parental" classes and class functions in the first order scripts and make "daughter" classes that inherit there parents in the second order scripts.

A good example could be seen in the script of (todo).py

(todo) use route like '.'

### Archived files under the "files" folder

This is the folder where I put all the acrhived files.

* "Project" folder : include all python scripts to initiate a new project with.


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

#### Restart a previous case

script: 
	bash_scripts/parse_case.sh
	(dependence) awk_states/parse_log_last_step: parse output at the last time step
        (dependence) awk_states/parse_snapshot: parse information of the last snapshot saved

Operations:
* a. restart a case if it doesn't reach the end time I set here (e.g. due to a maintainance of the server).
* b. loop in a directory and restart cases where "End time" is not reached

Note:

If I want run a case to 60e6, I'll use 59e6 here.

Example:

	Lib_parse_case check_time_restart `pwd` 59e6 high2


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

#### Explain the json file inputs:

Below is an example of json file used for creating cases, one can always run this command to see the description for each entry.

	Lib_TwoDSubduction0_Case -h

Every option is assigned as key and value.
Together they are related with either a variable in the prm file or a defined function to change multiple parameters at the same time.
For example, the "include fast first step" option would actually copy the prm file, change the option of "End time" to 0.0, solver scheme
to "No advection, No stokes" and dump this back to a file called "case_f.prm".
The goal of this approach is to reduce the process of model setup to what we would mental address the task.
Say, when we think of changing the geometry, normally we have to vary a lot of entries in the prm file,
e.g. the function you are using for your initial temperature field, entries in your world builder file.
But with this, we reduce it back to just change the option of geometry from box to chunk.

	{
		"base directory": "${ASPECT_LAB_DIR}/files/TwoDSubduction/210103", 
		"output directory": "${ASPECT_LAB_DIR}/.test/TwoDSubduction_cases",
		"name": "peierls_test",
		"include fast first step": 1,
		"geometry": "chunk",
		"potential temperature": 1573.0,
		"boundary condition": {"model": "all free slip"},
		"use world builder": 1, 
		"world builder": 
		{
			"plate age method": "adjust box width",
			"subducting plate": 
			{
				"age trench": 40e6, 
				"sp rate": 0.05
			},
			"overiding plate": 
			{
				"age": 40e6, 
				"transit": {"age": 20e6, "length": 700e3}
			}
		},
		"Include peierls creep": 0,
		"additional files": ["job_rome.sh", "job_high2.sh", "job_p-billen.sh"]
	}


### Create new group

Here, the task is to create mutiple cases at a time and hold them in a single folder.

I made use of the previous script I used for create a single case (Cases.py), and write a new script Group.py.
In this new scripts, I read in a json file with the configurations and write a bundle of json case-configuring files for interfaces defined in Cases.py.

#### workflow

1. test new feature with a "case.json".
   1. implement this option in the Cases.py script
   2. put option into the case.json file
2. Generate a group of cases
   1. save the previous "case.json", as well as the ".prm" file, ".wb" file to a separate folder under "files" (name with date)
   2. edit a group.json, manage both "base features" and "features". The first one is universal settings for the whole group, while the second one create differences among cases.
   3. run with "Lib_foo_group create_group ...", make sure you manage the pop-up options correctly (whether to delete, whether to update.)

#### Explain the json file inputs:

First, this json file is dependent on the previous json file in the section of "Create new case".

Here is an example:
	{
	"base name": "eba_cdpt_cart", 
	"base json": "${ASPECT_LAB_DIR}/files/TwoDSubduction/220109/case.json",
	"base directory": "/home/lochy/ASPECT_PROJECT/aspectLib/files/TwoDSubduction/220109",
	"output directory": "${TwoDSubduction_DIR}/EBA_CDPT_cart1",
	"base features":[
		{
		"name": "Geometry of the model",
		"key": ["geometry"],
		"unit": "",
		"values": ["box"],
		"abbreviating strings": [""]
		}
	],
	"features":[
		{
		"name": "Age of the subducting plate",
		"unit": "yr",
		"key": ["world builder", "subducting plate", "age trench"],
		"values": [40e6, 80e6],
		"abbreviation by value": 1,
		"abbreviating value options": ["SA", 1e-6]
		},		
		{
		"name": "Age of the overiding plate",
		"key": ["world builder", "overiding plate", "age"],
		"unit": "yr",
		"values": [20e6, 40e6],
		"abbreviation by value": 1,
		"abbreviating value options": ["OA", 1e-6]
		},
		{
		"name": "coupling the eclogite phase to shear zone viscosity",
		"key": ["coupling the eclogite phase to shear zone viscosity"],
		"unit": "",
		"values": [0, 1],
		"abbreviating strings": ["", "CpEcl"],
		"if abbreviating": [0, 1]
		},
		{
		"name": "Width of the box",
		"key": ["world builder", "box width before adjusting"],
		"unit": "km",
		"values": [6.783e6, 8.896e6],
		"abbreviating strings": ["width61", "width80"],
		"if abbreviating": [0, 1]
		}
	],
	"bindings": [[0, 0, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 0, 1]]
	}


Same as before, one could run with this command to get a full detail:
	
	Lib_TwoDSubduction0_Group -h

Note that the main difference from the json file for a single case is we add in the option of 
"base features" and "features" in this one.
They are nested options as of options in the case of "Create new case".
Both of these are lists to start with, so each entry in the list is a separate option.
For example, the 2nd entry in the "features": 
		
		{
		"name": "Age of the subducting plate",
		"unit": "yr",
		"key": ["world builder", "subducting plate", "age trench"],
		"values": [40e6, 80e6],
		"abbreviation by value": 1,
		"abbreviating value options": ["SA", 1e-6]
		},		

Here, keys are a listed, so this is related to a series of nested keys in a json file for "Create new case".
Note this part in the previous section:
		
		"world builder": 
		{
			"plate age method": "adjust box width",
			"subducting plate": 
			{
				"age trench": 40e6, 
				"sp rate": 0.05
			},
			"overiding plate": 
			{
				"age": 40e6, 
				"transit": {"age": 20e6, "length": 700e3}
			}
		}

"values" is also a list. This means we could assign different value to a variable,
in this case, the age of the overiding plate, and we will create a new case for each of 
these different values.

Then I'll discuss how the abbreviation options work.
To start with, it's no more than add a string in the case name, for example:

	"abbreviating value options": ["SA", 1e-6]

This means we will append a string of "SA" + "something" to our case name.
While "something" is the "value" * 1e-6, in the case, 40e6 * 1e-6 = 40 and 80e6 * 1e-6 = 80, specifically.
So you'll expect your case names to have "SA40" and "SA80" in them.

For the "base feature", it's a special type of feature that we want to apply on the whole group to start with.
That's why there need to be just one value in the list, and the abbreviating option is always like this:

	"abbreviating strings": [""]

Because we don't want anything to be appended to the case name in this special occasion.


## Cooperate with other models / analytic results / laws

### Cooperate with CDPT phase transition models

Script: Phase Transition.py

* print phase-transition outputs
* print entropy changes on phase transitions

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

The package is depending of the vtk package (version > 9.0.1), when compiling this package, one need to provide
the path to the vtk package.

	cmake -DVTK_DIR:PATH=/home/lochy/VTK-9.0.1/build ..

the modules in the package could be built separately, e.g

	make TwoDSubduction_SlabAnalysis

To run this module, one need also to provide an input file with parameters.
This file should have contents like:

	(contents of TwoDSubduction_SlabAnalysis.input_s000154)
	/home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT1/eba_cdpt_SA80.0_OA40.0/output/solution/solution-00166.pvtu
	/home/lochy/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT1/eba_cdpt_SA80.0_OA40.0/temp_output/depth_average_output_s1334
	0.6278449223041909

Run with:

	./TwoDSubduction_SlabAnalysis TwoDSubduction_SlabAnalysis.input_s000154


### Linear plots from fields in Aspect's outputs

For this, I defined a LINEARPLOT class to read in the default file format in aspect.
The format of the file looks like:
	# 1: Time step number
	# 2: Time (years)
	# 3: Time step size (years)
	0.0 0.0 0.0
	...
These also work as a standard interface to parse data from these files
shilofue/PhaseTransition.py

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

## Schedule

I have a class SCHEDULER in Scheduler.py to schedule task to run later.
This would put commands to run into a file and save that file(e.g. sync data, plot, generate animation).
I could then choose what I want to schedule based on the project I am working on.
See an example with:
	Lib_TwoDSubduction -h

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

## Archived Files

### Phase transitions configuration

* files/TwoDSubduction/phases_1_0.json - the CDPT model I use.

## a small sub-project - testing the overshoot of temperature in the latent heat bench mark

This is a sub-project I named "LatentHeatBK" where I tested the overshooting of temperature.
As of Feb 2022, I have assembled the "Cases.py", where I could change the mesh resolution, velocity, as well as which phase transition model to use.
