[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lhy11009/aspectLib/HEAD)

# Logistics

I store all my cases in a structured file systems and wrote my scripts in a hierachical package. 

## Structured file systems

## Hierachical layout of packages

On the first order, there are scripts below 'python_scripts(shilofue)' that should be applicable for all the models I run.

On the second order, there are subfolders under the 'python_scripts' folder (e.g. TwoDSubduction0). Those contains scripts for that specific program.

In the scripts, I define "parental" classes and class functions in the first order scripts and make "daughter" classes that inherit there parents in the second order scripts.
A good example could be seen in the script of (todo).py
