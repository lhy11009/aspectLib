## aspectLib

This is a repository of my scripts post-processing ASPECT

### Python scripts in shilofue folder

#### Input Options

Output file type, png or pdf:

    -ft, --file_type

#### shilofue/PlotDepthAverage.py

Plot the buoyancy and the buoyancy ratio (buoyancy / density in the second file)
with two depth_average files, the first will be taken as the reference profile.

    python -m shilofue.PlotDepthAverage plot_buoyancy\
    -i {case0}/output/depth_average.txt \
    -i1 {case0}/output/depth_average.txt

Plot the density change between two profiles and compared to HeFesto Output.
The path of the files are not hard coded into the "CompareHefestoBuoyancy" function.

    python -m shilofue.PlotDepthAverage compare_hefesto_buoyancy

#### shilofue/PostHefesto.py

Plot the buoyancy and density ratio from two HeFESTo profiles

    python -m shilofue.PostHefesto plot_hefesto_buoyancy -i hefesto_buoyancy_input

input file as an example:

    files/configure_files/hefesto_buoyancy_input

Create hefesto cases from json files

    python -m shilofue.PostHefesto create_case -i `pwd`/config_hefesto.json

json file as an example

    files/json_examples/config_hefesto_case.json

Assemble the parallel files. Do this after getting the results from the server

    python -m shilofue.PostHefesto assemble_parallel_files -i /home/lochy/ASPECT_PROJECT/HeFESTo_DIR/test_hefesto_parallel1

#### shilofue/Rheology.py

Plot summary of mantle rheology from multiple rheologies:

    python -m shilofue.Rheology compare_mantle_rheology