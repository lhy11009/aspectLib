## aspectLib

This is a repository of my scripts post-processing ASPECT

### Python scripts in shilofue folder

#### shilofue/PlotDepthAverage.py

Plot the buoyancy and the buoyancy ratio (buoyancy / density in the second file)
with two depth_average files, the first will be taken as the reference profile.

    python -m shilofue.PlotDepthAverage plot_buoyancy\
    -i {case0}/output/depth_average.txt \
    -i1 {case0}/output/depth_average.txt
