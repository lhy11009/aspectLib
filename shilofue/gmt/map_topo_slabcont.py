# to use this first: # > conda activate pygmt import os import pygmt from pathlib import Path
import os
import pygmt
from pathlib import Path

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
# This is a list of names of subduction zones
loc = ['alu','cal','cam','car','cas','cot','hal','hel','him','hin','izu',
	'ker','kur','mak','man','mue','pam',
    'phi','png','puy','ryu','sam','sco','sol','sul','sum','van'];
# Here, we use pygmt to create a figure
fig = pygmt.Figure()
pygmt.config(FONT='Times-Roman')
pygmt.config(FONT_LABEL='Times-Roman,12p')
# Initiate a basemap within that figure
fig.basemap(region="g", projection="W15c", frame=True)
grid = pygmt.datasets.load_earth_relief(resolution="10m")
fig.grdimage(grid=grid,cmap="gray")
fig.coast(shorelines="1/0.5p,black")
# make cpt(color) file
cptfile = os.path.join(ASPECT_LAB_DIR, "shilofue", "gmt", "slabdepth.cpt")
assert(os.path.isfile(cptfile))
# pygmt.makecpt(cmap="buda",series=(-700,0,100), color_model="+c",output=cptfile)
# pygmt.makecpt(cmap="batlow",series=(-700,0,100), color_model="+c-700-0",output=cptfile)
# plot the slab 2d data set
slab2dir = os.path.join(ASPECT_LAB_DIR, "large_data_files", "slab2", "Slab2Distribute_Mar2018")
for i in range(len(loc)):
    #print(loc[i])
    # use path to look for grd file, this is implemented because there is an appendix in the file name
    grdfile_name_base = loc[i] + '_slab2_dep'
    pathlist = Path(slab2dir).rglob(grdfile_name_base + '*'+ '.grd')
    for path in pathlist:
        grdfile = str(path)
        print(grdfile)
    # grdfile = os.path.join(slab2dir, loc[i] + '_slab2_dep.grd')
    # fig.grdimage(grid=grdfile, nan_transparent="gray", cmap='buda')
    fig.grdimage(grid=grdfile, nan_transparent="gray", cmap=cptfile)
# fig.colorbar(cmap=cptfile,position="JMR",box=False,frame=["x+lDepth", "y+lkm"])
fig.savefig('map_topo_slabs.png')
