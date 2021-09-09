#!/Users/billen/opt/miniconda3/bin/python

from cmcrameri import cm
import matplotlib.pyplot as plt
import numpy as np

datadir = '/Users/billen/Box-Sync/BA18-Mod4/model-files/'
filepre = 'mt17_25'
ts = 4000
tstep = str(ts)
print(tstep)

datafile = datadir + filepre + '.cap00.' + tstep + '_slicedata'

# lat lon radius vlat vlon vrad T eta
dat = np.genfromtxt(datafile)
# crust, harz, cont-crust, phase-transitions
#dat1 = np.genfromtxt("/Users/billen/Box-Sync/BA18-Mod4/model-files/mt17_25.opt00.3000_slicedata")


rad2m = (180./np.pi)*111.11*1000;  # radians to meters (at equator)
ms2cmyr = 100*3600*24*376.25;   #  m/s to cm/yr

# Constants
kappa = 1e-6;  # m^2/s  thermal diffusivity
Re = 6371137;  # m, Earth radius
eps0 = kappa/(Re*Re);  # 1/s  (strain-rate)
v0 = (kappa/Re);  # m/s, (velocity)
eta0 = 1e20;  # Pa-s, (viscosity)
To = 1400; # C (temperature)

nodey = 1153;
nodez = 513;

lon = (np.array([x[1] for x in dat]).reshape(-1, nodez)*rad2m).transpose()
rad = (np.array([x[2] for x in dat]).reshape(-1, nodez)*Re).transpose()
T = (np.array([x[6] for x in dat]).reshape(-1, nodez)*To).transpose()
visc = (np.array([x[7] for x in dat]).reshape(-1, nodez)*eta0).transpose()

xmin = 3750;
xmax = 4750;
ymin = 0;
ymax = 1000;

fig = plt.figure()
ax1 = plt.subplot(1,1,1)

#c1 = ax1.pcolormesh(lon/1000, (Re-rad)/1000, np.log10(visc), cmap = cm.grayC,shading='auto', rasterized=True)
c1 = ax1.pcolormesh(lon/1000, (Re-rad)/1000, np.log10(visc), cmap = cm.davos_r,shading='auto', rasterized=True)
#c1.set_clim(19,24)
#ax1.contour(lon/1000, (Re-rad)/1000, T, [1000], colors='black')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymax, ymin)
ax1.set_xlabel('Distance (km)');
ax1.set_ylabel('Depth (km)');
ax1.set_title('Viscosity (Pa s)');
ax1.set_aspect('equal', adjustable='box')
fig.colorbar(c1, ax=ax1) 
#ax1.set_clim(0, 1600)  (this not working 
plt.tight_layout()

pdffile = filepre + '_' + tstep + '_visc_davos.pdf'
fig.savefig(pdffile,bbox_inches='tight', dpi=150)
plt.close(fig)

print('here')
