#!/usr/local/bin/python
#
# Lacerda@Corrego - 24/Mar/2016
#
import sys
import numpy as np
import tables as tbl
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.objects import runstats
from CALIFAUtils.scripts import ma_mask_xyz
from CALIFAUtils.plots import plot_text_ax, density_contour

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

def plotmuZR(ax, mu, Z, z = None, cmap = 'viridis', mask = None, labels = True, zlabel = None, N = False):
    if mask is None: mask = np.zeros_like(mu, dtype = np.bool_)
    extent = [0, 5, -2, 0.4]
    plt.axis(extent)
    bins = [30, 30]
    if z is None:
        xm, ym = ma_mask_xyz(np.ma.log10(mu), Z, mask = mask)
        density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range = [extent[0:2], extent[2:4]], colors = [ 'b', 'y', 'r' ])
        sc = ax.scatter(xm, ym, marker = 'o', c = '0.5', s = 10, edgecolor = 'none', alpha = 0.3)
    else:
        xm, ym, zm = ma_mask_xyz(np.ma.log10(mu), Z, z, mask = mask)
        sc = ax.scatter(xm, ym, c = zm, cmap = cmap, marker = 'o', s = 30, edgecolor = 'none')#, alpha = 0.6)
        cb = plt.colorbar(sc)#, ticks = [0, .5, 1, 1.5, 2])
        cb.set_label(zlabel)
    if labels:
        ax.set_xlabel(r'$\langle\log\ \mu_\star\rangle_R$ [$M_\odot pc^{-2}$]')
        ax.set_ylabel(r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]')
    if not N:
        N = xm.count()
    plot_text_ax(ax, '%d' % N, 0.01, 0.99, 20, 'top', 'left', 'k')
    return ax

h5fname = sys.argv[1]
SFRgroup = '/SFR05050525'
Zgroup = '/Zstar05'
fpref = '.png'

iT = 1
h5file = tbl.open_file(h5fname, mode = 'r')
tbl_gals = h5file.root.pycasso.main
tbl_zones = h5file.root.pycasso.zones
tbl_integrated = h5file.root.pycasso.integrated
node_SFR = h5file.get_node(SFRgroup)
tSF__T = node_SFR.tSF.cols.age[:]
tbl_zones_neb = node_SFR.zones_neb
tbl_zones_SF = node_SFR.zones_SF.read_where('id_tSF == 1')
tbl_integrated_neb = node_SFR.integrated_neb
tbl_integrated_SF = node_SFR.integrated_SF.read_where('id_tSF == 1')
node_Zstar = h5file.get_node(Zgroup)
tbl_zones_Z = node_Zstar.zones_Z.read_where('id_tZ == 5')

Hb_obs__g = tbl_zones.cols.F_obs_Hb[:]
O3_obs__g = tbl_zones.cols.F_obs_O3[:]
Ha_obs__g = tbl_zones.cols.F_obs_Ha[:]
N2_obs__g = tbl_zones.cols.F_obs_N2[:]
SNHb = Hb_obs__g / tbl_zones.cols.eF_obs_Hb[:]
SNO3 = O3_obs__g / tbl_zones.cols.eF_obs_O3[:]
SNHa = Ha_obs__g / tbl_zones.cols.eF_obs_Ha[:]
SNN2 = N2_obs__g / tbl_zones.cols.eF_obs_N2[:]

x_Y = tbl_zones_SF['xY']
flag_residual = tbl_zones.cols.flag_residual[:]
flag_BPT = tbl_zones_neb.cols.flag_BPT[:]
flag_RGB = tbl_zones.cols.flag_RGB[:]
tau_V = tbl_zones.cols.tau_V[:]
tau_V_neb = tbl_zones.cols.tau_V_neb[:]
etau_V_neb = tbl_zones.cols.etau_V_neb[:]
zone_distance_HLR = tbl_zones.cols.distance_HLR[:]

mu = tbl_zones.cols.McorSD[:]
aZ_mass = tbl_zones_Z['alogZ_mass']

f = plt.figure()
NCols = 3
NRows = 2
page_size_inches = [15, 10]
f.set_size_inches(page_size_inches)
f.set_dpi(100)
grid_shape = (NRows, NCols)

ax = plt.subplot2grid(grid_shape, (0, 0))
ax = plotmuZR(ax, mu, aZ_mass, labels = False)#, N = 226176)
ax.yaxis.set_major_locator(MaxNLocator(4))
plt.setp(ax.get_xticklabels(), visible = False)

ax = plt.subplot2grid(grid_shape, (0, 1))
mask = (flag_residual > 0)
ax = plotmuZR(ax, mu, aZ_mass, mask = mask, labels = False)
plt.setp(ax.get_xticklabels(), visible = False)
plt.setp(ax.get_yticklabels(), visible = False)

ax = plt.subplot2grid(grid_shape, (0, 2))
mask |= (flag_BPT > 0)
mask |= (flag_RGB > 0)
mask = mask | (SNHb < 3) | (SNO3 < 3) | (SNHa < 3) | (SNN2 < 3)
mask = mask | np.isnan(SNHb) | np.isnan(SNO3) | np.isnan(SNHa) | np.isnan(SNN2)
ax = plotmuZR(ax, mu, aZ_mass, mask = mask, labels = False)
plt.setp(ax.get_xticklabels(), visible = False)
plt.setp(ax.get_yticklabels(), visible = False)

ax = plt.subplot2grid(grid_shape, (1, 0))
mask = (mask | np.isnan(x_Y) | (x_Y < 0.05))
ax = plotmuZR(ax, mu, aZ_mass, mask = mask)
ax.xaxis.set_major_locator(MaxNLocator(4))
ax.yaxis.set_major_locator(MaxNLocator(4))

ax = plt.subplot2grid(grid_shape, (1, 1))
mask = mask | (tau_V < 0.05) | (tau_V_neb < 0.05) | (etau_V_neb >= 0.25)
mask = mask | np.isnan(tau_V) | np.isnan(tau_V_neb) | np.isnan(etau_V_neb)
ax = plotmuZR(ax, mu, aZ_mass, mask = mask, labels = False)
ax.xaxis.set_major_locator(MaxNLocator(4))
plt.setp(ax.get_yticklabels(), visible = False)

ax = plt.subplot2grid(grid_shape, (1, 2))
mask = (mask | np.isnan(zone_distance_HLR) | (zone_distance_HLR < 0.7) | (zone_distance_HLR > 3))
ax = plotmuZR(ax, mu, aZ_mass, mask = mask, labels = False)#, N = 16479)
ax.xaxis.set_major_locator(MaxNLocator(4))
plt.setp(ax.get_yticklabels(), visible = False)

f.subplots_adjust(left = 0.1, right = 0.90, top = 0.95, bottom = 0.1, hspace = 0, wspace = 0)
f.savefig('muZ%s' % fpref)

h5file.close()
