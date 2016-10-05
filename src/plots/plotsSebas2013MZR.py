import os
import sys
import numpy as np
import matplotlib as mpl
from pytu.objects import runstats
from matplotlib import pyplot as plt
from pytu.functions import ma_mask_xyz
from matplotlib.ticker import MaxNLocator
from pytu.plots import plot_text_ax, density_contour

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
dflt_kw_scatter = dict(cmap='viridis', marker='o', s=5, edgecolor='none')

mSFRSD = np.ma.load('mSFRSD.pkl')
mZneb = np.ma.load('mZneb.pkl')
mMcorSD = np.ma.load('mMcorSD.pkl')
mZneb_res = np.ma.load('mZneb_res.pkl')
mtype = np.ma.load('m_type.pkl')
mtypeorig = np.ma.load('m_type_orig.pkl')

# m = np.bitwise_or(np.hstack(SFRSD_mask), np.bitwise_or(np.hstack(Zneb_mask), np.hstack(McorSD_mask)))
# m_res = np.bitwise_or(np.hstack(SFRSD_mask), np.bitwise_or(np.hstack(Zneb_res_mask), np.hstack(McorSD_mask)))

f = plt.figure(dpi=100, figsize=(10, 8))
N_cols, N_rows = 2, 2
ax = f.gca()
x, y = np.ma.log10(mMcorSD), mZneb
z = mtype
xm, ym, zm = ma_mask_xyz(x, y, z)
# z = np.ma.log10(mSFRSD.ravel())
sc = ax.scatter(xm, ym, c=zm, **dflt_kw_scatter)
bins = [50, 50]
extent = [-0.5, 4.5, 8, 8.8]
density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range=[extent[0:2], extent[2:4]], levels_confidence=[0.05, 0.25, 0.45, 0.55, 0.75, 0.95], colors='k')
plt.colorbar(sc)
ax.set_ylim(extent[2:4])
N_y_out = (ym.compressed() < 8).astype('int').sum()
print N_y_out
f.savefig('MZR.png')
