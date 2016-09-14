import sys
import numpy as np
import tables as tbl
import matplotlib as mpl
mpl.rcParams['font.size'] = 13
mpl.rcParams['axes.labelsize'] = 13
mpl.rcParams['axes.titlesize'] = 13
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from pytu.plots import plot_histo_ax
from pytu.plots import next_row_col

if __name__ == '__main__':
    h5fname = sys.argv[1]
    # h5fname = '/home/lacerda/dev/astro/PhD_dust2gas/hdf5/v20_q050.d15a_rawdata_SFR_Zstar.h5'

    fpref = '.png'

    h5file = tbl.open_file(h5fname, mode='r')
    tbl_gals = h5file.root.pycasso.main
    tbl_zones = h5file.root.pycasso.zones

    Hb_obs__g = tbl_zones.cols.F_obs_Hb[:]
    O3_obs__g = tbl_zones.cols.F_obs_O3[:]
    Ha_obs__g = tbl_zones.cols.F_obs_Ha[:]
    N2_obs__g = tbl_zones.cols.F_obs_N2[:]
    SNHb = Hb_obs__g / tbl_zones.cols.eF_obs_Hb[:]
    SNO3 = O3_obs__g / tbl_zones.cols.eF_obs_O3[:]
    SNHa = Ha_obs__g / tbl_zones.cols.eF_obs_Ha[:]
    SNN2 = N2_obs__g / tbl_zones.cols.eF_obs_N2[:]

    f = plt.figure()
    NRows = 2
    NCols = 2
    page_size_inches = [NCols * 4, NRows * 4]
    f.set_size_inches(page_size_inches)
    f.set_dpi(100)
    grid_shape = (NRows, NCols)
    plot_order = [r'SN(H${}_{\beta}$)', r'SN(H${}_{\alpha}$)', r'SN(O[III])', r'SN(N[II])']
    plot_dict = {r'SN(H${}_{\beta}$)': SNHb, r'SN(H${}_{\alpha}$)': SNHa, r'SN(O[III])': SNO3, r'SN(N[II])': SNN2}
    row, col = 0, 0
    for k in plot_order:
        ax = plt.subplot2grid(grid_shape, loc=(row, col))
        SN = plot_dict[k]
        m = np.bitwise_or(np.isinf(SN), np.isnan(SN))
        m = np.bitwise_or(m, np.greater(SN, 1e4))
        SNm = np.ma.masked_array(SN, mask=m)
        x = SNm.compressed()
        ax.set_title(r'N = %d ($\geq$ 3: %d)' % (len(x), (x >= 3).astype('int').sum()))
        ax.set_xlabel(k)
        ax.set_ylabel(r'frac')
        kwh = {'bins': 20, 'range': [0, 20], 'normed': True}
        ax = plot_histo_ax(ax, x, first=True, fs=10, kwargs_histo=kwh)
        ax.axvline(x=3, ls='--', c='k')
        row, col = next_row_col(row, col, NRows, NCols)
    f.subplots_adjust(hspace=0.35, wspace=0.4, left=0.1, right=0.95, top=0.95)
    f.savefig('teste%s' % fpref)
    h5file.close()
