#
# Lacerda@Corrego - 24/Mar/2016
#
import sys
import numpy as np
import tables as tbl
import matplotlib as mpl
mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plot_text_ax, density_contour
from CALIFAUtils.scripts import ma_mask_xyz


def plotBPT(ax, N2Ha, O3Hb, z=None, cmap='viridis', mask=None, labels=True, N=False):
    if mask is None:
        mask = np.zeros_like(O3Hb, dtype=np.bool_)
    from CALIFAUtils.lines import Lines
    extent = [-1.5, 1, -1.5, 1.5]
    if z is None:
        bins = [30, 30]
        xm, ym = ma_mask_xyz(N2Ha, O3Hb, mask=mask)
        density_contour(xm.compressed(), ym.compressed(), bins[0], bins[1], ax, range=[extent[0:2], extent[2:4]], colors=['b', 'y', 'r'])
        # counts, xe, ye, im = ax.hist2d(xm.compressed(), ym.compressed(), range=[extent[0:2], extent[2:4]], bins=bins, cmap='Blues')
        sc = ax.scatter(xm, ym, marker='o', c='0.5', s=10, edgecolor='none', alpha=0.4)
    else:
        xm, ym, z = ma_mask_xyz(N2Ha, O3Hb, z, mask=mask)
        sc = ax.scatter(xm, ym, c=z, cmap=cmap, marker='o', s=30, edgecolor='none')  # , alpha=0.6)
        cb = plt.colorbar(sc, ticks=[0, .5, 1, 1.5, 2])
        cb.set_label('R [HLR]')
    if labels:
        ax.set_xlabel(r'$\log\ [NII]/H\alpha$')
        ax.set_ylabel(r'$\log\ [OIII]/H\beta$')
    L = Lines()
    if not N:
        N = xm.count()
    plt.axis(extent)
    plot_text_ax(ax, '%d' % N, 0.01, 0.99, 20, 'top', 'left', 'k')
    plot_text_ax(ax, 'S06', 0.25, 0.02, 20, 'bottom', 'left', 'k')
    plot_text_ax(ax, 'K01', 0.95, 0.02, 20, 'bottom', 'right', 'k')
    plot_text_ax(ax, 'CF10', 0.92, 0.98, 20, 'top', 'right', 'k', rotation=38)  # 44.62)
    ax.plot(L.x['S06'], L.y['S06'], 'k--', label='S06')
    ax.plot(L.x['K01'], L.y['K01'], 'k--', label='K01')
    ax.plot(L.x['CF10'], L.y['CF10'], 'k--', label='CF10')
    return ax


if __name__ == '__main__':
    h5fname = sys.argv[1]
    # h5fname = '/home/lacerda/dev/astro/PhD_dust2gas/hdf5/v20_q050.d15a_rawdata_SFR_Zstar.h5'

    SFRgroup = '/SFR05050525'
    fpref = '.png'

    iT = 1
    h5file = tbl.open_file(h5fname, mode='r')
    tbl_gals = h5file.root.pycasso.main
    tbl_zones = h5file.root.pycasso.zones
    tbl_integrated = h5file.root.pycasso.integrated
    node_SFR = h5file.get_node(SFRgroup)
    tSF__T = node_SFR.tSF.cols.age[:]
    tbl_SFR_zones_neb = node_SFR.zones_neb
    tbl_SFR_zones_SF = node_SFR.zones_SF.read_where('id_tSF == 1')
    tbl_SFR_integrated_neb = node_SFR.integrated_neb
    tbl_SFR_integrated_SF = node_SFR.integrated_SF.read_where('id_tSF == 1')

    Hb_obs__g = tbl_zones.cols.F_obs_Hb[:]
    O3_obs__g = tbl_zones.cols.F_obs_O3[:]
    Ha_obs__g = tbl_zones.cols.F_obs_Ha[:]
    N2_obs__g = tbl_zones.cols.F_obs_N2[:]
    SNHb = Hb_obs__g / tbl_zones.cols.eF_obs_Hb[:]
    SNO3 = O3_obs__g / tbl_zones.cols.eF_obs_O3[:]
    SNHa = Ha_obs__g / tbl_zones.cols.eF_obs_Ha[:]
    SNN2 = N2_obs__g / tbl_zones.cols.eF_obs_N2[:]
    # Razao entre os fluxos de N2/Ha e O3/Hb
    N2Ha__g = np.log10(N2_obs__g) - np.log10(Ha_obs__g)
    O3Hb__g = np.log10(O3_obs__g) - np.log10(Hb_obs__g)

    x_Y = tbl_SFR_zones_SF['xY']
    flag_residual = tbl_zones.cols.flag_residual[:]
    flag_BPT = tbl_SFR_zones_neb.cols.flag_BPT[:]
    flag_RGB = tbl_zones.cols.flag_RGB[:]
    tau_V = tbl_zones.cols.tau_V[:]
    tau_V_neb = tbl_zones.cols.tau_V_neb[:]
    etau_V_neb = tbl_zones.cols.etau_V_neb[:]
    zone_distance_HLR = tbl_zones.cols.distance_HLR[:]

    f = plt.figure()
    NCols = 3
    NRows = 2
    page_size_inches = [15, 10]
    f.set_size_inches(page_size_inches)
    f.set_dpi(100)
    grid_shape = (NRows, NCols)

    ax = plt.subplot2grid(grid_shape, (0, 0))
    ax = plotBPT(ax, N2Ha__g, O3Hb__g, labels=False, N=226176)
    ax.yaxis.set_major_locator(MaxNLocator(4))
    plt.setp(ax.get_xticklabels(), visible=False)

    ax = plt.subplot2grid(grid_shape, (0, 1))
    mask = (flag_residual > 0)
    ax = plotBPT(ax, N2Ha__g, O3Hb__g, mask=mask, labels=False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)

    ax = plt.subplot2grid(grid_shape, (0, 2))
    mask |= (flag_BPT > 0)
    mask |= (flag_RGB > 0)
    mask = mask | (SNHb < 3) | (SNO3 < 3) | (SNHa < 3) | (SNN2 < 3)
    mask = mask | np.isnan(SNHb) | np.isnan(SNO3) | np.isnan(SNHa) | np.isnan(SNN2)
    mask = mask | np.isinf(SNHb) | np.isinf(SNO3) | np.isinf(SNHa) | np.isinf(SNN2)
    ax = plotBPT(ax, N2Ha__g, O3Hb__g, mask=mask, labels=False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)

    ax = plt.subplot2grid(grid_shape, (1, 0))
    mask = (mask | np.isnan(x_Y) | np.isinf(x_Y) | (x_Y < 0.05))
    ax = plotBPT(ax, N2Ha__g, O3Hb__g, mask=mask)
    ax.xaxis.set_major_locator(MaxNLocator(4, prune='upper'))
    ax.yaxis.set_major_locator(MaxNLocator(4))

    ax = plt.subplot2grid(grid_shape, (1, 1))
    mask = mask | (tau_V < 0.05) | (tau_V_neb < 0.05) | (etau_V_neb >= 0.25)
    mask = mask | np.isnan(tau_V) | np.isnan(tau_V_neb) | np.isnan(etau_V_neb)
    mask = mask | np.isinf(tau_V) | np.isinf(tau_V_neb) | np.isinf(etau_V_neb)
    ax = plotBPT(ax, N2Ha__g, O3Hb__g, mask=mask, labels=False)
    ax.xaxis.set_major_locator(MaxNLocator(4, prune='upper'))
    plt.setp(ax.get_yticklabels(), visible=False)

    ax = plt.subplot2grid(grid_shape, (1, 2))
    mask = mask | np.isnan(zone_distance_HLR) | np.isinf(zone_distance_HLR)
    mask = mask | (zone_distance_HLR < 0.7) | (zone_distance_HLR > 3)
    ax = plotBPT(ax, N2Ha__g, O3Hb__g, mask=mask, labels=False, N=16479)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    plt.setp(ax.get_yticklabels(), visible=False)

    f.subplots_adjust(left=0.1, right=0.90, top=0.95, bottom=0.1, hspace=0, wspace=0)
    f.savefig('BPT%s' % fpref)

    h5file.close()
