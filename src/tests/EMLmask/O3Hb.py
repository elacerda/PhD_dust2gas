import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pycasso import fitsQ3DataCube
from pytu.plots import plot_text_ax, plot_OLSbisector_ax
from pytu.functions import ma_mask_xyz
from CALIFAUtils.objects import GasProp
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.scripts import create_zones_masks_gal
# from pytu.functions import ma_mask_xyz
# from CALIFAUtils.scripts import calc_SFR
from pytu.objects import runstats
# from CALIFAUtils.plots import DrawHLRCircle, plot_gal_img_ax


mpl.use('agg')


def plotBPT(ax, O3Hb, N2Ha, z=None, cmap='viridis', mask=None, npts=True, clean=False, sc_kwargs=None):
    sc_kw = dict(marker='o', s=30, edgecolor='none')
    if isinstance(sc_kwargs, dict):
        sc_kw.update(sc_kwargs)
    if mask is None:
        mask = np.ones_like(O3Hb, dtype=np.bool_)
    from CALIFAUtils.lines import Lines
    L = Lines()
    # plt.axis([-1, 0.5, -1, 1])
    ax.axis([-1, 0.8, -1.5, 1])
    if z is None:
        sc = ax.scatter(N2Ha[mask], O3Hb[mask], **sc_kw)
    else:
        sc = ax.scatter(N2Ha[mask], O3Hb[mask], c=z[mask], cmap=cmap, marker='o', s=30, edgecolor='none', **sc_kw)
        cb = plt.colorbar(sc, ticks=[0, .5, 1, 1.5, 2])
        cb.set_label('R [HLR]', fontsize=15)
    if not clean:
        trans_angle = ax.transData.transform_angles(np.array((np.arctan(1.01) * 180 / np.pi,)), np.array((-1, 0.5)).reshape((1, 2)))[0]
        ax.set_xlabel(r'$\log\ [NII]/H\alpha$', fontsize=14)
        ax.set_ylabel(r'$\log\ [OIII]/H\beta$', fontsize=14)
        plot_text_ax(ax, 'S06', 0.25, 0.02, 14, 'bottom', 'left', 'k')
        plot_text_ax(ax, 'K01', 0.95, 0.02, 14, 'bottom', 'right', 'k')
        plot_text_ax(ax, 'CF10', 0.92, 0.98, 14, 'top', 'right', 'k', rotation=trans_angle)
        ax.plot(L.x['S06'], L.y['S06'], 'k--', label='S06')
        ax.plot(L.x['K01'], L.y['K01'], 'k--', label='K01')
        ax.plot(L.x['CF10'], L.y['CF10'], 'k--', label='CF10')
    if npts:
        plot_text_ax(ax, 'N: %d' % len(O3Hb[mask]), 0.01, 1.01, 14, 'bottom', 'left', 'k')
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    return ax


mpl.rcParams['font.size'] = 16
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14


def get_subplot(N_rows, N_cols, N_plot):
    return (N_rows * 100) + (N_cols * 10) + N_plot


if __name__ == '__main__':
    califaID = sys.argv[1]
    # Initial definitions
    pycasso_cube_dir = '%s/CALIFA/gal_fits/v20_q050.d15a/' % os.environ['HOME']
    eml_cube_dir = '%s/CALIFA/rgb-gas/v20_q050.d15a/' % os.environ['HOME']
    gasprop_cube_dir = '%s/CALIFA/rgb-gas/v20_q050.d15a/prop/' % os.environ['HOME']
    h5file = '%s/dev/astro/PhD_dust2gas/runs/v20_q050.d15a.h5' % os.environ['HOME']
    pycasso_cube_suffix = '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.fits'
    eml_cube_suffix = '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.fits'
    gasprop_cube_suffix = '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.GasProp.fits'
    img_dir = '%s/CALIFA/images/' % os.environ['HOME']
    pycasso_cube_file = pycasso_cube_dir + califaID + pycasso_cube_suffix
    eml_cube_file = eml_cube_dir + califaID + eml_cube_suffix
    gasprop_cube_file = gasprop_cube_dir + califaID + gasprop_cube_suffix
    dflt_kw_imshow = dict(origin='lower', interpolation='nearest', aspect='equal', cmap='viridis')
    dflt_kw_runstats = dict(smooth=True, sigma=1.2, debug=True, frac=0.08, gs_prc=True, poly1d=True)
    AVtoTauV = 1. / (np.log10(np.exp(1)) / 0.4)
    tSF__T = np.array([1, 3.2, 10, 50, 100]) * 1e7
    N_T = len(tSF__T)
    dflt_kw_runstats = dict(smooth = True, sigma = 1.2, debug = True, frac = 0.08, gs_prc = True, poly1d = True, OLS=True)

    # Load all data from fits
    K = fitsQ3DataCube(pycasso_cube_file)
    K.loadEmLinesDataCube(eml_cube_file)
    K.GP = GasProp(gasprop_cube_file)

    # Set galaxy geometry using ellipses
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)

    aux = create_zones_masks_gal(K, tSF__T,
                                 mintauv=0.05,
                                 mintauvneb=0.05,
                                 maxtauvneberr=0.25,
                                 minpopx=0.05,
                                 minEWHb=3,
                                 minSNR=3,
                                 minSNRHb=3,
                                 nolinecuts=False,
                                 rgbcuts=True,
                                 underS06=False,
                                 whanSF=False,
                                 filter_residual=True,
                                 summary=True)
    mask__Tz = aux[0]
    mask_syn__Tz = aux[1]
    mask_eml__z = aux[2]
    mask_popx__Tz = aux[3]
    mask_tau_V__z = aux[4]
    mask_residual__z = aux[5]
    mask_tau_V_neb__z = aux[6]
    mask_tau_V_neb_err__z = aux[7]
    mask_EW_Hb__z = aux[8]
    mask_whan__z = aux[9]
    mask_bpt__z = aux[10]
    mask_lines_dict__Lz = aux[11]

    f_obs__lz = K.f_obs
    f_syn__lz = K.f_syn

    wl_band = (K.l_obs > 4800) & (K.l_obs < 5100)
    mask_r = np.less(K.zoneDistance_HLR, 0.7)
    mask_Hb = mask_lines_dict__Lz['4861']
    mask_O3 = mask_lines_dict__Lz['5007']
    m_Hb_O3 = np.where(~mask_Hb & ~mask_O3 & ~mask_r)[0]
    m_Hb_notO3 = np.where(~mask_Hb & mask_O3 & ~mask_r)[0]
    m_notHb_O3 = np.where(mask_Hb & ~mask_O3 & ~mask_r)[0]
    m_notHb_notO3 = np.where(mask_Hb & mask_O3 & ~mask_r)[0]
    # print m_Hb_O3,
    # m_Hb_O3 = np.bitwise_and(m_Hb_O3, np.greater_equal(K.zoneDistance_HLR, 0.7))
    # m_Hb_notO3 = np.bitwise_and(m_Hb_notO3, np.greater_equal(K.zoneDistance_HLR, 0.7))
    # m_notHb_O3 = np.bitwise_and(m_notHb_O3, np.greater_equal(K.zoneDistance_HLR, 0.7))
    # m_notHb_notO3 = np.bitwise_and(m_notHb_notO3, np.greater_equal(K.zoneDistance_HLR, 0.7))

    zones = np.arange(K.N_zone)
    zones_Hb_O3 = zones[m_Hb_O3]
    zones_Hb_notO3 = zones[m_Hb_notO3]
    zones_notHb_O3 = zones[m_notHb_O3]
    zones_notHb_notO3 = zones[m_notHb_notO3]
    f_obs_Hb_O3 = f_obs__lz[:, m_Hb_O3]
    f_obs_Hb_notO3 = f_obs__lz[:, m_Hb_notO3]
    f_obs_notHb_O3 = f_obs__lz[:, m_notHb_O3]
    f_obs_notHb_notO3 = f_obs__lz[:, m_notHb_notO3]
    f_syn_Hb_O3 = f_syn__lz[:, m_Hb_O3]
    f_syn_Hb_notO3 = f_syn__lz[:, m_Hb_notO3]
    f_syn_notHb_O3 = f_syn__lz[:, m_notHb_O3]
    f_syn_notHb_notO3 = f_syn__lz[:, m_notHb_notO3]
    f_res_Hb_O3 = f_obs_Hb_O3 - f_syn_Hb_O3
    f_res_Hb_notO3 = f_obs_Hb_notO3 - f_syn_Hb_notO3
    f_res_notHb_O3 = f_obs_notHb_O3 - f_syn_notHb_O3
    f_res_notHb_notO3 = f_obs_notHb_notO3 - f_syn_notHb_notO3

    Hb_central_wl = '4861'
    O3_central_wl = '5007'
    Ha_central_wl = '6563'
    N2_central_wl = '6583'
    i_Hb = K.EL.lines.index(Hb_central_wl)
    i_O3 = K.EL.lines.index(O3_central_wl)
    i_Ha = K.EL.lines.index(Ha_central_wl)
    i_N2 = K.EL.lines.index(N2_central_wl)



    f = plt.figure(dpi=100, figsize=(15, 10))
    f.suptitle(r'H$\beta$ - pos:4861$\pm%.1f\AA$ - max$\sigma$:%.1f - minsnr:%.2f -- [OIII] - pos:5007$\pm%.1f\AA$ - max$\sigma$:%.1f - minsnr:%.2f' %
               (K.GP._dlcons['4861']['pos'], K.GP._dlcons['4861']['sigma'], K.GP._dlcons['4861']['SN'],
                K.GP._dlcons['5007']['pos'], K.GP._dlcons['5007']['sigma'], K.GP._dlcons['5007']['SN'])
               )
    multipl_res = 1.
    ax = plt.subplot(221)
    ax.set_title('both OK')
    ax.plot(K.l_obs[wl_band], f_obs_Hb_O3[wl_band, 0], 'k')
    ax.plot(K.l_obs[wl_band], f_syn_Hb_O3[wl_band, 0], 'g')
    ax.plot(K.l_obs[wl_band], multipl_res * f_res_Hb_O3[wl_band, 0], 'r')
    plot_text_ax(ax, txt='zone: %d' % zones_Hb_O3[0], xpos=0.01, ypos=0.99, va='top', ha='left')
    txt = r'$\ \ $Hb:%.1f$\AA\ \ $O3:%.1f$\AA$' % (K.EL.pos[i_Hb][m_Hb_O3][0], K.EL.pos[i_O3][m_Hb_O3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.55, va='top', ha='right')
    txt = r'$\sigma$:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.sigma[i_Hb][m_Hb_O3][0], K.EL.sigma[i_O3][m_Hb_O3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.50, va='top', ha='right')
    txt = r'SN:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_Hb_O3][0], K.EL.snr__Lz[i_O3][m_Hb_O3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.45, va='top', ha='right')
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')
    ax = plt.subplot(222)
    ax.set_title('none OK')
    ax.plot(K.l_obs[wl_band], f_obs_notHb_notO3[wl_band, 0], 'k')
    ax.plot(K.l_obs[wl_band], f_syn_notHb_notO3[wl_band, 0], 'g')
    ax.plot(K.l_obs[wl_band], multipl_res * f_res_notHb_notO3[wl_band, 0], 'r')
    plot_text_ax(ax, txt='zone: %d' % zones_notHb_notO3[0], xpos=0.01, ypos=0.99, va='top', ha='left')
    txt = r'$\ \ $Hb:%.1f$\AA\ \ $O3:%.1f$\AA$' % (K.EL.pos[i_Hb][m_notHb_notO3][0], K.EL.pos[i_O3][m_notHb_notO3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.55, va='top', ha='right')
    txt = r'$\sigma$:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.sigma[i_Hb][m_notHb_notO3][0], K.EL.sigma[i_O3][m_notHb_notO3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.50, va='top', ha='right')
    txt = r'SN:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_notHb_notO3][0], K.EL.snr__Lz[i_O3][m_notHb_notO3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.45, va='top', ha='right')
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')
    ax = plt.subplot(223)
    ax.set_title('O3 Ok')
    ax.plot(K.l_obs[wl_band], f_obs_notHb_O3[wl_band, 0], 'k')
    ax.plot(K.l_obs[wl_band], f_syn_notHb_O3[wl_band, 0], 'g')
    ax.plot(K.l_obs[wl_band], multipl_res * f_res_notHb_O3[wl_band, 0], 'k')
    plot_text_ax(ax, txt='zone: %d' % zones_notHb_O3[0], xpos=0.01, ypos=0.99, va='top', ha='left')
    txt = r'$\ \ $Hb:%.1f$\AA\ \ $O3:%.1f$\AA$' % (K.EL.pos[i_Hb][m_notHb_O3][0], K.EL.pos[i_O3][m_notHb_O3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.55, va='top', ha='right')
    txt = r'$\sigma$:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.sigma[i_Hb][m_notHb_O3][0], K.EL.sigma[i_O3][m_notHb_O3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.50, va='top', ha='right')
    txt = r'SN:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_notHb_O3][0], K.EL.snr__Lz[i_O3][m_notHb_O3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.45, va='top', ha='right')
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')
    ax = plt.subplot(224)
    ax.set_title('Hb Ok')
    ax.plot(K.l_obs[wl_band], f_obs_Hb_notO3[wl_band, 0], 'k')
    ax.plot(K.l_obs[wl_band], f_syn_Hb_notO3[wl_band, 0], 'g')
    ax.plot(K.l_obs[wl_band], multipl_res * f_res_Hb_notO3[wl_band, 0], 'r')
    plot_text_ax(ax, txt='zone: %d' % zones_Hb_notO3[0], xpos=0.01, ypos=0.99, va='top', ha='left')
    txt = r'$\ \ $Hb:%.1f$\AA\ \ $O3:%.1f$\AA$' % (K.EL.pos[i_Hb][m_Hb_notO3][0], K.EL.pos[i_O3][m_Hb_notO3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.55, va='top', ha='right')
    txt = r'$\sigma$:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.sigma[i_Hb][m_Hb_notO3][0], K.EL.sigma[i_O3][m_Hb_notO3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.50, va='top', ha='right')
    txt = r'SN:$\ \ $Hb:%.1f$\ \ $O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_Hb_notO3][0], K.EL.snr__Lz[i_O3][m_Hb_notO3][0])
    plot_text_ax(ax, txt=txt, xpos=0.99, ypos=0.45, va='top', ha='right')
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')
    f.subplots_adjust(wspace=0.1, hspace=0.25, left=0.05, bottom=0.05, right=0.95, top=0.90)
    f.savefig('%s-O3Hb.png' % califaID)
    plt.close(f)

    O3Hb = np.log10(K.EL.O3_obs__z) - np.log10(K.EL.Hb_obs__z)
    N2Ha = np.log10(K.EL.N2_obs__z) - np.log10(K.EL.Ha_obs__z)
    O3Hb_Hb_O3 = np.log10(K.EL.O3_obs__z[m_Hb_O3]) - np.log10(K.EL.Hb_obs__z[m_Hb_O3])
    N2Ha_Hb_O3 = np.log10(K.EL.N2_obs__z[m_Hb_O3]) - np.log10(K.EL.Ha_obs__z[m_Hb_O3])
    O3Hb_Hb_notO3 = np.log10(K.EL.O3_obs__z[m_Hb_notO3]) - np.log10(K.EL.Hb_obs__z[m_Hb_notO3])
    N2Ha_Hb_notO3 = np.log10(K.EL.N2_obs__z[m_Hb_notO3]) - np.log10(K.EL.Ha_obs__z[m_Hb_notO3])
    O3Hb_notHb_O3 = np.log10(K.EL.O3_obs__z[m_notHb_O3]) - np.log10(K.EL.Hb_obs__z[m_notHb_O3])
    N2Ha_notHb_O3 = np.log10(K.EL.N2_obs__z[m_notHb_O3]) - np.log10(K.EL.Ha_obs__z[m_notHb_O3])
    O3Hb_notHb_notO3 = np.log10(K.EL.O3_obs__z[m_notHb_notO3]) - np.log10(K.EL.Hb_obs__z[m_notHb_notO3])
    N2Ha_notHb_notO3 = np.log10(K.EL.N2_obs__z[m_notHb_notO3]) - np.log10(K.EL.Ha_obs__z[m_notHb_notO3])

    f = plt.figure(dpi=100, figsize=(15, 10))
    f.suptitle(r'%d zones (%d > 0.7HLR) H$\beta$ - pos:4861$\pm%.1f\AA$ - max$\sigma$:%.1f - minsnr:%.2f -- [OIII] - pos:5007$\pm%.1f\AA$ - max$\sigma$:%.1f - minsnr:%.2f' %
               (K.N_zone, (~mask_r).astype('int').sum(),
                K.GP._dlcons['4861']['pos'], K.GP._dlcons['4861']['sigma'], K.GP._dlcons['4861']['SN'],
                K.GP._dlcons['5007']['pos'], K.GP._dlcons['5007']['sigma'], K.GP._dlcons['5007']['SN']),
               fontsize=14)
    ax = plt.subplot(221)
    ax.set_title('both OK')
    ax = plotBPT(ax, O3Hb, N2Ha, npts=False, clean=True, sc_kwargs=dict(c='grey', alpha=0.4, ))
    plotBPT(ax, O3Hb_Hb_O3, N2Ha_Hb_O3, sc_kwargs=dict(c='k', edgecolor='w', s=10))
    ax = plt.subplot(222)
    ax.set_title('none OK')
    xm, ym = ma_mask_xyz(N2Ha_notHb_notO3, O3Hb_notHb_notO3)
    rs = runstats(xm.compressed(), ym.compressed(), **dflt_kw_runstats)
    # ax.plot(rs.xS, rs.yS, '.-', c = 'b', lw=2, label = 'smoothed median')
    ax = plot_OLSbisector_ax(ax, rs.OLS_median_slope, rs.OLS_median_intercept,
                        pos_x = 0.98, pos_y = 0.02, fs = 14, c = 'b', OLS = True, x_rms = xm, y_rms = ym)
    # ax.plot(rs.xMedian, rs.yMedian, '.-', c = 'b', lw=2, label = 'median')
    # ax.plot(rs.xMean, rs.yMean, '.-', c = 'r', lw=2, label = 'mean')
    ax.legend(loc='best')
    plotBPT(ax, O3Hb, N2Ha, npts=False, clean=True, sc_kwargs=dict(c='grey', alpha=0.4))
    plotBPT(ax, O3Hb_notHb_notO3, N2Ha_notHb_notO3, sc_kwargs=dict(c='k', edgecolor='w', s=10))
    ax = plt.subplot(223)
    ax.set_title('O3 Ok')
    plotBPT(ax, O3Hb, N2Ha, npts=False, clean=True, sc_kwargs=dict(c='grey', alpha=0.4))
    plotBPT(ax, O3Hb_notHb_O3, N2Ha_notHb_O3, sc_kwargs=dict(c='k', edgecolor='w', s=10))
    ax = plt.subplot(224)
    ax.set_title('Hb Ok')
    plotBPT(ax, O3Hb, N2Ha, npts=False, clean=True, sc_kwargs=dict(c='grey', alpha=0.4))
    plotBPT(ax, O3Hb_Hb_notO3, N2Ha_Hb_notO3, sc_kwargs=dict(c='k', edgecolor='w', s=10))
    f.subplots_adjust(wspace=0.25, hspace=0.25, left=0.1, bottom=0.1, right=0.95, top=0.90)
    f.savefig('%s-BPT.png' % califaID)
    plt.close(f)

    K.GP.close()
    K.EL.close()
    K.close()
