import os
import sys
import numpy as np
import tables as tbl
import matplotlib as mpl
import matplotlib.pyplot as plt
from pycasso import fitsQ3DataCube
from pytu.plots import plot_text_ax
from CALIFAUtils.objects import GasProp
from CALIFAUtils.scripts import create_zones_masks_gal
# from pytu.functions import ma_mask_xyz
# from CALIFAUtils.scripts import calc_SFR
# from pytu.objects import runstats
# from CALIFAUtils.plots import DrawHLRCircle, plot_gal_img_ax


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

    # truth table
    mask_l = (K.l_obs > 4800) & (K.l_obs < 5100)
    m_Hb_O3 = np.where((mask_lines_dict__Lz['4861']==False) & (mask_lines_dict__Lz['5007']==False))[0]
    m_Hb_notO3 = np.where((mask_lines_dict__Lz['4861']==False) & (mask_lines_dict__Lz['5007']==True))[0]
    m_notHb_O3 = np.where((mask_lines_dict__Lz['4861']==True) & (mask_lines_dict__Lz['5007']==False))[0]
    m_notHb_notO3 = np.where((mask_lines_dict__Lz['4861']==True) & (mask_lines_dict__Lz['5007']==True))[0]
    f_obs_Hb_O3 = f_obs__lz[:, m_Hb_O3]
    f_obs_Hb_notO3 = f_obs__lz[:, m_Hb_notO3]
    f_obs_notHb_O3 = f_obs__lz[:, m_notHb_O3]
    f_obs_notHb_notO3 = f_obs__lz[:, m_notHb_notO3]

    Hb_central_wl = '4861'
    O3_central_wl = '5007'
    i_Hb = K.EL.lines.index(Hb_central_wl)
    i_O3 = K.EL.lines.index(O3_central_wl)

    f = plt.figure(dpi=100, figsize=(15, 10))
    f.suptitle(r'H$\beta$ - pos:4861$\pm%.1f\AA$ - max$\sigma$:%.1f - minsnr:%.2f -- [OIII] - pos:5007$\pm%.1f\AA$ - max$\sigma$:%.1f - minsnr:%.2f' %
               (K.GP._dlcons['4861']['pos'], K.GP._dlcons['4861']['sigma'], K.GP._dlcons['4861']['SN'],
                K.GP._dlcons['5007']['pos'], K.GP._dlcons['5007']['sigma'], K.GP._dlcons['5007']['SN'])
               )

    ax = plt.subplot(221)
    ax.plot(K.l_obs[mask_l], f_obs_Hb_O3[mask_l, 0])
    ax.set_title('both OK')
    plot_text_ax(ax, txt='pos: Hb:%.1f O3:%.1f' % (K.EL.pos[i_Hb][m_Hb_O3][0], K.EL.pos[i_O3][m_Hb_O3][0]), xpos=0.99, ypos=0.99, va='top', ha='right')
    plot_text_ax(ax, txt=r'$\sigma$: Hb:%.1f O3:%.1f' % (K.EL.sigma[i_Hb][m_Hb_O3][0], K.EL.sigma[i_O3][m_Hb_O3][0]), xpos=0.99, ypos=0.94, va='top', ha='right')
    plot_text_ax(ax, txt='SN: Hb:%.1f O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_Hb_O3][0], K.EL.snr__Lz[i_O3][m_Hb_O3][0]), xpos=0.99, ypos=0.89, va='top', ha='right')
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')

    ax = plt.subplot(222)
    ax.set_title('none OK')
    plot_text_ax(ax, txt='pos: Hb:%.1f O3:%.1f' % (K.EL.pos[i_Hb][m_notHb_notO3][0], K.EL.pos[i_O3][m_notHb_notO3][0]), xpos=0.99, ypos=0.99, va='top', ha='right')
    plot_text_ax(ax, txt=r'$\sigma$: Hb:%.1f O3:%.1f' % (K.EL.sigma[i_Hb][m_notHb_notO3][0], K.EL.sigma[i_O3][m_notHb_notO3][0]), xpos=0.99, ypos=0.94, va='top', ha='right')
    plot_text_ax(ax, txt='SN: Hb:%.1f O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_notHb_notO3][0], K.EL.snr__Lz[i_O3][m_notHb_notO3][0]), xpos=0.99, ypos=0.89, va='top', ha='right')
    ax.plot(K.l_obs[mask_l], f_obs_notHb_notO3[mask_l, 0])
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')

    ax = plt.subplot(223)
    plot_text_ax(ax, txt='pos: Hb:%.1f O3:%.1f' % (K.EL.pos[i_Hb][m_notHb_O3][0], K.EL.pos[i_O3][m_notHb_O3][0]), xpos=0.99, ypos=0.99, va='top', ha='right')
    plot_text_ax(ax, txt=r'$\sigma$: Hb:%.1f O3:%.1f' % (K.EL.sigma[i_Hb][m_notHb_O3][0], K.EL.sigma[i_O3][m_notHb_O3][0]), xpos=0.99, ypos=0.94, va='top', ha='right')
    plot_text_ax(ax, txt='SN: Hb:%.1f O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_notHb_O3][0], K.EL.snr__Lz[i_O3][m_notHb_O3][0]), xpos=0.99, ypos=0.89, va='top', ha='right')
    ax.set_title('O3 Ok')
    ax.plot(K.l_obs[mask_l], f_obs_notHb_O3[mask_l, 0])
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')

    ax = plt.subplot(224)
    plot_text_ax(ax, txt='pos: Hb:%.1f O3:%.1f' % (K.EL.pos[i_Hb][m_Hb_notO3][0], K.EL.pos[i_O3][m_Hb_notO3][0]), xpos=0.99, ypos=0.99, va='top', ha='right')
    plot_text_ax(ax, txt=r'$\sigma$: Hb:%.1f O3:%.1f' % (K.EL.sigma[i_Hb][m_Hb_notO3][0], K.EL.sigma[i_O3][m_Hb_notO3][0]), xpos=0.99, ypos=0.94, va='top', ha='right')
    plot_text_ax(ax, txt='SN: Hb:%.1f O3:%.1f' % (K.EL.snr__Lz[i_Hb][m_Hb_notO3][0], K.EL.snr__Lz[i_O3][m_Hb_notO3][0]), xpos=0.99, ypos=0.89, va='top', ha='right')
    ax.set_title('Hb Ok')
    ax.plot(K.l_obs[mask_l], f_obs_Hb_notO3[mask_l, 0])
    ax.axvline(x=4861, ls='--')
    ax.axvline(x=5007, ls='--')

    f.subplots_adjust(wspace=0.1, hspace=0.25, left=0.05, bottom=0.05, right=0.95, top=0.90)
    f.savefig('O3Hb.png')
    plt.close(f)

    K.GP.close()
    K.EL.close()
    K.close()
