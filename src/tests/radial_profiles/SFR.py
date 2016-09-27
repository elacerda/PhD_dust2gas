import os
import sys
import numpy as np
import tables as tbl
import matplotlib as mpl
import matplotlib.pyplot as plt
from pycasso import fitsQ3DataCube
from CALIFAUtils.objects import GasProp
from CALIFAUtils.scripts import calc_SFR
from CALIFAUtils.plots import DrawHLRCircle, plot_gal_img_ax


mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
#mpl.rcParams['font.family'] = 'serif'
#mpl.rcParams['font.serif'] = 'Times New Roman'


def get_subplot(N_rows, N_cols, N_plot):
    return (N_rows * 100) + (N_cols * 10) + N_plot


if __name__ == '__main__':
    califaID = sys.argv[1]
    # Initial definitions
    rbinini = 0.
    rbinfin = 3.
    rbinstep = 0.1
    R_bin__r = np.arange(rbinini, rbinfin + rbinstep, rbinstep)
    R_bin_center__r = (R_bin__r[:-1] + R_bin__r[1:]) / 2.0
    N_R_bins = len(R_bin_center__r)
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
    default_kwargs_imshow = dict(origin='lower', interpolation='nearest', aspect='equal', cmap='viridis')
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

    SFR = calc_SFR(K, 32e6)[0]
    SFR__yx = K.zoneToYX(SFR, surface_density=False)
    SFRSD__yx = K.zoneToYX(SFR)
    aSFRSD__r = K.radialProfile(SFRSD__yx, R_bin__r, K.HLR_pix)
    alogSFRSD__r = K.radialProfile(np.ma.log10(SFRSD__yx), R_bin__r, K.HLR_pix)

    f = plt.figure(dpi=100, figsize=(15, 10))
    N_rows = 2
    N_cols = 3
    img_file = '%s%s.jpg' % (img_dir, califaID)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 1))
    plot_gal_img_ax(ax, img_file, califaID, 0.02, 0.98, 16, K)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 2))
    ax.set_title('log and de-zone')
    im = ax.imshow(K.zoneToYX(np.ma.log10(SFR), extensive=False, surface_density=False), vmin=-8, vmax=-1, **default_kwargs_imshow)
    DrawHLRCircle(ax, K, color='k')
    f.colorbar(im)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 3))
    ax.set_title('de-zone and log')
    im = ax.imshow(np.ma.log10(K.zoneToYX(SFR, extensive=False, surface_density=False)), vmin=-8, vmax=-1, **default_kwargs_imshow)
    DrawHLRCircle(ax, K, color='k')
    f.colorbar(im)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 4))
    ax.set_title('de-zone weights')
    im = ax.imshow(K._dezonificationWeight, origin='lower', interpolation='nearest', aspect='equal', cmap='viridis')
    DrawHLRCircle(ax, K, color='k')
    f.colorbar(im)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 5))
    ax.set_title('extensive - log and de-zone')
    im = ax.imshow(K.zoneToYX(np.ma.log10(SFR), extensive=True, surface_density=False), vmin=-8, vmax=-1, **default_kwargs_imshow)
    DrawHLRCircle(ax, K, color='k')
    f.colorbar(im)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 6))
    ax.set_title('extensive - de-zone and log')
    im = ax.imshow(np.ma.log10(K.zoneToYX(SFR, extensive=True, surface_density=True)), vmin=-8, vmax=-1, **default_kwargs_imshow)
    DrawHLRCircle(ax, K, color='k')
    f.colorbar(im)
    f.savefig('SFR.png')
    plt.close(f)
    K.close()

    f = plt.figure(dpi=100, figsize=(15, 10))
    N_rows = 2
    N_cols = 2
    img_file = '%s%s.jpg' % (img_dir, califaID)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 1))
    plot_gal_img_ax(ax, img_file, califaID, 0.02, 0.98, 14, K)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 2))
    ax.plot(R_bin_center__r, np.ma.log10(aSFRSD__r), label=r'$\log \langle \Sigma_{SFR} \rangle$')
    ax.plot(R_bin_center__r, alogSFRSD__r, label=r'$\langle \log \Sigma_{SFR} \rangle$')
    ax.set_xlabel('R [ HLR ]')
    ax.set_ylabel(r'$\Sigma_{SFR}\ [M_\odot yr^{-1} pc^{-2}]$')
    ax.legend(loc='best', fontsize=14)
    ax = plt.subplot(get_subplot(N_rows, N_cols, 3))
    ax.set_title(r'$\Sigma_{SFR}$')
    im = ax.imshow(SFRSD__yx*1e7, **default_kwargs_imshow)
    DrawHLRCircle(ax, K, color='k')
    cb = f.colorbar(im)
    cb.set_label(r'$\Sigma_{SFR}\ [10^{-7}\ M_\odot yr^{-1} pc^{-2}]$')
    ax = plt.subplot(get_subplot(N_rows, N_cols, 4))
    ax.set_title(r'$\log \Sigma_{SFR}$')
    im = ax.imshow(np.ma.log10(SFRSD__yx), vmax=-6, vmin=-9, **default_kwargs_imshow)
    DrawHLRCircle(ax, K, color='k')
    cb = f.colorbar(im)
    cb.set_label(r'$\log \Sigma_{SFR}\ [M_\odot yr^{-1} pc^{-2}]$')
    f.subplots_adjust(wspace=0.1, hspace=0.25, left=0.03, bottom=0.1, right=0.95, top=0.95)
    f.savefig('%s-SFRSD.png' % califaID)
    plt.close(f)
    K.close()
