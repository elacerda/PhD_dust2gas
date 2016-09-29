import os
import sys
import numpy as np
import tables as tbl
import matplotlib.pyplot as plt
from pycasso import fitsQ3DataCube
from CALIFAUtils.objects import GasProp
from CALIFAUtils.scripts import calc_xY, prop_Y
from CALIFAUtils.plots import DrawHLRCircle, plot_gal_img_ax


def radialProfileWeighted(v__yx, w__yx, **kwargs):
    r_func = kwargs.get('r_func', None)
    rad_scale = kwargs.get('rad_scale', None)
    bin_r = kwargs.get('bin_r', None)

    v__r = None

    if r_func:
        w__r = r_func(w__yx, bin_r=bin_r, mode='sum', rad_scale=rad_scale)
        v_w__r = r_func(v__yx * w__yx, bin_r=bin_r, mode='sum', rad_scale=rad_scale)
        v__r = v_w__r / w__r

    return v__r

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
    pycasso_cube_dir = eml_cube_dir = gasprop_cube_dir = '%s/LOCAL/data/' % os.environ['HOME']
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

    tau_V__r = np.ma.masked_all((N_R_bins))
    tau_V_npts__r = np.ma.masked_all((N_R_bins))
    tau_V_me__r = np.ma.masked_all((N_R_bins))
    tau_V_me_npts__r = np.ma.masked_all((N_R_bins))
    tau_V_L__r = np.ma.masked_all((N_R_bins))

    tau_V__yx = K.A_V__yx * AVtoTauV
    tau_V__r, tau_V_npts__r = K.radialProfile(prop=tau_V__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, return_npts=True)
    tau_V_me__r, tau_V_me_npts__r = K.radialProfile(prop=tau_V__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, mode='mean_exact', return_npts=True)
    fobs_norm__yx = K.zoneToYX(K.fobs_norm)
    fobs_norm_corr__yx = fobs_norm__yx * np.exp(tau_V__yx)
    aux0 = K.radialProfile(prop=fobs_norm__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, mode='sum')
    aux1 = K.radialProfile(prop=fobs_norm_corr__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, mode='sum')
    tau_V_L__r = np.log(aux1) - np.log(aux0)

    N_rows, N_cols = 1, 3
    f, axArr = plt.subplots(N_rows, N_cols)
    f.set_dpi(100)
    f.set_size_inches(15, 5)
    ax = axArr[0]
    img_file = '%s%s.jpg' % (img_dir, califaID)
    plot_gal_img_ax(ax, img_file, califaID, 0.02, 0.98, 16, K)
    ax = axArr[1]
    im = ax.imshow(tau_V__yx, **default_kwargs_imshow)
    DrawHLRCircle(ax, K)
    # print im.get_window_extent()
    print ax.__dict__.keys()
    print ax._current_image

    ax = axArr[2]

    ax.plot(R_bin_center__r, tau_V__r, ls='-', label='tau_V')
    ax.plot(R_bin_center__r, tau_V_me__r, ls='--', label='tau_V_me')
    ax.plot(R_bin_center__r, tau_V_L__r, ls=':', label='tau_V_L')
    ax.legend(loc='best')
    ax.set_xlabel('R')
    ax.set_ylabel(r'$\tau_V$')
    ax.set_aspect('equal')
    plt.subplots_adjust(wspace = 0.1, left = 0.05, bottom = 0.1, right = 0.9, top = 0.95)
    f.savefig('%s-tau_V.png' % califaID)
    plt.close(f)

    # x_Y__Tr = np.ma.masked_all((N_T, N_R_bins))
    # x_Y_npts__Tr = np.ma.masked_all((N_T, N_R_bins))
    # x_Y_me__Tr = np.ma.masked_all((N_T, N_R_bins))
    # x_Y_me_npts__Tr = np.ma.masked_all((N_T, N_R_bins))
    # x_Y_L__Tr = np.ma.masked_all((N_T, N_R_bins))
    #
    # for i_T, tSF in enumerate(tSF__T):
    #     x_Y__z, integrated_x_Y = calc_xY(K, tSF)
    #     x_Y__yx = K.zoneToYX(x_Y__z, extensive=False)
    #     x_Y__Tr[i_T, :], x_Y_npts__Tr[i_T, :] = K.radialProfile(prop=x_Y__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, return_npts=True)
    #     x_Y_me__Tr[i_T, :], x_Y_me_npts__Tr[i_T, :] = K.radialProfile(prop=x_Y__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, mode='mean_exact', return_npts=True)
    #
    #     Lobn_Y__z = prop_Y(K.Lobn__tZz, tSF, K.ageBase)
    #     LobnSD__yx = K.LobnSD__yx
    #     LobnSD_Y__yx = K.zoneToYX(Lobn_Y__z)
    #     LobnSD_sum__r = K.radialProfile(prop=LobnSD__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, mode='sum')
    #     LobnSD_Y_sum__r = K.radialProfile(prop=LobnSD_Y__yx, bin_r=R_bin__r, rad_scale=K.HLR_pix, mode='sum')
    #     x_Y_L__Tr[i_T, :] = LobnSD_Y_sum__r / LobnSD_sum__r
    #
    #     N_rows, N_cols = 2, 2
    #     f, axArr = plt.subplots(N_rows, N_cols)
    #     f.set_dpi(100)
    #     f.set_size_inches(10, 10)
    #     ax = axArr[0, 0]
    #     img_file = '%s%s.jpg' % (img_dir, califaID)
    #     plot_gal_img_ax(ax, img_file, califaID, 0.02, 0.98, 16, K)
    #     ax = axArr[0, 1]
    #     im = ax.imshow(x_Y__yx, **default_kwargs_imshow)
    #     DrawHLRCircle(ax, K)
    #     ax = axArr[1, 0]
    #     ax.plot(R_bin_center__r, x_Y__Tr[i_T, :], ls='-', label='x_Y')
    #     ax.plot(R_bin_center__r, x_Y_me__Tr[i_T, :], ls='--', label='x_Y_me')
    #     ax.plot(R_bin_center__r, x_Y_L__Tr[i_T, :], ls=':', label='x_Y_L')
    #     ax.legend(loc='best')
    #     ax.set_xlabel('R')
    #     ax.set_ylabel(r'$x_Y$')
    #     plt.subplots_adjust(wspace = 0.1, left = 0.1, bottom = 0.1, right = 0.9, top = 0.95)
    #     f.savefig('%s-x_Y_%.2fMyr.png' % (califaID, tSF/1e6))
    #     plt.close(f)
    #
    #     # Mcor_Y__z = prop_Y(K.Mcor__tZz, tSF, K.ageBase)
    #     # McorSD__yx = K.McorSD__yx
    #     # McorSD_Y__yx = K.zoneToYX(Mcor_Y__z)
