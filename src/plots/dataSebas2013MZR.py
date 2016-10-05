import os
import sys
import time
import numpy as np
from pycasso import fitsQ3DataCube
from matplotlib import pyplot as plt
from pytu.objects import tupperware_none
from CALIFAUtils.objects import GasProp
from CALIFAUtils.scripts import create_zones_masks_gal, calc_SFR, \
                                get_morfologia, sort_gals, my_morf


def verify_files(K, califaID, EL = True, GP = True):
    if K is None:
        print '<<< %s galaxy: miss files' % califaID
        return 0, False
    if EL and K.EL is None:
        print '<<< %s galaxy: miss EmLines files' % califaID
        return 1, False
        if K.EL.flux[0, :].sum() == 0.:
            print '<<< %s EmLines FITS problem' % califaID
            return 2, False
    if GP and K.GP._hdulist is None:
        print '<<< %s galaxy: miss gasprop file' % califaID
        return 2, False
    # Problem in FITS file
    return 0, True


def pickle_ma_dict(to_pickle=None, file_pref='out-'):
    if to_pickle is not None:
        for k in to_pickle.keys():
            fname = '%s%s.pkl' % (file_pref, k)
            print fname
            N = len(to_pickle[k])
            data = [to_pickle[k][i].data for i in xrange(N)]
            mask = [to_pickle[k][i].mask for i in xrange(N)]
            np.ma.masked_array(np.hstack(data), mask=np.hstack(mask)).dump(fname)


def pickle_ma_ravel_dict(to_pickle=None, file_pref='out-'):
    if to_pickle is not None:
        for k in to_pickle.keys():
            fname = '%s%s.pkl' % (file_pref, k)
            print fname
            N = len(to_pickle[k])
            data = [to_pickle[k][i].data.ravel() for i in xrange(N)]
            mask = [to_pickle[k][i].mask.ravel() for i in xrange(N)]
            np.ma.masked_array(np.hstack(data), mask=np.hstack(mask)).dump(fname)

debug = True
max_gals = 10
mintauv = 0.05,
mintauvneb = 0.05,
maxtauvneberr = 0.25,
minpopx = 0.05,
minEWHb = 3,
minSNR = 3,
minSNRHb = 3,
nolinecuts = False,
rgbcuts = True,
underS06 = False,
whanSF = False,
filter_residual = True,
fpref = '.png'
rbinini = 0.
rbinfin = 3.
rbinstep = 0.1
R_bin__r = np.arange(rbinini, rbinfin + rbinstep, rbinstep)
R_bin_center__r = (R_bin__r[:-1] + R_bin__r[1:]) / 2.0
N_R_bins = len(R_bin_center__r)
tSF__T = np.array([1, 3.2, 10, 50, 100]) * 1e7
N_T = len(tSF__T)
iT = 1  # 32 Myrs
dflt_kw_imshow = dict(origin='lower', interpolation='nearest', aspect='equal', cmap='viridis')
dflt_kw_runstats = dict(smooth=True, sigma=1.2, debug=True, frac=0.08, gs_prc=True, poly1d=True)
dflt_kw_scatter = dict(cmap='viridis', marker='o', s=5, edgecolor='none')
work_dir = '%s/CALIFA/' % os.environ['HOME']
pycasso_cube_dir = '%sgal_fits/v20_q050.d15a/' % work_dir
eml_cube_dir = '%srgb-gas/v20_q050.d15a/' % work_dir
gasprop_cube_dir = '%srgb-gas/v20_q050.d15a/prop/' % work_dir
pycasso_cube_suffix = '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.fits'
eml_cube_suffix = '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.fits'
gasprop_cube_suffix = '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.GasProp.fits'
img_dir = '%simages/' % work_dir
dflt_kw_create_zones_masks = dict(
    return_mask_lines_separated=True,
    mask_lines_snr_only=True,
    summary=True,
    mintauv=mintauv,
    mintauvneb=mintauvneb,
    maxtauvneberr=maxtauvneberr,
    minpopx=minpopx,
    minEWHb=minEWHb,
    minSNR=minSNR,
    minSNRHb=minSNRHb,
    nolinecuts=nolinecuts,
    rgbcuts=rgbcuts,
    underS06=underS06,
    whanSF=whanSF,
    filter_residual=filter_residual,
)

if __name__ == '__main__':
    list_gals = sys.argv[1]
    gals, _ = sort_gals(gals=list_gals, work_dir=work_dir)
    N_gals = len(gals)
    if debug:
        max_gals = max_gals
        if N_gals > max_gals:
            N_gals = max_gals

    masks = {}
    masks['mask__gTz'] = []
    masks['mask_syn__gTz'] = []
    masks['mask_eml__gz'] = []
    masks['mask_popx__gTz'] = []
    masks['mask_tau_V__gz'] = []
    masks['mask_residual__gz'] = []
    masks['mask_tau_V_neb__gz'] = []
    masks['mask_tau_V_neb_err__gz'] = []
    masks['mask_EW_Hb__gz'] = []
    masks['mask_whan__gz'] = []
    masks['mask_bpt__gz'] = []
    masks['mask_lines_dict__gLmz'] = []

    D = tupperware_none()
    D.pixels = {}
    D.pixels['SFRSD'] = []
    D.pixels['McorSD'] = []
    D.pixels['logOH_M13'] = []
    D.pixels['logOH_PP04'] = []
    D.pixels['m_type'] = []
    D.pixels['m_type_orig'] = []
    D.zones = {}
    D.zones['SFR'] = []
    D.zones['Mcor'] = []
    D.zones['SFRSD'] = []
    D.zones['McorSD'] = []
    D.zones['logOH_M13'] = []
    D.zones['logOH_PP04'] = []
    D.zones['m_type'] = []
    D.zones['m_type_orig'] = []
    D.radial_profiles = {}
    D.radial_profiles['aSFRSD'] = []
    D.radial_profiles['alogSFRSD'] = []
    D.radial_profiles['aSFRSD_std'] = []
    D.radial_profiles['alogSFRSD_std'] = []
    D.radial_profiles['aMcorSD'] = []
    D.radial_profiles['aMcorSD_std'] = []
    D.radial_profiles['alogMcorSD'] = []
    D.radial_profiles['alogMcorSD_std'] = []
    D.radial_profiles['alogOH_M13'] = []
    D.radial_profiles['alogOH_M13_std'] = []
    D.radial_profiles['alogOH_PP04'] = []
    D.radial_profiles['alogOH_PP04_std'] = []
    # D.radial_profiles['m_type'] = []
    # D.radial_profiles['m_type_orig'] = []

    for g in gals[:max_gals]:
        t_init_gal = time.clock()

        califaID = g

        pycasso_cube_file = pycasso_cube_dir + califaID + pycasso_cube_suffix
        eml_cube_file = eml_cube_dir + califaID + eml_cube_suffix
        gasprop_cube_file = gasprop_cube_dir + califaID + gasprop_cube_suffix

        # Load all data from fits
        K = fitsQ3DataCube(pycasso_cube_file)
        K.loadEmLinesDataCube(eml_cube_file)
        K.GP = GasProp(gasprop_cube_file)

        sit, verify = verify_files(K, califaID, True, True)
        if verify is not True:
            print '<<< ', califaID, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue

        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)

        mto = get_morfologia(K.califaID)[0]
        mt = my_morf(mto)
        aux = np.ma.masked_array([mto for i in range(K.N_zone)], mask=np.zeros((K.N_zone), dtype='bool'), dtype='object')
        D.zones['m_type_orig'].append(aux)
        aux = np.ma.masked_array(np.ones((K.N_zone)) * mt, mask=np.zeros((K.N_zone), dtype='bool'))
        D.zones['m_type'].append(aux)
        aux = np.ma.masked_array([mto for i in range(K.N_x * K.N_y)], mask=np.zeros((K.N_x * K.N_y), dtype='bool'), dtype='object')
        D.pixels['m_type_orig'].append(aux)
        aux = np.ma.masked_array(np.ones((K.N_x * K.N_y)) * mt, mask=np.zeros((K.N_x * K.N_y), dtype='bool'))
        D.pixels['m_type'].append(aux)

        aux = create_zones_masks_gal(K, tSF__T, **dflt_kw_create_zones_masks)
        masks['mask__gTz'].append(aux[0])
        masks['mask_syn__gTz'].append(aux[1])
        masks['mask_eml__gz'].append(aux[2])
        masks['mask_popx__gTz'].append(aux[3])
        masks['mask_tau_V__gz'].append(aux[4])
        masks['mask_residual__gz'].append(aux[5])
        masks['mask_tau_V_neb__gz'].append(aux[6])
        masks['mask_tau_V_neb_err__gz'].append(aux[7])
        masks['mask_EW_Hb__gz'].append(aux[8])
        masks['mask_whan__gz'].append(aux[9])
        masks['mask_bpt__gz'].append(aux[10])
        masks['mask_lines_dict__gLmz'].append(aux[11])

        # masks SNR bpt lines and synQC+min(popx,tauv) also...
        # means that only includes zones with:
        # synQC = ~(K.filterResidual(w2=4600))
        # frac_popx_young >= 0.05
        # tau_V >= 0.05
        m_aux = np.copy(masks['mask_syn__gTz'][-1][iT])
        m_aux = np.bitwise_or(m_aux, masks['mask_lines_dict__gLmz'][-1]['4861']['SNR'])
        m_aux = np.bitwise_or(m_aux, masks['mask_lines_dict__gLmz'][-1]['5007']['SNR'])
        m_aux = np.bitwise_or(m_aux, masks['mask_lines_dict__gLmz'][-1]['6563']['SNR'])
        m_aux = np.bitwise_or(m_aux, masks['mask_lines_dict__gLmz'][-1]['6583']['SNR'])

        SFR = np.ma.masked_array(calc_SFR(K, tSF__T[iT])[0])
        SFR[m_aux] = np.ma.masked
        D.zones['SFR'].append(SFR)
        D.zones['SFRSD'].append(SFR/K.zoneArea_pc2)
        SFRSD__yx = K.zoneToYX(SFR)
        D.pixels['SFRSD'].append(SFRSD__yx)
        D.radial_profiles['aSFRSD'].append(K.radialProfile(SFRSD__yx, R_bin__r, K.HLR_pix))
        D.radial_profiles['aSFRSD_std'].append(K.radialProfile(SFRSD__yx, R_bin__r, K.HLR_pix, mode='std'))
        D.radial_profiles['alogSFRSD'].append(K.radialProfile(np.ma.log10(SFRSD__yx), R_bin__r, K.HLR_pix))
        D.radial_profiles['alogSFRSD_std'].append(K.radialProfile(np.ma.log10(SFRSD__yx), R_bin__r, K.HLR_pix, mode='std'))

        Mcor = np.ma.masked_array(K.Mcor__z)
        Mcor[m_aux] = np.ma.masked
        D.zones['Mcor'].append(SFR)
        D.zones['McorSD'].append(SFR/K.zoneArea_pc2)
        McorSD__yx = K.zoneToYX(Mcor)
        D.pixels['McorSD'].append(McorSD__yx)
        D.radial_profiles['aMcorSD'].append(K.radialProfile(McorSD__yx, R_bin__r, K.HLR_pix))
        D.radial_profiles['aMcorSD_std'].append(K.radialProfile(McorSD__yx, R_bin__r, K.HLR_pix, mode='std'))
        D.radial_profiles['alogMcorSD'].append(K.radialProfile(np.ma.log10(McorSD__yx), R_bin__r, K.HLR_pix))
        D.radial_profiles['alogMcorSD_std'].append(K.radialProfile(np.ma.log10(McorSD__yx), R_bin__r, K.HLR_pix, mode='std'))

        Zneb_M13 = np.ma.masked_array(K.EL.Zneb_M13__z)
        Zneb_M13[m_aux] = np.ma.masked
        D.zones['logOH_M13'].append(Zneb_M13)
        Zneb_M13__yx = K.zoneToYX(Zneb_M13, extensive=False)
        D.pixels['logOH_M13'].append(Zneb_M13__yx)
        D.radial_profiles['alogOH_M13'].append(K.radialProfile(Zneb_M13__yx, R_bin__r, K.HLR_pix))
        D.radial_profiles['alogOH_M13_std'].append(K.radialProfile(Zneb_M13__yx, R_bin__r, K.HLR_pix, mode='std'))

        Zneb_PP04 = np.ma.masked_array(K.EL.Zneb_PP04__z)
        Zneb_PP04[m_aux] = np.ma.masked
        D.zones['logOH_PP04'].append(Zneb_PP04)
        Zneb_PP04__yx = K.zoneToYX(Zneb_PP04, extensive=False)
        D.pixels['logOH_PP04'].append(Zneb_PP04__yx)
        D.radial_profiles['alogOH_PP04'].append(K.radialProfile(Zneb_PP04__yx, R_bin__r, K.HLR_pix))
        D.radial_profiles['alogOH_PP04_std'].append(K.radialProfile(Zneb_PP04__yx, R_bin__r, K.HLR_pix, mode='std'))

        K.close()
        t_fin_gal = time.clock()
        print 'gal: ', califaID, ' time: ', t_fin_gal - t_init_gal, 's'

t_ini_pickle = time.clock()
pickle_ma_dict(to_pickle=D.zones, file_pref='zones-')
pickle_ma_ravel_dict(to_pickle=D.pixels, file_pref='pixels-')
pickle_ma_dict(to_pickle=D.radial_profiles, file_pref='radprof-')
t_fin_pickle = time.clock()
print 'pickles time: ', t_fin_gal - t_init_gal, 's'
