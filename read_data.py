#!/usr/bin/python
#
# Lacerda@Saco - 14/Mar/2016
#
from CALIFAUtils.scripts import calc_alogZ_Stuff
from CALIFAUtils.scripts import create_masks_gal 
from pystarlight.util.base import StarlightBase
from CALIFAUtils.scripts import read_gal_cubes
from CALIFAUtils.scripts import get_morfologia
from pystarlight.util.constants import L_sun
from pystarlight.util import redenninglaws
from CALIFAUtils.objects import stack_gals
from CALIFAUtils.scripts import sort_gals
from CALIFAUtils.scripts import debug_var
from CALIFAUtils.scripts import calc_SFR
from CALIFAUtils.scripts import calc_xY
from CALIFAUtils.scripts import my_morf
import argparse as ap
import numpy as np
import time
import h5py
import sys

class ALLGals(stack_gals):
    def __init__(self, N_gals, N_T, N_U):
        stack_gals.__init__(self)
        self.N_gals = N_gals
        self.N_T = N_T
        self.N_U = N_U
        self.range_T = xrange(self.N_T)
        self.range_U = xrange(self.N_U)
        self._init_arrays()
        
    def mask_gal(self, iGal):
        for v in self.__dict__.keys():
            if isinstance(self.__dict__[v], np.ma.core.MaskedArray):
                self.__dict__[v][..., iGal] = np.ma.masked
        
    def _init_arrays(self):
        N_gals = self.N_gals
        N_T = self.N_T
        N_U = self.N_U

        # Gal
        self.new1d('zone_area_pc2__g')
        self.new1d('zone_dist_HLR__g')
        self.califaIDs = np.ma.masked_all((N_gals), dtype = '|S5')
        self.N_zones__g = np.ma.masked_all((N_gals), dtype = int)
        self.morfType_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.ba_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.ba_PyCASSO_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.Mr_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.ur_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.parsecPerPixel__g = np.ma.masked_all((N_gals), dtype = np.float_)

        # Synthetic
        self.new1d('Mcor__g')
        self.new1d('McorSD__g')
        self.new1d('at_flux__g')
        self.new1d('at_mass__g')
        self.new2d('x_Y__Tg', N_T)
        self.new2d_masked('tau_V__Tg', N_T)
        self.new2d_masked('SFR__Tg', N_T)
        self.new2d_masked('SFRSD__Tg', N_T)
        self.new2d_masked('Mcor__Tg', N_T)
        self.new2d_masked('McorSD__Tg', N_T)
        self.new2d_masked('at_flux__Tg', N_T)
        self.new2d_masked('at_mass__Tg', N_T)
        self.new2d_masked('alogZ_mass__Ug', N_U)
        self.new2d_masked('alogZ_flux__Ug', N_U)
        self.at_flux_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.Mcor_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.McorSD_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.HLR_pix_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.HMR_pix_GAL__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.alogZ_mass_GAL__Ug = np.ma.masked_all((N_U, N_gals), dtype = np.float_)
        self.alogZ_flux_GAL__Ug = np.ma.masked_all((N_U, N_gals), dtype = np.float_)
        self.integrated_x_Y__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.float_)
        self.integrated_SFR__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.float_)
        self.integrated_SFRSD__Tg = np.ma.masked_all((N_T, N_gals), dtype = np.float_)
        self.integrated_tau_V__g = np.ma.masked_all((N_gals), dtype = np.float_)
        
        # Nebular
        self.new1d_masked('tau_V_neb__g')
        self.new1d_masked('tau_V_neb_err__g')
        self.new1d_masked('SFR_Ha__g')
        self.new1d_masked('SFRSD_Ha__g')
        self.new1d_masked('F_obs_Hb__g')
        self.new1d_masked('F_obs_O3__g')
        self.new1d_masked('F_obs_Ha__g')
        self.new1d_masked('F_obs_N2__g')
        self.new1d_masked('eF_obs_Hb__g')
        self.new1d_masked('eF_obs_O3__g')
        self.new1d_masked('eF_obs_Ha__g')
        self.new1d_masked('eF_obs_N2__g')
        self.new1d_masked('F_int_Ha__g')
        self.new1d_masked('F_int_Hb__g')
        self.new1d_masked('F_int_N2__g')
        self.new1d_masked('F_int_O3__g')
        self.new1d_masked('baseline_Hb__g')
        self.new1d_masked('baseline_O3__g')
        self.new1d_masked('baseline_Ha__g')
        self.new1d_masked('baseline_N2__g')
        self.new1d_masked('L_int_Ha__g')
        self.new1d_masked('L_int_Ha_err__g')
        self.new1d_masked('L_obs_Ha__g')
        self.new1d_masked('L_obs_Ha_err__g')
        self.new1d_masked('EW_Ha__g')
        self.new1d_masked('EW_Hb__g')
        self.new1d_masked('logOH_M13__g')
        self.integrated_EW_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_EW_Hb__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_tau_V_neb__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_tau_V_neb_err__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_obs_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_obs_Hb__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_obs_O3__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_obs_N2__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_eF_obs_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_eF_obs_Hb__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_eF_obs_O3__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_eF_obs_N2__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_int_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_int_Hb__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_int_O3__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_F_int_N2__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_baseline_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_baseline_Hb__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_baseline_O3__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_baseline_N2__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_L_int_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_L_obs_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_SFR_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_SFRSD_Ha__g = np.ma.masked_all((N_gals), dtype = np.float_)
        self.integrated_logOH_M13__g = np.ma.masked_all((N_gals), dtype = np.float_)
        
    def stack_zones_data(self):
        self.stack()
        
    def integrated_mask(self):
        aux = np.less(self.integrated_tau_V_neb__g, 0)
        self.integrated_tau_V_neb__g[aux] = np.ma.masked 

    def create_dict_h5(self):
        D = {}
        for v in self.__dict__.keys():
            if v[0] != '_':
                suffix = v.split('_')[-1]
                if isinstance(self.__dict__[v], np.ma.core.MaskedArray):
                    tmp_data = {'masked/data/%s' % v : self.__dict__[v].data}
                    tmp_mask = {'masked/mask/%s' % v : self.__dict__[v].mask}
                else:
                    if suffix == 'Tg':
                        tmp_data = {'masked/data/%s/%d' % (v, i) : self.__dict__[v][i].data for i in self.range_T}
                        tmp_mask = {'masked/mask/%s/%d' % (v, i) : self.__dict__[v][i].mask for i in self.range_T}
                    elif suffix == 'Ug':
                        tmp_data = {'masked/data/%s/%d' % (v, i) : self.__dict__[v][i].data for i in self.range_U}
                        tmp_mask = {'masked/mask/%s/%d' % (v, i) : self.__dict__[v][i].mask for i in self.range_U}
                    else:
                        tmp_data = {}
                        tmp_mask = {}
                D.update(tmp_data)
                D.update(tmp_mask)
        return D      

def parse_args(args_str):
    default_args = {
        'debug' : False,
        'spiral' : False,
        'hdf5' : 'output.h5',
        'whanSF' : None,
        'rgbcuts' : False,
        'underS06' : False,
        'nolinecuts' : False,
        'filter_residual' : False,
        'minpopx' : np.finfo(np.float_).min,
        'mintauv' : np.finfo(np.float_).min,
        'minEWHb' : np.finfo(np.float_).min,
        'mintauvneb' : np.finfo(np.float_).min,
        'maxtauvneberr' : np.finfo(np.float_).max,
        'v_run' : -1,
        'minSNR' : 3,
        'minSNRHb' : 3,
        'gals' : '/Users/lacerda/CALIFA/listv20_q050.d15a.txt',
        'pycasso_cube_dir' : '/Users/lacerda/CALIFA/gal_fits/v20_q050.d15a',
        'pycasso_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.fits',
        'eml_cube_dir' : '/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a',
        'eml_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.fits', 
        'gasprop_cube_dir' : '/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a/prop',
        'gasprop_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.GasProp.fits', 
    }
    
    parser = ap.ArgumentParser(description = '%s' % args_str)
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--nolinecuts' ,
                        action = 'store_true',
                        default = default_args['nolinecuts'])
    parser.add_argument('--spiral', '-S',
                        action = 'store_true',
                        default = default_args['spiral'])
    parser.add_argument('--filter_residual', '-R',
                        action = 'store_true',
                        default = default_args['filter_residual'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['hdf5'])
    parser.add_argument('--v_run',
                        metavar = 'INT',
                        type = int,
                        default = default_args['v_run'])
    parser.add_argument('--underS06',
                        action = 'store_true',
                        default = default_args['underS06'])
    parser.add_argument('--whanSF',
                        metavar = 'INT',
                        type = int,
                        #action = 'store_true',
                        default = default_args['whanSF'])
    parser.add_argument('--rgbcuts',
                        action = 'store_true',
                        default = default_args['rgbcuts'])
    parser.add_argument('--minpopx',
                        help = 'Negative to disable mask in popx',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minpopx'])
    parser.add_argument('--minEWHb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minEWHb'])    
    parser.add_argument('--minSNR',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minSNR'])    
    parser.add_argument('--minSNRHb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minSNRHb'])    
    parser.add_argument('--mintauv',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['mintauv'])
    parser.add_argument('--mintauvneb',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['mintauvneb'])
    parser.add_argument('--maxtauvneberr',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['maxtauvneberr'])
    parser.add_argument('--gals', '-G',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['gals'])
    parser.add_argument('--pycasso_cube_dir',
                        metavar = 'DIR',
                        type = str,
                        default = default_args['pycasso_cube_dir'])
    parser.add_argument('--pycasso_cube_suffix',
                        metavar = 'SUFFIX',
                        type = str,
                        default = default_args['pycasso_cube_suffix'])
    parser.add_argument('--eml_cube_dir',
                        metavar = 'DIR',
                        type = str,
                        default = default_args['eml_cube_dir'])
    parser.add_argument('--eml_cube_suffix',
                        metavar = 'SUFFIX',
                        type = str,
                        default = default_args['eml_cube_suffix'])
    parser.add_argument('--gasprop_cube_dir',
                        metavar = 'DIR',
                        type = str,
                        default = default_args['gasprop_cube_dir'])
    parser.add_argument('--gasprop_cube_suffix',
                        metavar = 'SUFFIX',
                        type = str,
                        default = default_args['gasprop_cube_suffix'])
    
    args = parser.parse_args()

    if args.nolinecuts:
        args.rgbcuts = False

    return args

def verify_files(K, califaID, EL = True, GP = True):
    if K is None:
        print '<<< %s galaxy: miss files' % califaID
        return 0, False
    if EL == True and K.EL is None:
        print '<<< %s galaxy: miss EmLines files' % califaID
        return 1, False
        if K.EL.flux[0, :].sum() == 0.:
            print '<<< %s EmLines FITS problem' % califaID
            return 2, False
    if GP is True and K.GP._hdulist is None:
        print '<<< %s galaxy: miss gasprop file' % califaID
        return 2, False
    # Problem in FITS file
    return 0, True       

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse arguments 
    args = parse_args(sys.argv[0])
    debug_var(True, args = args.__dict__)    
    
    Zsun = 0.019

    # Reading galaxies file,
    gals, _ = sort_gals(gals = args.gals, order = 1)
    N_gals = len(gals)
    maxGals = None
    if args.debug:
        maxGals = 10
        if N_gals > maxGals:
            N_gals = maxGals

    # SFR-time-scale array (index __T)
    base = StarlightBase('/Users/lacerda/LOCAL/data/BASE.CALIFA.gsd6.h5', 'gsd6e', hdf5 = True)
    #tSF__T = np.array([0.032 , 0.3 , 1.5, 14.2]) * 1.e9
    tSF__T = np.asarray(base.ageBase)
    #tSF__T = np.array([1, 3.2, 10, 100]) * 1e7
    N_T = len(tSF__T)

    # Z-time-scale array (index __U).
    tZ__U = np.array([1.0 , 2.0 , 5.0 , 8.0 , 11.3 , 14.2]) * 1.e9
    N_U = len(tZ__U)

    ALL = ALLGals(N_gals, N_T, N_U)

    ALL.new2d('msk__Tg', N_T)
    ALL.new2d('msk_syn__Tg', N_T)
    ALL.new2d('msk_popx__Tg', N_T)
    ALL.new2d('msk_lines_dict__Lg', 4)
    ALL.new1d('msk_eml__g')
    ALL.new1d('msk_tau_V__g')
    ALL.new1d('msk_residual__g')
    ALL.new1d('msk_tau_V_neb__g')
    ALL.new1d('msk_tau_V_neb_err__g')
    ALL.new1d('msk_EW_Hb__g')
    ALL.new1d('msk_whan__g')
    ALL.new1d('msk_bpt__g')

    for iGal, gal in enumerate(gals[0:maxGals]):
        K = read_gal_cubes(gal, debug = True, 
                           pycasso_cube_dir = args.pycasso_cube_dir, pycasso_cube_suffix = args.pycasso_cube_suffix,
                           eml_cube_dir = args.eml_cube_dir, eml_cube_suffix = args.eml_cube_suffix,
                           gasprop_cube_dir = args.gasprop_cube_dir, gasprop_cube_suffix = args.gasprop_cube_suffix,
                           )

        t_init_gal = time.clock()
        califaID = gals[iGal] 
        
        # Saving for later :D        
        ALL.califaIDs[iGal] = califaID
        
        sit, verify = verify_files(K, califaID, EL = True, GP = True)
        
        if verify is not True:
            ALL.mask_gal(iGal)
            print '<<< ', califaID, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue

        tipos, tipo, tipo_m, tipo_p = get_morfologia(califaID)
        my_type = my_morf(tipos)
        
        # Only spiral
        if args.spiral and my_type < 0:  # between Sa ... Sd
            ALL.mask_gal(iGal)
            K.GP.close()
            K.EL.close()
            K.close()
            print '<<< %s galaxy: is not a spiral (type: %f (%s))' % (califaID, my_type, tipos) 
            continue

        N_zone = K.N_zone

        # Setup elliptical-rings geometry
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        # Saving for later :D
        ALL.N_zones__g[iGal] = N_zone
        ALL.ba_PyCASSO_GAL__g[iGal] = ba
        
        print '>>> Doing %d %s %s (%d)' % (iGal, califaID, tipos, my_type)
        
        # zone distance in HLR
        #zoneDistHLR = np.sqrt((K.zonePos['x'] - K.x0) ** 2. + (K.zonePos['y'] - K.y0) ** 2.) / K.HLR_pix
        zone_distance_HLR = K.zoneDistance_HLR
        zone_area_pc2 = K.zoneArea_pc2
        ba_GAL = np.float(K.masterListData['ba'])
        Mr_GAL = np.float(K.masterListData['Mr'])
        ur_GAL = np.float(K.masterListData['u-r'])
        
        # Saving for later :D
        ALL.append1d('zone_dist_HLR__g', zone_distance_HLR)
        ALL.append1d('zone_area_pc2__g', zone_area_pc2)
        ALL.ba_GAL__g[iGal] = ba_GAL
        ALL.Mr_GAL__g[iGal] = Mr_GAL
        ALL.ur_GAL__g[iGal] = ur_GAL
        ALL.morfType_GAL__g[iGal] = my_type

        # masks: more info. in create_masks() 
        mask__Tz, mask_syn__Tz, mask_eml__z, mask_popx__Tz, \
        mask_tau_V__z, mask_residual__z, mask_tau_V_neb__z, \
        mask_tau_V_neb_err__z, mask_EW_Hb__z, mask_whan__z, \
        mask_bpt__z, mask_lines_dict__Lz = create_masks_gal(K = K, tSF__T = tSF__T, args = args, debug = False)
        
        ALL.append1d('msk_eml__g', mask_eml__z)
        ALL.append1d('msk_tau_V__g', mask_tau_V__z)
        ALL.append1d('msk_residual__g', mask_residual__z)
        ALL.append1d('msk_tau_V_neb__g', mask_tau_V_neb__z)
        ALL.append1d('msk_tau_V_neb_err__g', mask_tau_V_neb_err__z)
        ALL.append1d('msk_EW_Hb__g', mask_EW_Hb__z)
        ALL.append1d('msk_whan__g', mask_whan__z)
        ALL.append1d('msk_bpt__g', mask_bpt__z)
        ALL.append2d('msk_lines_dict__Lg', 0, mask_lines_dict__Lz['4861'])
        ALL.append2d('msk_lines_dict__Lg', 1, mask_lines_dict__Lz['5007'])
        ALL.append2d('msk_lines_dict__Lg', 2, mask_lines_dict__Lz['6563'])
        ALL.append2d('msk_lines_dict__Lg', 3, mask_lines_dict__Lz['6583'])
        for iT in xrange(N_T):
            ALL.append2d('msk__Tg', iT, mask__Tz[iT])
            ALL.append2d('msk_syn__Tg', iT, mask_syn__Tz[iT])
            ALL.append2d('msk_popx__Tg', iT, mask_popx__Tz[iT])

        ####################################################
        ####### STARLIGHT ##################################
        ####################################################
        AVtotauV = 1. / (np.log10(np.exp(1)) / 0.4)
        ALL.integrated_tau_V__g[iGal] = K.integrated_keywords['A_V'] * AVtotauV 
        
        # Compute galaxy-wide mu (cf eq 2 in GD14) - following Andre's tip.
        Mcor__z = K.Mcor__z
        McorSD__z = K.Mcor__z / K.zoneArea_pc2
        Mcor_GAL = K.Mcor_tot.sum()
        McorSD_GAL = K.McorSD__yx.mean()
        at_flux__z = K.at_flux__z
        at_mass__z = K.at_mass__z
        HLR_pix = K.HLR_pix
        HMR_pix = K.getHalfRadius(K.McorSD__yx)
        parsecPerPixel = K.parsecPerPixel
        
        # Compute & store galaxy-wide at_flux
        numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
        at_flux_GAL = numerator__z.sum() / denominator__z.sum()
        
        # Saving for later :D
        ALL.append1d('Mcor__g', Mcor__z)
        ALL.append1d('McorSD__g', McorSD__z)
        ALL.append1d('at_flux__g', at_flux__z)
        ALL.append1d('at_mass__g', at_mass__z)
        ALL.Mcor_GAL__g[iGal] = Mcor_GAL
        ALL.McorSD_GAL__g[iGal] = McorSD_GAL
        ALL.HLR_pix_GAL__g[iGal] = HLR_pix
        ALL.HMR_pix_GAL__g[iGal] = HMR_pix
        ALL.parsecPerPixel__g[iGal] = parsecPerPixel
        ALL.at_flux_GAL__g[iGal] = at_flux_GAL
        
        # Composition by StarForming time scale
        for iT, tSF in enumerate(tSF__T):
            Mcor__z = np.ma.masked_array(K.Mcor__z)
            McorSD__z = np.ma.masked_array(K.Mcor__z / K.zoneArea_pc2)
            tau_V__z = np.ma.masked_array(K.tau_V__z)
            at_flux__z = np.ma.masked_array(K.at_flux__z)
            at_mass__z = np.ma.masked_array(K.at_mass__z)

            x_Y__z, integrated_x_Y = calc_xY(K, tSF)
            aux = calc_SFR(K, tSF)
            SFR__z = np.ma.masked_array(aux[0])
            SFRSD__z = np.ma.masked_array(aux[1])

            tau_V__z[mask_syn__Tz[iT]] = np.ma.masked
            SFR__z[mask_syn__Tz[iT]] = np.ma.masked
            SFRSD__z[mask_syn__Tz[iT]] = np.ma.masked
            Mcor__z[mask_syn__Tz[iT]] = np.ma.masked
            McorSD__z[mask_syn__Tz[iT]] = np.ma.masked
            at_flux__z[mask_syn__Tz[iT]] = np.ma.masked
            at_mass__z[mask_syn__Tz[iT]] = np.ma.masked
            integrated_SFR = SFR__z.sum()

            # Saving for later :D                
            ALL.append2d('x_Y__Tg', iT, x_Y__z)
            ALL.append2d_masked('tau_V__Tg', iT, tau_V__z.data, tau_V__z.mask)
            ALL.append2d_masked('Mcor__Tg', iT, Mcor__z.data, Mcor__z.mask)
            ALL.append2d_masked('McorSD__Tg', iT, McorSD__z.data, McorSD__z.mask)
            ALL.append2d_masked('SFR__Tg', iT, SFR__z.data, SFR__z.mask)
            ALL.append2d_masked('SFRSD__Tg', iT, SFRSD__z.data, SFRSD__z.mask)
            ALL.append2d_masked('at_flux__Tg', iT, at_flux__z.data, at_flux__z.mask)
            ALL.append2d_masked('at_mass__Tg', iT, at_mass__z.data, at_mass__z.mask)
            ALL.integrated_x_Y__Tg[iT, iGal] = integrated_x_Y
            ALL.integrated_SFR__Tg[iT, iGal] = integrated_SFR
            ALL.integrated_SFRSD__Tg[iT, iGal] = integrated_SFR / K.zoneArea_pc2.sum()
                        
        for iU, tZ in enumerate(tZ__U):
            aux = calc_alogZ_Stuff(K, tZ, args.minpopx, None)
            alogZ_mass__z = aux[0]
            alogZ_flux__z = aux[1]
            alogZ_mass_GAL = aux[2]
            alogZ_flux_GAL = aux[3]
            
            # Saving for later :D
            ALL.append2d_masked('alogZ_mass__Ug', iU, alogZ_mass__z.data, alogZ_mass__z.mask)
            ALL.append2d_masked('alogZ_flux__Ug', iU, alogZ_flux__z.data, alogZ_flux__z.mask)
            ALL.alogZ_mass_GAL__Ug[iU, iGal] = alogZ_mass_GAL
            ALL.alogZ_flux_GAL__Ug[iU, iGal] = alogZ_flux_GAL
        ####################################################
        ####################################################
        ####################################################    
        
        ####################################################
        ######## EmLines ###################################
        ####################################################
        Hb_central_wl = '4861'
        O3_central_wl = '5007'
        Ha_central_wl = '6563'
        N2_central_wl = '6583'
        i_Hb = K.EL.lines.index(Hb_central_wl)
        i_O3 = K.EL.lines.index(O3_central_wl)
        i_Ha = K.EL.lines.index(Ha_central_wl)
        i_N2 = K.EL.lines.index(N2_central_wl)

        ########## tau_V #########
        mask_tau_V_neb_aux__z = np.zeros((K.N_zone), dtype = np.bool_)
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_lines_dict__Lz[Ha_central_wl])
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_lines_dict__Lz[Hb_central_wl])
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_tau_V_neb__z)
        mask_tau_V_neb_aux__z = np.bitwise_or(mask_tau_V_neb_aux__z, mask_tau_V_neb_err__z)
        
        print '# N_mask_tau_V_neb_aux (Ha+Hb+tauVNeb+tauVNebErr): ', mask_tau_V_neb_aux__z.astype(int).sum()

        tau_V_neb__z = np.ma.masked_array(K.EL.tau_V_neb__z, mask = mask_tau_V_neb_aux__z)
        tau_V_neb_err__z = np.ma.masked_array(K.EL.tau_V_neb_err__z, mask = mask_tau_V_neb_aux__z)

        # Saving for later :D
        ALL.append1d_masked('tau_V_neb__g', tau_V_neb__z.data, tau_V_neb__z.mask)
        ALL.append1d_masked('tau_V_neb_err__g', tau_V_neb_err__z.data, tau_V_neb_err__z.mask)
        ALL.integrated_tau_V_neb__g[iGal] = K.EL.integrated_tau_V_neb
        ALL.integrated_tau_V_neb_err__g[iGal] = K.EL.integrated_tau_V_neb_err
        ##########################

        ########### EW ###########
        EW_Ha__z = np.ma.masked_array(K.EL.EW[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        EW_Hb__z = np.ma.masked_array(K.EL.EW[i_Hb, :], mask = mask_lines_dict__Lz[Hb_central_wl])
        baseline_Hb__z = K.EL.baseline[i_Hb]
        baseline_Ha__z = K.EL.baseline[i_Ha]

        # Saving for later :D        
        ALL.append1d_masked('EW_Hb__g', EW_Hb__z.data, EW_Hb__z.mask)
        ALL.append1d_masked('EW_Ha__g', EW_Ha__z.data, EW_Ha__z.mask)
        ALL.append1d_masked('baseline_Hb__g', K.EL.baseline[i_Hb], EW_Hb__z.mask)
        ALL.append1d_masked('baseline_Ha__g', K.EL.baseline[i_Ha], EW_Ha__z.mask)
        ALL.integrated_EW_Ha__g[iGal] = K.EL.integrated_EW[i_Ha]
        ALL.integrated_EW_Hb__g[iGal] = K.EL.integrated_EW[i_Hb]
        ALL.integrated_baseline_Ha__g[iGal] = K.EL.integrated_baseline[i_Ha]
        ALL.integrated_baseline_Hb__g[iGal] = K.EL.integrated_baseline[i_Hb]
        ##########################
                
        #### intrinsic Ha Lum ####
        q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
        expqtau = [ np.ma.exp(qcard * tau_V_neb__z) for qcard in q ]
        integrated_expqtau = [ np.ma.exp(qcard * K.EL.integrated_tau_V_neb) for qcard in q ]
        F_obs_Ha__z = np.ma.masked_array(K.EL.flux[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        L_obs__Lz = K.EL._F_to_L(K.EL.flux) / L_sun
        L_obs_err__Lz = K.EL._F_to_L(K.EL.eflux) / L_sun        
        L_obs_Ha__z = np.ma.masked_array(L_obs__Lz[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        L_obs_Hb__z = np.ma.masked_array(L_obs__Lz[i_Hb, :], mask = mask_lines_dict__Lz[Hb_central_wl])
        L_obs_Ha_err__z = np.ma.masked_array(L_obs_err__Lz[i_Ha, :], mask = mask_lines_dict__Lz[Ha_central_wl])
        L_obs_Hb_err__z = np.ma.masked_array(L_obs_err__Lz[i_Hb, :], mask = mask_lines_dict__Lz[Hb_central_wl])
        L_obs_HaHb__z = L_obs_Ha__z / L_obs_Hb__z
        # L_int_Ha__Lz intrinsic Ha luminosity 
        # For the zones where I don't have values for tau_V_neb I don't correct the Lum_Ha
        L_int_Ha__z = np.where(~mask_tau_V_neb_aux__z, L_obs_Ha__z * expqtau[2], L_obs_Ha__z)
        L_int_Ha__z = np.ma.masked_array(L_int_Ha__z, mask = mask_tau_V_neb_aux__z)
        # L_int_Ha_err__Lz intrinsic Ha luminosity propagated error
        qq = q[2] / (q[0] - q[2])
        a = L_obs_Ha_err__z
        b = qq * L_obs_HaHb__z * L_obs_Hb_err__z
        L_int_Ha_err__z = np.where(~mask_tau_V_neb_aux__z, L_obs_Ha_err__z, expqtau[2] * np.sqrt(a ** 2.0 + b ** 2.0))
        L_int_Ha_err__z = np.ma.masked_array(L_int_Ha_err__z, mask = mask_tau_V_neb_aux__z)
        integrated_L_obs_Ha = K.EL._F_to_L(K.EL.integrated_flux[i_Ha]) / L_sun
        integrated_L_int_Ha = integrated_L_obs_Ha * integrated_expqtau[2]
        
        # Saving for later :D
        ALL.append1d_masked('F_obs_Ha__g', F_obs_Ha__z.data, F_obs_Ha__z.mask)
        ALL.append1d_masked('L_obs_Ha__g', L_obs_Ha__z.data, L_obs_Ha__z.mask)
        ALL.append1d_masked('L_obs_Ha_err__g', L_obs_Ha_err__z.data, L_obs_Ha_err__z.mask)
        ALL.append1d_masked('L_int_Ha__g', L_int_Ha__z.data, L_int_Ha__z.mask)
        ALL.append1d_masked('L_int_Ha_err__g', L_int_Ha_err__z.data, L_int_Ha_err__z.mask)
        ALL.integrated_L_obs_Ha__g[iGal] = integrated_L_obs_Ha
        ALL.integrated_L_int_Ha__g[iGal] = integrated_L_int_Ha
        ##########################
        
        ###### OTH BPT LINES #####
        F_obs_Hb__z = np.ma.masked_array(K.EL.Hb_obs__z, mask = mask_lines_dict__Lz[Hb_central_wl])
        F_obs_O3__z = np.ma.masked_array(K.EL.O3_obs__z, mask = mask_lines_dict__Lz[O3_central_wl])
        F_obs_N2__z = np.ma.masked_array(K.EL.N2_obs__z, mask = mask_lines_dict__Lz[N2_central_wl])
        F_int_Hb__z = np.where(~mask_tau_V_neb_aux__z, F_obs_Hb__z * expqtau[0], F_obs_Hb__z)
        F_int_O3__z = np.where(~mask_tau_V_neb_aux__z, F_obs_O3__z * expqtau[1], F_obs_O3__z)
        F_int_Ha__z = np.where(~mask_tau_V_neb_aux__z, F_obs_Ha__z * expqtau[2], F_obs_Ha__z)
        F_int_N2__z = np.where(~mask_tau_V_neb_aux__z, F_obs_N2__z * expqtau[3], F_obs_N2__z)
        baseline_O3__z = K.EL.baseline[i_O3]
        baseline_N2__z = K.EL.baseline[i_N2]
        maskNone = np.zeros((K.N_zone), dtype = np.bool_)
        eF_obs_Hb__z = np.ma.masked_array(K.EL.eflux[i_Hb], mask = maskNone)
        eF_obs_O3__z = np.ma.masked_array(K.EL.eflux[i_O3], mask = maskNone)
        eF_obs_Ha__z = np.ma.masked_array(K.EL.eflux[i_Ha], mask = maskNone)
        eF_obs_N2__z = np.ma.masked_array(K.EL.eflux[i_N2], mask = maskNone)
        integrated_F_obs_Hb = K.EL.integrated_flux[i_Hb]
        integrated_F_obs_O3 = K.EL.integrated_flux[i_O3]
        integrated_F_obs_Ha = K.EL.integrated_flux[i_Ha] 
        integrated_F_obs_N2 = K.EL.integrated_flux[i_N2]
        integrated_F_int_Hb = K.EL.integrated_flux[i_Hb] * integrated_expqtau[0]
        integrated_F_int_O3 = K.EL.integrated_flux[i_O3] * integrated_expqtau[1]
        integrated_F_int_Ha = K.EL.integrated_flux[i_Ha] * integrated_expqtau[2] 
        integrated_F_int_N2 = K.EL.integrated_flux[i_N2] * integrated_expqtau[3]
        
        # Saving for later :D
        ALL.append1d_masked('F_obs_Hb__g', F_obs_Hb__z.data, F_obs_Hb__z.mask)
        ALL.append1d_masked('F_obs_O3__g', F_obs_O3__z.data, F_obs_O3__z.mask)
        ALL.append1d_masked('F_obs_N2__g', F_obs_N2__z.data, F_obs_N2__z.mask)
        ALL.append1d_masked('eF_obs_Hb__g', eF_obs_Hb__z.data, eF_obs_Hb__z.mask)
        ALL.append1d_masked('eF_obs_O3__g', eF_obs_O3__z.data, eF_obs_O3__z.mask)
        ALL.append1d_masked('eF_obs_Ha__g', eF_obs_Ha__z.data, eF_obs_Ha__z.mask)
        ALL.append1d_masked('eF_obs_N2__g', eF_obs_N2__z.data, eF_obs_N2__z.mask)
        ALL.append1d_masked('F_int_Hb__g', F_int_Hb__z, F_obs_Hb__z.mask)
        ALL.append1d_masked('F_int_O3__g', F_int_O3__z, F_obs_O3__z.mask)
        ALL.append1d_masked('F_int_Ha__g', F_int_Ha__z, F_obs_Ha__z.mask)
        ALL.append1d_masked('F_int_N2__g', F_int_N2__z, F_obs_N2__z.mask)
        ALL.append1d_masked('baseline_O3__g', K.EL.baseline[i_O3])
        ALL.append1d_masked('baseline_N2__g', K.EL.baseline[i_N2])
        ALL.integrated_F_obs_Hb__g[iGal] = integrated_F_obs_Hb
        ALL.integrated_F_obs_O3__g[iGal] = integrated_F_obs_O3
        ALL.integrated_F_obs_Ha__g[iGal] = integrated_F_obs_Ha 
        ALL.integrated_F_obs_N2__g[iGal] = integrated_F_obs_N2
        ALL.integrated_eF_obs_Hb__g[iGal] = K.EL.integrated_eflux[i_Hb]
        ALL.integrated_eF_obs_O3__g[iGal] = K.EL.integrated_eflux[i_O3]
        ALL.integrated_eF_obs_Ha__g[iGal] = K.EL.integrated_eflux[i_Ha]
        ALL.integrated_eF_obs_N2__g[iGal] = K.EL.integrated_eflux[i_N2]
        ALL.integrated_F_int_Hb__g[iGal] = integrated_F_int_Hb
        ALL.integrated_F_int_O3__g[iGal] = integrated_F_int_O3
        ALL.integrated_F_int_Ha__g[iGal] = integrated_F_int_Ha 
        ALL.integrated_F_int_N2__g[iGal] = integrated_F_int_N2
        ##########################

        #### SFR and SigmaSFR ####
        # 3.13 M_sun/yr was calculated using BC03 + Padova1994 + Salpeter
        k_SFR = 3.13
        SFR_Ha__z = np.ma.masked_array(k_SFR * L_int_Ha__z.data / 1e8, mask = L_int_Ha__z.mask)
        SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
        integrated_SFR_Ha = k_SFR * integrated_L_int_Ha / 1e8

        # Saving for later :D
        ALL.append1d_masked('SFR_Ha__g', SFR_Ha__z.data, SFR_Ha__z.mask)
        ALL.append1d_masked('SFRSD_Ha__g', SFRSD_Ha__z.data, SFRSD_Ha__z.mask)
        ALL.integrated_SFR_Ha__g[iGal] = integrated_SFR_Ha 
        ALL.integrated_SFRSD_Ha__g[iGal] = integrated_SFR_Ha / K.zoneArea_pc2.sum()
        ####################################################
        ####################################################
        ####################################################
        
        # M13 Zneb calib.
        mask_Zneb_aux__z = np.zeros((K.N_zone), dtype = np.bool_)
        mask_Zneb_aux__z = np.bitwise_or(mask_Zneb_aux__z, mask_tau_V_neb_aux__z)
        mask_Zneb_aux__z = np.bitwise_or(mask_Zneb_aux__z, mask_lines_dict__Lz[O3_central_wl])
        mask_Zneb_aux__z = np.bitwise_or(mask_Zneb_aux__z, mask_lines_dict__Lz[N2_central_wl])
        print '# N_mask_Zneb_aux__z (Ha+Hb+tauVNeb+tauVNebErr+O3+N2): ', mask_Zneb_aux__z.astype(int).sum()
        logOH__z = K.EL.Zneb_M13__z.copy()
        logOH__z[mask_Zneb_aux__z] = np.ma.masked
        ALL.append1d_masked('logOH_M13__g', logOH__z.data, logOH__z.mask)
        ALL.integrated_logOH_M13__g[iGal] = K.EL.integrated_Zneb_M13 
        
        K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (califaID, time.clock() - t_init_gal)

    t_init_stack = time.clock()        
    ALL.stack_zones_data()
    ALL.integrated_mask()
    print 'time stack: %.2f' % (time.clock() - t_init_stack)
