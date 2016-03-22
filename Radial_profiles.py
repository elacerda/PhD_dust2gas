#!/usr/bin/python
#
# Lacerda@Saco - 21/Mar/2016
#
from pystarlight.util.constants import L_sun
from pystarlight.util import redenninglaws
from CALIFAUtils.scripts import debug_var
from tables_description import zone_neb
from CALIFAUtils.objects import GasProp
from tables_description import zone_SF
from pycasso import fitsQ3DataCube
from tables_description import tSF
import argparse as ap
import tables as tbl
import numpy as np
import time
import sys

def parser_args(args_str):
    default_args = {
        'debug' : False,
        'hdf5' : 'output.h5',
        'gals' : '/Users/lacerda/CALIFA/listv20_q050.d15a.txt',
        'pycasso_cube_dir' : '/Users/lacerda/CALIFA/gal_fits/v20_q050.d15a',
        'pycasso_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.fits',
        'eml_cube_dir' : '/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a',
        'eml_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.fits', 
        'gasprop_cube_dir' : '/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a/prop',
        'gasprop_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.GasProp.fits',
        'group' : 'radial', 
        'SFRgroup' : 'mask05050525',
    }
    
    parser = ap.ArgumentParser(description = '%s' % args_str)
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default_args['debug'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['hdf5'])
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
    parser.add_argument('--group',
                        metavar = 'GROUPNAME',
                        type = str,
                        default = default_args['group'])
    parser.add_argument('--SFRgroup',
                        metavar = 'GROUPNAME',
                        type = str,
                        default = default_args['SFRgroup'])
    
    args = parser.parse_args()

    args = parser.parse_args()
    args.EL = (args.eml_cube_dir is not None)
    args.GP = (args.gasprop_cube_dir is not None)
    
    if args.pycasso_cube_dir[-1] != '/':
        args.pycasso_cube_dir += '/'
    if args.EL and (args.eml_cube_dir[-1] != '/'):
        args.eml_cube_dir += '/'
    if args.GP and (args.gasprop_cube_dir[-1] != '/'):
        args.gasprop_cube_dir += '/'
    if args.SFRgroup[0] != '/':
        args.SFRgroup = '/%s' % args.SFRgroup

    return args

def load_gal_cubes(args, califaID):
    pycasso_cube_file = args.pycasso_cube_dir + califaID + args.pycasso_cube_suffix
    K = fitsQ3DataCube(pycasso_cube_file)
    if args.EL:
        eml_cube_file = args.eml_cube_dir + califaID + args.eml_cube_suffix
        K.loadEmLinesDataCube(eml_cube_file)
    if args.GP:
        gasprop_cube_file = args.gasprop_cube_dir + califaID + args.gasprop_cube_suffix
        K.GP = GasProp(gasprop_cube_file)
    return K
    
if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse arguments 
    args = parser_args(sys.argv[0])
    debug_var(True, args = args.__dict__)
    
    h5file = tbl.open_file(args.hdf5, mode = 'r+')
    tbl_gals = h5file.root.pycasso.main
    tbl_zones = h5file.root.pycasso.zones
    tbl_integrated = h5file.root.pycasso.integrated
    
    node_SFR = h5file.get_node(args.SFRgroup)
    tbl_SFR_zones_neb = node_SFR.zones_neb
    tbl_SFR_zones_SF = node_SFR.zones_SF
    tbl_SFR_integrated_neb = node_SFR.integrated_neb
    tbl_SFR_integrated_SF = node_SFR.integrated_SF

    tSF__T = node_SFR.tSF.cols.age[:]
    
    for g in tbl_gals:
        iGal = g['id'] 
        # SYN Radial Profiles:
        for iT, tSF in enumerate(tSF__T):
            g_props_SF =  tbl_SFR_zones_SF.read_where('(id_gal == gid) & (id_tSF == iT)', {'gid' : iGal, 'iT' : iT})
            g_props_neb =  tbl_SFR_zones_neb.read_where('(id_gal == gid) & (id_tSF == iT)', {'gid' : iGal, 'iT' : iT})
            g_integrated_props_SF = tbl_SFR_integrated_SF.read_where('(id_gal == gid) & (id_tSF == iT)', {'gid' : iGal, 'iT' : iT})
            g_integrated_props_neb = tbl_SFR_integrated_neb.read_where('(id_gal == gid) & (id_tSF == iT)', {'gid' : iGal, 'iT' : iT})
            
            
                        
            ##########################################
            if (mask_radial.astype(int).sum() < (K.N_zone - minzones)):
                C.debug_var(args.debug, 
                            gal = califaID, 
                            iT = iT, 
                            tSF = '%.3f Myr' % (tSF / 1e6), 
                            radial_profiles = 'computing...',
                )
                
                x_Y__z = tbl_SFR_zones_SF.cols.xY[:]
                Mcor__z = np.ma.masked_array(ALL._Mcor__Tg[iT][iG], mask = mask_radial)
                McorSD__z = np.ma.masked_array(ALL._McorSD__Tg[iT][iG], mask = mask_radial)
                SFR__z = np.ma.masked_array(ALL._SFR__Tg[iT][iG], mask = mask_radial)
                SFRSD__z = np.ma.masked_array(ALL._SFRSD__Tg[iT][iG], mask = mask_radial)
                tau_V__z = np.ma.masked_array(ALL._tau_V__Tg[iT][iG], mask = mask_radial)
                at_flux__z = np.ma.masked_array(ALL._at_flux__Tg[iT][iG], mask = mask_radial)
                at_mass__z = np.ma.masked_array(ALL._at_mass__Tg[iT][iG], mask = mask_radial)
                
                #Mcor__yx = K.zoneToYX(Mcor__z, extensive = False, surface_density = False)
                #v__r, npts = K.radialProfile(Mcor__yx, Rbin__r, return_npts = True, rad_scale = K.HLR_pix, mode = 'sum', mask = mask__yx)
                #McorSD__r = v__r / (npts * K.parsecPerPixel**2.)
                #SFR__yx = K.zoneToYX(SFR__z, extensive = False, surface_density = False)
                #v__r, npts = K.radialProfile(SFR__yx, Rbin__r, return_npts = True, rad_scale = K.HLR_pix, mode = 'sum', mask = mask__yx)
                #aSFRSD__r = v__r / (npts * K.parsecPerPixel**2.)
                x_Y__r = K.zoneToRad(x_Y__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                McorSD__r = K.zoneToRad(McorSD__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                aSFRSD__r = K.zoneToRad(SFRSD__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                tau_V__r = K.zoneToRad(tau_V__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                at_flux__r = K.zoneToRad(at_flux__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                at_mass__r = K.zoneToRad(at_mass__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                at_flux_dezon__r = K.zoneToRad(at_flux__z, Rbin__r, rad_scale = K.HLR_pix)
                at_mass_dezon__r = K.zoneToRad(at_mass__z, Rbin__r, rad_scale = K.HLR_pix)
                
                x_Y_oneHLR = K.zoneToRad(x_Y__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                McorSD_oneHLR = K.zoneToRad(McorSD__z, Rbin_oneHLR, rad_scale = K.HLR_pix)
                aSFRSD_oneHLR = K.zoneToRad(SFRSD__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                tau_V_oneHLR = K.zoneToRad(tau_V__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                at_flux_oneHLR = K.zoneToRad(at_flux__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                at_mass_oneHLR = K.zoneToRad(at_mass__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                at_flux_dezon_oneHLR = K.zoneToRad(at_flux__z, Rbin_oneHLR, rad_scale = K.HLR_pix)
                at_mass_dezon_oneHLR = K.zoneToRad(at_mass__z, Rbin_oneHLR, rad_scale = K.HLR_pix)
    
                Lobn__yx = K.zoneToYX(np.ma.masked_array(K.Lobn__z, mask = mask_radial), extensive = False)
                at_flux__yx = K.zoneToYX(at_flux__z, extensive = False)
                w__r = K.radialProfile(Lobn__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                v_w__r = K.radialProfile(at_flux__yx * Lobn__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                at_flux_wei__r = v_w__r / w__r
                
                Mcor__yx = K.zoneToYX(Mcor__z, extensive = False)
                at_mass__yx = K.zoneToYX(at_mass__z, extensive = False)
                w__r = K.radialProfile(Mcor__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                v_w__r = K.radialProfile(at_mass__yx * Mcor__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                at_mass_wei__r = v_w__r / w__r
        
                ALL.x_Y__Trg[iT, :, iGal] = x_Y__r
                ALL.McorSD__Trg[iT, :, iGal] = McorSD__r
                ALL.aSFRSD__Trg[iT, :, iGal] = aSFRSD__r
                ALL.tau_V__Trg[iT, :, iGal] = tau_V__r
                ALL.at_flux__Trg[iT, :, iGal] = at_flux__r
                ALL.at_mass__Trg[iT, :, iGal] = at_mass__r
                ALL.at_flux_dezon__Trg[iT, :, iGal] = at_flux_dezon__r
                ALL.at_mass_dezon__Trg[iT, :, iGal] = at_mass_dezon__r
                ALL.at_flux_wei__Trg[iT, :, iGal] = at_flux_wei__r
                ALL.at_mass_wei__Trg[iT, :, iGal] = at_mass_wei__r
                ALL.x_Y_oneHLR__Tg[iT, iGal] = x_Y_oneHLR
                ALL.McorSD_oneHLR__Tg[iT, iGal] = McorSD_oneHLR
                ALL.aSFRSD_oneHLR__Tg[iT, iGal] = aSFRSD_oneHLR
                ALL.tau_V_oneHLR__Tg[iT, iGal] = tau_V_oneHLR
                ALL.at_flux_oneHLR__Tg[iT, iGal] = at_flux_oneHLR
                ALL.at_mass_oneHLR__Tg[iT, iGal] = at_mass_oneHLR
                ALL.at_flux_dezon_oneHLR__Tg[iT, iGal] = at_flux_dezon_oneHLR
                ALL.at_mass_dezon_oneHLR__Tg[iT, iGal] = at_mass_dezon_oneHLR
                
                #NEB
                SFRSD_Ha__z = np.ma.masked_array(ALL._SFRSD_Ha__g[iG], mask = mask_radial)
                tau_V_neb__z = np.ma.masked_array(ALL._tau_V_neb__g[iG], mask = mask_radial)
                logO3N2_M13__z = np.ma.masked_array(ALL._logO3N2_M13__g[iG], mask = mask_radial)
                EW_Ha__z = np.ma.masked_array(ALL._EW_Ha__g[iG], mask = mask_radial)
                EW_Hb__z = np.ma.masked_array(ALL._EW_Hb__g[iG], mask = mask_radial)
                F_obs_Hb__z = np.ma.masked_array(ALL._F_obs_Hb__g[iG], mask = mask_radial)
                F_obs_O3__z = np.ma.masked_array(ALL._F_obs_O3__g[iG], mask = mask_radial)
                F_obs_Ha__z = np.ma.masked_array(ALL._F_obs_Ha__g[iG], mask = mask_radial)
                F_obs_N2__z = np.ma.masked_array(ALL._F_obs_N2__g[iG], mask = mask_radial)
                F_int_Hb__z = np.ma.masked_array(ALL._F_int_Hb__g[iG], mask = mask_radial)
                F_int_O3__z = np.ma.masked_array(ALL._F_int_O3__g[iG], mask = mask_radial)
                F_int_Ha__z = np.ma.masked_array(ALL._F_int_Ha__g[iG], mask = mask_radial)
                F_int_N2__z = np.ma.masked_array(ALL._F_int_N2__g[iG], mask = mask_radial)
                baseline_Ha__z = np.ma.masked_array(ALL._baseline_Ha__g[iG], mask = mask_radial)
                baseline_Hb__z = np.ma.masked_array(ALL._baseline_Hb__g[iG], mask = mask_radial)
                
                tau_V_neb__r = K.zoneToRad(tau_V_neb__z, Rbin__r,  rad_scale = K.HLR_pix, extensive = False)
                tau_V_neb_oneHLR = K.zoneToRad(tau_V_neb__z, Rbin_oneHLR,  rad_scale = K.HLR_pix, extensive = False)
                logO3N2_M13__r = K.zoneToRad(logO3N2_M13__z, Rbin__r,  rad_scale = K.HLR_pix, extensive = False)
                logO3N2_M13_oneHLR = K.zoneToRad(logO3N2_M13__z, Rbin_oneHLR,  rad_scale = K.HLR_pix, extensive = False)
                aSFRSD_Ha__r = K.zoneToRad(SFRSD_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                aSFRSD_Ha_oneHLR = K.zoneToRad(SFRSD_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
    
                EW_Ha__r = K.zoneToRad(EW_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                EW_Hb__r = K.zoneToRad(EW_Hb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                EW_Ha_oneHLR = K.zoneToRad(EW_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                EW_Hb_oneHLR = K.zoneToRad(EW_Hb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
    
                F_Ha__yx = K.zoneToYX(F_obs_Ha__z, extensive = False)
                baseline_Ha__yx = K.zoneToYX(baseline_Ha__z, extensive = False)
                v__r = K.radialProfile(F_Ha__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                w__r = K.radialProfile(baseline_Ha__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                EW_Ha_wei__r = v__r / w__r
    
                F_Hb__yx = K.zoneToYX(F_obs_Hb__z, extensive = False)
                baseline_Hb__yx = K.zoneToYX(baseline_Hb__z, extensive = False)
                v__r = K.radialProfile(F_Hb__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                w__r = K.radialProfile(baseline_Hb__yx, bin_r = Rbin__r, mode = 'sum', rad_scale = K.HLR_pix)
                EW_Hb_wei__r = v__r / w__r
    
                F_obs_Ha__r = K.zoneToRad(F_obs_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_Hb__r = K.zoneToRad(F_obs_Hb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_O3__r = K.zoneToRad(F_obs_O3__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_N2__r = K.zoneToRad(F_obs_N2__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_Ha__r = K.zoneToRad(F_int_Ha__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_Hb__r = K.zoneToRad(F_int_Hb__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_O3__r = K.zoneToRad(F_int_O3__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_int_N2__r = K.zoneToRad(F_int_N2__z, Rbin__r, rad_scale = K.HLR_pix, extensive = False)
                F_obs_Ha_oneHLR = K.zoneToRad(F_obs_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_obs_Hb_oneHLR = K.zoneToRad(F_obs_Hb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_obs_O3_oneHLR = K.zoneToRad(F_obs_O3__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_obs_N2_oneHLR = K.zoneToRad(F_obs_N2__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_Ha_oneHLR = K.zoneToRad(F_int_Ha__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_Hb_oneHLR = K.zoneToRad(F_int_Hb__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_O3_oneHLR = K.zoneToRad(F_int_O3__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                F_int_N2_oneHLR = K.zoneToRad(F_int_N2__z, Rbin_oneHLR, rad_scale = K.HLR_pix, extensive = False)
                
                ALL.tau_V_neb__Trg[iT, :, iGal] = tau_V_neb__r
                ALL.logO3N2_M13__Trg[iT, :, iGal] = logO3N2_M13__r
                ALL.EW_Ha__Trg[iT, :, iGal] = EW_Ha__r
                ALL.EW_Hb__Trg[iT, :, iGal] = EW_Hb__r
                ALL.EW_Ha_wei__Trg[iT, :, iGal] = EW_Ha_wei__r
                ALL.EW_Hb_wei__Trg[iT, :, iGal] = EW_Hb_wei__r
                ALL.F_obs_Ha__Trg[iT, :, iGal] = F_obs_Ha__r
                ALL.F_obs_Hb__Trg[iT, :, iGal] = F_obs_Hb__r
                ALL.F_obs_O3__Trg[iT, :, iGal] = F_obs_O3__r
                ALL.F_obs_N2__Trg[iT, :, iGal] = F_obs_N2__r
                ALL.F_int_Ha__Trg[iT, :, iGal] = F_int_Ha__r
                ALL.F_int_Hb__Trg[iT, :, iGal] = F_int_Hb__r
                ALL.F_int_O3__Trg[iT, :, iGal] = F_int_O3__r
                ALL.F_int_N2__Trg[iT, :, iGal] = F_int_N2__r
                ALL.aSFRSD_Ha__Trg[iT, :, iGal] = aSFRSD_Ha__r
                ALL.tau_V_neb_oneHLR__Tg[iT, iGal] = tau_V_neb_oneHLR
                ALL.logO3N2_M13_oneHLR__Tg[iT, iGal] = logO3N2_M13_oneHLR
                ALL.EW_Ha_oneHLR__Tg[iT, iGal] = EW_Ha_oneHLR
                ALL.EW_Hb_oneHLR__Tg[iT, iGal] = EW_Hb_oneHLR
                ALL.F_obs_Ha_oneHLR__Tg[iT, iGal] = F_obs_Ha_oneHLR
                ALL.F_obs_Hb_oneHLR__Tg[iT, iGal] = F_obs_Hb_oneHLR
                ALL.F_obs_O3_oneHLR__Tg[iT, iGal] = F_obs_O3_oneHLR
                ALL.F_obs_N2_oneHLR__Tg[iT, iGal] = F_obs_N2_oneHLR
                ALL.F_int_Ha_oneHLR__Tg[iT, iGal] = F_int_Ha_oneHLR
                ALL.F_int_Hb_oneHLR__Tg[iT, iGal] = F_int_Hb_oneHLR
                ALL.F_int_O3_oneHLR__Tg[iT, iGal] = F_int_O3_oneHLR
                ALL.F_int_N2_oneHLR__Tg[iT, iGal] = F_int_N2_oneHLR
                ALL.aSFRSD_Ha_oneHLR__Tg[iT, iGal] = aSFRSD_Ha_oneHLR
