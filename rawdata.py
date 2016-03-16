#!/usr/bin/python
#
# Lacerda@Saco - 15/Mar/2016
#
from CALIFAUtils.scripts import read_gal_cubes
from CALIFAUtils.scripts import get_morfologia
from CALIFAUtils.scripts import sort_gals
from CALIFAUtils.scripts import debug_var
from CALIFAUtils.scripts import my_morf
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
    
    args = parser.parse_args()

    return args

class galaxy(tbl.IsDescription):
    id                          = tbl.UInt16Col(pos = 1)
    name                        = tbl.StringCol(20, pos = 2)
    califaID                    = tbl.StringCol(5, pos = 3)
    N_zones                     = tbl.UInt16Col(pos = 4)
    distance_Mpc                = tbl.Float64Col(pos = 5)
    redshift                    = tbl.Float64Col(pos = 6)
    m_type                      = tbl.Int32Col(pos = 7) 
    ba                          = tbl.Float64Col(pos = 8)
    ba_PyCASSO                  = tbl.Float64Col(pos = 9)
    ParsecPerPixel              = tbl.Float64Col(pos = 10)
    Mr                          = tbl.Float64Col(pos = 11)
    ur                          = tbl.Float64Col(pos = 12)
    
    # Syn
    HLR_pix                     = tbl.Float64Col(pos = 13)
    HMR_pix                     = tbl.Float64Col(pos = 14)
    integrated_tau_V            = tbl.Float64Col(pos = 15)
    at_flux_GAL                 = tbl.Float64Col(pos = 16)
    
    # Neb
    integrated_tau_V_neb        = tbl.Float64Col(pos = 17)
    integrated_etau_V_neb       = tbl.Float64Col(pos = 18)
    integrated_F_obs_Hb         = tbl.Float64Col(pos = 19)
    integrated_F_obs_O3         = tbl.Float64Col(pos = 20)
    integrated_F_obs_Ha         = tbl.Float64Col(pos = 21)
    integrated_F_obs_N2         = tbl.Float64Col(pos = 22)
    integrated_eF_obs_Hb        = tbl.Float64Col(pos = 23)
    integrated_eF_obs_O3        = tbl.Float64Col(pos = 24)
    integrated_eF_obs_Ha        = tbl.Float64Col(pos = 25)
    integrated_eF_obs_N2        = tbl.Float64Col(pos = 26)
    integrated_baseline_Hb      = tbl.Float64Col(pos = 27)
    integrated_baseline_O3      = tbl.Float64Col(pos = 28)
    integrated_baseline_Ha      = tbl.Float64Col(pos = 29)
    integrated_baseline_N2      = tbl.Float64Col(pos = 30)
    integrated_EW_Hb            = tbl.Float64Col(pos = 31)
    integrated_EW_O3            = tbl.Float64Col(pos = 32)
    integrated_EW_Ha            = tbl.Float64Col(pos = 33)
    integrated_EW_N2            = tbl.Float64Col(pos = 34)
    
class zone(tbl.IsDescription):
    id_gal                  = tbl.UInt16Col(pos = 1)
    id_zone                 = tbl.UInt16Col(pos = 2)
    distance_HLR            = tbl.Float64Col(pos = 3)
    area_pc2                = tbl.Float64Col(pos = 4)
    
    # Syn
    Mcor                    = tbl.Float64Col(pos = 5)
    tau_V                   = tbl.Float64Col(pos = 6)
    McorSD                  = tbl.Float64Col(pos = 7)
    at_flux                 = tbl.Float64Col(pos = 8)
    at_mass                 = tbl.Float64Col(pos = 9)
    
    # Neb
    tau_V_neb               = tbl.Float64Col(pos = 10)
    etau_V_neb              = tbl.Float64Col(pos = 11)
    F_obs__Hb               = tbl.Float64Col(pos = 12)
    F_obs__O3               = tbl.Float64Col(pos = 13)
    F_obs__Ha               = tbl.Float64Col(pos = 14)
    F_obs__N2               = tbl.Float64Col(pos = 15)
    eF_obs__Hb              = tbl.Float64Col(pos = 16)
    eF_obs__O3              = tbl.Float64Col(pos = 17)
    eF_obs__Ha              = tbl.Float64Col(pos = 18)
    eF_obs__N2              = tbl.Float64Col(pos = 19)
    baseline_Hb             = tbl.Float64Col(pos = 20)
    baseline_O3             = tbl.Float64Col(pos = 21)
    baseline_Ha             = tbl.Float64Col(pos = 22)
    baseline_N2             = tbl.Float64Col(pos = 23)
    EW_Hb                   = tbl.Float64Col(pos = 24)
    EW_O3                   = tbl.Float64Col(pos = 25)
    EW_Ha                   = tbl.Float64Col(pos = 26)
    EW_N2                   = tbl.Float64Col(pos = 27)
    logOH                   = tbl.Float64Col(pos = 28)
    
if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse arguments 
    args = parser_args(sys.argv[0])
    debug_var(True, args = args.__dict__)    

    h5file = tbl.open_file('teste.h5', mode =  'w', title = 'Test file')
    group = h5file.create_group('/', 'pycasso', 'Galaxy PyCASSO data', filters = tbl.Filters(1))
    table_main = h5file.create_table(group, 'main', galaxy, 'Main data')
    table_zone = h5file.create_table(group, 'zones', zone, 'Zone data')
    
    gals, _ = sort_gals(args.gals, order = 1)
    
    r = table_main.row
    for iGal, gal in enumerate(gals):
        t_init_gal = time.clock()

        K = read_gal_cubes(gal, 
                           debug = args.debug, 
                           pycasso_cube_dir = args.pycasso_cube_dir, 
                           pycasso_cube_suffix = args.pycasso_cube_suffix,
                           eml_cube_dir = args.eml_cube_dir, 
                           eml_cube_suffix = args.eml_cube_suffix,
                           gasprop_cube_dir = args.gasprop_cube_dir, 
                           gasprop_cube_suffix = args.gasprop_cube_suffix)
        
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
    
        Hb_central_wl = '4861'
        O3_central_wl = '5007'
        Ha_central_wl = '6563'
        N2_central_wl = '6583'
        i_Hb = K.EL.lines.index(Hb_central_wl)
        i_O3 = K.EL.lines.index(O3_central_wl)
        i_Ha = K.EL.lines.index(Ha_central_wl)
        i_N2 = K.EL.lines.index(N2_central_wl)
        
        numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
        at_flux_GAL = numerator__z.sum() / denominator__z.sum()
        AVtoTauV = 1. / (np.log10(np.exp(1)) / 0.4)
        m_type = my_morf(get_morfologia(K.califaID)[0])
    
        # Must follow pos order in galaxy class
        main_data = [(
            iGal, 
            K.galaxyName, 
            K.califaID, 
            K.N_zone, 
            K.distance_Mpc, 
            K.redshift, 
            m_type, 
            np.float(K.masterListData['ba']), 
            ba, 
            K.parsecPerPixel, 
            np.float(K.masterListData['Mr']), 
            np.float(K.masterListData['u-r']), 
            K.HLR_pix, 
            K.getHalfRadius(K.McorSD__yx), 
            K.integrated_keywords['A_V'] * AVtoTauV, 
            at_flux_GAL,
            K.EL.integrated_tau_V_neb, 
            K.EL.integrated_tau_V_neb_err,
            K.EL.integrated_flux[i_Hb], 
            K.EL.integrated_flux[i_O3],
            K.EL.integrated_flux[i_Ha], 
            K.EL.integrated_flux[i_N2],
            K.EL.integrated_eflux[i_Hb], 
            K.EL.integrated_eflux[i_O3],
            K.EL.integrated_eflux[i_Ha], 
            K.EL.integrated_eflux[i_N2],
            K.EL.integrated_baseline[i_Hb], 
            K.EL.integrated_baseline[i_O3],
            K.EL.integrated_baseline[i_Ha], 
            K.EL.integrated_baseline[i_N2],
            K.EL.integrated_EW[i_Hb], 
            K.EL.integrated_EW[i_O3],
            K.EL.integrated_EW[i_Ha], 
            K.EL.integrated_EW[i_N2],
        )]
        table_main.append(main_data)
        table_main.flush()
    
        # zip to transpose arrays to col_arrays
        zone_data = zip(
            np.zeros((K.N_zone), dtype = int) + iGal,
            np.arange(K.N_zone),
            K.zoneDistance_HLR,
            K.zoneArea_pc2,
            # Syn
            K.Mcor__z,
            K.tau_V__z,
            K.Mcor__z / K.zoneArea_pc2,
            K.at_flux__z,
            K.at_mass__z,
            # Neb
            K.EL.tau_V_neb__z,
            K.EL.tau_V_neb_err__z,
            K.EL.flux[i_Hb], 
            K.EL.flux[i_O3], 
            K.EL.flux[i_Ha], 
            K.EL.flux[i_N2], 
            K.EL.eflux[i_Hb], 
            K.EL.eflux[i_O3], 
            K.EL.eflux[i_Ha], 
            K.EL.eflux[i_N2], 
            K.EL.baseline[i_Hb], 
            K.EL.baseline[i_O3], 
            K.EL.baseline[i_Ha], 
            K.EL.baseline[i_N2], 
            K.EL.EW[i_Hb], 
            K.EL.EW[i_O3], 
            K.EL.EW[i_Ha], 
            K.EL.EW[i_N2], 
            K.EL.Zneb_M13__z,
        )
        table_zone.append(zone_data)
        table_zone.flush()

        K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (gal, time.clock() - t_init_gal)
    h5file.close()
    print 'total time: %.2f' % (time.clock() - t_init_prog)    