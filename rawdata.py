#!/usr/bin/python
#
# Lacerda@Saco - 15/Mar/2016
#
from CALIFAUtils.scripts import read_gal_cubes
from CALIFAUtils.scripts import get_morfologia
from CALIFAUtils.scripts import sort_gals
from CALIFAUtils.scripts import debug_var
from CALIFAUtils.objects import GasProp
from CALIFAUtils.scripts import my_morf
from tables_description import galaxy
from tables_description import zone
from pycasso import fitsQ3DataCube
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
        'morph_file' : '/Users/lacerda/CALIFA/morph_eye_class.csv',  
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
    parser.add_argument('--morph_file', '-M',
                        metavar = 'FILE',
                        type = str,
                        default = default_args['morph_file'])
    
    args = parser.parse_args()
    args.EL = (args.eml_cube_dir is not None)
    args.GP = (args.gasprop_cube_dir is not None)
    
    if args.pycasso_cube_dir[-1] != '/':
        args.pycasso_cube_dir += '/'
    if args.EL and (args.eml_cube_dir[-1] != '/'):
        args.eml_cube_dir += '/'
    if args.GP and (args.gasprop_cube_dir[-1] != '/'):
        args.gasprop_cube_dir += '/'

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
    
if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse arguments 
    args = parser_args(sys.argv[0])
    debug_var(True, args = args.__dict__)    

    h5file = tbl.open_file(args.hdf5, mode =  'w', title = 'Test file')
    group = h5file.create_group('/', 'pycasso', 'Galaxy PyCASSO data', filters = tbl.Filters(1))
    tbl_main = h5file.create_table(group, 'main', galaxy, 'Main data')
    tbl_zone = h5file.create_table(group, 'zones', zone, 'Zone data')
    tbl_integrated = h5file.create_table(group, 'integrated', zone, 'Integrated data')
        
    gals, _ = sort_gals(args.gals, order = 1)
    N_gals = len(gals)
    max_gals = N_gals
    if args.debug:
        max_gals = 10
        
    id_zone_ini = id_zone_fin = 0
    
    for iGal, gal in enumerate(gals[0:max_gals]):
        t_init_gal = time.clock()
        
        K = load_gal_cubes(args, gal)
        sit, verify = verify_files(K, gal, EL = args.EL, GP = args.GP)
        
        if verify is not True:
            print '<<< ', gal, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue

        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
    
        id_zone_ini = id_zone_fin
        id_zone_fin = id_zone_ini + K.N_zone
        
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
        m_type = my_morf(get_morfologia(K.califaID, morph_file = args.morph_file)[0])
    
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
        )]
        tbl_main.append(main_data)
        tbl_main.flush()
        
        # integrated data stored as a zone with id -1
        r = tbl_integrated.row
        r['id_gal'] = iGal
        r['tau_V']  =  K.integrated_keywords['A_V'] * AVtoTauV 
        r['at_flux'] = at_flux_GAL
        r['tau_V_neb'] = K.EL.integrated_tau_V_neb 
        r['etau_V_neb'] = K.EL.integrated_tau_V_neb_err
        r['F_obs_Hb'] = K.EL.integrated_flux[i_Hb]
        r['F_obs_O3'] = K.EL.integrated_flux[i_O3]
        r['F_obs_Ha'] = K.EL.integrated_flux[i_Ha] 
        r['F_obs_N2'] = K.EL.integrated_flux[i_N2]
        r['eF_obs_Hb'] = K.EL.integrated_eflux[i_Hb] 
        r['eF_obs_O3'] = K.EL.integrated_eflux[i_O3]
        r['eF_obs_Ha'] = K.EL.integrated_eflux[i_Ha] 
        r['eF_obs_N2'] = K.EL.integrated_eflux[i_N2]
        r['baseline_Hb'] = K.EL.integrated_baseline[i_Hb] 
        r['baseline_O3'] = K.EL.integrated_baseline[i_O3]
        r['baseline_Ha'] = K.EL.integrated_baseline[i_Ha] 
        r['baseline_N2'] = K.EL.integrated_baseline[i_N2]
        r['EW_Hb'] = K.EL.integrated_EW[i_Hb] 
        r['EW_O3'] = K.EL.integrated_EW[i_O3]
        r['EW_Ha'] = K.EL.integrated_EW[i_Ha] 
        r['EW_N2'] = K.EL.integrated_EW[i_N2]
        r['sigma_Hb'] = K.EL.integrated_sigma[i_Hb]
        r['sigma_O3'] = K.EL.integrated_sigma[i_O3]
        r['sigma_Ha'] = K.EL.integrated_sigma[i_Ha]
        r['sigma_N2'] = K.EL.integrated_sigma[i_N2]
        r['esigma_Hb'] = K.EL.integrated_esigma[i_Hb]
        r['esigma_O3'] = K.EL.integrated_esigma[i_O3]
        r['esigma_Ha'] = K.EL.integrated_esigma[i_Ha]
        r['esigma_N2'] = K.EL.integrated_esigma[i_N2]
        r['pos_Hb'] = K.EL.integrated_pos[i_Hb]
        r['pos_O3'] = K.EL.integrated_pos[i_O3]
        r['pos_Ha'] = K.EL.integrated_pos[i_Ha]
        r['pos_N2'] = K.EL.integrated_pos[i_N2]
        r['epos_Hb'] = K.EL.integrated_epos[i_Hb]
        r['epos_O3'] = K.EL.integrated_epos[i_O3]
        r['epos_Ha'] = K.EL.integrated_epos[i_Ha]
        r['epos_N2'] = K.EL.integrated_epos[i_N2]
        r.append()
        tbl_integrated.flush()
        del r
    
        # zip to transpose arrays to col_arrays
        zone_data = zip(
            np.zeros((K.N_zone), dtype = int) + iGal,
            np.arange(id_zone_ini, id_zone_fin),
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
            K.EL.sigma[i_Hb],
            K.EL.sigma[i_O3],
            K.EL.sigma[i_Ha],
            K.EL.sigma[i_N2],
            K.EL.esigma[i_Hb],
            K.EL.esigma[i_O3],
            K.EL.esigma[i_Ha],
            K.EL.esigma[i_N2],
            K.EL.pos[i_Hb],
            K.EL.pos[i_O3],
            K.EL.pos[i_Ha],
            K.EL.pos[i_N2],
            K.EL.epos[i_Hb],
            K.EL.epos[i_O3],
            K.EL.epos[i_Ha],
            K.EL.epos[i_N2],
            K.EL.Zneb_M13__z,
        )
        tbl_zone.append(zone_data)
        tbl_zone.flush()

        K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (gal, time.clock() - t_init_gal)
        
    tbl_main.cols.id.create_csindex()
    tbl_main.cols.califaID.create_csindex()
    tbl_main.cols.m_type.create_csindex()
    tbl_zone.cols.id_gal.create_csindex()
    tbl_zone.cols.id.create_csindex()
    tbl_zone.cols.i_zone.create_index()
    tbl_zone.cols.tau_V.create_index()
    tbl_integrated.cols.id_gal.create_csindex()

    tbl_main.flush()
    tbl_zone.flush()
    tbl_integrated.flush()

    h5file.close()
    print 'total time: %.2f' % (time.clock() - t_init_prog)    