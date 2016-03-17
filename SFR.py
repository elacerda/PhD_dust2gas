#!/usr/bin/python
#
# Lacerda@Saco - 16/Mar/2016
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

def calc_agebins(ages, age = None):
    # Define ranges for age-bins
    # ToDo: This age-bin-edges thing could be made more elegant & general.
    aCen__t = ages
    aLow__t = np.empty_like(ages)
    aUpp__t = np.empty_like(ages)
    aLow__t[0] = 0.
    aLow__t[1:] = (aCen__t[1:] + aCen__t[:-1]) / 2
    aUpp__t[:-1] = aLow__t[1:]
    aUpp__t[-1] = aCen__t[-1]
    # Find index of age-bin corresponding to the last bin fully within < tSF
    age_index = -1
    if age is not None:
        age_index = np.where(aLow__t < age)[0][-1]
    return aCen__t, aLow__t, aUpp__t, age_index 

def calc_xY(K, tY):
    _, aLow__t, aUpp__t, indY = calc_agebins(K.ageBase, tY)
    # Compute xY__z
    x__tZz = K.popx / K.popx.sum(axis = 1).sum(axis = 0)
    integrated_x__tZ = K.integrated_popx / K.integrated_popx.sum()
    aux1__z = x__tZz[:indY, :, :].sum(axis = 1).sum(axis = 0)
    aux2__z = x__tZz[indY, :, :].sum(axis = 0) * (tY - aLow__t[indY]) / (aUpp__t[indY] - aLow__t[indY])
    integrated_aux1 = integrated_x__tZ[:indY, :].sum()
    integrated_aux2 = integrated_x__tZ[indY, :].sum(axis = 0) * (tY - aLow__t[indY]) / (aUpp__t[indY] - aLow__t[indY])
    return (aux1__z + aux2__z), (integrated_aux1 + integrated_aux2)

def calc_SFR(K, tSF):
    '''
    Add up (in Mini and x) populations younger than tSF to compute SFR's and xY.
    First for zones (__z), and then images (__yx).

    tSF is a arbitrary/continuous number; it'll cover full age-bins plus a last one which will be
    only partially covered. The mass (light) within this last age-bin is scaled by the ratio
    (tSF - bin-lower-age) / bin-size. (P ex, if tSF is such that only 34% of the last-bin width is
    covered, then only 34% of its mass is added to the total.)

    Since our bases span so many ages, this "exact" calculation is just a refinement over the simpler
    method of just adding upp all age-bins satisfying K.agebase < tSF.

    OBS: Note that we are NOT dezonifying SFR__z. (Among other reasons, it'll be compared to the un-dezonifiable tauV!)

    Cid@IAA - 27/Jan/2015
    '''
    _, aLow__t, aUpp__t, indSF = calc_agebins(K.ageBase, tSF)
    
    # Compute SFR__z
    aux1__z = K.Mini__tZz[:indSF, :, :].sum(axis = 1).sum(axis = 0)
    aux2__z = K.Mini__tZz[indSF, :, :].sum(axis = 0) * (tSF - aLow__t[indSF]) / (aUpp__t[indSF] - aLow__t[indSF])
    SFR__z = (aux1__z + aux2__z) / tSF
    SFRSD__z = SFR__z / K.zoneArea_pc2
    
    #aux1__z = K.MiniSD__tZz[:indSF, :, :].sum(axis = 1).sum(axis = 0)
    #aux2__z = K.MiniSD__tZz[indSF, :, :].sum(axis = 0) * (tSF - aLow__t[indSF]) / (aUpp__t[indSF] - aLow__t[indSF])

    return SFR__z, SFRSD__z

def parser_args(args_str):
    default_args = {
        'debug' : False,
        'hdf5' : 'output.h5',
        'minpopx' : np.finfo(np.float_).min,
        'mintauv' : np.finfo(np.float_).min,
        'mintauvneb' : np.finfo(np.float_).min,
        'maxtauvneberr' : np.finfo(np.float_).max,
        'gals' : '/Users/lacerda/CALIFA/listv20_q050.d15a.txt',
        'pycasso_cube_dir' : '/Users/lacerda/CALIFA/gal_fits/v20_q050.d15a',
        'pycasso_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.fits',
        'eml_cube_dir' : '/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a',
        'eml_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.fits', 
        'gasprop_cube_dir' : '/Users/lacerda/CALIFA/rgb-gas/v20_q050.d15a/prop',
        'gasprop_cube_suffix' : '_synthesis_eBR_v20_q050.d15a512.ps03.k1.mE.CCM.Bgsd6e.EML.MC100.GasProp.fits',
        'group' : 'young' 
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
    parser.add_argument('--minpopx',
                        help = 'Negative to disable mask in popx',
                        metavar = 'FRAC',
                        type = float,
                        default = default_args['minpopx'])
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
    parser.add_argument('--group',
                        metavar = 'GROUPNAME',
                        type = str,
                        default = default_args['group'])
    
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

    tSF__T = np.array([1, 3.2, 10, 50, 100]) * 1e7
    N_T = len(tSF__T)
    
    h5file = tbl.open_file(args.hdf5, mode = 'r+')
    tbl_gals = h5file.root.pycasso.main
    tbl_zone = h5file.root.pycasso.zones
    tbl_integrated = h5file.root.pycasso.integrated
    
    group_description = 'minpopx:%.2f/mintauV:%.2f/mintauVneb:%.2f/maxtauVneberr:%.2f - zones calculation' % (args.minpopx, args.mintauv, args.mintauvneb, args.maxtauvneberr)
    group = h5file.create_group('/', args.group, group_description, 
                                filters = tbl.Filters(1))
    tbl_tSF = h5file.create_table(group, 'tSF', tSF, 'tSF data')
    tbl_zone_SF = h5file.create_table(group, 'zones_SF', zone_SF, 'tSF data')
    tbl_integrated_SF = h5file.create_table(group, 'integrated_SF', zone_SF, 'tSF data')
    tbl_zone_neb = h5file.create_table(group, 'zones_neb', zone_neb, 'Zone SF data')
    tbl_integrated_neb = h5file.create_table(group, 'integrated_neb', zone_neb, 'Zone SF data')
    
    tbl_tSF.append([x for x in enumerate(tSF__T)])
    tbl_tSF.cols.id.create_index()
    tbl_tSF.flush()
    
    for g in tbl_gals:
        t_init_gal = time.clock()
        
        K = load_gal_cubes(args, g['califaID'])
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)
        
        for iT, tSF in enumerate(tSF__T):
            x_Y__z, integrated_x_Y = calc_xY(K, tSF)
            SFR__z, SFRSD__z = calc_SFR(K, tSF)
            
            tmp = np.zeros((g['N_zone']), dtype = np.int)
            id_tSF = tmp + iT
            id_gal = tmp + g['id']
            zone_SF_data = zip(
                id_gal, 
                np.arange(g['N_zone']),
                id_tSF, 
                x_Y__z, 
                SFR__z, 
                SFRSD__z,
            )
            tbl_zone_SF.append(zone_SF_data)
            tbl_zone_SF.flush()
            del tmp

            integrated_SFR = SFR__z.sum()
            integrated_SFRSD = integrated_SFR / K.zoneArea_pc2.sum()
            integrated_SF_data = [(
                g['id'], 
                -1,
                iT, 
                integrated_x_Y,
                integrated_SFR,
                integrated_SFRSD, 
            )]
            tbl_integrated_SF.append(integrated_SF_data)
            tbl_integrated_SF.flush()
            
        g_props = tbl_zone.read_where('id_gal == gid', {'gid' : g['id']})
        
        # zone index sorted by id
        _izS = np.argsort(g_props['id'])
        
        tau_V_neb = g_props['tau_V_neb'][_izS]
        etau_V_neb = g_props['etau_V_neb'][_izS]
        snr_Hb = g_props['F_obs_Hb'][_izS] / g_props['eF_obs_Hb'][_izS]
        snr_O3 = g_props['F_obs_O3'][_izS] / g_props['eF_obs_O3'][_izS]
        snr_Ha = g_props['F_obs_Ha'][_izS] / g_props['eF_obs_Ha'][_izS]
        snr_N2 = g_props['F_obs_N2'][_izS] / g_props['eF_obs_N2'][_izS]
        flag_BPT = np.bitwise_or(np.less(snr_Hb, 3), np.less(snr_O3, 3))
        flag_BPT = np.bitwise_or(flag_BPT, np.less(snr_Ha, 3))
        flag_BPT = np.bitwise_or(flag_BPT, np.less(snr_N2, 3))
        flag_WHAN = np.bitwise_or(np.less(snr_Ha, 3), np.less(snr_N2, 3))
        mask_HaHb = np.bitwise_or(np.less(snr_Hb, 3), np.less(snr_Ha, 3))
        mask_tau_V_neb = np.bitwise_or(np.less(tau_V_neb, args.mintauvneb),
                                       np.greater_equal(etau_V_neb, args.maxtauvneberr))
        mask_neb = np.bitwise_or(mask_tau_V_neb, mask_HaHb)
        
        q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])
        expqtau = [ np.ma.exp(qcard * tau_V_neb) for qcard in q ]
        L_obs_Ha__z = K.EL._F_to_L(g_props['F_obs_Ha'][_izS], g['distance_Mpc']) / L_sun
        L_int_Ha__z = np.where(~(mask_neb), L_obs_Ha__z * expqtau[2], L_obs_Ha__z)
        SFR_Ha__z = 3.13 * L_int_Ha__z 
        SFRSD_Ha__z = SFR_Ha__z / K.zoneArea_pc2
        
        tmp = np.ones((g['N_zone']), dtype = np.int)
        zone_neb_data = zip(
            tmp * g['id'], np.arange(g['N_zone']),
            L_obs_Ha__z, L_int_Ha__z,
            SFR_Ha__z, SFRSD_Ha__z,
            flag_BPT, flag_WHAN,
        )
        tbl_zone_neb.append(zone_neb_data)
        tbl_zone_neb.flush()

        g_int_props = tbl_integrated.read_where('id_gal == gid', {'gid' : g['id']})
        integrated_tau_V_neb = g_int_props['tau_V_neb']
        integrated_expqtau = [ np.ma.exp(qcard * integrated_tau_V_neb) for qcard in q ]
        integrated_tau_V_neb = g_int_props['tau_V_neb']
        integrated_etau_V_neb = g_int_props['etau_V_neb']
        integrated_snr_Hb = g_int_props['F_obs_Hb'] / g_int_props['eF_obs_Hb']
        integrated_snr_O3 = g_int_props['F_obs_O3'] / g_int_props['eF_obs_O3']
        integrated_snr_Ha = g_int_props['F_obs_Ha'] / g_int_props['eF_obs_Ha']
        integrated_snr_N2 = g_int_props['F_obs_N2'] / g_int_props['eF_obs_N2']
        integrated_flag_BPT = np.bitwise_or(np.less(integrated_snr_Hb, 3), np.less(integrated_snr_O3, 3))
        integrated_flag_BPT = np.bitwise_or(integrated_flag_BPT, np.less(integrated_snr_Ha, 3))
        integrated_flag_BPT = np.bitwise_or(integrated_flag_BPT, np.less(integrated_snr_N2, 3))
        integrated_flag_WHAN = np.bitwise_or(np.less(integrated_snr_Ha, 3), np.less(integrated_snr_N2, 3))
        integrated_mask_HaHb = np.bitwise_or(np.less(integrated_snr_Hb, 3), np.less(integrated_snr_Ha, 3))
        integrated_mask_tau_V_neb = np.bitwise_or(np.less(integrated_tau_V_neb, args.mintauvneb),
                                                  np.greater_equal(integrated_etau_V_neb, args.maxtauvneberr))
        integrated_mask_neb = np.bitwise_or(integrated_mask_tau_V_neb, integrated_mask_HaHb)
        integrated_L_obs_Ha = K.EL._F_to_L(g_int_props['F_obs_Ha'], g['distance_Mpc']) / L_sun
        if integrated_mask_neb:
            integrated_L_int_Ha = integrated_L_obs_Ha
        else:
            integrated_L_int_Ha = integrated_L_obs_Ha * integrated_expqtau[2] 
        integrated_SFR_Ha = SFR_Ha__z.sum()
        integrated_SFRSD_Ha = integrated_SFR_Ha / K.zoneArea_pc2.sum()
        integrated_zone_neb_data = [(
            g['id'], -1,
            integrated_L_obs_Ha, integrated_L_int_Ha,
            integrated_SFR_Ha, integrated_SFRSD_Ha,
            integrated_flag_BPT, integrated_flag_WHAN
        )]
        tbl_integrated_neb.append(integrated_zone_neb_data)
        tbl_integrated_neb.flush()
        
        K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (g['califaID'], time.clock() - t_init_gal)

    tbl_zone_SF.cols.id_gal.create_index()
    tbl_zone_SF.cols.id_zone.create_index()
    tbl_zone_SF.cols.id_tSF.create_index()
    tbl_zone_neb.cols.id_gal.create_index()
    tbl_zone_neb.cols.id_zone.create_index()
    tbl_zone_neb.cols.flag_BPT.create_index()
    tbl_zone_neb.cols.flag_WHAN.create_index()
    tbl_integrated_SF.cols.id_gal.create_index()
    tbl_integrated_SF.cols.id_tSF.create_index()
    tbl_integrated_neb.cols.id_gal.create_index()
    tbl_integrated_neb.cols.flag_BPT.create_index()
    tbl_integrated_neb.cols.flag_WHAN.create_index()
    
    tbl_zone_SF.flush()
    tbl_zone_neb.flush()
    tbl_integrated_SF.flush()
    tbl_integrated_neb.flush()

    h5file.close()
    print 'total time: %.2f' % (time.clock() - t_init_prog)    