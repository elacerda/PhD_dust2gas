#
# Lacerda@Saco - 16/Mar/2016
#
import os
import sys
import time
import numpy as np
import tables as tbl
from tables_description import tSF as tStarForm
from tables_description import zone_SF
from tables_description import zone_neb
from pytu.functions import debug_var
from pystarlight.util import redenninglaws
from pystarlight.util.constants import L_sun
from rawdata import fix_dir_args
from rawdata import verify_files
from rawdata import load_gal_cubes
from pytu.objects import CustomArgumentParser


def parser_args(default_args_file='default.args'):
    default_args = {
        'debug': False,
        'group': 'young',
        'hdf5': 'output.h5',
        'minpopx': np.finfo(np.float_).min,
        'mintauv': np.finfo(np.float_).min,
        'mintauvneb': np.finfo(np.float_).min,
        'maxtauvneberr': np.finfo(np.float_).max,
    }

    parser = CustomArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--debug', '-D', action='store_true',
                        default=default_args['debug'])
    parser.add_argument('--hdf5', '-H', metavar='FILE', type=str,
                        default=default_args['hdf5'])
    parser.add_argument('--pycasso_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--eml_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--gasprop_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--minpopx', metavar='FRAC', type=float,
                        help='Negative to disable mask in popx',
                        default=default_args['minpopx'])
    parser.add_argument('--mintauv', metavar='FRAC', type=float,
                        default=default_args['mintauv'])
    parser.add_argument('--mintauvneb', metavar='FRAC', type=float,
                        default=default_args['mintauvneb'])
    parser.add_argument('--maxtauvneberr', metavar='FRAC', type=float,
                        default=default_args['maxtauvneberr'])
    parser.add_argument('--group', metavar='GROUPNAME', type=str,
                        default=default_args['group'])

    args_list = sys.argv[1:]
    # if exists file default.args, load default args
    if os.path.isfile(default_args_file):
        args_list.insert(0, '@%s' % default_args_file)
    debug_var(True, args_list=args_list)
    return fix_dir_args(parser.parse_args(args=args_list))


def calc_agebins(ages, age=None):
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
    x__tZz = K.popx / K.popx.sum(axis=1).sum(axis=0)
    integrated_x__tZ = K.integrated_popx / K.integrated_popx.sum()
    aux1__z = x__tZz[:indY, :, :].sum(axis=1).sum(axis=0)
    aux2__z = x__tZz[indY, :, :].sum(axis=0) * (tY - aLow__t[indY]) / (aUpp__t[indY] - aLow__t[indY])
    integrated_aux1 = integrated_x__tZ[:indY, :].sum()
    integrated_aux2 = integrated_x__tZ[indY, :].sum(axis=0) * (tY - aLow__t[indY]) / (aUpp__t[indY] - aLow__t[indY])
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
    aux1__z = K.Mini__tZz[:indSF, :, :].sum(axis=1).sum(axis=0)
    aux2__z = K.Mini__tZz[indSF, :, :].sum(axis=0) * (tSF - aLow__t[indSF]) / (aUpp__t[indSF] - aLow__t[indSF])
    SFR__z = (aux1__z + aux2__z) / tSF
    SFRSD__z = SFR__z / K.zoneArea_pc2

    # aux1__z = K.MiniSD__tZz[:indSF, :, :].sum(axis=1).sum(axis=0)
    # aux2__z = K.MiniSD__tZz[indSF, :, :].sum(axis=0) * (tSF - aLow__t[indSF]) / (aUpp__t[indSF] - aLow__t[indSF])

    return SFR__z, SFRSD__z


if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse arguments
    args = parser_args(default_args_file='SFR_default.args')
    debug_var(True, args=args.__dict__)

    tSF__T = np.array([1, 3.2, 10, 50, 100]) * 1e7
    N_T = len(tSF__T)

    q = redenninglaws.Cardelli_RedLaw([4861, 5007, 6563, 6583])

    h5file = tbl.open_file(args.hdf5, mode='r+')
    tbl_gals = h5file.root.pycasso.main
    tbl_zones = h5file.root.pycasso.zones
    tbl_integrated = h5file.root.pycasso.integrated

    group_description = 'minpopx:%.2f' % args.minpopx
    group_description += '/mintauV:%.2f' % args.mintauv
    group_description += '/mintauVneb:%.2f/' % args.mintauvneb
    group_description += 'maxtauVneberr:%.2f' % args.maxtauvneberr
    group_description += ' - SFR calculation'
    group = h5file.create_group('/', args.group, group_description,
                                filters=tbl.Filters(1))
    tbl_tSF = h5file.create_table(group, 'tSF', tStarForm, 'tSF data')
    tbl_zone_SF = h5file.create_table(group, 'zones_SF', zone_SF, 'tSF data')
    tbl_integrated_SF = h5file.create_table(group, 'integrated_SF', zone_SF, 'tSF data')
    tbl_zone_neb = h5file.create_table(group, 'zones_neb', zone_neb, 'Zone SF data')
    tbl_integrated_neb = h5file.create_table(group, 'integrated_neb', zone_neb, 'Zone SF data')

    tbl_tSF.append([t for t in enumerate(tSF__T)])
    tbl_tSF.cols.id.create_csindex()
    tbl_tSF.flush()

    for g in tbl_gals:
        t_init_gal = time.clock()

        K, _ = load_gal_cubes(args, g['califaID'],
                              pycasso_cube_file=g['pycasso_cube_filename'],
                              eml_cube_file=g['eml_cube_filename'],
                              gasprop_cube_file=g['gasprop_cube_filename'])
        sit, verify = verify_files(K, g['califaID'], EL=args.EL, GP=args.GP)
        # set files situation
        if verify is not True:
            print '<<< ', g['califaID'], sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue

        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)

        g_props__z = tbl_zones.read_where('id_gal == gid', {'gid': g['id']})
        g_int_props = tbl_integrated.read_where('id_gal == gid', {'gid': g['id']})
        id_zones = g_props__z['id']
        _izS = np.argsort(g_props__z['id'])  # zone index sorted by id

        # NEB #
        tau_V_neb__z = g_props__z['tau_V_neb'][_izS]
        etau_V_neb__z = g_props__z['etau_V_neb'][_izS]
        snr_Hb__z = g_props__z['F_obs_Hb'][_izS] / g_props__z['eF_obs_Hb'][_izS]
        snr_O3__z = g_props__z['F_obs_O3'][_izS] / g_props__z['eF_obs_O3'][_izS]
        snr_Ha__z = g_props__z['F_obs_Ha'][_izS] / g_props__z['eF_obs_Ha'][_izS]
        snr_N2__z = g_props__z['F_obs_N2'][_izS] / g_props__z['eF_obs_N2'][_izS]

        flag_BPT__z = np.bitwise_or(np.less(snr_Hb__z, 3), np.less(snr_O3__z, 3))
        flag_BPT__z = np.bitwise_or(flag_BPT__z, np.less(snr_Ha__z, 3))
        flag_BPT__z = np.bitwise_or(flag_BPT__z, np.less(snr_N2__z, 3))
        flag_WHAN__z = np.bitwise_or(np.less(snr_Ha__z, 3), np.less(snr_N2__z, 3))
        flag_tau_V_neb__z = np.less(tau_V_neb__z, args.mintauvneb)
        flag_etau_V_neb__z = np.greater(etau_V_neb__z, args.maxtauvneberr)
        mask_HaHb__z = np.bitwise_or(np.less(snr_Hb__z, 3), np.less(snr_Ha__z, 3))
        mask_tau_V_neb__z = np.bitwise_or(flag_tau_V_neb__z, flag_etau_V_neb__z)
        mask_neb = np.bitwise_or(mask_tau_V_neb__z, mask_HaHb__z)
        expqtau = [np.ma.exp(qcard * tau_V_neb__z) for qcard in q]
        L_obs_Ha__z = K.EL._F_to_L(g_props__z['F_obs_Ha'][_izS], g['distance_Mpc']) / L_sun
        L_int_Ha__z = np.where(~(mask_neb), L_obs_Ha__z * expqtau[2], L_obs_Ha__z)
        SFR_Ha__z = 3.13 * L_int_Ha__z
        SFRSD_Ha__z = SFR_Ha__z / g_props__z['area_pc2']

        tmp = np.ones((g['N_zone']), dtype=np.int)
        zone_neb_data = zip(
            tmp * g['id'], id_zones, np.arange(g['N_zone']),
            L_obs_Ha__z, L_int_Ha__z,
            SFR_Ha__z, SFRSD_Ha__z,
            flag_BPT__z, flag_WHAN__z,
            flag_tau_V_neb__z, flag_etau_V_neb__z
        )
        tbl_zone_neb.append(zone_neb_data)
        tbl_zone_neb.flush()

        integrated_tau_V_neb = g_int_props['tau_V_neb']
        integrated_etau_V_neb = g_int_props['etau_V_neb']
        integrated_expqtau = [ np.ma.exp(qcard * integrated_tau_V_neb) for qcard in q ]
        integrated_snr_Hb = g_int_props['F_obs_Hb'] / g_int_props['eF_obs_Hb']
        integrated_snr_O3 = g_int_props['F_obs_O3'] / g_int_props['eF_obs_O3']
        integrated_snr_Ha = g_int_props['F_obs_Ha'] / g_int_props['eF_obs_Ha']
        integrated_snr_N2 = g_int_props['F_obs_N2'] / g_int_props['eF_obs_N2']
        integrated_flag_BPT = np.bitwise_or(np.less(integrated_snr_Hb, 3), np.less(integrated_snr_O3, 3))
        integrated_flag_BPT = np.bitwise_or(integrated_flag_BPT, np.less(integrated_snr_Ha, 3))
        integrated_flag_BPT = np.bitwise_or(integrated_flag_BPT, np.less(integrated_snr_N2, 3))
        integrated_flag_WHAN = np.bitwise_or(np.less(integrated_snr_Ha, 3), np.less(integrated_snr_N2, 3))
        integrated_flag_tau_V_neb = np.less(integrated_tau_V_neb, args.mintauvneb)
        integrated_flag_etau_V_neb = np.greater(integrated_etau_V_neb, args.maxtauvneberr)
        integrated_mask_HaHb = np.bitwise_or(np.less(integrated_snr_Hb, 3), np.less(integrated_snr_Ha, 3))
        integrated_mask_tau_V_neb = np.bitwise_or(integrated_flag_tau_V_neb, integrated_flag_etau_V_neb)
        integrated_mask_neb = np.bitwise_or(integrated_mask_tau_V_neb, integrated_mask_HaHb)
        integrated_L_obs_Ha = K.EL._F_to_L(g_int_props['F_obs_Ha'], g['distance_Mpc']) / L_sun
        if integrated_mask_neb:
            integrated_L_int_Ha = integrated_L_obs_Ha
        else:
            integrated_L_int_Ha = integrated_L_obs_Ha * integrated_expqtau[2]
        integrated_SFR_Ha = SFR_Ha__z.sum()
        integrated_SFRSD_Ha = integrated_SFR_Ha / g_props__z['area_pc2'].sum()
        integrated_zone_neb_data = [(
            g['id'],
            -1,
            -1,
            integrated_L_obs_Ha, integrated_L_int_Ha,
            integrated_SFR_Ha, integrated_SFRSD_Ha,
            integrated_flag_BPT, integrated_flag_WHAN,
            integrated_flag_tau_V_neb, integrated_flag_etau_V_neb
        )]
        tbl_integrated_neb.append(integrated_zone_neb_data)
        tbl_integrated_neb.flush()

        # SYN #
        tau_V__z = g_props__z['tau_V'][_izS]
        integrated_tau_V = g_int_props['tau_V']

        for iT, tSF in enumerate(tSF__T):
            x_Y__z, integrated_x_Y = calc_xY(K, tSF)
            SFR__z, SFRSD__z = calc_SFR(K, tSF)

            flag_xY = np.less(x_Y__z, args.minpopx)
            flag_tau_V = np.less(tau_V__z, args.mintauv)

            tmp = np.zeros((g['N_zone']), dtype=np.int)
            id_tSF = tmp + iT
            id_gal = tmp + g['id']
            zone_SF_data = zip(
                id_gal,
                id_zones,
                np.arange(g['N_zone']),
                id_tSF,
                x_Y__z,
                SFR__z,
                SFRSD__z,
                flag_xY,
                flag_tau_V
            )
            tbl_zone_SF.append(zone_SF_data)
            tbl_zone_SF.flush()
            del tmp

            integrated_flag_xY = np.less(integrated_x_Y, args.minpopx)
            integrated_flag_tau_V = np.less(integrated_tau_V, args.mintauv)

            integrated_SFR = SFR__z.sum()
            integrated_SFRSD = integrated_SFR / g_props__z['area_pc2'].sum()
            integrated_SF_data = [(
                g['id'],
                -1,
                -1,
                iT,
                integrated_x_Y,
                integrated_SFR,
                integrated_SFRSD,
                integrated_flag_xY,
                integrated_flag_tau_V
            )]
            tbl_integrated_SF.append(integrated_SF_data)
            tbl_integrated_SF.flush()

        K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (g['califaID'], time.clock() - t_init_gal)

    tbl_zone_SF.cols.id_gal.create_csindex()
    tbl_zone_SF.cols.id_zone.create_csindex()
    tbl_zone_SF.cols.id_tSF.create_csindex()
    tbl_zone_SF.cols.i_zone.create_index()
    tbl_zone_SF.cols.flag_xY.create_index()
    tbl_zone_SF.cols.flag_tau_V.create_index()

    tbl_zone_neb.cols.id_gal.create_csindex()
    tbl_zone_neb.cols.id_zone.create_csindex()
    tbl_zone_neb.cols.i_zone.create_index()
    tbl_zone_neb.cols.flag_BPT.create_index()
    tbl_zone_neb.cols.flag_WHAN.create_index()
    tbl_zone_neb.cols.flag_tau_V_neb.create_index()
    tbl_zone_neb.cols.flag_etau_V_neb.create_index()

    tbl_integrated_SF.cols.id_gal.create_csindex()
    tbl_integrated_SF.cols.id_tSF.create_csindex()
    tbl_integrated_SF.cols.flag_xY.create_index()
    tbl_integrated_SF.cols.flag_tau_V.create_index()

    tbl_integrated_neb.cols.id_gal.create_csindex()
    tbl_integrated_neb.cols.flag_BPT.create_index()
    tbl_integrated_neb.cols.flag_WHAN.create_index()
    tbl_integrated_neb.cols.flag_tau_V_neb.create_index()
    tbl_integrated_neb.cols.flag_etau_V_neb.create_index()

    tbl_zone_SF.flush()
    tbl_zone_neb.flush()
    tbl_integrated_SF.flush()
    tbl_integrated_neb.flush()

    h5file.close()
    print 'total time: %.2f' % (time.clock() - t_init_prog)
