#
# Lacerda@Saco - 25/Mar/2016
#
import os
import sys
import time
import numpy as np
import tables as tbl
from PhD_dust2gas import *
from pytu.functions import debug_var
from pystarlight.util import redenninglaws
from PhD_dust2gas.tables_description import zone_Z
from PhD_dust2gas.tables_description import tZ as tMet
from pytu.objects import CustomArgumentParser


def parser_args(default_args_file='default.args'):
    default_args = {
        'debug': False,
        'hdf5': 'output.h5',
        'group': 'Zstar',
        'group_SF': 'SFR05050525'
    }

    parser = CustomArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--debug', '-D', action='store_true',
                        default=default_args['debug'])
    parser.add_argument('--hdf5', '-H', metavar='FILE', type=str,
                        default=default_args['hdf5'])
    parser.add_argument('--pycasso_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--pycasso_cube_suffix', metavar='SUFFIX', type=str)
    parser.add_argument('--eml_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--eml_cube_suffix', metavar='SUFFIX', type=str)
    parser.add_argument('--gasprop_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--gasprop_cube_suffix', metavar='SUFFIX', type=str)
    parser.add_argument('--group', metavar='GROUPNAME', type=str,
                        default=default_args['group'])
    parser.add_argument('--group_SF', metavar='GROUPNAME', type=str,
                        default=default_args['group_SF'])

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


def calc_alogZ_Stuff(K, tZ, xOkMin):
    '''
    Compute average logZ (alogZ_*) both for each zone (*__z) and the galaxy-wide
    average (*_GAL, computed a la GD14).

    Only st-pops satisfying the input flag__t (ageBase-related) mask are considered!
    This allows us to compute alogZ_* for, say, < 2 Gyr or 1--7 Gyr populations,
    as well as the whole-age range (using flag__t=True for all base ages)
    with a same function and saving the trouble of keeping separate variables for the same thing:-)

    ==> return alogZ_mass_GAL, alogZ_flux_GAL, isOkFrac_GAL

    Cid@Lagoa - 05/Jun/2014

    !!HELP!! ATT: My way of computing alogZ_*__z produces nan', which are ugly but harmless.
    I tried to fix it using masked arrays:

    alogZ_mass__z  = np.ma.masked_array( numerator__z/(denominator__z+0e-30) , mask=(denominator__z == 0))

    but this did not work!

    Cid@Lagoa - 20/Jun/2014

    Correct nan problems using:
    alogZ_mass__z[np.isnan(alogZ_mass__z)] = np.ma.masked
    Lacerda@Granada - 19/Feb/2015

    removed radial profiles inside this func.
    Lacerda@Granada - 23/Feb/2015
    '''
    # --------------------------------------------------------------------------
    # Initialization
    Zsun = 0.019
    # Define log of base metallicities **in solar units** for convenience
    logZBase__Z = np.log10(K.metBase / Zsun)
    # --------------------------------------------------------------------------
    flag__t = K.ageBase <= tZ

    # --------------------------------------------------------------------------
    # Define alogZ_****__z: flux & mass weighted average logZ for each zone
    # ==> alogZ_mass__z - ATT: There may be nan's here depending on flag__t!
    numerator__z = np.tensordot(K.Mcor__tZz[flag__t, :, :], logZBase__Z, (1, 0)).sum(axis=0)
    denominator__z = K.Mcor__tZz[flag__t, :, :].sum(axis=1).sum(axis=0)
    alogZ_mass__z = np.ma.masked_array(numerator__z / denominator__z)
    alogZ_mass__z[np.isnan(alogZ_mass__z)] = np.ma.masked

    # ==> alogZ_flux__z - ATT: There may be nan's here depending on flag__t!
    numerator__z = np.tensordot(K.Lobn__tZz[flag__t, :, :], logZBase__Z, (1, 0)).sum(axis=0)
    denominator__z = K.Lobn__tZz[flag__t, :, :].sum(axis=1).sum(axis=0)
    alogZ_flux__z = np.ma.masked_array(numerator__z / denominator__z)
    alogZ_flux__z[np.isnan(alogZ_mass__z)] = np.ma.masked
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Def galaxy-wide averages of alogZ in light & mass, but **discards** zones having
    # too little light fractions in the ages given by flag__t
    isOk__z = np.ones_like(K.Mcor__z, dtype=np.bool)

    # Define Ok flag: Zones with light fraction x < xOkMin are not reliable for alogZ (& etc) estimation!
    x_Y__z, int_x_Y__z = calc_xY(K, tZ)

    if xOkMin >= 0.:
        # EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # x__tZz = K.popx / K.popx.sum(axis=1).sum(axis=0)
        # xOk__z = x__tZz[flag__t, :, :].sum(axis=1).sum(axis=0)
        # EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        isOk__z = (x_Y__z >= xOkMin)
        isOk_int = (int_x_Y__z >= xOkMin)

    # Fraction of all zones which are Ok in the isOk__z sense. Useful to censor galaxies whose
    # galaxy-wide averages are based on too few zones (hence not representative)
    # OBS: isOkFrac_GAL is not used in this function, but it's returned to be used by the caller
    isOkFrac_GAL = (1.0 * isOk__z.sum()) / (1.0 * K.N_zone)

    # Galaxy wide averages of logZ - ATT: Only isOk__z zones are considered in this averaging!
    numerator__z = K.Mcor__tZz[flag__t, :, :].sum(axis=1).sum(axis=0) * alogZ_mass__z
    denominator__z = K.Mcor__tZz[flag__t, :, :].sum(axis=1).sum(axis=0)
    alogZ_mass_GAL = numerator__z[isOk__z].sum() / denominator__z[isOk__z].sum()

    numerator__z = K.Lobn__tZz[flag__t, :, :].sum(axis=1).sum(axis=0) * alogZ_flux__z
    denominator__z = K.Lobn__tZz[flag__t, :, :].sum(axis=1).sum(axis=0)
    alogZ_flux_GAL = numerator__z[isOk__z].sum() / denominator__z[isOk__z].sum()

    return alogZ_mass__z, alogZ_flux__z, alogZ_mass_GAL, alogZ_flux_GAL, isOk__z, isOk_int, isOkFrac_GAL


if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()
    lines_bpt = [4861, 5007, 6563, 6583]

    # Parse arguments
    args = parser_args(default_args_file='Zstar_default.args')
    debug_var(True, args=args.__dict__)

    tZ__U = np.array([1.0, 2.0, 5.0, 8.0, 11.3, 14.2]) * 1.e9
    N_U = len(tZ__U)

    q = redenninglaws.Cardelli_RedLaw(lines_bpt)

    h5file = tbl.open_file(args.hdf5, mode='r+')
    tbl_gals = h5file.root.pycasso.main
    tbl_zones = h5file.root.pycasso.zones
    tbl_integrated = h5file.root.pycasso.integrated

    minpopx = np.float(h5file.get_node('/'+args.group_SF)._g_gettitle().split('/')[0].split(':')[-1])

    group_description = 'minpopx:%.2f - zones calculation' % minpopx
    group = h5file.create_group('/', args.group, group_description,
                                filters=tbl.Filters(1))
    tbl_tZ = h5file.create_table(group, 'tZ', tMet, 'tZ data')
    tbl_zone_Z = h5file.create_table(group, 'zones_Z', zone_Z, 'Z data')
    tbl_integrated_Z = h5file.create_table(group, 'integrated_Z', zone_Z, 'Z data')

    tbl_tZ.append([tZ for tZ in enumerate(tZ__U)])
    tbl_tZ.cols.id.create_csindex()
    tbl_tZ.flush()

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

        g_props = tbl_zones.read_where('id_gal == gid', {'gid': g['id']})
        g_int_props = tbl_integrated.read_where('id_gal == gid', {'gid': g['id']})
        id_zones = g_props['id']
        _izS = np.argsort(g_props['id'])  # zone index sorted by id

        # SYN #
        for iU, tZ in enumerate(tZ__U):
            # All calculations inside calc_alogZ_Stuff().
            aux = calc_alogZ_Stuff(K, tZ, minpopx)
            alogZ_mass__z = aux[0]
            alogZ_flux__z = aux[1]
            alogZ_mass_GAL = aux[2]
            alogZ_flux_GAL = aux[3]
            flagOk_xY = aux[4]
            integrated_flagOk_xY = aux[5]
            isOkFrac_GAL = aux[6]

            flag_xY = ~(flagOk_xY)

            tmp = np.zeros((g['N_zone']), dtype=np.int)
            id_tZ = tmp + iU
            id_gal = tmp + g['id']
            zone_Z_data = zip(
                id_gal,
                id_zones,
                np.arange(g['N_zone']),
                id_tZ,
                alogZ_flux__z,
                alogZ_mass__z,
                flag_xY
            )
            tbl_zone_Z.append(zone_Z_data)
            tbl_zone_Z.flush()
            del tmp

            integrated_flag_xY = ~(integrated_flagOk_xY)
            integrated_Z_data = [(
                g['id'],
                -1,
                -1,
                iU,
                alogZ_flux_GAL,
                alogZ_mass_GAL,
                integrated_flag_xY
            )]
            tbl_integrated_Z.append(integrated_Z_data)
            tbl_integrated_Z.flush()

        K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (g['califaID'], time.clock() - t_init_gal)

    tbl_zone_Z.cols.id_gal.create_csindex()
    tbl_zone_Z.cols.id_zone.create_csindex()
    tbl_zone_Z.cols.id_tZ.create_csindex()
    tbl_zone_Z.cols.i_zone.create_index()
    tbl_zone_Z.cols.flag_xY.create_index()

    tbl_integrated_Z.cols.id_gal.create_csindex()
    tbl_integrated_Z.cols.id_tZ.create_csindex()
    tbl_integrated_Z.cols.flag_xY.create_index()

    tbl_zone_Z.flush()
    tbl_integrated_Z.flush()

    h5file.close()
    print 'total time: %.2f' % (time.clock() - t_init_prog)
