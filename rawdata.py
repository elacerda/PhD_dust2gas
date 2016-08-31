
# Lacerda@Saco - 15/Mar/2016
#
import os
import sys
import time
import tables as tbl
import numpy as np
from pycasso import fitsQ3DataCube
from tables_description import zone
from tables_description import galaxy
from CALIFAUtils.scripts import my_morf
from CALIFAUtils.objects import GasProp
from pytu.functions import debug_var
from CALIFAUtils.scripts import sort_gals
from CALIFAUtils.scripts import get_morfologia
from pytu.objects import CustomArgumentParser


def parser_args(default_args_file='default.args'):
    '''
        Parse the command line args
        With fromfile_prefix_chars=@ we can read and parse command line args
        inside a file with @file.txt.
        default args inside default_args_file
    '''
    default_args = {
        'debug': False,
        'hdf5': 'output.h5',
        'gals': 'listv20_q050.d15a.txt',
    }

    parser = CustomArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--debug', '-D', action='store_true',
                        default=default_args['debug'])
    parser.add_argument('--hdf5', '-H', metavar='FILE', type=str,
                        default=default_args['hdf5'])
    parser.add_argument('--gals', '-G', metavar='FILE', type=str,
                        default=default_args['gals'])
    parser.add_argument('--pycasso_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--pycasso_cube_suffix', metavar='SUFFIX', type=str,
                        required=True)
    parser.add_argument('--eml_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--eml_cube_suffix', metavar='SUFFIX', type=str,
                        required=True)
    parser.add_argument('--gasprop_cube_dir', metavar='DIR', type=str,
                        required=True)
    parser.add_argument('--gasprop_cube_suffix', metavar='SUFFIX', type=str,
                        required=True)
    parser.add_argument('--morph_file', '-M', metavar='FILE', type=str)

    args_list = sys.argv[1:]
    # if exists file default.args, load default args
    if os.path.isfile(default_args_file):
        args_list.insert(0, '@%s' % default_args_file)
    debug_var(True, args_list=args_list)
    return fix_dir_args(parser.parse_args(args=args_list))


def fix_dir_args(args):
    args.EL = (args.eml_cube_dir is not None)
    args.GP = (args.gasprop_cube_dir is not None)

    if args.pycasso_cube_dir[-1] != '/':
        args.pycasso_cube_dir += '/'
    if args.EL and (args.eml_cube_dir[-1] != '/'):
        args.eml_cube_dir += '/'
    if args.GP and (args.gasprop_cube_dir[-1] != '/'):
        args.gasprop_cube_dir += '/'

    return args


def load_gal_cubes(args, califaID,
                   pycasso_cube_file=None,
                   eml_cube_file=None,
                   gasprop_cube_file=None):
    '''
        Open PyCASSO SUPERFITS (K), EmissionLines FITS (K.EL) and
        GasProp FITS (K.GP)

        The directories and suffixes are in:
        args.[pycasso_cube_[dir,suffix],
              eml_cube_[dir,suffix],
              gasprop_cube_[dir,suffix]]
    '''
    if pycasso_cube_file is None:
        pycasso_cube_file = califaID + args.pycasso_cube_suffix
    K_cube = fitsQ3DataCube(args.pycasso_cube_dir + pycasso_cube_file)
    filenames = [pycasso_cube_file, '', '']
    if args.EL:
        if eml_cube_file is None:
            eml_cube_file = califaID + args.eml_cube_suffix
        K_cube.loadEmLinesDataCube(args.eml_cube_dir + eml_cube_file)
        filenames.insert(1, eml_cube_file)
    if args.GP:
        # GasProp in CALIFAUtils.objects
        if gasprop_cube_file is None:
            gasprop_cube_file = califaID + args.gasprop_cube_suffix
        K_cube.GP = GasProp(args.gasprop_cube_dir + gasprop_cube_file)
        filenames.insert(2, gasprop_cube_file)
    return K_cube, filenames


def verify_files(K, califaID, EL=True, GP=True):
    if K is None:
        print '<<< %s galaxy: miss files' % califaID
        return 0, False
    if EL:
        if K.EL is None:
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


if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()

    # Parse args
    args = parser_args(default_args_file='rawdata_default.args')
    debug_var(True, args=args.__dict__)

    # open h5file
    h5file = tbl.open_file(args.hdf5, mode='w', title='SFR data')

    # create h5 groups and tables within
    group = h5file.create_group('/', 'pycasso', 'Galaxy PyCASSO data',
                                filters=tbl.Filters(1))
    tbl_main = h5file.create_table(group, 'main', galaxy, 'Main data')
    tbl_zone = h5file.create_table(group, 'zones', zone, 'Zone data')
    tbl_integrated = h5file.create_table(group, 'integrated', zone,
                                         'Integrated data')

    # read gals from args.gals file
    gals, _ = sort_gals(args.gals, order=1)
    N_gals = len(gals)
    max_gals = N_gals
    if args.debug:
        max_gals = 10

    # set ids as zero
    id_zone_ini = id_zone_fin = id_gal = 0

    # main loop thru galaxies
    for iGal, gal in enumerate(gals[0:max_gals]):
        t_init_gal = time.clock()

        # load PyCASSO Superfits, EMLines and GasProps
        K, cube_filenames = load_gal_cubes(args, gal)
        sit, verify = verify_files(K, gal, EL=args.EL, GP=args.GP)

        # set files situation
        if verify is not True:
            print '<<< ', gal, sit
            if sit == 1:
                K.close()
            elif sit == 2:
                K.EL.close()
                K.close()
            continue

        # Set galaxy geometry using ellipses
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)

        # setting galaxy inital and final zone ids
        id_zone_ini = id_zone_fin
        id_zone_fin = id_zone_ini + K.N_zone

        # get lines index in K.GP.lines list
        Hb_central_wl = '4861'
        O3_central_wl = '5007'
        Ha_central_wl = '6563'
        N2_central_wl = '6583'
        i_Hb = K.EL.lines.index(Hb_central_wl)
        i_O3 = K.EL.lines.index(O3_central_wl)
        i_Ha = K.EL.lines.index(Ha_central_wl)
        i_N2 = K.EL.lines.index(N2_central_wl)

        # Flag to checkup all RGB EMLines flags:
        # central wl, sigma, S/N, flux >= 0
        # flag_RGB True means that the zone do not fullfill all RGB quality
        # requirements
        flag_RGB_OK__z = np.zeros((K.N_zone), dtype=np.bool_)
        for l in [Hb_central_wl, O3_central_wl, Ha_central_wl, N2_central_wl]:
            pos = K.GP._dlcons[l]['pos']
            sigma = K.GP._dlcons[l]['sigma']
            snr = K.GP._dlcons[l]['SN']
            tmp = K.EL._setMaskLineFluxNeg(l)
            tmp |= K.EL._setMaskLineDisplacement(l, pos)
            tmp |= K.EL._setMaskLineSigma(l, sigma)
            tmp |= K.EL._setMaskLineSNR(l, snr)
            flag_RGB__z = np.bitwise_or(flag_RGB__z, tmp)

        # calculate at_flux_GAL (whole galaxy average at_flux)
        numerator__z = K.Lobn__tZz.sum(axis=1).sum(axis=0) * K.at_flux__z
        denominator__z = K.Lobn__tZz.sum(axis=1).sum(axis=0)
        at_flux_GAL = numerator__z.sum() / denominator__z.sum()

        # AV to tauV
        AVtoTauV = 1. / (np.log10(np.exp(1)) / 0.4)

        '''
            morphological type:
            'S0' : -1, 'S0a' : -1
            'Sa' : 0, 'Sab' : 1, 'Sb' : 2, 'Sbc' : 3, 'Sc' : 4,
            'Scd' : 5, 'Sd' : 6, 'Sdm' : 7, 'Sm' : 7, 'Ir' : 7,
            'E0' : -2, 'E1' : -2, 'E2' : -2, 'E3' : -2, 'E4' : -2, 'E5' : -2,
            'E6' : -2, 'E7' : -2,
        '''
        m_type_orig = get_morfologia(K.califaID, morph_file=args.morph_file)[0]
        m_type = my_morf(m_type_orig)

        # Prepare galaxy main data to fill main table
        # the galaxy ID will be iGal, so, gaps between ids could exists
        main_data = [(
            iGal,
            K.galaxyName,
            K.califaID,
            K.N_zone,
            K.distance_Mpc,
            K.redshift,
            m_type,
            m_type_orig,
            np.float(K.masterListData['ba']),
            ba,
            K.parsecPerPixel,
            np.float(K.masterListData['Mr']),
            np.float(K.masterListData['u-r']),
            K.HLR_pix,
            K.getHalfRadius(K.McorSD__yx),
            cube_filenames[0],
            cube_filenames[1],
            cube_filenames[2],
        )]

        # append data to main table
        tbl_main.append(main_data)
        tbl_main.flush()

        # preparare integrated data.
        # integrated data stored as a zone with id -1
        row_zone = tbl_integrated.row
        row_zone['id_gal'] = iGal
        row_zone['tau_V'] = K.integrated_keywords['A_V']*AVtoTauV
        row_zone['at_flux'] = at_flux_GAL
        row_zone['tau_V_neb'] = K.EL.integrated_tau_V_neb
        row_zone['etau_V_neb'] = K.EL.integrated_tau_V_neb_err
        row_zone['F_obs_Hb'] = K.EL.integrated_flux[i_Hb]
        row_zone['F_obs_O3'] = K.EL.integrated_flux[i_O3]
        row_zone['F_obs_Ha'] = K.EL.integrated_flux[i_Ha]
        row_zone['F_obs_N2'] = K.EL.integrated_flux[i_N2]
        row_zone['eF_obs_Hb'] = K.EL.integrated_eflux[i_Hb]
        row_zone['eF_obs_O3'] = K.EL.integrated_eflux[i_O3]
        row_zone['eF_obs_Ha'] = K.EL.integrated_eflux[i_Ha]
        row_zone['eF_obs_N2'] = K.EL.integrated_eflux[i_N2]
        row_zone['baseline_Hb'] = K.EL.integrated_baseline[i_Hb]
        row_zone['baseline_O3'] = K.EL.integrated_baseline[i_O3]
        row_zone['baseline_Ha'] = K.EL.integrated_baseline[i_Ha]
        row_zone['baseline_N2'] = K.EL.integrated_baseline[i_N2]
        row_zone['EW_Hb'] = K.EL.integrated_EW[i_Hb]
        row_zone['EW_O3'] = K.EL.integrated_EW[i_O3]
        row_zone['EW_Ha'] = K.EL.integrated_EW[i_Ha]
        row_zone['EW_N2'] = K.EL.integrated_EW[i_N2]
        row_zone['sigma_Hb'] = K.EL.integrated_sigma[i_Hb]
        row_zone['sigma_O3'] = K.EL.integrated_sigma[i_O3]
        row_zone['sigma_Ha'] = K.EL.integrated_sigma[i_Ha]
        row_zone['sigma_N2'] = K.EL.integrated_sigma[i_N2]
        row_zone['esigma_Hb'] = K.EL.integrated_esigma[i_Hb]
        row_zone['esigma_O3'] = K.EL.integrated_esigma[i_O3]
        row_zone['esigma_Ha'] = K.EL.integrated_esigma[i_Ha]
        row_zone['esigma_N2'] = K.EL.integrated_esigma[i_N2]
        row_zone['pos_Hb'] = K.EL.integrated_pos[i_Hb]
        row_zone['pos_O3'] = K.EL.integrated_pos[i_O3]
        row_zone['pos_Ha'] = K.EL.integrated_pos[i_Ha]
        row_zone['pos_N2'] = K.EL.integrated_pos[i_N2]
        row_zone['epos_Hb'] = K.EL.integrated_epos[i_Hb]
        row_zone['epos_O3'] = K.EL.integrated_epos[i_O3]
        row_zone['epos_Ha'] = K.EL.integrated_epos[i_Ha]
        row_zone['epos_N2'] = K.EL.integrated_epos[i_N2]
        # append integrated data
        row_zone.append()
        tbl_integrated.flush()

        del row_zone

        # prepare zone data that will be stored in tbl_zone
        # zip to transpose arrays to col_arrays
        zone_data = zip(
            np.zeros((K.N_zone), dtype=int) + id_gal,
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
            flag_RGB__z,
            ~(K.filterResidual(w2=4600)),
        )
        # apppend zone data
        tbl_zone.append(zone_data)
        tbl_zone.flush()

        # close opened fits
        K.GP.close()
        K.EL.close()
        K.close()
        del K
        print 'time per galaxy: %s %.2f' % (gal, time.clock() - t_init_gal)

        id_gal += 1

    # creating primary indexes to relational searches between tables
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

    # close h5files
    h5file.close()

    print 'total time: %.2f' % (time.clock() - t_init_prog)
