#
# Lacerda@Saco - 25/Mar/2016
#
import os
import sys
import time
import numpy as np
import tables as tbl
from PhD_dust2gas import *
from CALIFAUtils.scripts import create_zones_masks_gal
from pytu.objects import CustomArgumentParser, tupperware_none


def parser_args(default_args_file='default.args'):
    dflt = {
        'debug': False,
        'hdf5': 'output.h5',
        # 'group': 'RadialProfiles',
        'group_Zstar': 'Zstar',
        'group_SF': 'SFR05050525',
        'group_description': None,
        'rbinini': 0.,
        'rbinfin': 3.,
        'rbinstep': 0.1,
    }

    parser = CustomArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('--debug', '-D', action='store_true', default=dflt['debug'])
    parser.add_argument('--hdf5', '-H', metavar='FILE', type=str, default=dflt['hdf5'])
    parser.add_argument('--group_description', metavar='STR', type=str, default=dflt['group_description'])
    parser.add_argument('--group_Zstar', metavar='STR', type=str, default=dflt['group_Zstar'])
    parser.add_argument('--group_SF', metavar='STR', type=str, default=dflt['group_SF'])
    parser.add_argument('--rbinini', metavar='HLR', type=float, default=dflt['rbinini'])
    parser.add_argument('--rbinfin', metavar='HLR', type=float, default=dflt['rbinfin'])
    parser.add_argument('--rbinstep', metavar='HLR', type=float, default=dflt['rbinstep'])
    parser.add_argument('--pycasso_cube_dir', metavar='DIR', type=str, required=True)
    parser.add_argument('--eml_cube_dir', metavar='DIR', type=str, required=True)
    parser.add_argument('--gasprop_cube_dir', metavar='DIR', type=str, required=True)
    args_list = sys.argv[1:]
    # if exists file default.args, load default args
    if os.path.isfile(default_args_file):
        args_list.insert(0, '@%s' % default_args_file)
    debug_var(True, args_list=args_list)
    # return parser.parse_args(args=args_list)
    return fix_dir_args(parser.parse_args(args=args_list))


if __name__ == '__main__':
    # Saving the initial time
    t_init_prog = time.clock()
    lines_bpt = [4861, 5007, 6563, 6583]

    # Parse arguments
    args = parser_args(default_args_file='radial_profiles_default.args')

    R_bin__r = np.arange(args.rbinini, args.rbinfin + args.rbinstep, args.rbinstep)
    R_bin_center__r = (R_bin__r[:-1] + R_bin__r[1:]) / 2.0
    N_R_bins = len(R_bin_center__r)

    # read hdf5 tables
    h5file = tbl.open_file(args.hdf5, mode='r+')
    tbl_gals = h5file.root.pycasso.main
    tbl_zones = h5file.root.pycasso.zones
    group_SF = h5file.get_node('/'+args.group_SF)
    tbl_tSF = group_SF.tSF
    tbl_zones_SF = group_SF.zones_SF
    group_Zstar = h5file.get_node('/'+args.group_Zstar)
    tbl_tZ = group_Zstar.tZ
    tbl_zones_Z = group_Zstar.zones_Z

    # create radial_profiles hdf5 table
    # group_description = ''
    # if args.group_description is not None:
    #    group_description = args.group_description
    # group = h5file.create_group('/', args.group, group_description,
    #                             filters=tbl.Filters(1))
    # table_description = 'Light fraction in young populations'
    # tbl_x_Y = h5file.create_table(group, 'xY__T', radprof_td, table_description)

    gTr_shape = (len(tbl_gals), len(tbl_tSF), N_R_bins)
    gUr_shape = (len(tbl_gals), len(tbl_tZ), N_R_bins)
    R = tupperware_none()
    R.x_Y__gTr = np.ma.masked_all(gTr_shape)
    R.x_Y_std__gTr = np.ma.masked_all(gTr_shape)
    R.McorSD__gTr = np.ma.masked_all(gTr_shape)
    R.McorSD_std__gTr = np.ma.masked_all(gTr_shape)
    R.SFRSD__gTr = np.ma.masked_all(gTr_shape)
    R.SFRSD_std__gTr = np.ma.masked_all(gTr_shape)
    R.SFRSD_neb__gTr = np.ma.masked_all(gTr_shape)
    R.SFRSD_neb_std__gTr = np.ma.masked_all(gTr_shape)
    R.tau_V__gTr = np.ma.masked_all(gTr_shape)
    R.tau_V_std__gTr = np.ma.masked_all(gTr_shape)
    R.tau_V_neb__gTr = np.ma.masked_all(gTr_shape)
    R.tau_V_neb_std__gTr = np.ma.masked_all(gTr_shape)
    R.logOH__gTr = np.ma.masked_all(gTr_shape)
    R.logOH_std__gTr = np.ma.masked_all(gTr_shape)
    # R.alogZ_flux__gUr = np.ma.masked_all(gUr_shape)
    # R.alogZ_mass__gUr = np.ma.masked_all(gUr_shape)

    for i_gal, g in enumerate(tbl_gals):
        t_init_gal = time.clock()

        # Load all data from fits
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

        # Set galaxy geometry using ellipses
        pa, ba = K.getEllipseParams()
        K.setGeometry(pa, ba)

        aux = create_zones_masks_gal(K, tbl_tSF.cols.age[:],
                                     return_mask_lines_separated=True,
                                     mask_lines_snr_only=True,
                                     mintauv=0.05,
                                     mintauvneb=0.05,
                                     maxtauvneberr=0.25,
                                     minpopx=0.05,
                                     minEWHb=3,
                                     minSNR=3,
                                     minSNRHb=3,
                                     nolinecuts=False,
                                     rgbcuts=True,
                                     underS06=False,
                                     whanSF=False,
                                     filter_residual=True,
                                     summary=True)
        mask__Tz = aux[0]
        mask_syn__Tz = aux[1]
        mask_eml__z = aux[2]
        mask_popx__Tz = aux[3]
        mask_tau_V__z = aux[4]
        mask_residual__z = aux[5]
        mask_tau_V_neb__z = aux[6]
        mask_tau_V_neb_err__z = aux[7]
        mask_EW_Hb__z = aux[8]
        mask_whan__z = aux[9]
        mask_bpt__z = aux[10]
        mask_lines_dict__Lmz = aux[11]

        debug_var(True, galaxy=g['califaID'], N_zone=g['N_zone'])
        g_props__z = tbl_zones.read_where('id_gal == gid', {'gid': g['id']})
        g_neb_props__z = group_SF.zones_neb.read_where('id_gal == gid', {'gid': g['id']})

        Mcor__z = g_props__z['Mcor']
        McorSD__z = g_props__z['McorSD']
        tau_V__z = g_props__z['tau_V']
        at_mass__z = g_props__z['at_mass']
        at_flux__z = g_props__z['at_flux']
        SFR_neb__z = g_neb_props__z['SFR']
        SFRSD_neb__z = g_neb_props__z['SFRSD']
        tau_V_neb__z = g_props__z['tau_V_neb']
        logOH__z = g_props__z['logOH']

        # # Mask: zones to create radial profile.
        # mask_residual = np.copy(g_props__z['flag_residual'])
        # mask_neb = np.bitwise_or(g_props__z['flag_RGB'], g_neb_props__z['flag_tau_V_neb'])
        # mask_neb = np.bitwise_or(mask_neb, g_neb_props__z['flag_etau_V_neb'])

        # for i_U, tZ in enumerate(tbl_tZ):
        #     # reading galaxies data
        #     g_props_Z__z = tbl_zones_Z.read_where('(id_gal == gid) & (id_tZ == iU)', {'gid': g['id'], 'iU': tZ['id']})
        #     id_zones = g_props_Z__z['id_zone']
        #     _izS = np.argsort(id_zones)  # zone index sorted by id
        #
        #     alogZ_flux__z = g_props_Z__z['alogZ_flux']
        #     alogZ_mass__z = g_props_Z__z['alogZ_mass']
        #
        #     alogZ_flux__gUr = np.ma.masked_all(gUr_shape)
        #     alogZ_mass__gUr = np.ma.masked_all(gUr_shape)

        for i_T, tSF in enumerate(tbl_tSF):
            # reading galaxies data
            g_props_SF__z = tbl_zones_SF.read_where('(id_gal == gid) & (id_tSF == iT)', {'gid': g['id'], 'iT': tSF['id']})
            id_zones = g_props_SF__z['id_zone']
            _izS = np.argsort(id_zones)  # zone index sorted by id

            x_Y__z = g_props_SF__z['xY']
            SFR__z = g_props_SF__z['SFR']
            SFRSD__z = g_props_SF__z['SFRSD']

            mask_tmp = np.bitwise_or(mask__Tz[i_T], ~(np.isfinite(Mcor__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(McorSD__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(at_mass__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(at_flux__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(tau_V__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(tau_V_neb__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(SFRSD__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(SFRSD_neb__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(logOH__z)))
            mask_tmp = np.bitwise_or(mask_tmp, ~(np.isfinite(x_Y__z)))

            N_zones_masked = mask_tmp.astype('int').sum()

            debug_var(True, pref='>>>', Galaxy=g['califaID'], tSF=tSF['age'], zones=g['N_zone'], masked=N_zones_masked)

            x_Y_m__z = np.ma.masked_array(x_Y__z, mask=mask_tmp)
            Mcor_m__z = np.ma.masked_array(Mcor__z, mask=mask_tmp)
            McorSD_m__z = np.ma.masked_array(McorSD__z, mask=mask_tmp)
            at_mass_m__z = np.ma.masked_array(at_mass__z, mask=mask_tmp)
            at_flux_m__z = np.ma.masked_array(at_flux__z, mask=mask_tmp)
            tau_V_m__z = np.ma.masked_array(tau_V__z, mask=mask_tmp)
            tau_V_neb_m__z = np.ma.masked_array(tau_V_neb__z, mask=mask_tmp)
            SFR_m__z = np.ma.masked_array(SFR__z, mask=mask_tmp)
            SFR_neb_m__z = np.ma.masked_array(SFR_neb__z, mask=mask_tmp)
            logOH_m__z = np.ma.masked_array(logOH__z, mask=mask_tmp)

            R.x_Y__gTr[i_gal, i_T, :] = K.zoneToRad(x_Y_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False)
            R.x_Y_std__gTr[i_gal, i_T, :] = K.zoneToRad(x_Y_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False, mode='std')
            R.McorSD__gTr[i_gal, i_T, :] = K.zoneToRad(McorSD_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False)
            R.McorSD_std__gTr[i_gal, i_T, :] = K.zoneToRad(McorSD_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False, mode='std')
            R.SFRSD__gTr[i_gal, i_T, :] = K.zoneToRad(SFR_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=True, surface_density=True)
            R.SFRSD_std__gTr[i_gal, i_T, :] = K.zoneToRad(SFR_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=True, surface_density=True, mode='std')
            R.SFRSD_neb__gTr[i_gal, i_T, :] = K.zoneToRad(SFR_neb_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=True, surface_density=True)
            R.SFRSD_neb_std__gTr[i_gal, i_T, :] = K.zoneToRad(SFR_neb_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=True, surface_density=True, mode='std')
            R.tau_V__gTr[i_gal, i_T, :] = K.zoneToRad(tau_V_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False)
            R.tau_V_std__gTr[i_gal, i_T, :] = K.zoneToRad(tau_V_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False, mode='std')
            R.tau_V_neb__gTr[i_gal, i_T, :] = K.zoneToRad(tau_V_neb_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False)
            R.tau_V_neb_std__gTr[i_gal, i_T, :] = K.zoneToRad(tau_V_neb_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False, mode='std')
            R.logOH__gTr[i_gal, i_T, :] = K.zoneToRad(logOH_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False)
            R.logOH_std__gTr[i_gal, i_T, :] = K.zoneToRad(logOH_m__z, R_bin__r, rad_scale=K.HLR_pix, extensive=False, mode='std')

        K.EL.close()
        K.GP.close()
        K.close()

    h5file.close()

    pickle_dir = '%s/dev/astro/PhD_dust2gas/runs/pickles' % os.environ['HOME']
    pickle_ext = '.pkl'

    print 'dumping pickles...'
    for k in R.__dict__:
        attr = getattr(R, k)
        if isinstance(attr, np.ma.MaskedArray):
            debug_var(True, attr=k)
            attr.dump('%s/%s%s' % (pickle_dir, k, pickle_ext))

    # x_Y__gTr.dump('%s/x_Y__gTr%s' % (pickle_dir, pickle_ext))
    # McorSD__gTr.dump('%s/McorSD__gTr%s' % (pickle_dir, pickle_ext))
    # SFRSD__gTr.dump('%s/SFRSD__gTr%s' % (pickle_dir, pickle_ext))
    # SFRSD_neb__gTr.dump('%s/SFRSD_neb__gTr%s' % (pickle_dir, pickle_ext))
    # tau_V__gTr.dump('%s/tau_V__gTr%s' % (pickle_dir, pickle_ext))
    # tau_V_neb__gTr.dump('%s/tau_V_neb__gTr%s' % (pickle_dir, pickle_ext))
    # logOH__gTr.dump('%s/logOH__gTr%s' % (pickle_dir, pickle_ext))
    # x_Y_std__gTr.dump('%s/x_Y_std__gTr%s' % (pickle_dir, pickle_ext))
    # McorSD_std__gTr.dump('%s/McorSD_std__gTr%s' % (pickle_dir, pickle_ext))
    # SFRSD_std__gTr.dump('%s/SFRSD_std__gTr%s' % (pickle_dir, pickle_ext))
    # SFRSD_neb_std__gTr.dump('%s/SFRSD_neb_std__gTr%s' % (pickle_dir, pickle_ext))
    # tau_V_std__gTr.dump('%s/tau_V_std__gTr%s' % (pickle_dir, pickle_ext))
    # tau_V_neb_std__gTr.dump('%s/tau_V_neb_std__gTr%s' % (pickle_dir, pickle_ext))
    # logOH_std__gTr.dump('%s/logOH_std__gTr%s' % (pickle_dir, pickle_ext))
