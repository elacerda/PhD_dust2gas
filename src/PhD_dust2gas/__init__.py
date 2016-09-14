from CALIFAUtils.scripts import debug_var


__author__ = 'lacerda@astro.ufsc.br'


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
    from pycasso import fitsQ3DataCube
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
        from CALIFAUtils.objects import GasProp
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
