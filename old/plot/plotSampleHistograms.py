#!/usr/bin/python
#
# Lacerda@Granada - 02/Dec/2015
#
#from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from os.path import basename
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

mpl.rcParams['font.size'] = 14
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10 
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'

default_sc_kwargs = dict(marker = 'o', s = 10, edgecolor = 'none', cmap = 'Spectral')
default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.02)
default_ols_plot_kwargs = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')
default_ols_kwargs = dict(c = 'r', pos_x = 0.98, pos_y = 0.01, fs = 15, rms = True, text = True, kwargs_plot = default_ols_plot_kwargs)

A4Size_inches = [ 8.267, 11.692 ]
LetterSize_inches = [ 8.5, 11 ]

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'hdf5' : None,
        'output' : None,
        'itSF' : 11,
        'maskradius' : None,
        'slice_gals' : None,
        'dryrun' : False,
        'bamin' : 0,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--dryrun',
                        action = 'store_true',
                        default = default['dryrun'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--itSF', '-T',
                        help = 'SF age index.',
                        metavar = '',
                        type = int,
                        default = default['itSF'])
    parser.add_argument('--maskradius', '-R',
                        help = 'initial RDisc value in HLR',
                        metavar = 'NUM',
                        type = float,
                        default = default['maskradius'])
    parser.add_argument('--output', '-o',
                        help = 'Name of the output PDF file.',
                        metavar = 'FILENAME',
                        type = str,
                        default = default['output'])
    parser.add_argument('--bamin', '-B',
                        help = 'min b/a',
                        metavar = '',
                        type = float,
                        default = default['bamin'])

    return parser.parse_args()

if __name__ == '__main__':
    args = parser_args()
    
    C.debug_var(args.debug, args = args)
    
    H = C.H5SFRData(args.hdf5)
    iT = args.itSF
    iU = -1

    fnamesuffix = '.pdf'
    minR = args.maskradius

    if minR is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        maxR = H.Rbin__r[-1]
        maskRadiusOk__g = (H.zone_dist_HLR__g >= minR) & (H.zone_dist_HLR__g <= maxR) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * ((H.RbinCenter__r >= minR) & (H.RbinCenter__r <= maxR))).T
        fnamesuffix = '_maskradius%s' % fnamesuffix
        
    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
        gals_slice__integr = np.ones(H.califaIDs_all.shape, dtype = np.bool)
        gals_txt = ''
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)
        gals_txt = basename(args.slice_gals)
        fnamesuffix = '_%s%s' % ('.'.join(gals_txt.split('.')[:-1]), fnamesuffix)
        #integrated
        l_gals, _ = C.sort_gals(args.slice_gals)
        l_gals = sorted(l_gals)
        gals_slice__integr = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
        for g in l_gals:
            i = H.califaIDs_all.tolist().index(g)
            gals_slice__integr[i] = True    
    
    ba_max = args.bamin
    mask_GAL__g = np.bitwise_or(np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool), np.less(H.ba_GAL__g, ba_max))
    
    mask__g = np.bitwise_or(np.ma.log10(H.SFRSD__Tg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Tg[iT]).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.SFRSD_Ha__g * 1e6).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    mask__g = np.bitwise_or(mask__g, H.logO3N2_M13__g.mask)
    #mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), ba_max))
    mask__g = np.bitwise_or(mask__g, ~maskRadiusOk__g)
    mask__g = np.bitwise_or(mask__g, ~gals_slice__g)

    mask__rg = mask__rg = np.bitwise_or(~maskRadiusOk__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), ba_max))
    mask__rg = np.bitwise_or(mask__rg, ~gals_slice__rg)
    
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    
    xm, ym, zm = C.ma_mask_xyz(x = np.ma.log10(H.tau_V__Tg[iT]), y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), z = np.ma.log10(H.x_Y__Tg[iT]), mask = mask__g)
    C.OLS_bisector(xm, ym, debug = True)
    
    print ((~mask__g).sum()),((~mask__rg).sum()), NgalsOkZones, NgalsOkRbins

    count = 0
    fraclim = 5
    
    for gal in np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]):
        zoneDistance_HLR = H.get_prop_gal(H.zone_dist_HLR__g, gal)
        igal_all = H.califaIDs_all.tolist().index(gal)
        igal = H.califaIDs.tolist().index(gal)
        Nz = len(zoneDistance_HLR[zoneDistance_HLR >= minR])
        mask__z = H.get_prop_gal(mask__g, gal)
        mask__r = H.get_prop_gal(mask__rg, gal)
        NzOk = (~(mask__z)).astype(int).sum()
        NrOk = (~(mask__r)).astype(int).sum()
        NzOkfrac = 100*NzOk / Nz 
        NrOkfrac = 100*NrOk / H.NRbins 
        
        if (NzOkfrac < fraclim) or (NzOkfrac < fraclim):
            count += 1
            pref = '>>>'
        else:
            pref = ''
        print pref, gal, len(zoneDistance_HLR), Nz, mask__z.sum(), mask__r.sum(), NzOk, NrOk, NzOkfrac, NrOkfrac
        
    print count
    D_gals = {
        gal : {
            'morfType' : H.morfType_GAL__g[H.califaIDs_all.tolist().index(gal)], 
            'Mcor' : H.Mcor_GAL__g[H.califaIDs_all.tolist().index(gal)], 
            'McorSD' : H.McorSD_GAL__g[H.califaIDs_all.tolist().index(gal)], 
            'atflux' : H.at_flux_GAL__g[H.califaIDs_all.tolist().index(gal)]
        } for gal in np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g])
    }

    ######################### Morf #########################
    ######################### Morf #########################
    ######################### Morf #########################
    
    maskSaSab = (H.morfType_GAL__g < 2) & (H.morfType_GAL__g >= 0)
    maskSb = (H.morfType_GAL__g == 2.)
    maskSbc = (H.morfType_GAL__g == 3.)
    maskScScd = (H.morfType_GAL__g < 6.) & (H.morfType_GAL__g >= 4.)
    maskSdSmIrr = (H.morfType_GAL__g >= 6) 
    mask_morf = [ maskSaSab, maskSb, maskSbc, maskScScd, maskSdSmIrr  ]
    color_morf = [ 'orange', 'green', '#00D0C9', '#0076C9', 'blue' ]
    label_morf = [ 'Sa + Sab', 'Sb', 'Sbc', 'Sc + Scd', 'Sd + Sm + Irr' ]
    tagname_morf = [ l.replace(' + ', '') for l in label_morf ]

    if args.dryrun: sys.exit()
    ####################################################



    with PdfPages('sample%s' % fnamesuffix) as pdf:
        f = plt.figure()
        sc_kwargs = default_sc_kwargs.copy()
        rs_kwargs = default_rs_kwargs.copy()
        ols_kwargs = default_ols_kwargs.copy()
        suptitle_R = r'Nzones:%d(%d gals) $t_{SF}$:%.2fMyr ' % ((~mask__g).sum(), NgalsOkZones, (H.tSF__T[iT] / 1e6))
        NRows = 2
        NCols = 3
        page_size_inches = [11, 11]
        f.set_size_inches(page_size_inches)
        grid_shape = (NRows, NCols)
        f.suptitle(suptitle_R, fontsize = 15)

        
    
        f.subplots_adjust(bottom = 0.15, hspace = 0.25, wspace = 0.25, right = 0.95, left = 0.07)
        plt.close(f)
        pdf.savefig(f)
        
        
