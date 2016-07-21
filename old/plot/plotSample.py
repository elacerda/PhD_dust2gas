#!/usr/bin/python
#
# Lacerda@Granada - 02/Dec/2015
#
#from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import plotOLSbisectorAxis
from CALIFAUtils.plots import plot_text_ax
from matplotlib import pyplot as plt
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
default_rs_kwargs = dict(smooth = True, sigma = 1.2, frac = 0.02, inverse = True)
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

    minR = 0
    fnamesuffix = '.pdf'
    
    if args.maskradius is None:
        maskRadiusOk__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        maskRadiusOk__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        minR = args.maskradius
        maskRadiusOk__g = (H.zone_dist_HLR__g >= args.maskradius) & (H.zone_dist_HLR__g <= H.Rbin__r[-1]) 
        maskRadiusOk__rg = (np.ones((H.NRbins, H.N_gals_all), dtype = np.bool).T * (H.RbinCenter__r >= args.maskradius)).T
        
    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)

    if args.slice_gals is None:
        N_gals = H.N_gals
        gals_slice__g = np.ones_like(H.zone_dist_HLR__g, dtype = np.bool)
        gals_slice__rg = np.ones((H.NRbins, H.N_gals_all), dtype = np.bool)
        gals_slice__integr = np.ones(H.califaIDs_all.shape, dtype = np.bool)
        gals_txt = ''
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)
        gals_txt = (args.slice_gals).split('/')[-1]
        fnamesuffix = '_%s%s' % (gals_txt, fnamesuffix)
        #integrated
        l_gals, _ = C.sort_gals(args.slice_gals)
        l_gals = sorted(l_gals)
        gals_slice__integr = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
        for g in l_gals:
            i = H.califaIDs_all.tolist().index(g)
            gals_slice__integr[i] = True

    ##########################
    ######### MASKS ##########
    ##########################
    #mask_GAL__g = np.bitwise_or(np.less(H.integrated_EW_Ha__g, 3.), np.less(H.ba_GAL__g, args.bamin))
    mask_GAL__g = np.bitwise_or(np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool), np.less(H.ba_GAL__g, args.bamin))
    #mask_GAL__g = np.bitwise_or(np.less(H.integrated_EW_Ha__g, 3.), np.zeros_like(H.ba_GAL__g, dtype = np.bool))
    #mask_GAL__g = np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool)
    mask_GAL__g = np.bitwise_or(mask_GAL__g, ~gals_slice__integr)
    
    mask__g = np.bitwise_or(np.ma.log10(H.SFRSD__Tg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Tg[iT]).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.SFRSD_Ha__g * 1e6).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    mask__g = np.bitwise_or(mask__g, H.O_O3N2_M13__g.mask)
    mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), args.bamin))
    mask__g = np.bitwise_or(mask__g, ~maskRadiusOk__g)
    mask__g = np.bitwise_or(mask__g, ~gals_slice__g)
    #mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    #mask__g = ~maskRadiusOk__g
    
    mask__rg = np.bitwise_or(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Trg[iT]).mask)
    mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.aSFRSD_Ha__rg * 1e6).mask)
    mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.tau_V_neb__rg).mask)
    mask__rg = np.bitwise_or(mask__rg, H.O_O3N2_M13__rg.mask)
    mask__rg = np.bitwise_or(mask__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), args.bamin))
    mask__rg = np.bitwise_or(mask__rg, ~maskRadiusOk__rg)
    mask__rg = np.bitwise_or(mask__rg, ~gals_slice__rg)
    #mask__rg = np.bitwise_or(mask__rg, np.less(H.EW_Ha__rg, 3.))
    #mask__rg = ~maskRadiusOk__rg
    
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    
    print NgalsOkZones, NgalsOkRbins
    
    if args.dryrun: sys.exit()
    ####################################################

    with PdfPages('sample.pdf') as pdf:
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
        xm, ym, zm = C.ma_mask_xyz(x = np.ma.log10(H.tau_V__Tg[iT]), y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), z = np.ma.log10(H.x_Y__Tg[iT]), mask = mask__g)
        #xm, ym, zm = C.ma_mask_xyz(x = np.ma.log10(H.tau_V__Tg[iT]), y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), z = np.ma.log10(H.x_Y__Tg[iT]))
        xlabel = r'$\log\ \tau_V$'
        ylabel = r'$\log\ \Sigma_{SFR} (t_\star = %.1f)\ [M_\odot\ yr^{-1}\ kpc^{-2}$]' % (H.tSF__T[iT] / 1e6)
        zlabel = r'$\log\ x_Y$ [ frac ]'
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xran = [-2, 0.5]
        # yran = [-3, 0]
        # zran = [-3, 0]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        xran = [xm.compressed().min(), xm.compressed().max()]
        yran = [ym.compressed().min(), ym.compressed().max()]
        zran = [zm.compressed().min(), zm.compressed().max()]
        hxyzw, xyzwxedges, xyzwyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), weights = 10 ** zm.compressed(), range = [xran, yran])
        hxy, xyxedges, xyyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), range = [xran, yran])
        hxz, xzxedges, xzyedges = np.histogram2d(xm.compressed(), zm.compressed(), bins = (50, 50), range = [xran, zran])
        hyz, yzxedges, yzyedges = np.histogram2d(ym.compressed(), zm.compressed(), bins = (50, 50), range = [yran, zran])
        #colors=['black','green', 'blue','red']
        colors = ['black', 'black', 'black', 'black']
        #cmap = 'viridis'
        #cmap = 'cubehelix'
        #cmap = 'Spectral'
        cmap = 'Blues'
        linewidths = (2, 1.66, 1.34, 1)
    
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        X, Y = np.meshgrid(xyxedges, xyyedges)
        im = ax.pcolormesh(X, Y, hxy.T, cmap = cmap)
        #ax.set_aspect('equal')
        levels = (500, 100, 25, 5)
        extent = xran + yran
        cset = ax.contour(hxy.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
        ax.set_xlim(xran)
        ax.set_ylim(yran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        X, Y = np.meshgrid(xzxedges, xzyedges)
        im = ax.pcolormesh(X, Y, hxz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = xran + zran
        cset = ax.contour(hxz.T, levels, origin = 'lower', colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, zm, **ols_kwargs)
        ax.axhline(y = np.log10(0.03), c = 'r', ls = '--')
        ax.text(-2.7, np.log10(0.03), '0.03', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.05), c = 'r', ls = '--')
        ax.text(-2.7, np.log10(0.05), '0.05', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.15), c = 'r', ls = '--')
        ax.text(-2.7, np.log10(0.15), '0.15', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.set_xlim(xran)
        ax.set_ylim(zran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        X, Y = np.meshgrid(yzxedges, yzyedges)
        im = ax.pcolormesh(X, Y, hyz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = yran + zran
        cset = ax.contour(hyz.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(ym.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, ym, zm, **ols_kwargs)
        ax.axhline(y = np.log10(0.03), c = 'r', ls = '--')
        ax.text(-5, np.log10(0.03), '0.03', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.05), c = 'r', ls = '--')
        ax.text(-5, np.log10(0.05), '0.05', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.15), c = 'r', ls = '--')
        ax.text(-5, np.log10(0.15), '0.15', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.set_xlim(yran)
        ax.set_ylim(zran)
        ax.set_xlabel(ylabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
    
        mask_tmp = np.bitwise_or(mask__g, H.x_Y__Tg[iT] < 0.03)
        xm, ym, zm = C.ma_mask_xyz(x = np.ma.log10(H.tau_V__Tg[iT]), y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), z = np.ma.log10(H.x_Y__Tg[iT]), mask = mask_tmp)
        xlabel = r'$\log\ \tau_V$'
        ylabel = r'$\log\ \Sigma_{SFR} (t_\star = %.1f)\ [M_\odot\ yr^{-1}\ kpc^{-2}$]' % (H.tSF__T[iT] / 1e6)
        zlabel = r'$\log\ x_Y$ [ frac ]'
        xran = [-4, 0.5]
        yran = [-3.5, 1]
        zran = [np.log10(0.03), 0]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xran = [xm.compressed().min(), xm.compressed().max()]
        # yran = [ym.compressed().min(), ym.compressed().max()]
        # zran = [zm.compressed().min(), zm.compressed().max()]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        hxyzw, xyzwxedges, xyzwyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), weights = 10 ** zm.compressed(), range = [xran, yran])
        hxy, xyxedges, xyyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), range = [xran, yran])
        hxz, xzxedges, xzyedges = np.histogram2d(xm.compressed(), zm.compressed(), bins = (50, 50), range = [xran, zran])
        hyz, yzxedges, yzyedges = np.histogram2d(ym.compressed(), zm.compressed(), bins = (50, 50), range = [yran, zran])
        colors = ['black', 'black', 'black', 'black']
        cmap = 'Blues'
        linewidths = (2, 1.66, 1.34, 1)
        ols_kwargs.update(dict(pos_x = 0.98, pos_y = 0.93))
    
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        X, Y = np.meshgrid(xyxedges, xyyedges)
        im = ax.pcolormesh(X, Y, hxy.T, cmap = cmap)
        #ax.set_aspect('equal')
        levels = (500, 100, 25, 5)
        extent = xran + yran
        cset = ax.contour(hxy.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
        ax.set_xlim(xran)
        ax.set_ylim(yran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        X, Y = np.meshgrid(xzxedges, xzyedges)
        im = ax.pcolormesh(X, Y, hxz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = xran + zran
        cset = ax.contour(hxz.T, levels, origin = 'lower', colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, zm, **ols_kwargs)
        ax.set_xlim(xran)
        ax.set_ylim(zran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        X, Y = np.meshgrid(yzxedges, yzyedges)
        im = ax.pcolormesh(X, Y, hyz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = yran + zran
        cset = ax.contour(hyz.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(ym.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, ym, zm, **ols_kwargs)
        ax.set_xlim(yran)
        ax.set_ylim(zran)
        ax.set_xlabel(ylabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
    
        f.subplots_adjust(bottom = 0.15, hspace = 0.25, wspace = 0.25, right = 0.95, left = 0.07)
        f.savefig('sample.pdf')
        plt.close(f)
        pdf.savefig(f)
        
        #############################
        #############################
        #############################
        
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
        xm, ym, zm = C.ma_mask_xyz(x = np.ma.log10(H.tau_V_neb__g), y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), z = np.ma.log10(H.x_Y__Tg[iT]), mask = mask__g)
        #xm, ym, zm = C.ma_mask_xyz(x = np.ma.log10(H.tau_V_neb__g), y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), z = np.ma.log10(H.x_Y__Tg[iT]))
        xlabel = r'$\log\ \tau_V^{neb}$'
        ylabel = r'$\log\ \Sigma_{SFR} (t_\star = %.1f)\ [M_\odot\ yr^{-1}\ kpc^{-2}$]' % (H.tSF__T[iT] / 1e6)
        zlabel = r'$\log\ x_Y$ [ frac ]'
        xran = [xm.compressed().min(), xm.compressed().max()]
        yran = [ym.compressed().min(), ym.compressed().max()]
        zran = [zm.compressed().min(), zm.compressed().max()]
        hxyzw, xyzwxedges, xyzwyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), weights = 10 ** zm.compressed(), range = [xran, yran])
        hxy, xyxedges, xyyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), range = [xran, yran])
        hxz, xzxedges, xzyedges = np.histogram2d(xm.compressed(), zm.compressed(), bins = (50, 50), range = [xran, zran])
        hyz, yzxedges, yzyedges = np.histogram2d(ym.compressed(), zm.compressed(), bins = (50, 50), range = [yran, zran])
        #colors=['black','green', 'blue','red']
        colors = ['black', 'black', 'black', 'black']
        #cmap = 'viridis'
        #cmap = 'cubehelix'
        #cmap = 'Spectral'
        cmap = 'Blues'
        linewidths = (2, 1.66, 1.34, 1)
    
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        X, Y = np.meshgrid(xyxedges, xyyedges)
        im = ax.pcolormesh(X, Y, hxy.T, cmap = cmap)
        #ax.set_aspect('equal')
        levels = (500, 100, 25, 5)
        extent = xran + yran
        cset = ax.contour(hxy.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
        ax.set_xlim(xran)
        ax.set_ylim(yran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        X, Y = np.meshgrid(xzxedges, xzyedges)
        im = ax.pcolormesh(X, Y, hxz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = xran + zran
        cset = ax.contour(hxz.T, levels, origin = 'lower', colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, zm, **ols_kwargs)
        ax.axhline(y = np.log10(0.03), c = 'r', ls = '--')
        ax.text(-2.7, np.log10(0.03), '0.03', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.05), c = 'r', ls = '--')
        ax.text(-2.7, np.log10(0.05), '0.05', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.15), c = 'r', ls = '--')
        ax.text(-2.7, np.log10(0.15), '0.15', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.set_xlim(xran)
        ax.set_ylim(zran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        X, Y = np.meshgrid(yzxedges, yzyedges)
        im = ax.pcolormesh(X, Y, hyz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = yran + zran
        cset = ax.contour(hyz.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(ym.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, ym, zm, **ols_kwargs)
        ax.axhline(y = np.log10(0.03), c = 'r', ls = '--')
        ax.text(-5, np.log10(0.03), '0.03', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.05), c = 'r', ls = '--')
        ax.text(-5, np.log10(0.05), '0.05', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.axhline(y = np.log10(0.15), c = 'r', ls = '--')
        ax.text(-5, np.log10(0.15), '0.15', fontsize = 6, color = 'k', bbox = dict(facecolor = 'white', edgecolor = 'none'))
        ax.set_xlim(yran)
        ax.set_ylim(zran)
        ax.set_xlabel(ylabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
    
        mask_tmp = np.bitwise_or(mask__g, H.x_Y__Tg[iT] < 0.03)
        xm, ym, zm = C.ma_mask_xyz(x = np.ma.log10(H.tau_V_neb__g), y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6), z = np.ma.log10(H.x_Y__Tg[iT]), mask = mask_tmp)
        xlabel = r'$\log\ \tau_V^{neb}$'
        ylabel = r'$\log\ \Sigma_{SFR} (t_\star = %.1f)\ [M_\odot\ yr^{-1}\ kpc^{-2}$]' % (H.tSF__T[iT] / 1e6)
        zlabel = r'$\log\ x_Y$ [ frac ]'
        xran = [-3, 0.5]
        yran = [-3.5, 1]
        zran = [np.log10(0.03), 0]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        # xran = [xm.compressed().min(), xm.compressed().max()]
        # yran = [ym.compressed().min(), ym.compressed().max()]
        # zran = [zm.compressed().min(), zm.compressed().max()]
        #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        hxyzw, xyzwxedges, xyzwyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), weights = 10 ** zm.compressed(), range = [xran, yran])
        hxy, xyxedges, xyyedges = np.histogram2d(xm.compressed(), ym.compressed(), bins = (50, 50), range = [xran, yran])
        hxz, xzxedges, xzyedges = np.histogram2d(xm.compressed(), zm.compressed(), bins = (50, 50), range = [xran, zran])
        hyz, yzxedges, yzyedges = np.histogram2d(ym.compressed(), zm.compressed(), bins = (50, 50), range = [yran, zran])
        colors = ['black', 'black', 'black', 'black']
        cmap = 'Blues'
        linewidths = (2, 1.66, 1.34, 1)
        ols_kwargs.update(dict(pos_x = 0.98, pos_y = 0.93))
    
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        X, Y = np.meshgrid(xyxedges, xyyedges)
        im = ax.pcolormesh(X, Y, hxy.T, cmap = cmap)
        #ax.set_aspect('equal')
        levels = (500, 100, 25, 5)
        extent = xran + yran
        cset = ax.contour(hxy.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, ym, **ols_kwargs)
        ax.set_xlim(xran)
        ax.set_ylim(yran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        X, Y = np.meshgrid(xzxedges, xzyedges)
        im = ax.pcolormesh(X, Y, hxz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = xran + zran
        cset = ax.contour(hxz.T, levels, origin = 'lower', colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(xm.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, xm, zm, **ols_kwargs)
        ax.set_xlim(xran)
        ax.set_ylim(zran)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
        
        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        X, Y = np.meshgrid(yzxedges, yzyedges)
        im = ax.pcolormesh(X, Y, hyz.T, cmap = cmap)
        levels = (500, 100, 25, 5)
        extent = yran + zran
        cset = ax.contour(hyz.T, levels, origin = 'lower' , colors = colors, linewidths = linewidths, extent = extent)
        for c in cset.collections:
            c.set_linestyle('solid')
        ax.clabel(cset, inline = 1, fontsize = 10, fmt = '%.0f')
        rs = C.runstats(ym.compressed(), zm.compressed(), **rs_kwargs)
        ax.plot(rs.xMedian, rs.yMedian, 'k--', lw = 2)
        ax.plot(rs.inv_yMedian, rs.inv_xMedian, 'k--', lw = 2)
        plotOLSbisectorAxis(ax, ym, zm, **ols_kwargs)
        ax.set_xlim(yran)
        ax.set_ylim(zran)
        ax.set_xlabel(ylabel)
        ax.set_ylabel(zlabel)
        f.colorbar(mappable = im, ax = ax)
    
        f.subplots_adjust(bottom = 0.15, hspace = 0.25, wspace = 0.25, right = 0.95, left = 0.07)
        plt.close(f)
        pdf.savefig(f)
        
        
