#!/usr/bin/python
#
# Lacerda@Corrego - 16/Nov/2015
#
from matplotlib.backends.backend_pdf import PdfPages
from CALIFAUtils.plots import plot_text_ax
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt
import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

#RNuc = 0.5

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# mpl.rcParams['font.size'] = 12
# mpl.rcParams['axes.labelsize'] = 10
# mpl.rcParams['axes.titlesize'] = 12
# mpl.rcParams['xtick.labelsize'] = 10
# mpl.rcParams['ytick.labelsize'] = 10
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

mpl.rcParams['font.size'] = 20
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
A4Size_inches = [ 8.267, 11.692 ]
LetterSize_inches = [ 8.5, 11 ]

default_sc_kwargs = dict(marker = 'o', s = 10, alpha = 0.9, edgecolor = 'none', cmap = 'Spectral')
default_rs_kwargs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30)
default_ols_kwargs = dict(c = 'k', pos_x = 0.98, pos_y = 0.01, fs = 12, rms = True, text = True)
default_ols_plot_kwargs = dict(c = 'r', ls = '--', lw = 2, label = 'OLS')

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
        'plot_SKmet' : False,
        'plot_closedbox' : False,
        'plot_SKdev' : False,
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
    parser.add_argument('--plot_SKmet',
                        action = 'store_true',
                        default = default['plot_SKmet'])
    parser.add_argument('--plot_SKdev',
                        action = 'store_true',
                        default = default['plot_SKdev'])
    parser.add_argument('--plot_closedbox',
                        action = 'store_true',
                        default = default['plot_closedbox'])


    return parser.parse_args()

def givethatrealmtfkingmorf(notrealmtfkingmorf = None):
    mtype = {
        'E0' : 0,
        'E1' : 0,
        'E2' : 0,
        'E3' : 0,
        'E4' : 0,
        'E5' : 0,
        'E6' : 0,
        'E7' : 0,
        'S0' : 1,
        'S0a' : 1,
        'Sa' : 2,
        'Sab' : 2,
        'Sb' : 3,
        'Sbc' : 4,
        'Sc' : 5,
        'Scd' : 6,
        'Sd' : 6,
        'Sdm' : -1,
        'Sm' : -1,
        'Ir' : -1,
    }
    realmtfkingmof = mtype[notrealmtfkingmorf]
    return realmtfkingmof   
    
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# def plot_radial_SigmaSF(f, H, mask__g, mask__rg, mask_GAL__g):
#     sc_kwargs = default_sc_kwargs.copy()
#     rs_kwargs = default_rs_kwargs.copy()
#     ols_kwargs = default_ols_kwargs.copy()
#     ols_plot_kwargs = default_ols_plot_kwargs.copy()
#     NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
#     NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
#     #suptitle_R = r'NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
#     NRows = 2
#     NCols = 8
#     page_size_inches = [12, 8]
#     f.set_size_inches(page_size_inches)
#     grid_shape = (NRows, NCols)
#     #f.suptitle(suptitle_R, fontsize = 15)
# 
#     colortipo = ['brown', 'red', 'orange', 'green', '#00D0C9', '#0076C9', 'blue']
#     Ntype = len(colortipo)
#     mtypes = [ 0, 1, 2, 3, 4, 5, 6 ]
#     mtype_labels = [ 'E', 'S0', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd' ]
#     halfbinstep = np.diff(mtypes)[0]/2.
#     tickpos = np.linspace(mtypes[0] + halfbinstep[0], mtypes[-1] - halfbinstep[0], Ntype)
#      
#     cmap = mpl.colors.ListedColormap(colortipo)
# 
#     sc = ax.scatter(x, y, c = z_mtype, cmap = cmap, marker = 'o', s = 10, edgecolor = 'none', label = '')
# 
#     cb = f.colorbar(sc, ticks = tickpos)
#     cb.ax.set_yticklabels(mtype_labels)
#     cb.ax.yaxis.set_ticks_position('none')    
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

def plot_gas_profile_zones(f, H, gas_dict, mask__g, mask__rg, mask_GAL__g):
    #gasnames = [ 'GRVTauVStar', 'GRVTauVNeb', 'BRTauVStar', 'BRTauVNeb', 'SKSFRSDStar', 'SKSFRSDHa' ]
    gasnames = [ 'GRVTauVStar', 'GRVTauVNeb', 'BRTauVStar', 'BRTauVNeb']
    rs_kwargs = default_rs_kwargs.copy()
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    suptitle_R = r'NZones:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%% ' % ((~mask__g).sum(), NgalsOkZones, (H.tSF__T[iT] / 1e6), H.xOkMin * 100.)
    NRows = 2
    NCols = len(gasnames)
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [11, 8]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 15)
    
    col = 0
    for gasname in gasnames:
        gas = gas_dict[gasname]
        
        ax = plt.subplot2grid(grid_shape, loc = (0, col))
        ylim = [0, 2]
        xlim = [minR, 3]
        plot_text_ax(ax, '%s' % gasname , 0.95, 0.95, 12, 'top', 'right', 'white')
        xm, ym = C.ma_mask_xyz(H.zone_dist_HLR__g, gas['SigmaGas__g'], mask = mask__g)
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xS, np.ma.log10(rs.yS), c = 'white', label = gasname)
        #ax.plot(H.zone_dist_HLR__g, np.ma.log10(gas['SigmaGas__g']), c = 'white', ls = '--')
        xm, ym = C.ma_mask_xyz(xm, np.ma.log10(ym))
        ax.hist2d(xm.compressed(), ym.compressed(), bins = 100, cmap = 'viridis', range = [xlim, ylim])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([minR, 1, 2, 3])
        plt.setp(ax.get_xticklabels(), visible = False)
        if col == 0:
            ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [$M_\odot pc^{-2}$]')
        else:
            plt.setp(ax.get_yticklabels(), visible = False)
        ax = plt.subplot2grid(grid_shape, loc = (1, col))
        ylim = [0, 0.3]
        xlim = [minR, 3]
        plot_text_ax(ax, '%s' % gasname , 0.95, 0.95, 12, 'top', 'right', 'white')
        xm, ym = C.ma_mask_xyz(H.zone_dist_HLR__g, gas['f_gas__g'], mask = mask__g)
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xS, rs.yS, c = 'white', label = gasname)
        #ax.plot(H.zone_dist_HLR__g, gas['f_gas__g'], c = 'white', ls = '--')
        ax.hist2d(xm.compressed(), ym.compressed(), bins = 100, cmap = 'viridis', range = [xlim, ylim])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([minR, 1, 2, 3])
        if col == 0:
            ax.set_xlabel(r'R [HLR]')
            ax.set_ylabel(r'$f_{gas}$')
        else:
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
        col += 1
        
    f.subplots_adjust(bottom = 0.15, hspace = 0.15, wspace = 0., right = 0.95, left = 0.10)

def plot_gas_profile(f, H, gas_dict, mask__g, mask__rg, mask_GAL__g):
    #gasnames = [ 'GRVTauVStar', 'GRVTauVNeb', 'BRTauVStar', 'BRTauVNeb', 'SKSFRSDStar', 'SKSFRSDHa' ]
    gasnames = [ 'GRVTauVStar', 'GRVTauVNeb', 'BRTauVStar', 'BRTauVNeb']
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%' % ((~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100.)
    NRows = 2
    NCols = len(gasnames)
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [11, 8]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 15)
    
    col = 0
    for gasname in gasnames:
        gas = gas_dict[gasname]
        
        # SigmaGas radial mean
        SigmaGas__rg = gas['SigmaGas__rg']
        f_gas__rg = gas['f_gas__rg']
        morf_mask = np.bitwise_or(np.less(H.morfType_GAL__g, 9), np.greater(H.morfType_GAL__g, 12))
        aSigmaGas__r = np.zeros((H.NRbins), dtype = np.float_)
        af_gas__r = np.zeros((H.NRbins), dtype = np.float_)
        for ir in xrange(H.NRbins):
            test = False
            #test = H.RbinCenter__r[ir] > (2 + H.RbinStep)
            if test:
                SigmaGas_aux = SigmaGas__rg[ir, ~morf_mask]
                f_gas_aux = f_gas__rg[ir, ~morf_mask]
            else:
                SigmaGas_aux = np.where(SigmaGas__rg.mask[ir, ~morf_mask], 0., SigmaGas__rg[ir, ~morf_mask])
                f_gas_aux = np.where(f_gas__rg.mask[ir, ~morf_mask], 0., f_gas__rg[ir, ~morf_mask])
            aSigmaGas__r[ir] = SigmaGas_aux.mean()
            af_gas__r[ir] = f_gas_aux.mean()
            
        ax = plt.subplot2grid(grid_shape, loc = (0, col))
        ylim = [0, 2]
        xlim = [minR, 3]
        plot_text_ax(ax, '%s' % gasname , 0.95, 0.95, 12, 'top', 'right', 'white')
        xm, ym = C.ma_mask_xyz(H.Rtoplot(), SigmaGas__rg, mask = mask__rg)
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xS, np.ma.log10(rs.yS), c = 'white', label = gasname)
        ax.plot(H.RbinCenter__r, np.ma.log10(aSigmaGas__r), c = 'white', ls = '--')
        xm, ym = C.ma_mask_xyz(xm, np.ma.log10(ym))
        ax.hist2d(xm.compressed(), ym.compressed(), bins = 100, cmap = 'viridis', range = [[minR, 3], [0, 2]])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([minR, 1, 2, 3])
        plt.setp(ax.get_xticklabels(), visible = False)
        if col == 0:
            ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [$M_\odot pc^{-2}$]')
        else:
            plt.setp(ax.get_yticklabels(), visible = False)
        ax = plt.subplot2grid(grid_shape, loc = (1, col))
        ylim = [0, 0.3]
        xlim = [minR, 3]
        plot_text_ax(ax, '%s' % gasname , 0.95, 0.95, 12, 'top', 'right', 'white')
        xm, ym = C.ma_mask_xyz(H.Rtoplot(), f_gas__rg, mask = mask__rg)
        rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
        ax.plot(rs.xS, rs.yS, c = 'white', label = gasname)
        ax.plot(H.RbinCenter__r, af_gas__r, c = 'white', ls = '--')
        ax.hist2d(xm.compressed(), ym.compressed(), bins = 100, cmap = 'viridis', range = [xlim, ylim])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([minR, 1, 2, 3])
        if col == 0:
            ax.set_xlabel(r'R [HLR]')
            ax.set_ylabel(r'$f_{gas}$')
        else:
            plt.setp(ax.get_xticklabels(), visible = False)
            plt.setp(ax.get_yticklabels(), visible = False)
        col += 1
        
    f.subplots_adjust(bottom = 0.15, hspace = 0.15, wspace = 0., right = 0.95, left = 0.10)

def plot_depletion(f, H, gas_dict, gasname, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    gas = gas_dict[gasname]
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'%s NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (gasname, (~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    NRows = 1
    NCols = 2
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [10, 6]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 15)
    
    tdep = gas['SigmaGas__rg'] / (H.aSFRSD__Trg[iT] * 1e6)

    input_kw = dict(
                f = f,
                debug = True,
                #xlim = (-10, 2),
                #ylim = (8., 8.8),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 4,
                x_minor_locator = .8,
                y_major_locator = .2,
                y_minor_locator = .05,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
    )
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = H.at_flux__Trg[iT],
                    y = np.ma.log10(tdep),
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'$\langle \log\ t \rangle_L (R)$ [yr]',
                    ylabel = r'$t_{dep} [yr]$',
                    **input_kw
    )
    ####
    tdep = gas['SigmaGas__rg'] / (H.aSFRSD_Ha__rg * 1e6)

    input_kw = dict(
                f = f,
                debug = True,
                #xlim = (-10, 2),
                #ylim = (8., 8.8),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 4,
                x_minor_locator = .8,
                y_major_locator = .2,
                y_minor_locator = .05,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
    )
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = H.at_flux__Trg[iT],
                    y = np.ma.log10(tdep),
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'$\langle \log\ t \rangle_L (R)$ [yr]',
                    ylabel = r'$t_{dep} [yr]$',
                    **input_kw
    )
    
    f.subplots_adjust(bottom = 0.2, hspace = 0.4)    

def plot_cover(f, H, gas_dict, gasname, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    gas = gas_dict[gasname]
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'%s NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (gasname, (~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [7, 6]
    f.suptitle(suptitle_R, fontsize = 15)

    input_kw = dict(
                f = f,
                debug = True,
                xlim = (-7, 0),
                ylim = (8., 8.8),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 2,
                x_minor_locator = .4,
                y_major_locator = .2,
                y_minor_locator = .05,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
    )

    ### RBINS ln F_GAS vs ZNeb M13 ###
    ax = f.gca()        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['f_gas__rg']),
                    y = H.O_O3N2_M13__rg,
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'$\ln\ f_{gas}$ (zones) ',
                    ylabel = r'$12 + \log(O/H)$ (zones)',
                    **input_kw
    )
    ax = r_kw['ax']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####
    
def plot_Z_Mratio(f, H, gas_dict, gasname, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    gas = gas_dict[gasname]
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'%s NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (gasname, (~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    NRows = 1
    NCols = 2
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [10, 6]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 15)
    
    input_kw = dict(
                f = f,
                debug = True,
                ylim = (8.0, 8.8),
                #xlim = (0, 100),
                zlim = (minR, 3),
                x_major_locator = 1,
                x_minor_locator = .2,
                y_major_locator = 0.2,
                y_minor_locator = .05,
                kwargs_scatter = sc_kwargs,
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                running_stats = True,
                kwargs_rs = rs_kwargs,
                #ols = True,
                #kwargs_ols = ols_kwargs,
                #kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
                #write_N = True,
                #kwargs_zbins_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30),
                spearmanr = True,
                #legend = True,
    )
    
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    r_kw = C.plot_zbins(
                    ax = f.gca(),
                    x = np.ma.log10(H.McorSD__Trg[iT]/gas['SigmaGas__rg']),
                    y = H.O_O3N2_M13__rg,
                    z = H.Rtoplot(),
                    xlabel = r'$\frac{\mu_\star}{\Sigma_{gas}}$',
                    ylabel = r'12 + $\log$ O/H',
                    zlabel = r'R [HLR]',
                    add_mask = mask__rg,
                    colorbar = False,
                    **input_kw
    )
    rs = r_kw['rs']    
    mask = rs.xS > 1.4
    yieldNeb = rs.yS[mask].min()
    

    input_kw['ylim'] = (-1.8, 0.4)
    input_kw['y_major_locator'] = 1.
    input_kw['y_minor_locator'] = .2

    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    r_kw = C.plot_zbins(
                    ax = f.gca(),
                    x = np.ma.log10(H.McorSD__Trg[iT]/gas['SigmaGas__rg']),
                    y = H.alogZ_mass__Urg[-1],
                    z = H.Rtoplot(),
                    xlabel = r'$\frac{\mu_\star}{\Sigma_{gas}}$',
                    ylabel = r'$\langle \log\ Z_\star \rangle_M (R)$ (t < %.2f Gyr) [$Z_\odot$]' % (H.tZ__U[iU] / 1e9),
                    zlabel = r'R [HLR]',
                    add_mask = mask__rg,
                    colorbar = True,
                    **input_kw
    )    
    rs = r_kw['rs']    
    mask = rs.xS > 1.4
    yieldStar = rs.yS[mask].min()
    print gasname, yieldNeb, yieldStar
    f.subplots_adjust(bottom = 0.2, hspace = 0.4)

def plot_SKGas(f, H, gas_dict, gasname, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    gas = gas_dict[gasname]
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'%s NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (gasname, (~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [8, 6]
    f.set_size_inches(page_size_inches)
    #grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 10)

    input_kw = dict(
                f = f,
                debug = True,
                xlim = (-1, 3),
                ylim = (-3, 0.5),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 1,
                x_minor_locator = .2,
                y_major_locator = 1,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
                spearmanr = True, 
    )

    #ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    ax = f.gca()        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['SigmaGas__rg']),
                    y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'$\log\ \Sigma_{gas}\ [M_\odot pc^{-2}]$',
                    ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$',
                    colorbar = True, 
                    **input_kw
    )
    f.subplots_adjust(bottom = 0.2)
    #ax = r_kw['ax']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)

def plot_Z(f, H, iT, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    #page_size_inches = [8, 6]
    #f.set_size_inches(page_size_inches)
    f.suptitle(suptitle_R, fontsize = 10)
    
    input_kw = dict(
                f = f,
                debug = True,
                ylim = (8.0, 8.8),
                xlim = (-1.8, 0.4),
                zlim = (minR, 3),
                y_major_locator = 0.2,
                y_minor_locator = .05,
                x_major_locator = 1,
                x_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
                #write_N = True,
                #kwargs_zbins_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30),
                spearmanr = True,
                #legend = True,
    )

    r_kw = C.plot_zbins(
                    ax = f.gca(),
                    x = H.alogZ_mass__Urg[-1],
                    y = H.O_O3N2_M13__rg,
                    z = H.Rtoplot(),
                    xlabel = r'$\langle \log\ Z_\star \rangle_M (R)$ [$Z_\odot$]',
                    ylabel = r'12 + $\log$ O/H',
                    zlabel = r'R [HLR]',
                    add_mask = mask__rg,
                    colorbar = True,
                    **input_kw
    )
    ax = r_kw['ax']
    f = r_kw['f']
    ax.plot(ax.get_xlim(), np.asarray(ax.get_xlim()) - 8.69, 'k--', lw = 1)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    #f.subplots_adjust(bottom = 0.15, right = 0.95)
    
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # for i in range(H.N_U):
    #     xm, ym = C.ma_mask_xyz(H.O_O3N2_M13__rg, H.alogZ_mass__Urg[i], mask = mask__rg)
    #     rs = C.runstats(xm.compressed(), ym.compressed(), **rs_kwargs)
    #     ax.plot(rs.xS, rs.yS, '--k', lw = 3)#,  label = r'($t_Z$ < %.2f Gyr)' % (self.tZ__U[iU] / 1e9))
    #     ax.text(rs.xS[0] - 0.03 * (xmax - xmin), rs.yS[0], '%.2f Gyr' % (H.tZ__U[i] / 1e9), fontsize = 8, color = 'k',
    #             #transform = ax.transAxes,
    #             verticalalignment = 'center', horizontalalignment = 'right',
    #             bbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.))
    #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
         
def plot_SFR(f, H, iT, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    NRows = 1
    NCols = 2
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [10, 6]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 15)
    
    input_kw = dict(
                f = f,
                debug = True,
                xlim = (-3.5, 1),
                ylim = (-3.5, 1),
                zlim = (minR, 3),
                x_major_locator = 1,
                x_minor_locator = .2,
                y_major_locator = 1,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
                #write_N = True,
                #kwargs_zbins_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30),
                spearmanr = True,
                #legend = True,
    )

    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    input_kw.update        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                    y = np.ma.log10(H.aSFRSD_Ha__rg * 1e6),
                    z = H.Rtoplot(),
                    xlabel = r'$\log\ \Sigma_{SFR}^\star(t_\star, R)\ [M_\odot yr^{-1} kpc^{-2}]$',
                    ylabel = r'$\log\ \Sigma_{SFR}^{neb}(R)\ [M_\odot yr^{-1} kpc^{-2}]$',
                    zlabel = r'R [HLR]',
                    add_mask = mask__rg,
                    colorbar = False,
                    **input_kw
    )
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--', lw = 1)

    input_kw['xlim'] = (-2, 2)
    input_kw['ylim'] = (-2, 2)
    
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    input_kw.update        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.aSFRSD__Trg[iT] / H.aSFRSD_oneHLR__Tg[iT]),
                    y = np.ma.log10(H.aSFRSD_Ha__rg / H.aSFRSD_Ha_oneHLR__g),
                    z = H.Rtoplot(),
                    xlabel = r'$\log\ \frac{\Sigma_{SFR}^\star}{\Sigma_{SFR}^\star(@1HLR)}$',
                    ylabel = r'$\log\ \frac{\Sigma_{SFR}^{neb}}{\Sigma_{SFR}^{neb}(@1HLR)}$',
                    zlabel = r'R [HLR]',
                    add_mask = mask__rg,
                    **input_kw
    )
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k--', lw = 1)
    
    f.subplots_adjust(bottom = 0.15, hspace = 0.4, wspace = 0.45, right = 0.95, left = 0.07)
    
def plot_tauV(f, H, iT, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    NRows = 1
    NCols = 2
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [10, 6]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 15)
    
    input_kw = dict(
                f = f,
                debug = True,
                xlim = (0, 1.5),
                ylim = (0, 3),
                zlim = (minR, 3),
                x_major_locator = .5,
                x_minor_locator = .1,
                y_major_locator = 1,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
                #write_N = True,
                #kwargs_zbins_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30),
                #spearmanr = True,
                #legend = True,
    )

    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = H.tau_V__Trg[iT],
                    y = H.tau_V_neb__rg,
                    z = H.Rtoplot(),
                    xlabel = r'$\tau_V^\star$',
                    ylabel = r'$\tau_V^{neb}$',
                    zlabel = r'R [HLR]',
                    add_mask = mask__rg,
                    **input_kw
    )

    input_kw = dict(
                f = f,
                debug = True,
                xlim = (0, 50),
                ylim = (0, 6),
                zlim = (minR, 3),
                x_major_locator = 10,
                x_minor_locator = 2,
                y_major_locator = 1,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                #ols = True,
                #kwargs_ols = ols_kwargs,
                #kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
                #write_N = True,
                running_stats = True,
                kwargs_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30),
                #spearmanr = True,
                #legend = True,
    )
    
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = H.x_Y__Trg[iT] * 100.,
                    y = H.tau_V_neb__rg / H.tau_V__Trg[iT],
                    z = H.Rtoplot(),
                    xlabel = r'$x_Y$ [%]',
                    ylabel = r'$\frac{\tau_V^{neb}}{\tau_V^\star}$',
                    zlabel = r'R [HLR]',
                    add_mask = mask__rg,
                    **input_kw
    )
    
    f.subplots_adjust(bottom = 0.15, hspace = 0.4, wspace = 0.45, right = 0.95, left = 0.07)

def plot_SK_dev(f, H, iT, z_dict, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'Nzones:%d(%d gals) NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((~mask__g).sum(), NgalsOkZones, (~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    NRows = 2
    NCols = 3
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [11, 8]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 10)

    z = z_dict['v']

    
    if z_dict['label'] == 'morph. type':
        mtypes_colors = ['brown', 'red', 'orange', 'green', '#00D0C9', '#0076C9', 'blue']
        mtypes_labels = [ 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd' ]
        mtypes = [ 9., 9.5, 10., 10.5, 11. , 11.5, 12. ] 
        mtypes_num = len(mtypes)
        halfbinstep = np.diff(mtypes)[0]/2.,
        tickpos = np.linspace(mtypes[0] + halfbinstep[0], mtypes[-1] - halfbinstep[0], mtypes_num)        

        input_kw = dict(
                    f = f,
                    debug = True,
                    xlim = (-1.5, 1),
                    ylim = (-3.5, 1),
                    #zlim = (0, 4),
                    x_major_locator = 1,
                    x_minor_locator = .2,
                    y_major_locator = 1,
                    y_minor_locator = .2,
                    kwargs_scatter = dict(marker = 'o', s = 10, alpha = 0.9, edgecolor = 'none', cmap = mpl.colors.ListedColormap(mtypes_colors)),
                    ols = True,
                    kwargs_ols = ols_kwargs,
                    kwargs_ols_plot = ols_plot_kwargs,
                    return_kwargs = True,
                    kwargs_sc_cb = dict(ticks = tickpos),
                    zticklabels = mtypes_labels, 
        )

        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        r_kw = C.plot_zbins(
                        ax = ax,
                        x = np.ma.log10(H.tau_V__Trg[iT]),
                        y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                        z = z,
                        xlabel = r'$\log\ \tau_V^\star$ (Rbins)',
                        ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$ (Rbins)',
                        zlabel = z_dict['label'],
                        add_mask = mask__rg,
                        **input_kw
        )
        ax = r_kw['ax']
        f = r_kw['f']
        xm = r_kw['xm']
        ym = r_kw['ym']
        zm = r_kw['zm']
        p_ols = [ r_kw['ols'][0], r_kw['ols'][1] ]    
        p_polyfit = np.polyfit(xm.compressed(), ym.compressed(), 1)
        ax.plot(xm.compressed(), np.polyval(p_polyfit, xm.compressed()), 'k:') 
        #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        R = ym - np.polyval(p_polyfit, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.scatter(propm.compressed(), Rm.compressed(), **sc_kwargs)
        for mt in np.unique(propm.compressed()):
            ax.plot(mt, np.median(R[propm == mt]), 'k*', lw = 10)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('polyfit o=1 rms: %.4f' % R.std())
        ax.set_xticks(mtypes)
        ax.set_xticklabels(mtypes_labels)
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        R = ym - np.polyval(p_ols, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.scatter(propm.compressed(), Rm.compressed(), **sc_kwargs)
        for mt in np.unique(propm.compressed()):
            ax.plot(mt, np.median(R[propm == mt]), 'k*', lw = 10) 
        ax.set_xticks(mtypes)
        ax.set_xticklabels(mtypes_labels)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('OLS bisect rms: %.4f' % R.std())
        ########################################
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        r_kw = C.plot_zbins(
                        ax = ax,
                        x = np.ma.log10(H.tau_V_neb__rg),
                        y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                        z = z,
                        xlabel = r'$\log\ \tau_V^{neb}$ (Rbins)',
                        ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$ (Rbins)',
                        zlabel = z_dict['label'],
                        add_mask = mask__rg,
                        **input_kw
        )
        ax = r_kw['ax']
        f = r_kw['f']
        xm = r_kw['xm']
        ym = r_kw['ym']
        zm = r_kw['zm']
        p_polyfit = np.polyfit(xm.compressed(), ym.compressed(), 1)
        p_ols = [ r_kw['ols'][0], r_kw['ols'][1] ]
        ax.plot(xm.compressed(), np.polyval(p_polyfit, xm.compressed()), 'k:') 
        #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        R = ym - np.polyval(p_polyfit, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.scatter(propm.compressed(), Rm.compressed(), **sc_kwargs)
        for mt in np.unique(propm.compressed()):
            ax.plot(mt, np.median(R[propm == mt]), 'k*', lw = 10) 
        ax.set_xticks(mtypes)
        ax.set_xticklabels(mtypes_labels)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('polyfit o=1 rms: %.4f' % R.std())
        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        R = ym - np.polyval(p_ols, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.scatter(propm.compressed(), Rm.compressed(), **sc_kwargs)
        for mt in np.unique(propm.compressed()):
            ax.plot(mt, np.median(R[propm == mt]), 'k*', lw = 10) 
        ax.set_xticks(mtypes)
        ax.set_xticklabels(mtypes_labels)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('OLS bisect rms: %.4f' % R.std())
        ########################################
    else:
        input_kw = dict(
                    f = f,
                    debug = True,
                    xlim = (-1.5, 1),
                    ylim = (-3.5, 1),
                    #zlim = (0, 4),
                    x_major_locator = 1,
                    x_minor_locator = .2,
                    y_major_locator = 1,
                    y_minor_locator = .2,
                    kwargs_scatter = sc_kwargs,
                    ols = True,
                    kwargs_ols = ols_kwargs,
                    kwargs_ols_plot = ols_plot_kwargs,
                    return_kwargs = True,
                    #write_N = True,
                    #kwargs_zbins_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30),
                    #spearmanr = True,
                    #legend = True,
        )
        ########################################
        ax = plt.subplot2grid(grid_shape, loc = (0, 0))
        r_kw = C.plot_zbins(
                        ax = ax,
                        x = np.ma.log10(H.tau_V__Trg[iT]),
                        y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                        z = z,
                        xlabel = r'$\log\ \tau_V^\star$ (Rbins)',
                        ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$ (Rbins)',
                        zlabel = z_dict['label'],
                        add_mask = mask__rg,
                        **input_kw
        )
        ax = r_kw['ax']
        f = r_kw['f']
        xm = r_kw['xm']
        ym = r_kw['ym']
        zm = r_kw['zm']
        p_ols = [ r_kw['ols'][0], r_kw['ols'][1] ]    
        p_polyfit = np.polyfit(xm.compressed(), ym.compressed(), 1)
        ax.plot(xm.compressed(), np.polyval(p_polyfit, xm.compressed()), 'k:') 
        #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
        ax = plt.subplot2grid(grid_shape, loc = (0, 1))
        R = ym - np.polyval(p_polyfit, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.hist2d(propm.compressed(), Rm.compressed(), bins=20, cmap = 'Blues')
        rs = C.runstats(propm.compressed(), Rm.compressed(), **rs_kwargs)
        ax.plot(rs.xS, rs.yS, 'k-', lw = 1)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('polyfit o=1 rms: %.4f' % R.std())
        ax = plt.subplot2grid(grid_shape, loc = (0, 2))
        R = ym - np.polyval(p_ols, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.hist2d(propm.compressed(), Rm.compressed(), bins=20, cmap = 'Blues')
        rs = C.runstats(propm.compressed(), Rm.compressed(), **rs_kwargs)
        ax.plot(rs.xS, rs.yS, 'k-', lw = 1)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('OLS bisect rms: %.4f' % R.std())
        ########################################
        ax = plt.subplot2grid(grid_shape, loc = (1, 0))
        r_kw = C.plot_zbins(
                        ax = ax,
                        x = np.ma.log10(H.tau_V_neb__rg),
                        y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                        z = z,
                        xlabel = r'$\log\ \tau_V^{neb}$ (Rbins)',
                        ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$ (Rbins)',
                        zlabel = z_dict['label'],
                        add_mask = mask__rg,
                        **input_kw
        )
        ax = r_kw['ax']
        f = r_kw['f']
        xm = r_kw['xm']
        ym = r_kw['ym']
        zm = r_kw['zm']
        p_polyfit = np.polyfit(xm.compressed(), ym.compressed(), 1)
        p_ols = [ r_kw['ols'][0], r_kw['ols'][1] ]
        ax.plot(xm.compressed(), np.polyval(p_polyfit, xm.compressed()), 'k:') 
        #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
        ax = plt.subplot2grid(grid_shape, loc = (1, 1))
        R = ym - np.polyval(p_polyfit, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.hist2d(propm.compressed(), Rm.compressed(), bins=20, cmap = 'Blues')
        rs = C.runstats(propm.compressed(), Rm.compressed(), **rs_kwargs)
        ax.plot(rs.xS, rs.yS, 'k-', lw = 1)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('polyfit o=1 rms: %.4f' % R.std())
        ax = plt.subplot2grid(grid_shape, loc = (1, 2))
        R = ym - np.polyval(p_ols, xm)
        propm, Rm = C.ma_mask_xyz(z, R)
        ax.hist2d(propm.compressed(), Rm.compressed(), bins=20, cmap = 'Blues')
        rs = C.runstats(propm.compressed(), Rm.compressed(), **rs_kwargs)
        ax.plot(rs.xS, rs.yS, 'k-', lw = 1)
        ax.set_ylim(-1.5,2.)  
        ax.set_title('OLS bisect rms: %.4f' % R.std())
        ########################################

    f.subplots_adjust(bottom = 0.15, hspace = 0.4, wspace = 0.45, right = 0.95, left = 0.07)

def plot_SK_Zneb_bins(f, H, iT, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    tSF = H.tSF__T[iT]
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'Nzones:%d(%d gals) NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % ((~mask__g).sum(), NgalsOkZones, (~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    NRows = 4
    NCols = 4
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [11, 11]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 10)

    input_kw = dict(
                f = f,
                debug = True,
                xlim = (-1.5, 1),
                ylim = (-3.5, 1),
                #zlim = (0, 4),
                x_major_locator = 1,
                x_minor_locator = .2,
                y_major_locator = 1,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
                write_N = True,
                kwargs_zbins_rs = dict(smooth = True, sigma = 1.2, overlap = 0.4, nBox = 30),
                spearmanr = True,
                #legend = True,
    )

    ### Zones logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    input_kw.update        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V__Tg[iT]),
                    y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6),
                    z = H.O_O3N2_M13__g,
                    #xlabel = r'$\log\ \tau_V^\star$ (zones)',
                    ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$ (zones)',
                    #zlabel = r'12 + $\log$ O/H',
                    add_mask = mask__g,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####
    
    ### RBINS logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V__Trg[iT]),
                    y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                    z = H.O_O3N2_M13__rg,
                    #xlabel = r'(Rbins)',
                    #ylabel = r'(Rbins)',
                    #add_mask = mask__rg,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### oneHLR logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 2))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_oneHLR__Tg[iT]),
                    y = np.ma.log10(H.aSFRSD_oneHLR__Tg[iT] * 1e6),
                    z = H.O_O3N2_M13_oneHLR__g,
                    #xlabel = r'(@1HLR)',
                    #ylabel = r'(@1HLR)',
                    add_mask = mask_GAL__g,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)

    ### integrated logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 3))
    input_kw.update        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.integrated_tau_V__g),
                    y = np.ma.log10(H.integrated_SFRSD__Tg[iT] * 1e6),
                    z = H.integrated_O_O3N2_M13__g,
                    #xlabel = r'(integrated)',
                    #ylabel = r'(integrated)',
                    add_mask = mask_GAL__g,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    
    ################################
    ################################
    ################################

    ### Zones logtauV vs logSFRSDHa vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 0))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_neb__g),
                    y = np.ma.log10(H.SFRSD_Ha__g * 1e6),
                    z = H.O_O3N2_M13__g,
                    #xlabel = r'$\log\ \tau_V^{neb}$ (zones)',
                    ylabel = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$ (zones)',
                    #zlabel = r'12 + $\log$ O/H',
                    add_mask = mask__g,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### RBINS logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 1))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_neb__rg),
                    y = np.ma.log10(H.aSFRSD_Ha__rg * 1e6),
                    z = H.O_O3N2_M13__rg,
                    #xlabel = r'(Rbins)',
                    #ylabel = r'(Rbins)',
                    add_mask = mask__rg,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### oneHLR logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 2))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_neb_oneHLR__g),
                    y = np.ma.log10(H.aSFRSD_Ha_oneHLR__g * 1e6),
                    z = H.O_O3N2_M13_oneHLR__g,
                    #xlabel = r'(@1HLR)',
                    #ylabel = r'(@1HLR)',
                    add_mask = mask_GAL__g,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)

    ### integrated logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 3))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.integrated_tau_V_neb__g),
                    y = np.ma.log10(H.integrated_SFRSD_Ha__g * 1e6),
                    z = H.integrated_O_O3N2_M13__g,
                    #xlabel = r'(integrated)',
                    #ylabel = r'(integrated)',
                    add_mask = mask_GAL__g,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ################################
    ################################
    ################################

    ### Zones logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 0))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V__Tg[iT]),
                    y = np.ma.log10(H.SFRSD_Ha__g * 1e6),
                    z = H.O_O3N2_M13__g,
                    #xlabel = r'$\log\ \tau_V^\star$ (zones)',
                    ylabel = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$ (zones)',
                    #zlabel = r'12 + $\log$ O/H',
                    add_mask = mask__g,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### RBINS logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 1))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V__Trg[iT]),
                    y = np.ma.log10(H.aSFRSD_Ha__rg * 1e6),
                    z = H.O_O3N2_M13__rg,
                    #xlabel = r'(Rbins)',
                    #ylabel = r'(Rbins)',
                    add_mask = mask__rg,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### oneHLR logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 2))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_oneHLR__Tg[iT]),
                    y = np.ma.log10(H.aSFRSD_Ha_oneHLR__g * 1e6),
                    z = H.O_O3N2_M13_oneHLR__g,
                    #xlabel = r'(@1HLR)',
                    #ylabel = r'(@1HLR)',
                    add_mask = mask_GAL__g,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)

    ### integrated logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 3))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.integrated_tau_V__g),
                    y = np.ma.log10(H.integrated_SFRSD_Ha__g * 1e6),
                    z = H.integrated_O_O3N2_M13__g,
                    #xlabel = r'(integrated)',
                    #ylabel = r'(integrated)',
                    add_mask = mask_GAL__g,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_xticklabels(), visible = False)
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)

    ################################
    ################################
    ################################

    ### Zones logtauV vs logSFRSDHa vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 0))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_neb__g),
                    y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6),
                    z = H.O_O3N2_M13__g,
                    xlabel = r'$\log\ \tau_V^{neb}$ (zones)',
                    ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$ (zones)',
                    #zlabel = r'12 + $\log$ O/H',
                    add_mask = mask__g,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### RBINS logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 1))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_neb__rg),
                    y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                    z = H.O_O3N2_M13__rg,
                    #xlabel = r'(Rbins)',
                    #ylabel = r'(Rbins)',
                    add_mask = mask__rg,
                    zbins = 3,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### oneHLR logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 2))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.tau_V_neb_oneHLR__g),
                    y = np.ma.log10(H.aSFRSD_oneHLR__Tg[iT] * 1e6),
                    z = H.O_O3N2_M13_oneHLR__g,
                    #xlabel = r'(@1HLR)',
                    #ylabel = r'(@1HLR)',
                    add_mask = mask_GAL__g,
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)

    ### integrated logtauV vs logSFRSD vs 12 + log (O/H) M13 [3 bins] ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 3))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(H.integrated_tau_V_neb__g),
                    y = np.ma.log10(H.integrated_SFRSD__Tg[iT] * 1e6),
                    z = H.integrated_O_O3N2_M13__g,
                    #xlabel = r'(integrated)',
                    #ylabel = r'(integrated)',
                    add_mask = mask_GAL__g,
                    **input_kw
    )
    ax = r_kw['ax']
    plt.setp(ax.get_yticklabels(), visible = False)    
    f = r_kw['f']
    #ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####
    #f.subplots_adjust(bottom = 0.15, hspace = 0.4, wspace = 0.45, right = 0.95, left = 0.07)
    f.subplots_adjust(bottom = 0.15, hspace = 0, wspace = 0, right = 0.95, left = 0.07)

def plot_closed_box_SK(f, H, gas_dict, gasname, mask__g, mask__rg, mask_GAL__g):
    sc_kwargs = default_sc_kwargs.copy()
    rs_kwargs = default_rs_kwargs.copy()
    ols_kwargs = default_ols_kwargs.copy()
    ols_plot_kwargs = default_ols_plot_kwargs.copy()
    gas = gas_dict[gasname]
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    suptitle_R = r'%s Nzones:%d(%d gals) NRbins:%d(%d gals)  $t_{SF}$:%.2fMyr  $x_Y$(min):%.0f%%  $\tau_V^\star$(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f' % (gasname, (~mask__g).sum(), NgalsOkZones, (~mask__rg).sum(), NgalsOkRbins, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
    NRows = 4
    NCols = 4
    #page_size_inches = LetterSize_inches[::-1]
    page_size_inches = [11, 11]
    f.set_size_inches(page_size_inches)
    grid_shape = (NRows, NCols)
    f.suptitle(suptitle_R, fontsize = 10)

    input_kw = dict(
                f = f,
                debug = True,
                xlim = (-10, 2),
                ylim = (8., 8.8),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 4,
                x_minor_locator = .8,
                y_major_locator = .2,
                y_minor_locator = .05,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
    )
    
    ### Zones ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 0))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['f_gas__g']),
                    y = H.O_O3N2_M13__g,
                    z = H.zone_dist_HLR__g,
                    add_mask = mask__g,
                    xlabel = r'$\ln\ f_{gas}$ (zones) ',
                    ylabel = r'$12 + \log(O/H)$ (zones)',
                    zlabel = r'R [HLR]',
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    f = r_kw['f']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### RBINS ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 1))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['f_gas__rg']),
                    y = H.O_O3N2_M13__rg,
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'(Rbins)',
                    ylabel = r'(Rbins)',
                    **input_kw
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### oneHLR ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 2))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['f_gas_oneHLR']),
                    y = H.O_O3N2_M13_oneHLR__g,
                    add_mask = mask_GAL__g,
                    xlabel = r'(@1HLR)',
                    ylabel = r'(@1HLR)',
                    **input_kw
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ####

    ### integrated ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (0, 3))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['integrated_f_gas']),
                    y = H.integrated_O_O3N2_M13__g,
                    add_mask = mask_GAL__g,
                    xlabel = r'(integrated)',
                    ylabel = r'(integrated)',
                    **input_kw
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    input_kw = dict(
                f = f,
                debug = True,
                xlim = (-10, 2),
                ylim = (-2.5, 1),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 4,
                x_minor_locator = .8,
                y_major_locator = 1.,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
    )
    
    ### Zones ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 0))
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['f_gas__g']),
                    y = H.alogZ_mass__Ug[-1],
                    z = H.zone_dist_HLR__g,
                    add_mask = mask__g,
                    xlabel = r'$\ln\ f_{gas}$ (zones) ',
                    ylabel = r'$\langle \log\ Z_\star \rangle_M$ (t < %.2f Gyr) [$Z_\odot$] (zones)' % (H.tZ__U[-1] / 1e9),
                    zlabel = r'R [HLR]',
                    colorbar = False,
                    **input_kw
    )
    ax = r_kw['ax']
    f = r_kw['f']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### RBINS ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 1))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['f_gas__rg']),
                    y = H.alogZ_mass__Urg[-1],
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'(Rbins)',
                    ylabel = r'(Rbins)',
                    **input_kw
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    ### oneHLR ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 2))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['f_gas_oneHLR']),
                    y = H.alogZ_mass_oneHLR__Ug[-1],
                    add_mask = mask_GAL__g,
                    xlabel = r'(@1HLR)',
                    ylabel = r'(@1HLR)',
                    **input_kw
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ####

    ### integrated ln F_GAS vs ZNeb M13 ###
    ax = plt.subplot2grid(grid_shape, loc = (1, 3))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log(gas['integrated_f_gas']),
                    y = H.alogZ_mass_GAL__Ug[-1],
                    add_mask = mask_GAL__g,
                    xlabel = r'(integrated)',
                    ylabel = r'(integrated)',
                    **input_kw
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    ####

    kw_input = dict(
                f = f,
                debug = True,
                xlim = (-.5, 3.),
                ylim = (-3, 0.5),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 1.,
                x_minor_locator = .2,
                y_major_locator = 1.,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
    )
    
    ### Zones SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 0))
            
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['SigmaGas__g']),
                    y = np.ma.log10(H.SFRSD__Tg[iT] * 1e6),
                    z = H.zone_dist_HLR__g,
                    add_mask = mask__g,
                    xlabel = r'$\log\ \Sigma_{gas}\ [M_\odot pc^{-2}]$ (zones)',
                    ylabel = r'$\log\ \Sigma_{SFR}^\star(t_\star)\ [M_\odot yr^{-1} kpc^{-2}]$ (zones)',
                    zlabel = r'R [HLR]',
                    colorbar = False,
                    **kw_input 
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    p = [ SK_slope, np.log10(SK_zero) ]
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####

    ### Rbins SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 1))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['SigmaGas__rg']),
                    y = np.ma.log10(H.aSFRSD__Trg[iT] * 1e6),
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'(Rbins)',
                    ylabel = r'(Rbins)',
                    **kw_input 
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####
    
    ### oneHLR SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 2))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['SigmaGas_oneHLR']),
                    y = np.ma.log10(H.aSFRSD_oneHLR__Tg[iT] * 1e6),
                    add_mask = mask_GAL__g,
                    xlabel = r'(@1HLR)',
                    ylabel = r'(@1HLR)',
                    **kw_input
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####
    
    ### integrated SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (2, 3))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['integrated_SigmaGas']),
                    y = np.ma.log10(H.integrated_SFRSD__Tg[iT] * 1e6),
                    add_mask = mask_GAL__g,
                    xlabel = r'(integrated)',
                    ylabel = r'(integrated)',
                    **kw_input
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####        

    kw_input = dict(
                f = f,
                debug = True,
                xlim = (-.5, 3.),
                ylim = (-3, 0.5),
                zlim = (minR, 3),
                kwargs_sc_cb = dict(ticks = [ minR, 1, 2, 3 ]),
                x_major_locator = 1.,
                x_minor_locator = .2,
                y_major_locator = 1.,
                y_minor_locator = .2,
                kwargs_scatter = sc_kwargs,
                ols = True,
                kwargs_ols = ols_kwargs,
                kwargs_ols_plot = ols_plot_kwargs,
                return_kwargs = True,
    )
    

    ### Zones SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 0))
            
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['SigmaGas__g']),
                    y = np.ma.log10(H.SFRSD_Ha__g * 1e6),
                    z = H.zone_dist_HLR__g,
                    add_mask = mask__g,
                    xlabel = r'$\log\ \Sigma_{gas}\ [M_\odot pc^{-2}]$ (zones)',
                    ylabel = r'$\log\ \Sigma_{SFR}^{H\alpha}\ [M_\odot yr^{-1} kpc^{-2}]$ (zones)',
                    zlabel = r'R [HLR]',
                    colorbar = False,
                    **kw_input 
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    p = [ SK_slope, np.log10(SK_zero) ]
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####

    ### Rbins SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 1))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['SigmaGas__rg']),
                    y = np.ma.log10(H.aSFRSD_Ha__rg * 1e6),
                    z = H.Rtoplot(),
                    add_mask = mask__rg,
                    xlabel = r'(Rbins)',
                    ylabel = r'(Rbins)',
                    **kw_input 
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####
    
    ### oneHLR SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 2))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['SigmaGas_oneHLR']),
                    y = np.ma.log10(H.aSFRSD_Ha_oneHLR__g * 1e6),
                    add_mask = mask_GAL__g,
                    xlabel = r'(@1HLR)',
                    ylabel = r'(@1HLR)',
                    **kw_input
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####
    
    ### integrated SigmaGas vs SFRSD syn ###
    ax = plt.subplot2grid(grid_shape, loc = (3, 3))        
    r_kw = C.plot_zbins(
                    ax = ax,
                    x = np.ma.log10(gas['integrated_SigmaGas']),
                    y = np.ma.log10(H.integrated_SFRSD_Ha__g * 1e6),
                    add_mask = mask_GAL__g,
                    xlabel = r'(integrated)',
                    ylabel = r'(integrated)',
                    **kw_input
    )
    ax = r_kw['ax']
    ax.set_title('N:%d Rs:%.4f (p: %.1f) ' % (r_kw['xm'].count(), r_kw['Rs'], 100 * r_kw['Pvals']), fontsize = 8)
    f = r_kw['f']
    ax.plot(ax.get_xlim(), np.polyval(p, ax.get_xlim()), 'k:', label = 'SK')
    ####        
    f.subplots_adjust(bottom = 0.15, hspace = 0.4, wspace = 0.45, right = 0.95, left = 0.07)


if __name__ == '__main__':
    args = parser_args()
    
    C.debug_var(args.debug, args = args)
    
    H = C.H5SFRData(args.hdf5)
    iT = args.itSF
    iU = -1

    minR = 0
    
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
    
    Area_GAL__g = (H.Mcor_GAL__g / H.McorSD_GAL__g)
    dustdim = 0.2 # md / rhod

    gasnames = [ 'GRVTauVStar', 'GRVTauVNeb', 'BRTauVStar', 'BRTauVNeb', 'SKSFRSDStar', 'SKSFRSDHa' ]
    gas_dict = {}
    
    ######################
    # SK Law 
    ######################
    SK_zero = 1.6e-4
    SK_slope = 1.4
    aux = 1e6 ** (1. / SK_slope)
    SK_SigmaGas__g = aux * (H.SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__g = aux * (H.SFRSD_Ha__g / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas__rg = aux * (H.aSFRSD__Trg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__rg = aux * (H.aSFRSD_Ha__rg / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_oneHLR__g = aux * (H.aSFRSD_oneHLR__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha_oneHLR__g = aux * (H.aSFRSD_Ha_oneHLR__g / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas = aux * (H.integrated_SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas_Ha = aux * (H.integrated_SFRSD_Ha__g / SK_zero) ** (1. / SK_slope) 
    SK_DGR__g = dustdim * H.tau_V__Tg[iT] / SK_SigmaGas__g
    SK_DGR_Ha__g = dustdim * H.tau_V_neb__g / SK_SigmaGas_Ha__g
    SK_DGR__rg = dustdim * H.tau_V__Trg[iT] / SK_SigmaGas__rg
    SK_DGR_Ha__rg = dustdim * H.tau_V_neb__rg / SK_SigmaGas_Ha__rg
    SK_DGR_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / SK_SigmaGas_oneHLR__g
    SK_DGR_Ha_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__g / SK_SigmaGas_Ha_oneHLR__g
    SK_integrated_DGR = dustdim * H.integrated_tau_V__g / SK_integrated_SigmaGas
    SK_integrated_DGR_Ha = dustdim * H.integrated_tau_V_neb__g / SK_integrated_SigmaGas_Ha
    SK_GSR__g = SK_SigmaGas__g / H.McorSD__Tg[iT]
    SK_GSR_Ha__g = SK_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    SK_GSR__rg = SK_SigmaGas__rg / H.McorSD__Trg[iT]
    SK_GSR_Ha__rg = SK_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    SK_GSR_oneHLR__g = SK_SigmaGas_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    SK_GSR_Ha_oneHLR__g = SK_SigmaGas_Ha_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    SK_integrated_GSR = SK_integrated_SigmaGas / H.McorSD_GAL__g
    SK_integrated_GSR_Ha = SK_integrated_SigmaGas_Ha / H.McorSD_GAL__g
    SK_f_gas__g = 1. / (1. + 1. / SK_GSR__g)
    SK_f_gas_Ha__g = 1. / (1. + 1. / SK_GSR_Ha__g)
    SK_f_gas__rg = 1. / (1. + 1. / SK_GSR__rg)
    SK_f_gas_Ha__rg = 1. / (1. + 1. / SK_GSR_Ha__rg)
    SK_f_gas_oneHLR__g = 1. / (1. + 1. / SK_GSR_oneHLR__g)
    SK_f_gas_Ha_oneHLR__g = 1. / (1. + 1. / SK_GSR_Ha_oneHLR__g)
    SK_integrated_f_gas = 1. / (1. + 1. / SK_integrated_GSR)
    SK_integrated_f_gas_Ha = 1. / (1. + 1. / SK_integrated_GSR_Ha)
    gas_dict['SKSFRSDStar'] = dict(
        c = 'g',
        SigmaGas__g = SK_SigmaGas__g,
        SigmaGas__rg = SK_SigmaGas__rg,
        integrated_SigmaGas = SK_integrated_SigmaGas,
        SigmaGas_oneHLR = SK_SigmaGas_oneHLR__g,
        DGR__g = SK_DGR__g,
        DGR__rg = SK_DGR__rg,
        integrated_DGR = SK_integrated_DGR,
        DGR_oneHLR = SK_DGR_oneHLR__g,
        GSR__g = SK_GSR__g,
        GSR__rg = SK_GSR__rg,
        integrated_GSR = SK_integrated_GSR,
        GSR_oneHLR = SK_GSR_oneHLR__g,
        f_gas__g = SK_f_gas__g,
        f_gas__rg = SK_f_gas__rg,
        integrated_f_gas = SK_integrated_f_gas,
        f_gas_oneHLR = SK_f_gas_oneHLR__g,
    )
    gas_dict['SKSFRSDHa'] = dict(
        c = 'k',
        SigmaGas__g = SK_SigmaGas_Ha__g,
        SigmaGas__rg = SK_SigmaGas_Ha__rg,
        integrated_SigmaGas = SK_integrated_SigmaGas_Ha,
        SigmaGas_oneHLR = SK_SigmaGas_Ha_oneHLR__g,
        DGR__g = SK_DGR_Ha__g,
        DGR__rg = SK_DGR_Ha__rg,
        integrated_DGR = SK_integrated_DGR_Ha,
        DGR_oneHLR = SK_DGR_Ha_oneHLR__g,
        GSR__g = SK_GSR_Ha__g,
        GSR__rg = SK_GSR_Ha__rg,
        integrated_GSR = SK_integrated_GSR_Ha,
        GSR_oneHLR = SK_GSR_Ha_oneHLR__g,
        f_gas__g = SK_f_gas_Ha__g,
        f_gas__rg = SK_f_gas_Ha__rg,
        integrated_f_gas = SK_integrated_f_gas_Ha,
        f_gas_oneHLR = SK_f_gas_Ha_oneHLR__g,
    ) 
    ######################
    
    ######################
    # Guiderdoni & Rocca-Volmerange (1987)
    # DGR Calculado para Zsun 
    ######################
    GRV_DGR = 0.059
    GRV_SigmaGas__g = H.tau_V__Tg[iT] / GRV_DGR
    GRV_SigmaGas_Ha__g = H.tau_V_neb__g / GRV_DGR
    GRV_SigmaGas__rg = H.tau_V__Trg[iT] / GRV_DGR
    GRV_SigmaGas_Ha__rg = H.tau_V_neb__rg / GRV_DGR
    GRV_SigmaGas_oneHLR__g = H.tau_V_oneHLR__Tg[iT] / GRV_DGR
    GRV_SigmaGas_Ha_oneHLR__g = H.tau_V_neb_oneHLR__g / GRV_DGR
    GRV_integrated_SigmaGas = H.integrated_tau_V__g / GRV_DGR
    GRV_integrated_SigmaGas_Ha = H.integrated_tau_V_neb__g / GRV_DGR
    GRV_GSR__g = GRV_SigmaGas__g / H.McorSD__Tg[iT]
    GRV_GSR_Ha__g = GRV_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    GRV_GSR__rg = GRV_SigmaGas__rg / H.McorSD__Trg[iT]
    GRV_GSR_Ha__rg = GRV_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    GRV_GSR_oneHLR__g = GRV_SigmaGas_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    GRV_GSR_Ha_oneHLR__g = GRV_SigmaGas_Ha_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    GRV_integrated_GSR = GRV_integrated_SigmaGas / H.McorSD_GAL__g
    GRV_integrated_GSR_Ha = GRV_integrated_SigmaGas_Ha / H.McorSD_GAL__g 
    GRV_f_gas__g = 1. / (1. + 1. / GRV_GSR__g)
    GRV_f_gas_Ha__g = 1. / (1. + 1. / GRV_GSR_Ha__g)
    GRV_f_gas__rg = 1. / (1. + 1. / GRV_GSR__rg)
    GRV_f_gas_Ha__rg = 1. / (1. + 1. / GRV_GSR_Ha__rg)
    GRV_f_gas_oneHLR__g = 1. / (1. + 1. / GRV_GSR_oneHLR__g)
    GRV_f_gas_Ha_oneHLR__g = 1. / (1. + 1. / GRV_GSR_Ha_oneHLR__g)
    GRV_integrated_f_gas = 1. / (1. + 1. / GRV_integrated_GSR)
    GRV_integrated_f_gas_Ha = 1. / (1. + 1. / GRV_integrated_GSR_Ha)
    gas_dict['GRVTauVStar'] = dict(
        c = 'y',
        SigmaGas__g = GRV_SigmaGas__g,
        SigmaGas__rg = GRV_SigmaGas__rg,
        integrated_SigmaGas = GRV_integrated_SigmaGas,
        SigmaGas_oneHLR = GRV_SigmaGas_oneHLR__g,
        DGR__g = GRV_DGR,
        DGR__rg = GRV_DGR,
        integrated_DGR = GRV_DGR,
        DGR_oneHLR = GRV_DGR,
        GSR__g = GRV_GSR__g,
        GSR__rg = GRV_GSR__rg,
        integrated_GSR = GRV_integrated_GSR,
        GSR_oneHLR = GRV_GSR_oneHLR__g,
        f_gas__g = GRV_f_gas__g,
        f_gas__rg = GRV_f_gas__rg,
        integrated_f_gas = GRV_integrated_f_gas,
        f_gas_oneHLR = GRV_f_gas_oneHLR__g,
    )
    gas_dict['GRVTauVNeb'] = dict(
        c = 'orange',
        SigmaGas__g = GRV_SigmaGas_Ha__g,
        SigmaGas__rg = GRV_SigmaGas_Ha__rg,
        integrated_SigmaGas = GRV_integrated_SigmaGas_Ha,
        SigmaGas_oneHLR = GRV_SigmaGas_Ha_oneHLR__g,
        DGR__g = GRV_DGR,
        DGR__rg = GRV_DGR,
        integrated_DGR = GRV_DGR,
        DGR_oneHLR = GRV_DGR,
        GSR__g = GRV_GSR_Ha__g,
        GSR__rg = GRV_GSR_Ha__rg,
        integrated_GSR = GRV_integrated_GSR_Ha,
        GSR_oneHLR = GRV_GSR_Ha_oneHLR__g,
        f_gas__g = GRV_f_gas_Ha__g,
        f_gas__rg = GRV_f_gas_Ha__rg,
        integrated_f_gas = GRV_integrated_f_gas_Ha,
        f_gas_oneHLR = GRV_f_gas_Ha_oneHLR__g,
    )

    ######################
    
    ######################
    #Brinchmann
    ######################
    DGR_conv_lim_sup = 1.1e-2
    DGR_conv_lim_inf = 5.3e-3
    DGR_interval = np.array([DGR_conv_lim_inf, DGR_conv_lim_sup])
    DGR_cte = DGR_interval.mean()
    OHSunBrinch_inv = 1 / (10.**(8.82 - 12))
    ######################
    # from OLS Bisector
    p_ols = np.array([1.87, -6.98])
    BR_OHBrinch_ols__g = np.ma.masked_all((H.O_O3N2_M13__g.shape))
    BR_OHBrinch_ols__g[~H.O_O3N2_M13__g.mask] = np.polyval(p_ols, H.O_O3N2_M13__g.compressed()) 
    BR_OHBrinch_ols__rg = np.ma.masked_all((H.O_O3N2_M13__rg.shape))
    BR_OHBrinch_ols__rg[~H.O_O3N2_M13__rg.mask] = np.polyval(p_ols, H.O_O3N2_M13__rg.compressed())
    BR_OHBrinch_ols_oneHLR__g = np.ma.masked_all((H.O_O3N2_M13_oneHLR__g.shape))
    BR_OHBrinch_ols_oneHLR__g[~H.O_O3N2_M13_oneHLR__g.mask] = np.polyval(p_ols, H.O_O3N2_M13_oneHLR__g.compressed())
    BR_integrated_OHBrinch_ols = np.ma.masked_all((H.integrated_O_O3N2_M13__g.shape))
    BR_integrated_OHBrinch_ols[~H.integrated_O_O3N2_M13__g.mask] = np.polyval(p_ols, H.integrated_O_O3N2_M13__g.compressed())
    BR_DGR_up_ols__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_up_ols__rg = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__rg = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_DGR_ols__g = DGR_cte * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_ols__rg = DGR_cte * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_DGR_ols_oneHLR__g = DGR_cte * (10 ** (BR_OHBrinch_ols_oneHLR__g - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols_oneHLR__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols_oneHLR__g - 12) * OHSunBrinch_inv)
    BR_DGR_up_ols_oneHLR__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols_oneHLR__g - 12) * OHSunBrinch_inv)
    BR_integrated_DGR_ols = DGR_cte * (10 ** (BR_integrated_OHBrinch_ols - 12) * OHSunBrinch_inv)
    BR_integrated_DGR_up_ols = DGR_conv_lim_sup * (10 ** (BR_integrated_OHBrinch_ols - 12) * OHSunBrinch_inv)
    BR_integrated_DGR_down_ols = DGR_conv_lim_inf * (10 ** (BR_integrated_OHBrinch_ols - 12) * OHSunBrinch_inv)
    BR_SigmaGas_up_ols__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_up_ols__g
    BR_SigmaGas_up_ols__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_up_ols__rg
    BR_SigmaGas_Ha_up_ols__g = dustdim * H.tau_V_neb__g / BR_DGR_up_ols__g
    BR_SigmaGas_Ha_up_ols__rg = dustdim * H.tau_V_neb__rg / BR_DGR_up_ols__rg
    BR_SigmaGas_down_ols__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_down_ols__g
    BR_SigmaGas_down_ols__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_down_ols__rg
    BR_SigmaGas_Ha_down_ols__g = dustdim * H.tau_V_neb__g / BR_DGR_down_ols__g
    BR_SigmaGas_Ha_down_ols__rg = dustdim * H.tau_V_neb__rg / BR_DGR_down_ols__rg
    BR_SigmaGas_ols__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_ols__g
    BR_SigmaGas_ols__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_ols__rg
    BR_SigmaGas_Ha_ols__g = dustdim * H.tau_V_neb__g / BR_DGR_ols__g
    BR_SigmaGas_Ha_ols__rg = dustdim * H.tau_V_neb__rg / BR_DGR_ols__rg
    BR_SigmaGas_ols_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / BR_DGR_ols_oneHLR__g
    BR_SigmaGas_up_ols_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / BR_DGR_up_ols_oneHLR__g
    BR_SigmaGas_down_ols_oneHLR__g = dustdim * H.tau_V_oneHLR__Tg[iT] / BR_DGR_down_ols_oneHLR__g
    BR_SigmaGas_Ha_ols_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__g / BR_DGR_ols_oneHLR__g
    BR_SigmaGas_Ha_up_ols_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__g / BR_DGR_up_ols_oneHLR__g
    BR_SigmaGas_Ha_down_ols_oneHLR__g = dustdim * H.tau_V_neb_oneHLR__g / BR_DGR_down_ols_oneHLR__g
    BR_integrated_SigmaGas_ols = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_ols
    BR_integrated_SigmaGas_up_ols = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_up_ols
    BR_integrated_SigmaGas_down_ols = dustdim * H.integrated_tau_V__g / BR_integrated_DGR_down_ols
    BR_integrated_SigmaGas_Ha_ols = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_ols
    BR_integrated_SigmaGas_Ha_up_ols = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_up_ols
    BR_integrated_SigmaGas_Ha_down_ols = dustdim * H.integrated_tau_V_neb__g / BR_integrated_DGR_down_ols
    BR_GSR_up_ols__g = BR_SigmaGas_up_ols__g / H.McorSD__Tg[iT]
    BR_GSR_up_ols__rg = BR_SigmaGas_up_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_up_ols__g = BR_SigmaGas_Ha_up_ols__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_up_ols__rg = BR_SigmaGas_Ha_up_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_down_ols__g = BR_SigmaGas_down_ols__g / H.McorSD__Tg[iT]
    BR_GSR_down_ols__rg = BR_SigmaGas_down_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_down_ols__g = BR_SigmaGas_Ha_down_ols__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_down_ols__rg = BR_SigmaGas_Ha_down_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_ols__g = BR_SigmaGas_ols__g / H.McorSD__Tg[iT]
    BR_GSR_ols__rg = BR_SigmaGas_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_ols__g = BR_SigmaGas_Ha_ols__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_ols__rg = BR_SigmaGas_Ha_ols__rg / H.McorSD__Trg[iT]
    BR_GSR_ols_oneHLR__g = BR_SigmaGas_ols_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_up_ols_oneHLR__g = BR_SigmaGas_up_ols_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_down_ols_oneHLR__g = BR_SigmaGas_down_ols_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_Ha_ols_oneHLR__g = BR_SigmaGas_Ha_ols_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_Ha_up_ols_oneHLR__g = BR_SigmaGas_Ha_up_ols_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_GSR_Ha_down_ols_oneHLR__g = BR_SigmaGas_Ha_down_ols_oneHLR__g / H.McorSD_oneHLR__Tg[iT]
    BR_integrated_GSR_ols = BR_integrated_SigmaGas_ols / H.McorSD_GAL__g
    BR_integrated_GSR_up_ols = BR_integrated_SigmaGas_up_ols / H.McorSD_GAL__g
    BR_integrated_GSR_down_ols = BR_integrated_SigmaGas_down_ols / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_ols = BR_integrated_SigmaGas_Ha_ols / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_up_ols = BR_integrated_SigmaGas_Ha_up_ols / H.McorSD_GAL__g
    BR_integrated_GSR_Ha_down_ols = BR_integrated_SigmaGas_Ha_down_ols / H.McorSD_GAL__g
    BR_f_gas_up_ols__g = 1. / (1. + 1. / BR_GSR_up_ols__g)
    BR_f_gas_up_ols__rg = 1. / (1. + 1. / BR_GSR_up_ols__rg)
    BR_f_gas_Ha_up_ols__g = 1. / (1. + 1. / BR_GSR_Ha_up_ols__g)
    BR_f_gas_Ha_up_ols__rg = 1. / (1. + 1. / BR_GSR_Ha_up_ols__rg)
    BR_f_gas_down_ols__g = 1. / (1. + 1. / BR_GSR_down_ols__g)
    BR_f_gas_down_ols__rg = 1. / (1. + 1. / BR_GSR_down_ols__rg)
    BR_f_gas_Ha_down_ols__g = 1. / (1. + 1. / BR_GSR_Ha_down_ols__g)
    BR_f_gas_Ha_down_ols__rg = 1. / (1. + 1. / BR_GSR_Ha_down_ols__rg)
    BR_f_gas_ols__g = 1. / (1. + 1. / BR_GSR_ols__g)
    BR_f_gas_ols__rg = 1. / (1. + 1. / BR_GSR_ols__rg)
    BR_f_gas_Ha_ols__g = 1. / (1. + 1. / BR_GSR_Ha_ols__g)
    BR_f_gas_Ha_ols__rg = 1. / (1. + 1. / BR_GSR_Ha_ols__rg)
    BR_f_gas_ols_oneHLR__g = 1. / (1. + 1. / BR_GSR_ols_oneHLR__g)
    BR_f_gas_up_ols_oneHLR__g = 1. / (1. + 1. / BR_GSR_up_ols_oneHLR__g)
    BR_f_gas_down_ols_oneHLR__g = 1. / (1. + 1. / BR_GSR_down_ols_oneHLR__g)
    BR_f_gas_Ha_ols_oneHLR__g = 1. / (1. + 1. / BR_GSR_Ha_ols_oneHLR__g)
    BR_f_gas_Ha_up_ols_oneHLR__g = 1. / (1. + 1. / BR_GSR_Ha_up_ols_oneHLR__g)
    BR_f_gas_Ha_down_ols_oneHLR__g = 1. / (1. + 1. / BR_GSR_Ha_down_ols_oneHLR__g)
    BR_integrated_f_gas_ols = 1. / (1. + 1. / BR_integrated_GSR_ols)
    BR_integrated_f_gas_up_ols = 1. / (1. + 1. / BR_integrated_GSR_up_ols)
    BR_integrated_f_gas_down_ols = 1. / (1. + 1. / BR_integrated_GSR_down_ols)
    BR_integrated_f_gas_Ha_ols = 1. / (1. + 1. / BR_integrated_GSR_Ha_ols)
    BR_integrated_f_gas_Ha_up_ols = 1. / (1. + 1. / BR_integrated_GSR_Ha_up_ols)
    BR_integrated_f_gas_Ha_down_ols = 1. / (1. + 1. / BR_integrated_GSR_Ha_down_ols)
    gas_dict['BRTauVStar'] = dict(
        c = 'b',
        SigmaGas__g = BR_SigmaGas_ols__g,
        SigmaGas__rg = BR_SigmaGas_ols__rg,
        integrated_SigmaGas = BR_integrated_SigmaGas_ols,
        SigmaGas_oneHLR = BR_SigmaGas_ols_oneHLR__g,
        DGR__g = BR_DGR_ols__g,
        DGR__rg = BR_DGR_ols__rg,
        integrated_DGR = BR_integrated_DGR_ols,
        DGR_oneHLR = BR_DGR_ols_oneHLR__g,
        GSR__g = BR_GSR_ols__g,
        GSR__rg = BR_GSR_ols__rg,
        integrated_GSR = BR_integrated_GSR_ols,
        GSR_oneHLR = BR_GSR_ols_oneHLR__g,
        f_gas__g = BR_f_gas_ols__g,
        f_gas__rg = BR_f_gas_ols__rg,
        integrated_f_gas = BR_integrated_f_gas_ols,
        f_gas_oneHLR = BR_f_gas_ols_oneHLR__g,
        up = dict(
            SigmaGas__g = BR_SigmaGas_up_ols__g,
            SigmaGas__rg = BR_SigmaGas_up_ols__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_up_ols,
            SigmaGas_oneHLR = BR_SigmaGas_up_ols_oneHLR__g,
            DGR__g = BR_DGR_up_ols__g,
            DGR__rg = BR_DGR_up_ols__rg,
            integrated_DGR = BR_integrated_DGR_up_ols,
            DGR_oneHLR = BR_DGR_up_ols_oneHLR__g,
            GSR__g = BR_GSR_up_ols__g,
            GSR__rg = BR_GSR_up_ols__rg,
            integrated_GSR = BR_integrated_GSR_up_ols,
            GSR_oneHLR = BR_GSR_up_ols_oneHLR__g,
            f_gas__g = BR_f_gas_up_ols__g,
            f_gas__rg = BR_f_gas_up_ols__rg,
            integrated_f_gas = BR_integrated_f_gas_up_ols,
            f_gas_oneHLR = BR_f_gas_up_ols_oneHLR__g,
        ),
        down = dict(
            SigmaGas__g = BR_SigmaGas_down_ols__g,
            SigmaGas__rg = BR_SigmaGas_down_ols__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_down_ols,
            SigmaGas_oneHLR = BR_SigmaGas_down_ols_oneHLR__g,
            DGR__g = BR_DGR_down_ols__g,
            DGR__rg = BR_DGR_down_ols__rg,
            integrated_DGR = BR_integrated_DGR_down_ols,
            DGR_oneHLR = BR_DGR_down_ols_oneHLR__g,
            GSR__g = BR_GSR_down_ols__g,
            GSR__rg = BR_GSR_down_ols__rg,
            integrated_GSR = BR_integrated_GSR_down_ols,
            GSR_oneHLR = BR_GSR_down_ols_oneHLR__g,
            f_gas__g = BR_f_gas_down_ols__g,
            f_gas__rg = BR_f_gas_down_ols__rg,
            integrated_f_gas = BR_integrated_f_gas_down_ols,
            f_gas_oneHLR = BR_f_gas_down_ols_oneHLR__g,
        ),
    )
    gas_dict['BRTauVNeb'] = dict(
        c = 'r',
        SigmaGas__g = BR_SigmaGas_Ha_ols__g,
        SigmaGas__rg = BR_SigmaGas_Ha_ols__rg,
        integrated_SigmaGas = BR_integrated_SigmaGas_Ha_ols,
        SigmaGas_oneHLR = BR_SigmaGas_Ha_ols_oneHLR__g,
        DGR__g = BR_DGR_ols__g,
        DGR__rg = BR_DGR_ols__rg,
        integrated_DGR = BR_integrated_DGR_ols,
        DGR_oneHLR = BR_DGR_ols_oneHLR__g,
        GSR__g = BR_GSR_Ha_ols__g,
        GSR__rg = BR_GSR_Ha_ols__rg,
        integrated_GSR = BR_integrated_GSR_Ha_ols,
        GSR_oneHLR = BR_GSR_Ha_ols_oneHLR__g,
        f_gas__g = BR_f_gas_Ha_ols__g,
        f_gas__rg = BR_f_gas_Ha_ols__rg,
        integrated_f_gas = BR_integrated_f_gas_Ha_ols,
        f_gas_oneHLR = BR_f_gas_Ha_ols_oneHLR__g,
        up = dict(
            SigmaGas__g = BR_SigmaGas_Ha_up_ols__g,
            SigmaGas__rg = BR_SigmaGas_Ha_up_ols__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_Ha_up_ols,
            SigmaGas_oneHLR = BR_SigmaGas_Ha_up_ols_oneHLR__g,
            DGR__g = BR_DGR_up_ols__g,
            DGR__rg = BR_DGR_up_ols__rg,
            integrated_DGR = BR_integrated_DGR_up_ols,
            DGR_oneHLR = BR_DGR_up_ols_oneHLR__g,
            GSR__g = BR_GSR_Ha_up_ols__g,
            GSR__rg = BR_GSR_Ha_up_ols__rg,
            integrated_GSR = BR_integrated_GSR_Ha_up_ols,
            GSR_oneHLR = BR_GSR_Ha_up_ols_oneHLR__g,
            f_gas__g = BR_f_gas_Ha_up_ols__g,
            f_gas__rg = BR_f_gas_Ha_up_ols__rg,
            integrated_f_gas = BR_integrated_f_gas_Ha_up_ols,
            f_gas_oneHLR = BR_f_gas_Ha_up_ols_oneHLR__g,
        ),
        down = dict(
            SigmaGas__g = BR_SigmaGas_Ha_down_ols__g,
            SigmaGas__rg = BR_SigmaGas_Ha_down_ols__rg,
            integrated_SigmaGas = BR_integrated_SigmaGas_Ha_down_ols,
            SigmaGas_oneHLR = BR_SigmaGas_Ha_down_ols_oneHLR__g,
            DGR__g = BR_DGR_down_ols__g,
            DGR__rg = BR_DGR_down_ols__rg,
            integrated_DGR = BR_integrated_DGR_down_ols,
            DGR_oneHLR = BR_DGR_down_ols_oneHLR__g,
            GSR__g = BR_GSR_Ha_down_ols__g,
            GSR__rg = BR_GSR_Ha_down_ols__rg,
            integrated_GSR = BR_integrated_GSR_Ha_down_ols,
            GSR_oneHLR = BR_GSR_Ha_down_ols_oneHLR__g,
            f_gas__g = BR_f_gas_Ha_down_ols__g,
            f_gas__rg = BR_f_gas_Ha_down_ols__rg,
            integrated_f_gas = BR_integrated_f_gas_Ha_down_ols,
            f_gas_oneHLR = BR_f_gas_Ha_down_ols_oneHLR__g,
        ),
    )

    ######################
    # from cubic polynomial fit
    #p_cubic = np.array([-4.91783872, 122.48149162, -1014.51941088, 2803.24285985])
    ######################
    
    ####################################
    #### Mass from CALIFA Datafiles ####
    ####################################
    ####################################
    '''
    FILE: M_H1_CALIFA.csv - Miguel A. Perez
    delimiter = ','
    comment = '#'
    columns:
        1 - CALIFA No
        2 - NED Name
        3 - Distance
        4 - RA(J2000.0)
        5 - DEC(J2000.0)
        6 - Sobs
        7 - Scor
        8 - Sabs
        9 - e_Sabs,
        10 - M(HI)abs
        11 - e_M(HI)abs
        
    FILE: CALIFA_HI_angel.dat - Angel R. Lopez-Sanchez
    delimiter = ','
    comments = ';'
    columns:
        1 - NED Name
        2 - redshift
        3 - log(Ms)
        4 - e_log(Ms)
        5 - 12+log(O/H)
        6 - e_12+log(O/H)
        7 - log(SFR)
        8 - e_log(SFR)
        9 - log(Mg = 1.32 MHI)
        10 - e_log(Mg = 1.32 MHI)
    '''
    dirs = C.CALIFAPaths()
    file_miguel = '%sM_HI_CALIFA.csv' % dirs.califa_work_dir
    dtype_miguel = np.dtype([('califaID', '|S5'), ('M_HI', np.float)])
    read_miguel = np.loadtxt(file_miguel, delimiter = ',', usecols = [0, 9], dtype = dtype_miguel)
    map_miguel = {}
    for i, g in enumerate(read_miguel['califaID']):
        map_miguel[g] = i
    aux = set(H.califaIDs.tolist())
    gals_miguel_intersect = sorted([g for g in map_miguel.keys() if g in aux])
    gals_miguel_slice = np.zeros(H.califaIDs_all.shape, dtype = np.bool)
    integrated_M_HI_miguel__g = np.ma.masked_all(H.califaIDs_all.shape, dtype = np.float_)
    for g in gals_miguel_intersect:
        i = H.califaIDs_all.tolist().index(g)
        i_r = map_miguel[g]
        #print g, i, i_r
        gals_miguel_slice[i] = True
        integrated_M_HI_miguel__g[i] = read_miguel['M_HI'][i_r]
    integrated_f_gas_miguel = 1. / (1. + 1 / ((integrated_M_HI_miguel__g) / H.Mcor_GAL__g))
    #integrated_f_gas_miguel = 1./(1. + 1/((1.32 * integrated_M_HI_miguel__g) / H.Mcor_GAL__g))
    #integrated_f_gas_angel = 1./(1. + 1/((10. ** integrated_log_M_gas_angel__g) / H.Mcor_GAL__g))

    Area_GAL__g = (H.Mcor_GAL__g / H.McorSD_GAL__g)
    
    ##########################
    ######### MASKS ##########
    ##########################
    
    ba_max = 0
    #ba_max = 0.7
    #mask_GAL__g = np.bitwise_or(np.less(H.integrated_EW_Ha__g, 3.), np.less(H.ba_GAL__g, ba_max))
    mask_GAL__g = np.bitwise_or(np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool), np.less(H.ba_GAL__g, ba_max))
    #mask_GAL__g = np.bitwise_or(np.less(H.integrated_EW_Ha__g, 3.), np.zeros_like(H.ba_GAL__g, dtype = np.bool))
    #mask_GAL__g = np.zeros_like(H.integrated_EW_Ha__g, dtype = np.bool)
    
    mask__g = np.bitwise_or(np.ma.log10(H.SFRSD__Tg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Tg[iT]).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.SFRSD_Ha__g * 1e6).mask)
    mask__g = np.bitwise_or(mask__g, np.ma.log10(H.tau_V_neb__g).mask)
    mask__g = np.bitwise_or(mask__g, H.O_O3N2_M13__g.mask)
    #mask__g = np.bitwise_or(mask__g, np.less(H.EW_Ha__g, 3.))
    mask__g = np.bitwise_or(mask__g, np.less(H.reply_arr_by_zones(H.ba_GAL__g), ba_max))
    mask__g = np.bitwise_or(mask__g, ~maskRadiusOk__g)
    mask__g = np.bitwise_or(mask__g, ~gals_slice__g)

    #mask__g = ~maskRadiusOk__g
    
    mask__rg = np.bitwise_or(np.ma.log10(H.aSFRSD__Trg[iT] * 1e6).mask, np.ma.log10(H.tau_V__Trg[iT]).mask)
    mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.aSFRSD_Ha__rg * 1e6).mask)
    mask__rg = np.bitwise_or(mask__rg, np.ma.log10(H.tau_V_neb__rg).mask)
    mask__rg = np.bitwise_or(mask__rg, H.O_O3N2_M13__rg.mask)
    #mask__rg = np.bitwise_or(mask__rg, np.less(H.EW_Ha__rg, 3.))
    mask__rg = np.bitwise_or(mask__rg, np.less(H.reply_arr_by_radius(H.ba_GAL__g), ba_max))
    mask__rg = np.bitwise_or(mask__rg, ~maskRadiusOk__rg)
    mask__rg = np.bitwise_or(mask__rg, ~gals_slice__rg)
    #mask__rg = ~maskRadiusOk__rg
    
    NgalsOkZones = len(np.unique(H.reply_arr_by_zones(H.califaIDs)[~mask__g]))  
    NgalsOkRbins = len(np.unique(H.reply_arr_by_radius(H.califaIDs_all)[~mask__rg]))
    
    if args.dryrun: sys.exit()
    ####################################################
    if args.plot_SKdev:
        props = [ 
            'atfluxR',
            'xYR',
            'logO3N2M13R',
            'alogZmassR',
            'logMcorSDR' ,
            'morfTypeR',
        ]
        if args.output is None:
            if args.maskradius is not None:
                output = 'SKdev_R%.1f.pdf' % args.maskradius
            else:
                output = 'SKdev.pdf'
        else:
            output = args.output
        with PdfPages(output) as pdf:
            #with PdfPages('SK_dev_%s.pdf' % prop_key) as pdf:
            for prop_key in props:
                _, prop_dict = H.get_plot_dict(iT, -1, prop_key)
                f = plt.figure()
                plot_SK_dev(f, H, iT, prop_dict, mask__g, mask__rg, mask_GAL__g)
                plt.close(f)
                pdf.savefig(f)
    ####################################################
    if args.plot_SKmet:
        if args.output is None:
            if args.maskradius is not None:
                output = 'SKmet_R%.1f.pdf' % args.maskradius
            else:
                output = 'SKmet.pdf'
        else:
            output = args.output
     
        with PdfPages(output) as pdf:
            f = plt.figure()
            plot_SK_Zneb_bins(f, H, iT, mask__g, mask__rg, mask_GAL__g)
            plt.close(f)
            pdf.savefig(f)
    ####################################################
    if args.plot_closedbox:
        if args.output is None:
            if args.maskradius is not None:
                output = 'closedbox_R%.1f.pdf' % args.maskradius
            else:
                output = 'closedbox.pdf'
        else:
            output = args.output

        with PdfPages(output) as pdf:
            for gasname in gasnames:
                if args.output is None:
                    if args.maskradius is not None:
                        output = '%s_R%.1f.pdf' % (gasname, args.maskradius)
                    else:
                        output = '%s.pdf' % gasname
                else:
                    output = args.output
                     
                f = plt.figure()
                plot_closed_box_SK(f, H, gas_dict, gasname, mask__g, mask__rg, mask_GAL__g)
                plt.close(f)
                pdf.savefig(f)
    ####################################################
    ####################################################
    
