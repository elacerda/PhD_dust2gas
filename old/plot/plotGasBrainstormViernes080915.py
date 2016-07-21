#!/usr/bin/python
#
# Lacerda@Granada - 13/Oct/2014
#
from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import MaxNLocator
from CALIFAUtils.plots import plot_text_ax
#from CALIFAUtils.plots import plot_zbins
from CALIFAUtils.objects import runstats
from matplotlib import pyplot as plt
#import matplotlib as mpl
import CALIFAUtils as C
import argparse as ap
import numpy as np
import sys

RNuc = 0.5

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# mpl.rcParams['font.size'] = 14
# mpl.rcParams['axes.labelsize'] = 12
# mpl.rcParams['axes.titlesize'] = 16
# mpl.rcParams['xtick.labelsize'] = 10
# mpl.rcParams['ytick.labelsize'] = 10 
# mpl.rcParams['font.family'] = 'serif'
# mpl.rcParams['font.serif'] = 'Times New Roman'
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

def plot_text_topright_color(ax, txt, c, fs = 11):
    kw_text = dict(pos_x = 0.99, pos_y = 0.99, c = c, fs = fs, va = 'top', ha = 'right',)
    plot_text_ax(ax, txt, **kw_text)

def parser_args():        
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])

    default = {
        'debug' : False,
        'scatter' : False,
        'hdf5' : None,
        'output' : None,
        'itSF' : 11,
        'maskradius' : None,
        'slice_gals' : None,
        'dryrun' : False,
    }
    
    parser.add_argument('--debug', '-D',
                        action = 'store_true',
                        default = default['debug'])
    parser.add_argument('--dryrun',
                        action = 'store_true',
                        default = default['dryrun'])
    parser.add_argument('--scatter', '-s',
                        action = 'store_true',
                        default = default['scatter'])
    parser.add_argument('--hdf5', '-H',
                        metavar = 'FILE',
                        type = str,
                        default = default['hdf5'])
    parser.add_argument('--slice_gals', '-S',
                        metavar = 'FILE',
                        type = str,
                        default = default['slice_gals'])
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

    return parser.parse_args()

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
        gals_txt = ''
    else:
        gals_slice__g, N_gals = H.get_mask_zones_list(args.slice_gals, return_ngals = True)
        gals_slice__rg, N_gals = H.get_mask_radius_list(args.slice_gals, return_ngals = True)
        gals_txt = (args.slice_gals).split('/')[-1]
    
    dustdim = 0.2 # md / rhod
    
    ######################
    # SK Law 
    ######################
    SK_zero = 1.6e-4
    #SK_slope = 1.4
    SK_slope = 1.4
    aux = 1e6 ** (1. / SK_slope)
    SK_SigmaGas__g = aux * (H.SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__g = aux * (H.SFRSD_Ha__g / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas__rg = aux * (H.aSFRSD__Trg[iT] / SK_zero) ** (1. / SK_slope)
    SK_SigmaGas_Ha__rg = aux * (H.aSFRSD_Ha__rg / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas = aux * (H.integrated_SFRSD__Tg[iT] / SK_zero) ** (1. / SK_slope)
    SK_integrated_SigmaGas_Ha = aux * (H.integrated_SFRSD_Ha__g / SK_zero) ** (1. / SK_slope) 
    SK_DGR__g = dustdim * H.tau_V__Tg[iT] / SK_SigmaGas__g
    SK_DGR_Ha__g = dustdim * H.tau_V_neb__g / SK_SigmaGas_Ha__g
    SK_DGR__rg = dustdim * H.tau_V__Trg[iT] / SK_SigmaGas__rg
    SK_DGR_Ha__rg = dustdim * H.tau_V_neb__rg / SK_SigmaGas_Ha__rg
    SK_integrated_DGR = dustdim * H.integrated_tau_V__g / SK_integrated_SigmaGas
    SK_integrated_DGR_Ha = dustdim * H.integrated_tau_V_neb__g / SK_integrated_SigmaGas_Ha
    SK_GSR__g = SK_SigmaGas__g / H.McorSD__Tg[iT]
    SK_GSR_Ha__g = SK_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    SK_GSR__rg = SK_SigmaGas__rg / H.McorSD__Trg[iT]
    SK_GSR_Ha__rg = SK_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    SK_integrated_GSR = SK_integrated_SigmaGas / H.McorSD_GAL__g
    SK_integrated_GSR_Ha = SK_integrated_SigmaGas_Ha / H.McorSD_GAL__g
    SK_f_gas__g = 1. / (1. + 1. / SK_GSR__g)
    SK_f_gas_Ha__g = 1. / (1. + 1. / SK_GSR_Ha__g)
    SK_f_gas__rg = 1. / (1. + 1. / SK_GSR__rg)
    SK_f_gas_Ha__rg = 1. / (1. + 1. / SK_GSR_Ha__rg)
    SK_integrated_f_gas = 1. / (1. + 1. / SK_integrated_GSR)
    SK_integrated_f_gas_Ha = 1. / (1. + 1. / SK_integrated_GSR_Ha)
    ######################
    
    ######################
    # Remy-Ruyer
    ######################
    RR_DGR = 10. ** (-2.21) # 0.006166
    RR_SigmaGas__g = dustdim * H.tau_V__Tg[iT] / RR_DGR
    RR_SigmaGas_Ha__g = dustdim * H.tau_V_neb__g / RR_DGR
    RR_SigmaGas__rg = dustdim * H.tau_V__Trg[iT] / RR_DGR
    RR_SigmaGas_Ha__rg = dustdim * H.tau_V_neb__rg / RR_DGR
    RR_integrated_SigmaGas = dustdim * H.integrated_tau_V__g / RR_DGR
    RR_integrated_SigmaGas_Ha = dustdim * H.integrated_tau_V_neb__g / RR_DGR
    RR_GSR__g = RR_SigmaGas__g / H.McorSD__Tg[iT]
    RR_GSR_Ha__g = RR_SigmaGas_Ha__g / H.McorSD__Tg[iT]
    RR_GSR__rg = RR_SigmaGas__rg / H.McorSD__Trg[iT]
    RR_GSR_Ha__rg = RR_SigmaGas_Ha__rg / H.McorSD__Trg[iT]
    RR_integrated_GSR = RR_integrated_SigmaGas / H.McorSD_GAL__g
    RR_integrated_GSR_Ha = RR_integrated_SigmaGas_Ha / H.McorSD_GAL__g 
    RR_f_gas__g = 1. / (1. + 1. / RR_GSR__g)
    RR_f_gas_Ha__g = 1. / (1. + 1. / RR_GSR_Ha__g)
    RR_f_gas__rg = 1. / (1. + 1. / RR_GSR__rg)
    RR_f_gas_Ha__rg = 1. / (1. + 1. / RR_GSR_Ha__rg)
    RR_integrated_f_gas = 1. / (1. + 1. / RR_integrated_GSR)
    RR_integrated_f_gas_Ha = 1. / (1. + 1. / RR_integrated_GSR_Ha)
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
    BR_DGR_up_ols__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_up_ols__rg = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_down_ols__rg = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
    BR_DGR_ols__g = DGR_cte * (10 ** (BR_OHBrinch_ols__g - 12) * OHSunBrinch_inv)
    BR_DGR_ols__rg = DGR_cte * (10 ** (BR_OHBrinch_ols__rg - 12) * OHSunBrinch_inv)
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
    ######################
    # from cubic polynomial fit
    p_cubic = np.array([-4.91783872, 122.48149162, -1014.51941088, 2803.24285985])
    BR_OHBrinch_cubic__g = np.ma.masked_all((H.O_O3N2_M13__g.shape))
    BR_OHBrinch_cubic__g[~H.O_O3N2_M13__g.mask] = np.polyval(p_cubic, H.O_O3N2_M13__g.compressed()) 
    BR_OHBrinch_cubic__rg = np.ma.masked_all((H.O_O3N2_M13__rg.shape))
    BR_OHBrinch_cubic__rg[~H.O_O3N2_M13__rg.mask] = np.polyval(p_cubic, H.O_O3N2_M13__rg.compressed()) 
    BR_DGR_up_cubic__g = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_cubic__g - 12) * OHSunBrinch_inv)
    BR_DGR_up_cubic__rg = DGR_conv_lim_sup * (10 ** (BR_OHBrinch_cubic__rg - 12) * OHSunBrinch_inv)
    BR_DGR_down_cubic__g = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_cubic__g - 12) * OHSunBrinch_inv)
    BR_DGR_down_cubic__rg = DGR_conv_lim_inf * (10 ** (BR_OHBrinch_cubic__rg - 12) * OHSunBrinch_inv)
    BR_DGR_cubic__g = DGR_cte * (10 ** (BR_OHBrinch_cubic__g - 12) * OHSunBrinch_inv)
    BR_DGR_cubic__rg = DGR_cte * (10 ** (BR_OHBrinch_cubic__rg - 12) * OHSunBrinch_inv)
    BR_SigmaGas_up_cubic__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_up_cubic__g
    BR_SigmaGas_up_cubic__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_up_cubic__rg
    BR_SigmaGas_Ha_up_cubic__g = dustdim * H.tau_V_neb__g / BR_DGR_up_cubic__g
    BR_SigmaGas_Ha_up_cubic__rg = dustdim * H.tau_V_neb__rg / BR_DGR_up_cubic__rg
    BR_SigmaGas_down_cubic__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_down_cubic__g
    BR_SigmaGas_down_cubic__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_down_cubic__rg
    BR_SigmaGas_Ha_down_cubic__g = dustdim * H.tau_V_neb__g / BR_DGR_down_cubic__g
    BR_SigmaGas_Ha_down_cubic__rg = dustdim * H.tau_V_neb__rg / BR_DGR_down_cubic__rg
    BR_SigmaGas_cubic__g = dustdim * H.tau_V__Tg[iT] / BR_DGR_cubic__g
    BR_SigmaGas_cubic__rg = dustdim * H.tau_V__Trg[iT] / BR_DGR_cubic__rg
    BR_SigmaGas_Ha_cubic__g = dustdim * H.tau_V_neb__g / BR_DGR_cubic__g
    BR_SigmaGas_Ha_cubic__rg = dustdim * H.tau_V_neb__rg / BR_DGR_cubic__rg
    BR_GSR_up_cubic__g = BR_SigmaGas_up_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_up_cubic__rg = BR_SigmaGas_up_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_up_cubic__g = BR_SigmaGas_Ha_up_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_up_cubic__rg = BR_SigmaGas_Ha_up_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_down_cubic__g = BR_SigmaGas_down_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_down_cubic__rg = BR_SigmaGas_down_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_down_cubic__g = BR_SigmaGas_Ha_down_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_down_cubic__rg = BR_SigmaGas_Ha_down_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_cubic__g = BR_SigmaGas_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_cubic__rg = BR_SigmaGas_cubic__rg / H.McorSD__Trg[iT]
    BR_GSR_Ha_cubic__g = BR_SigmaGas_Ha_cubic__g / H.McorSD__Tg[iT]
    BR_GSR_Ha_cubic__rg = BR_SigmaGas_Ha_cubic__rg / H.McorSD__Trg[iT]
    BR_f_gas_up_cubic__g = 1. / (1. + 1. / BR_GSR_up_cubic__g)
    BR_f_gas_up_cubic__rg = 1. / (1. + 1. / BR_GSR_up_cubic__rg)
    BR_f_gas_Ha_up_cubic__g = 1. / (1. + 1. / BR_GSR_Ha_up_cubic__g)
    BR_f_gas_Ha_up_cubic__rg = 1. / (1. + 1. / BR_GSR_Ha_up_cubic__rg)
    BR_f_gas_down_cubic__g = 1. / (1. + 1. / BR_GSR_down_cubic__g)
    BR_f_gas_down_cubic__rg = 1. / (1. + 1. / BR_GSR_down_cubic__rg)
    BR_f_gas_Ha_down_cubic__g = 1. / (1. + 1. / BR_GSR_Ha_down_cubic__g)
    BR_f_gas_Ha_down_cubic__rg = 1. / (1. + 1. / BR_GSR_Ha_down_cubic__rg)
    BR_f_gas_cubic__g = 1. / (1. + 1. / BR_GSR_cubic__g)
    BR_f_gas_cubic__rg = 1. / (1. + 1. / BR_GSR_cubic__rg)
    BR_f_gas_Ha_cubic__g = 1. / (1. + 1. / BR_GSR_Ha_cubic__g)
    BR_f_gas_Ha_cubic__rg = 1. / (1. + 1. / BR_GSR_Ha_cubic__rg)
    
    if args.dryrun:
        sys.exit('dry run...')

#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#     if args.debug is True:
#         import pandas as pd
#         import seaborn as sns
#         from scipy.stats import spearmanr
#         #sns.set(style = 'white')
#         sns.set_context('talk', font_scale = 1.)
#         
#         # Generate a random correlated bivariate dataset
#         xm, ym = C.ma_mask_xyz(x = H.alogZ_mass__Ug[-1], y = np.ma.log10(BR_f_gas_Ha_cubic__g))
#         x1 = pd.Series(xm.compressed(), name = r'$\langle \log\ Z_\star \rangle_M$ [$Z_\odot$]')
#         x2 = pd.Series(ym.compressed(), name = r'$\log\ f_{gas}$')
#         
#         # Show the joint distribution using kernel density estimation
#         NRows = 1
#         NCols = 4
#         page_size_inches = (NCols * 3.5, NRows * 4)
#         grid_shape = (NRows, NCols)
# 
#         f = plt.figure()
#         f.set_size_inches(page_size_inches)
#         ax = plt.subplot2grid(grid_shape, loc = (0, 0))
#         ax.set_axis_on()
#         sns.kdeplot(x1, x2, cmap = 'Blues', shade = True, shade_lowest = False)
#         sns.rugplot(x1, color = 'k', ax = ax)
#         #sns.rugplot(x2, vertical = True, color = 'k', ax = ax_joint)
#         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         # g = sns.jointplot(
#         #     x1,
#         #     x2,
#         #     space = 0,
#         #     stat_func = spearmanr,
#         #     kind = 'kde',
#         #     joint_kws = dict(s = 5, edgecolor = 'none'),
#         #     annot_kws = dict(stat = r'$R_s$', fontsize = 10),
#         #     marginal_kws = dict(color = ".5"),
#         # )
#         #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
#         rs = runstats(xm.compressed(), ym.compressed(), nBox = 20, smooth = True, sigma = 1.2, overlap = 0.4)
#         ax.plot(rs.xS, rs.yS, '.-', c = 'w')
#         #g = g.plot_joint(sns.kdeplot, n_levels=10, cmap = 'Blues_r')
#         # f = plt.gcf()
#         # f.set_dpi(100)
#         # f.set_size_inches(10, 8)
#         #g.fig.savefig('jointplot.png')
#         f.savefig('jointplot.png')
#         plt.close(f)
#         sys.exit()
#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    ######################
    ######################
    ######################
    
    if args.output is None:
        if args.maskradius is not None:
            output = 'gas_maskRadius%.1f.pdf' % args.maskradius
        else:
            output = 'gas.pdf'
    else:
        output = args.output
          
    default_sc_kwargs = dict(marker = 'o', s = 3, alpha = 0.9, edgecolor = 'none', label = '', vmax = H.Rbinfin, vmin = minR)
    #default_sc_kwargs = dict(marker = 'o', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
    default_rs_kwargs = dict(smooth = True, sigma = 1.2, overlap = 0.4)
    
    with PdfPages(output) as pdf:
        ##########################
        ######## deltaDGR ########
        ##########################
        props = [ 
            'alogSFRSDkpcR',
            'alogSFRSDHakpcR',
            'atfluxR',
            'xYR',
            'logO3N2M13R',
            'alogZmassR',
            'logMcorSDR' ,
            'logtauVR',
            'logtauVNebR' 
        ]
        txt_suptitle = r'$\Longrightarrow$ %s  NGals:%d  tSF:%.2f Myr  $x_Y$(min):%.0f%%  $\tau_V^\star $(min):%.2f  $\tau_V^{neb}$(min):%.2f  $\epsilon\tau_V^{neb}$(max):%.2f ' % (gals_txt, N_gals, (H.tSF__T[iT] / 1e6), H.xOkMin * 100., H.tauVOkMin, H.tauVNebOkMin, H.tauVNebErrMax)
        
        for prop_key in props:
            _, prop_dict = H.get_plot_dict(iT, -1, prop_key)
            x = prop_dict['v']
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            # norm = .5 * (x[10, :] + x[9, :])
            # x_norm = x/norm
            #EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
            x_norm = x 
    
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
    
            NRows = 1
            NCols = 4
            page_size_inches = (NCols * 3.5, NRows * 4)
            grid_shape = (NRows, NCols)
    
            f = plt.figure()
            f.set_size_inches(page_size_inches)
            f.suptitle(txt_suptitle, fontsize = 10)
            ax = plt.subplot2grid(grid_shape, loc = (0, 0))
            ax.set_axis_on()
            
            tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
            
            y = BR_DGR_ols__rg
            #y_norm = y / (0.5 * (y[10, :] + y[9, :]))
            y_norm = y
            xm, ym_BR_DGR_ols = C.ma_mask_xyz(x = x_norm,  y = np.ma.log10(y_norm), mask = tmp_mask)
            rs_BR_DGR_ols = runstats(xm.compressed(), ym_BR_DGR_ols.compressed(), nBox = 20, **rs_kwargs)

            y = BR_DGR_up_ols__rg
            #y_norm = y / (0.5 * (y[10, :] + y[9, :]))
            y_norm = y
            xm, ym_BR_DGR_up_ols = C.ma_mask_xyz(x = x_norm,  y = np.ma.log10(y_norm), mask = tmp_mask)
            rs_BR_DGR_up_ols = runstats(xm.compressed(), ym_BR_DGR_up_ols.compressed(), nBox = 20, **rs_kwargs)
            
            y = BR_DGR_down_ols__rg
            #y_norm = y / (0.5 * (y[10, :] + y[9, :]))
            y_norm = y
            xm, ym_BR_DGR_down_ols = C.ma_mask_xyz(x = x_norm,  y = np.ma.log10(y_norm), mask = tmp_mask)
            rs_BR_DGR_down_ols = runstats(xm.compressed(), ym_BR_DGR_down_ols.compressed(), nBox = 20, **rs_kwargs)
            
            y = SK_DGR__rg
            #y_norm = y / (0.5 * (y[10, :] + y[9, :]))
            y_norm = y
            xm, ym_SK_DGR = C.ma_mask_xyz(x = x_norm,  y = np.ma.log10(y_norm), mask = tmp_mask)
            rs_SK_DGR = runstats(xm.compressed(), ym_SK_DGR.compressed(), nBox = 20, **rs_kwargs)
            
            y = SK_DGR_Ha__rg
            #y_norm = y / (0.5 * (y[10, :] + y[9, :]))
            y_norm = y
            xm, ym_SK_DGR_Ha = C.ma_mask_xyz(x = x_norm,  y = np.ma.log10(y_norm), mask = tmp_mask)
            rs_SK_DGR_Ha = runstats(xm.compressed(), ym_SK_DGR_Ha.compressed(), nBox = 20, **rs_kwargs)
            
            ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'b')
            ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'g')
            ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'black')
            ax.fill_between(rs_BR_DGR_up_ols.xS, rs_BR_DGR_up_ols.yS, rs_BR_DGR_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            ax.set_xlim(prop_dict['lim'])
            #ax.set_ylim(0,15e-3)
            ax.set_ylim(-3.2, -1.8)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.grid()
            #ax.set_ylabel(r'$\frac{\log\ \delta_{DGR}}{\log\ \delta_{DGR}(@1HLR)}$')
            ax.set_ylabel(r'$\log\ \delta_{DGR}$')
            ax.set_xlabel(prop_dict['label'])
            
            ax = plt.subplot2grid(grid_shape, loc = (0, 1))
            ax.set_axis_on()
            ax.plot(rs_BR_DGR_ols.xS, rs_BR_DGR_ols.yS, '.-', c = 'k')
            plot_text_topright_color(ax, 'BR', c = 'b', fs = 11)
            #ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, marker = 'o', c = 'b', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_BR_DGR_ols.mask).compressed()
            ax.scatter(rs_BR_DGR_ols.x, rs_BR_DGR_ols.y, c = c, **sc_kwargs)
            ax.set_xlim(prop_dict['lim'])
            #ax.set_ylim(0,15e-3)
            ax.set_ylim(-3.2, -1.8)
            ax.xaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (0, 2))
            ax.set_axis_on()
            ax.plot(rs_SK_DGR.xS, rs_SK_DGR.yS, '.-', c = 'k')
            txt = r'SK syn'
            plot_text_topright_color(ax, txt, c = 'g', fs = 11)
            #ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, marker = 'o', c = 'g', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_SK_DGR.mask).compressed()
            ax.scatter(rs_SK_DGR.x, rs_SK_DGR.y, c = c, **sc_kwargs)
            ax.set_xlim(prop_dict['lim'])
            #ax.set_ylim(0,15e-3)
            ax.set_ylim(-3.2, -1.8)
            ax.xaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])
    
            ax = plt.subplot2grid(grid_shape, loc = (0, 3))
            ax.set_axis_on()
            ax.plot(rs_SK_DGR_Ha.xS, rs_SK_DGR_Ha.yS, '.-', c = 'k')
            txt = r'SK Neb'
            plot_text_topright_color(ax, txt, c = 'k', fs = 11)
            #ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, marker = 'o', c = 'black', s = 1, edgecolor = 'none', alpha = 0.8, label = '')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_SK_DGR_Ha.mask).compressed()
            sc = ax.scatter(rs_SK_DGR_Ha.x, rs_SK_DGR_Ha.y, c = c, **sc_kwargs)
            ax.set_xlim(prop_dict['lim'])
            #ax.set_ylim(0,15e-3)
            ax.set_ylim(-3.2, -1.8)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])
            cb = f.colorbar(sc)
            cb.set_label('R [HLR]')
            
            f.subplots_adjust(bottom = 0.15, wspace = 0.0, right = 0.9)
            pdf.savefig(f)
            #f.savefig('DGR_%s_%s' % (prop_key, output))
            plt.close(f)

        ##########################
        ######## SigmaGas ########
        ##########################
        for prop_key in props:
            _, prop_dict = H.get_plot_dict(iT, -1, prop_key)
            x = prop_dict['v']
    
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
    
            NRows = 2
            NCols = 3
            page_size_inches = (NCols * 3.5, NRows * 4)
            grid_shape = (NRows, NCols)
    
            f = plt.figure()
            f.set_size_inches(page_size_inches)
            f.suptitle(txt_suptitle, fontsize = 10)
            ax = plt.subplot2grid(grid_shape, loc = (0, 0))
            ax.set_axis_on()
            
            tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
            
            xm, ym_BR_SigmaGas_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_ols = runstats(xm.compressed(), ym_BR_SigmaGas_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_up_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_up_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_up_ols = runstats(xm.compressed(), ym_BR_SigmaGas_up_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_down_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_down_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_down_ols = runstats(xm.compressed(), ym_BR_SigmaGas_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_Ha_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_Ha_ols = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_SigmaGas_Ha_up_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_up_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_Ha_up_ols = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_up_ols.compressed(), nBox = 20, **rs_kwargs)            
            xm, ym_BR_SigmaGas_Ha_down_ols = C.ma_mask_xyz(x = x, y = np.log10(BR_SigmaGas_Ha_down_ols__rg), mask = tmp_mask)
            rs_BR_SigmaGas_Ha_down_ols = runstats(xm.compressed(), ym_BR_SigmaGas_Ha_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_SigmaGas = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas__rg), mask = tmp_mask)
            rs_SK_SigmaGas = runstats(xm.compressed(), ym_SK_SigmaGas.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_SigmaGas_Ha = C.ma_mask_xyz(x = x, y = np.log10(SK_SigmaGas_Ha__rg), mask = tmp_mask)
            rs_SK_SigmaGas_Ha = runstats(xm.compressed(), ym_SK_SigmaGas_Ha.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_RR_SigmaGas = C.ma_mask_xyz(x = x, y = np.log10(RR_SigmaGas__rg), mask = tmp_mask)
            rs_RR_SigmaGas = runstats(xm.compressed(), ym_RR_SigmaGas.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_BR_SigmaGas_ols.xS, rs_BR_SigmaGas_ols.yS, '.-', c = 'b')
            ax.plot(rs_BR_SigmaGas_Ha_ols.xS, rs_BR_SigmaGas_Ha_ols.yS, '.-', c = 'r')
            ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'g')
            ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'black')
            ax.plot(rs_RR_SigmaGas.xS, rs_RR_SigmaGas.yS, '.-', c = 'y')
            ax.fill_between(rs_BR_SigmaGas_up_ols.xS, rs_BR_SigmaGas_up_ols.yS, rs_BR_SigmaGas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            ax.fill_between(rs_BR_SigmaGas_Ha_up_ols.xS, rs_BR_SigmaGas_Ha_up_ols.yS, rs_BR_SigmaGas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0.4, 2)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.grid()
            ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (0, 1))
            ax.plot(rs_BR_SigmaGas_ols.xS, rs_BR_SigmaGas_ols.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_BR_SigmaGas_ols.mask).compressed()
            ax.scatter(rs_BR_SigmaGas_ols.x, rs_BR_SigmaGas_ols.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'BR Syn', c = 'b', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0.4, 2)
            ax.xaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (0, 2))
            ax.plot(rs_BR_SigmaGas_Ha_ols.xS, rs_BR_SigmaGas_Ha_ols.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_BR_SigmaGas_Ha_ols.mask).compressed()
            ax.scatter(rs_BR_SigmaGas_Ha_ols.x, rs_BR_SigmaGas_Ha_ols.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'BR Neb', c = 'r', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0.4, 2)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (1, 0))
            ax.plot(rs_SK_SigmaGas.xS, rs_SK_SigmaGas.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_SK_SigmaGas.mask).compressed()
            ax.scatter(rs_SK_SigmaGas.x, rs_SK_SigmaGas.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'SK Syn', c = 'g', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0.4, 2)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.grid()
            ax.set_xlabel(prop_dict['label'])
            ax.set_ylabel(r'$\log\ \Sigma_{gas}$ [M${}_\odot$ pc${}^{-2}$]')

            ax = plt.subplot2grid(grid_shape, loc = (1, 1))
            ax.plot(rs_SK_SigmaGas_Ha.xS, rs_SK_SigmaGas_Ha.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_SK_SigmaGas_Ha.mask).compressed()
            ax.scatter(rs_SK_SigmaGas_Ha.x, rs_SK_SigmaGas_Ha.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'SK Neb', c = 'k', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0.4, 2)
            ax.xaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (1, 2))
            ax.plot(rs_RR_SigmaGas.xS, rs_RR_SigmaGas.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_RR_SigmaGas.mask).compressed()
            sc = ax.scatter(rs_RR_SigmaGas.x, rs_RR_SigmaGas.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'RR', c = 'y', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0.4, 2)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])
            cb = f.colorbar(sc)
            cb.set_label('R [HLR]')
            
            f.subplots_adjust(bottom = 0.15, wspace = 0.0, right = 0.9)
            pdf.savefig(f)
            #f.savefig('DGR_%s_%s' % (prop_key, output))
            plt.close(f)

        ##########################
        ########## fGas ##########
        ##########################
        for prop_key in props:
            _, prop_dict = H.get_plot_dict(iT, -1, prop_key)
            x = prop_dict['v']
    
            rs_kwargs = default_rs_kwargs.copy()
            sc_kwargs = default_sc_kwargs.copy()
    
            NRows = 2
            NCols = 3
            page_size_inches = (NCols * 3.5, NRows * 4)
            grid_shape = (NRows, NCols)
    
            f = plt.figure()
            f.set_size_inches(page_size_inches)
            f.suptitle(txt_suptitle, fontsize = 10)
            ax = plt.subplot2grid(grid_shape, loc = (0, 0))
            ax.set_axis_on()
            
            tmp_mask = ~(maskRadiusOk__rg & gals_slice__rg)
            
            xm, ym_BR_f_gas_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_ols = runstats(xm.compressed(), ym_BR_f_gas_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_up_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_up_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_up_ols = runstats(xm.compressed(), ym_BR_f_gas_up_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_down_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_down_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_down_ols = runstats(xm.compressed(), ym_BR_f_gas_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_Ha_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_Ha_ols = runstats(xm.compressed(), ym_BR_f_gas_Ha_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_BR_f_gas_Ha_up_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_up_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_Ha_up_ols = runstats(xm.compressed(), ym_BR_f_gas_Ha_up_ols.compressed(), nBox = 20, **rs_kwargs)            
            xm, ym_BR_f_gas_Ha_down_ols = C.ma_mask_xyz(x = x, y = BR_f_gas_Ha_down_ols__rg, mask = tmp_mask)
            rs_BR_f_gas_Ha_down_ols = runstats(xm.compressed(), ym_BR_f_gas_Ha_down_ols.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_f_gas = C.ma_mask_xyz(x = x, y = SK_f_gas__rg, mask = tmp_mask)
            rs_SK_f_gas = runstats(xm.compressed(), ym_SK_f_gas.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_SK_f_gas_Ha = C.ma_mask_xyz(x = x, y = SK_f_gas_Ha__rg, mask = tmp_mask)
            rs_SK_f_gas_Ha = runstats(xm.compressed(), ym_SK_f_gas_Ha.compressed(), nBox = 20, **rs_kwargs)
            xm, ym_RR_f_gas = C.ma_mask_xyz(x = x, y = RR_f_gas__rg, mask = tmp_mask)
            rs_RR_f_gas = runstats(xm.compressed(), ym_RR_f_gas.compressed(), nBox = 20, **rs_kwargs)
            ax.plot(rs_BR_f_gas_ols.xS, rs_BR_f_gas_ols.yS, '.-', c = 'b')
            ax.plot(rs_BR_f_gas_Ha_ols.xS, rs_BR_f_gas_Ha_ols.yS, '.-', c = 'r')
            ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'g')
            ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'black')
            ax.plot(rs_RR_f_gas.xS, rs_RR_f_gas.yS, '.-', c = 'y')
            ax.fill_between(rs_BR_f_gas_up_ols.xS, rs_BR_f_gas_up_ols.yS, rs_BR_f_gas_down_ols.yS, edgecolor = 'k', facecolor = 'b', alpha = 0.4)            
            ax.fill_between(rs_BR_f_gas_Ha_up_ols.xS, rs_BR_f_gas_Ha_up_ols.yS, rs_BR_f_gas_Ha_down_ols.yS, edgecolor = 'k', facecolor = 'r', alpha = 0.4)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0, 0.3)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.grid()
            ax.set_ylabel(r'$f_{gas}$')
            ax.set_xlabel(prop_dict['label'])
            
            ax = plt.subplot2grid(grid_shape, loc = (0, 1))
            ax.plot(rs_BR_f_gas_ols.xS, rs_BR_f_gas_ols.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_BR_f_gas_ols.mask).compressed()
            ax.scatter(rs_BR_f_gas_ols.x, rs_BR_f_gas_ols.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'BR Syn', c = 'b', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0, 0.3)
            ax.xaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (0, 2))
            ax.plot(rs_BR_f_gas_Ha_ols.xS, rs_BR_f_gas_Ha_ols.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_BR_f_gas_Ha_ols.mask).compressed()
            ax.scatter(rs_BR_f_gas_Ha_ols.x, rs_BR_f_gas_Ha_ols.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'BR Neb', c = 'r', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0, 0.3)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (1, 0))
            ax.plot(rs_SK_f_gas.xS, rs_SK_f_gas.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_SK_f_gas.mask).compressed()
            ax.scatter(rs_SK_f_gas.x, rs_SK_f_gas.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'SK Syn', c = 'g', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0, 0.3)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            ax.grid()
            ax.set_xlabel(prop_dict['label'])
            ax.set_ylabel(r'$f_{gas}$')

            ax = plt.subplot2grid(grid_shape, loc = (1, 1))
            ax.plot(rs_SK_f_gas_Ha.xS, rs_SK_f_gas_Ha.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_SK_f_gas_Ha.mask).compressed()
            ax.scatter(rs_SK_f_gas_Ha.x, rs_SK_f_gas_Ha.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'SK Neb', c = 'k', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0, 0.3)
            ax.xaxis.set_major_locator(MaxNLocator(4, prune = 'both'))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])

            ax = plt.subplot2grid(grid_shape, loc = (1, 2))
            ax.plot(rs_RR_f_gas.xS, rs_RR_f_gas.yS, '.-', c = 'k')
            c = np.ma.masked_array(H.Rtoplot(), mask = ym_RR_f_gas.mask).compressed()
            sc = ax.scatter(rs_RR_f_gas.x, rs_RR_f_gas.y, c = c, **sc_kwargs)
            plot_text_topright_color(ax, 'RR', c = 'y', fs = 11)
            ax.set_xlim(prop_dict['lim'])
            ax.set_ylim(0, 0.3)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            plt.setp(ax.get_yticklabels(), visible = False)
            ax.grid()
            ax.set_xlabel(prop_dict['label'])
            cb = f.colorbar(sc)
            cb.set_label('R [HLR]')
            
            f.subplots_adjust(bottom = 0.15, wspace = 0.0, right = 0.9)
            pdf.savefig(f)
            #f.savefig('DGR_%s_%s' % (prop_key, output))
            plt.close(f)
