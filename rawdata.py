import CALIFAUtils as C
from tables import *
import numpy as np
import sys

class galaxy(IsDescription):
    id                          = UInt16Col(pos = 1)
    name                        = StringCol(20, pos = 2)
    califaID                    = StringCol(5, pos = 3)
    N_zones                     = UInt16Col(pos = 4)
    distance_Mpc                = Float64Col(pos = 5)
    redshift                    = Float64Col(pos = 6)
    m_type                      = Int32Col(pos = 7) 
    ba                          = Float64Col(pos = 8)
    ba_PyCASSO                  = Float64Col(pos = 9)
    ParsecPerPixel              = Float64Col(pos = 10)
    Mr                          = Float64Col(pos = 11)
    ur                          = Float64Col(pos = 12)
    
    # Syn
    HLR_pix                     = Float64Col(pos = 13)
    HMR_pix                     = Float64Col(pos = 14)
    integrated_tau_V            = Float64Col(pos = 15)
    at_flux_GAL                 = Float64Col(pos = 16)
    
    # Neb
    integrated_tau_V_neb        = Float64Col(pos = 17)
    integrated_etau_V_neb       = Float64Col(pos = 18)
    integrated_F_obs_Hb         = Float64Col(pos = 19)
    integrated_F_obs_O3         = Float64Col(pos = 20)
    integrated_F_obs_Ha         = Float64Col(pos = 21)
    integrated_F_obs_N2         = Float64Col(pos = 22)
    integrated_eF_obs_Hb        = Float64Col(pos = 23)
    integrated_eF_obs_O3        = Float64Col(pos = 24)
    integrated_eF_obs_Ha        = Float64Col(pos = 25)
    integrated_eF_obs_N2        = Float64Col(pos = 26)
    integrated_baseline_Hb      = Float64Col(pos = 27)
    integrated_baseline_O3      = Float64Col(pos = 28)
    integrated_baseline_Ha      = Float64Col(pos = 29)
    integrated_baseline_N2      = Float64Col(pos = 30)
    integrated_EW_Hb            = Float64Col(pos = 31)
    integrated_EW_O3            = Float64Col(pos = 32)
    integrated_EW_Ha            = Float64Col(pos = 33)
    integrated_EW_N2            = Float64Col(pos = 34)

class zone(IsDescription):
    id_gal                  = UInt16Col(pos = 1)
    id_zone                 = UInt16Col(pos = 2)
    distance_HLR            = Float64Col(pos = 3)
    area_pc2                = Float64Col(pos = 4)
    
    # Syn
    Mcor                    = Float64Col(pos = 5)
    tau_V                   = Float64Col(pos = 6)
    McorSD                  = Float64Col(pos = 7)
    at_flux                 = Float64Col(pos = 8)
    at_mass                 = Float64Col(pos = 9)
    
    # Neb
    tau_V_neb               = Float64Col(pos = 10)
    etau_V_neb              = Float64Col(pos = 11)
    F_obs__Hb               = Float64Col(pos = 12)
    F_obs__O3               = Float64Col(pos = 13)
    F_obs__Ha               = Float64Col(pos = 14)
    F_obs__N2               = Float64Col(pos = 15)
    eF_obs__Hb              = Float64Col(pos = 16)
    eF_obs__O3              = Float64Col(pos = 17)
    eF_obs__Ha              = Float64Col(pos = 18)
    eF_obs__N2              = Float64Col(pos = 19)
    baseline_Hb             = Float64Col(pos = 20)
    baseline_O3             = Float64Col(pos = 21)
    baseline_Ha             = Float64Col(pos = 22)
    baseline_N2             = Float64Col(pos = 23)
    EW_Hb                   = Float64Col(pos = 24)
    EW_O3                   = Float64Col(pos = 25)
    EW_Ha                   = Float64Col(pos = 26)
    EW_N2                   = Float64Col(pos = 27)
    logOH                   = Float64Col(pos = 28)
    
h5file = open_file('teste.h5', mode =  'w', title = 'Test file')
group = h5file.create_group('/', 'gal', 'Galaxy data')
table_main = h5file.create_table(group, 'main', galaxy, 'Main galaxy data')
table_zone = h5file.create_table(group, 'zones', zone, 'Zone data')

gals, _ = C.sort_gals('/Users/lacerda/CALIFA/listv20_q050.d15a.txt', order = 1)

r = table_main.row
for iGal, gal in enumerate(gals):
    K = C.read_one_cube(gal, EL = True, GP = True, v_run = 'last', debug = True)
    pa, ba = K.getEllipseParams()
    K.setGeometry(pa, ba)

    Hb_central_wl = '4861'
    O3_central_wl = '5007'
    Ha_central_wl = '6563'
    N2_central_wl = '6583'
    i_Hb = K.EL.lines.index(Hb_central_wl)
    i_O3 = K.EL.lines.index(O3_central_wl)
    i_Ha = K.EL.lines.index(Ha_central_wl)
    i_N2 = K.EL.lines.index(N2_central_wl)
    
    numerator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0) * K.at_flux__z
    denominator__z = K.Lobn__tZz.sum(axis = 1).sum(axis = 0)
    at_flux_GAL = numerator__z.sum() / denominator__z.sum()
    AVtoTauV = 1. / (np.log10(np.exp(1)) / 0.4)

    main_data = [(
        iGal, K.galaxyName, K.califaID, K.N_zone, K.distance_Mpc, K.redshift, 
        C.my_morf(C.get_morfologia(K.califaID)[0]), np.float(K.masterListData['ba']), 
        ba, K.parsecPerPixel, np.float(K.masterListData['Mr']), 
        np.float(K.masterListData['u-r']), 
        K.HLR_pix, K.getHalfRadius(K.McorSD__yx), 
        K.integrated_keywords['A_V'] * AVtoTauV, at_flux_GAL,
        K.EL.integrated_tau_V_neb, 
        K.EL.integrated_tau_V_neb_err,
        K.EL.integrated_flux[i_Hb],
        K.EL.integrated_flux[i_O3],
        K.EL.integrated_flux[i_Ha],
        K.EL.integrated_flux[i_N2],
        K.EL.integrated_eflux[i_Hb],
        K.EL.integrated_eflux[i_O3],
        K.EL.integrated_eflux[i_Ha],
        K.EL.integrated_eflux[i_N2],
        K.EL.integrated_baseline[i_Hb],
        K.EL.integrated_baseline[i_O3],
        K.EL.integrated_baseline[i_Ha],
        K.EL.integrated_baseline[i_N2],
        K.EL.integrated_EW[i_Hb],
        K.EL.integrated_EW[i_O3],
        K.EL.integrated_EW[i_Ha],
        K.EL.integrated_EW[i_N2],
    )]
    table_main.append(main_data)

    zone_data = zip(
        np.zeros((K.N_zone), dtype = int) + iGal,
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
        K.EL.Zneb_M13__z,
    )
    
    table_zone.append(zone_data)
    table_zone.flush()
    table_main.flush()
    
h5file.close()