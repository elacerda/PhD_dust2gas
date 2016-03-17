import tables as tbl
import numpy as np

class galaxy(tbl.IsDescription):
    id                          = tbl.UInt16Col(pos = 1)
    name                        = tbl.StringCol(20, pos = 2)
    califaID                    = tbl.StringCol(5, pos = 3)
    N_zone                      = tbl.UInt16Col(pos = 4)
    distance_Mpc                = tbl.Float64Col(pos = 5)
    redshift                    = tbl.Float64Col(pos = 6)
    m_type                      = tbl.Int32Col(pos = 7) 
    ba                          = tbl.Float64Col(pos = 8)
    ba_PyCASSO                  = tbl.Float64Col(pos = 9)
    ParsecPerPixel              = tbl.Float64Col(pos = 10)
    Mr                          = tbl.Float64Col(pos = 11)
    ur                          = tbl.Float64Col(pos = 12)
    
    # Syn
    HLR_pix                     = tbl.Float64Col(pos = 13)
    HMR_pix                     = tbl.Float64Col(pos = 14)
    
class zone(tbl.IsDescription):
    id_gal                      = tbl.UInt16Col(pos = 1)
    id                          = tbl.Int64Col(pos = 2, dflt = -1)
    distance_HLR                = tbl.Float64Col(pos = 3, dflt = np.nan)
    area_pc2                    = tbl.Float64Col(pos = 4, dflt = np.nan)
    
    # Syn
    Mcor                        = tbl.Float64Col(pos = 5, dflt = np.nan)
    tau_V                       = tbl.Float64Col(pos = 6, dflt = np.nan)
    McorSD                      = tbl.Float64Col(pos = 7, dflt = np.nan)
    at_flux                     = tbl.Float64Col(pos = 8, dflt = np.nan)
    at_mass                     = tbl.Float64Col(pos = 9, dflt = np.nan)
    
    # Neb
    tau_V_neb                   = tbl.Float64Col(pos = 10, dflt = np.nan)
    etau_V_neb                  = tbl.Float64Col(pos = 11, dflt = np.nan)
    F_obs_Hb                    = tbl.Float64Col(pos = 12, dflt = np.nan)
    F_obs_O3                    = tbl.Float64Col(pos = 13, dflt = np.nan)
    F_obs_Ha                    = tbl.Float64Col(pos = 14, dflt = np.nan)
    F_obs_N2                    = tbl.Float64Col(pos = 15, dflt = np.nan)
    eF_obs_Hb                   = tbl.Float64Col(pos = 16, dflt = np.nan)
    eF_obs_O3                   = tbl.Float64Col(pos = 17, dflt = np.nan)
    eF_obs_Ha                   = tbl.Float64Col(pos = 18, dflt = np.nan)
    eF_obs_N2                   = tbl.Float64Col(pos = 19, dflt = np.nan)
    baseline_Hb                 = tbl.Float64Col(pos = 20, dflt = np.nan)
    baseline_O3                 = tbl.Float64Col(pos = 21, dflt = np.nan)
    baseline_Ha                 = tbl.Float64Col(pos = 22, dflt = np.nan)
    baseline_N2                 = tbl.Float64Col(pos = 23, dflt = np.nan)
    EW_Hb                       = tbl.Float64Col(pos = 24, dflt = np.nan)
    EW_O3                       = tbl.Float64Col(pos = 25, dflt = np.nan)
    EW_Ha                       = tbl.Float64Col(pos = 26, dflt = np.nan)
    EW_N2                       = tbl.Float64Col(pos = 27, dflt = np.nan)
    sigma_Hb                    = tbl.Float64Col(pos = 28, dflt = np.nan)
    sigma_O3                    = tbl.Float64Col(pos = 29, dflt = np.nan)
    sigma_Ha                    = tbl.Float64Col(pos = 30, dflt = np.nan)
    sigma_N2                    = tbl.Float64Col(pos = 31, dflt = np.nan)
    esigma_Hb                   = tbl.Float64Col(pos = 32, dflt = np.nan)
    esigma_O3                   = tbl.Float64Col(pos = 33, dflt = np.nan)
    esigma_Ha                   = tbl.Float64Col(pos = 34, dflt = np.nan)
    esigma_N2                   = tbl.Float64Col(pos = 35, dflt = np.nan)
    pos_Hb                      = tbl.Float64Col(pos = 36, dflt = np.nan)
    pos_O3                      = tbl.Float64Col(pos = 37, dflt = np.nan)
    pos_Ha                      = tbl.Float64Col(pos = 38, dflt = np.nan)
    pos_N2                      = tbl.Float64Col(pos = 39, dflt = np.nan)
    epos_Hb                     = tbl.Float64Col(pos = 40, dflt = np.nan)
    epos_O3                     = tbl.Float64Col(pos = 41, dflt = np.nan)
    epos_Ha                     = tbl.Float64Col(pos = 42, dflt = np.nan)
    epos_N2                     = tbl.Float64Col(pos = 43, dflt = np.nan)
    logOH                       = tbl.Float64Col(pos = 44, dflt = np.nan)

class tSF(tbl.IsDescription):
    id                          = tbl.UInt8Col(pos = 1)
    age                         = tbl.Float64Col(pos = 2)
    
class zone_neb(tbl.IsDescription):
    id_gal                      = tbl.UInt16Col(pos = 1)
    id_zone                     = tbl.Int16Col(pos = 2, dflt = -1)
    L_obs_Ha                    = tbl.Float64Col(pos = 3, dflt = np.nan)
    L_int_Ha                    = tbl.Float64Col(pos = 4, dflt = np.nan)
    SFR                         = tbl.Float64Col(pos = 5, dflt = np.nan)
    SFRSD                       = tbl.Float64Col(pos = 6, dflt = np.nan)
    flag_BPT                    = tbl.UInt8Col(pos = 7, dflt = 0)
    flag_WHAN                   = tbl.UInt8Col(pos = 8, dflt = 0)
    
class zone_SF(tbl.IsDescription):
    id_gal                      = tbl.UInt16Col(pos = 1)
    id_zone                     = tbl.Int16Col(pos = 2, dflt = -1)
    id_tSF                      = tbl.UInt8Col(pos = 3)
    xY                          = tbl.Float64Col(pos = 4, dflt = np.nan)
    SFR                         = tbl.Float64Col(pos = 5, dflt = np.nan)
    SFRSD                       = tbl.Float64Col(pos = 6, dflt = np.nan)