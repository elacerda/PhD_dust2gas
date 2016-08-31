#
# Lacerda@Saco - 15/Mar/2016
#
import numpy as np
import tables as tbl


class galaxy(tbl.IsDescription):
    id = tbl.UInt16Col(pos=1)
    name = tbl.StringCol(20, pos=2)
    califaID = tbl.StringCol(5, pos=3)
    N_zone = tbl.UInt16Col(pos=4)
    distance_Mpc = tbl.Float64Col(pos=5)
    redshift = tbl.Float64Col(pos=6)
    m_type = tbl.Int32Col(pos=7)
    m_type_orig = tbl.StringCol(4, pos=8)
    ba = tbl.Float64Col(pos=9)
    ba_PyCASSO = tbl.Float64Col(pos=10)
    ParsecPerPixel = tbl.Float64Col(pos=11)
    Mr = tbl.Float64Col(pos=12)
    ur = tbl.Float64Col(pos=13)

    # Syn
    HLR_pix = tbl.Float64Col(pos=14)
    HMR_pix = tbl.Float64Col(pos=15)

    pycasso_cube_filename = tbl.StringCol(150, pos=16)
    eml_cube_filename = tbl.StringCol(150, pos=17)
    gasprop_cube_filename = tbl.StringCol(150, pos=18)


class zone(tbl.IsDescription):
    id_gal = tbl.UInt16Col(pos=1)
    id = tbl.Int64Col(pos=2, dflt=-1)
    i_zone = tbl.Int64Col(pos=3, dflt=-1)
    distance_HLR = tbl.Float64Col(pos=4, dflt=np.nan)
    area_pc2 = tbl.Float64Col(pos=5, dflt=np.nan)

    # Syn
    Mcor = tbl.Float64Col(pos=6, dflt=np.nan)
    tau_V = tbl.Float64Col(pos=7, dflt=np.nan)
    McorSD = tbl.Float64Col(pos=8, dflt=np.nan)
    at_flux = tbl.Float64Col(pos=9, dflt=np.nan)
    at_mass = tbl.Float64Col(pos=10, dflt=np.nan)

    # Neb
    tau_V_neb = tbl.Float64Col(pos=11, dflt=np.nan)
    etau_V_neb = tbl.Float64Col(pos=12, dflt=np.nan)
    F_obs_Hb = tbl.Float64Col(pos=13, dflt=np.nan)
    F_obs_O3 = tbl.Float64Col(pos=14, dflt=np.nan)
    F_obs_Ha = tbl.Float64Col(pos=15, dflt=np.nan)
    F_obs_N2 = tbl.Float64Col(pos=16, dflt=np.nan)
    eF_obs_Hb = tbl.Float64Col(pos=17, dflt=np.nan)
    eF_obs_O3 = tbl.Float64Col(pos=18, dflt=np.nan)
    eF_obs_Ha = tbl.Float64Col(pos=19, dflt=np.nan)
    eF_obs_N2 = tbl.Float64Col(pos=20, dflt=np.nan)
    baseline_Hb = tbl.Float64Col(pos=21, dflt=np.nan)
    baseline_O3 = tbl.Float64Col(pos=22, dflt=np.nan)
    baseline_Ha = tbl.Float64Col(pos=23, dflt=np.nan)
    baseline_N2 = tbl.Float64Col(pos=24, dflt=np.nan)
    EW_Hb = tbl.Float64Col(pos=25, dflt=np.nan)
    EW_O3 = tbl.Float64Col(pos=26, dflt=np.nan)
    EW_Ha = tbl.Float64Col(pos=27, dflt=np.nan)
    EW_N2 = tbl.Float64Col(pos=28, dflt=np.nan)
    sigma_Hb = tbl.Float64Col(pos=29, dflt=np.nan)
    sigma_O3 = tbl.Float64Col(pos=30, dflt=np.nan)
    sigma_Ha = tbl.Float64Col(pos=31, dflt=np.nan)
    sigma_N2 = tbl.Float64Col(pos=32, dflt=np.nan)
    esigma_Hb = tbl.Float64Col(pos=33, dflt=np.nan)
    esigma_O3 = tbl.Float64Col(pos=34, dflt=np.nan)
    esigma_Ha = tbl.Float64Col(pos=35, dflt=np.nan)
    esigma_N2 = tbl.Float64Col(pos=36, dflt=np.nan)
    pos_O3 = tbl.Float64Col(pos=38, dflt=np.nan)
    pos_Hb = tbl.Float64Col(pos=37, dflt=np.nan)
    pos_Ha = tbl.Float64Col(pos=39, dflt=np.nan)
    pos_N2 = tbl.Float64Col(pos=40, dflt=np.nan)
    epos_Hb = tbl.Float64Col(pos=41, dflt=np.nan)
    epos_O3 = tbl.Float64Col(pos=42, dflt=np.nan)
    epos_Ha = tbl.Float64Col(pos=43, dflt=np.nan)
    epos_N2 = tbl.Float64Col(pos=44, dflt=np.nan)
    logOH = tbl.Float64Col(pos=45, dflt=np.nan)

    flag_RGB = tbl.UInt8Col(pos=46, dflt=0)
    flag_residual = tbl.UInt8Col(pos=47, dflt=0)


# XXX TODO: XXX
# XXX TODO: XXX
# XXX TODO: XXX
# class radial_profile(tbl.IsDescription):
#     id_gal = tbl.UInt8Col(pos=1)
#     id = tbl.UInt8Col(pos=2)
#     N_bins = tbl.UInt8Col(pos=3, dflt=0)
#     zones_str = tbl.UInt8Col(pos=4)
#     mean_str = tbl.StringCol(pos=5)
#     dtype = tbl.StringCol(8, pos=6)
#     bin_r__str = tbl.StringCol(pos=7)
#     bin_center_r__str = tbl.StringCol(pos=8)


class tSF(tbl.IsDescription):
    id = tbl.UInt8Col(pos=1)
    age = tbl.Float64Col(pos=2)


class tZ(tbl.IsDescription):
    id = tbl.UInt8Col(pos=1)
    age = tbl.Float64Col(pos=2)


class zone_neb(tbl.IsDescription):
    id_gal = tbl.UInt16Col(pos=1)
    id_zone = tbl.Int64Col(pos=2, dflt=-1)
    i_zone = tbl.Int64Col(pos=3, dflt=-1)
    L_obs_Ha = tbl.Float64Col(pos=4, dflt=np.nan)
    L_int_Ha = tbl.Float64Col(pos=5, dflt=np.nan)
    SFR = tbl.Float64Col(pos=6, dflt=np.nan)
    SFRSD = tbl.Float64Col(pos=7, dflt=np.nan)
    flag_BPT = tbl.UInt8Col(pos=8, dflt=0)
    flag_WHAN = tbl.UInt8Col(pos=9, dflt=0)
    flag_tau_V_neb = tbl.UInt8Col(pos=10, dflt=0)
    flag_etau_V_neb = tbl.UInt8Col(pos=11, dflt=0)


class zone_SF(tbl.IsDescription):
    id_gal = tbl.UInt16Col(pos=1)
    id_zone = tbl.Int64Col(pos=2, dflt=-1)
    i_zone = tbl.Int64Col(pos=3, dflt=-1)
    id_tSF = tbl.UInt8Col(pos=4)
    xY = tbl.Float64Col(pos=5, dflt=np.nan)
    SFR = tbl.Float64Col(pos=6, dflt=np.nan)
    SFRSD = tbl.Float64Col(pos=7, dflt=np.nan)
    flag_xY = tbl.UInt8Col(pos=8, dflt=0)
    flag_tau_V = tbl.UInt8Col(pos=9, dflt=0)


class zone_Z(tbl.IsDescription):
    id_gal = tbl.UInt16Col(pos=1)
    id_zone = tbl.Int64Col(pos=2, dflt=-1)
    i_zone = tbl.Int64Col(pos=3, dflt=-1)
    id_tZ = tbl.UInt8Col(pos=4)
    alogZ_flux = tbl.Float64Col(pos=5, dflt=np.nan)
    alogZ_mass = tbl.Float64Col(pos=6, dflt=np.nan)
    flag_xY = tbl.UInt8Col(pos=7, dflt=0)
