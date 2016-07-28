#
# Lacerda@Macarrao - 27/Jul/2016
#
import os
import sys
import numpy as np
import tables as tbl
# import matplotlib as mpl
# from matplotlib import pyplot as plt
# local_path = os.getcwd()
# sys.path.append(local_path.rpartition('/')[0])
# import tables_description
from CALIFAUtils.scripts import my_morf_seq


if __name__ == '__main__':
    h5fname = sys.argv[1]
    SFRgroup = '/SFR05050525'
    fpref = '.png'

    iT = 1
    h5file = tbl.open_file(h5fname, mode='r')
    tbl_gals = h5file.root.pycasso.main
    tbl_zones = h5file.root.pycasso.zones
    tbl_integrated = h5file.root.pycasso.integrated

    balm_threshold = 2.85
    SNHa_threshold = 3.
    SNHb_threshold = 3.

    balm_low_frac_threshold = 0.5  # 50%

    print '### Best view if tabspace = 20 (in vim :set ts=20)'
    print '### %s report ### ' % sys.argv[0]
    print '### theoretical j_Ha/Hb: %.2f' % balm_threshold
    print '### SNHb: %.2f - SNHa: %.2f' % (SNHb_threshold, SNHa_threshold)
    print '### Balmer decrement zones fraction threshold: %.2f' % balm_low_frac_threshold
    print '### file: %s' % os.path.basename(h5fname)
    print '### N_gals: %d' % tbl_gals.nrows
    print '###'
    print "### MTYPE - morf type from Enrique's get_morfologia() (morph_eye_class.csv)"
    print "### BALMDECUNDER_FRAC - NZONES_BALMDECUNDER / NZONES_OK"
    print "### BALMDEC_INT - Obs. Ha/Hb integrated"
    print "### NZONES_OK - Number of zones OK"
    print "### NZONES_BALMDECUNDER - Number of zones OK with Ha/Hb < theoretical j_Ha/Hb"
    print "### NZONES - Galaxy number of zones"
    print '### Lines starting with * means that the galaxy has BALMDECUNDER_FRAC'
    print '###       bigger than zones fraction threshold'
    print '###'
    print "#CALIFAID\tMTYPE\tBALMDECUNDER_FRAC\tBALMDEC_INT\tNZONES_OK\tNZONES_BALMDECUNDER\tNZONES"

    balm = []
    m_type = []

    for g in tbl_gals:
        g_props__z = tbl_zones.read_where('id_gal == gid', {'gid': g['id']})
        g_int_props = tbl_integrated.read_where('id_gal == gid', {'gid': g['id']})

        Hb_obs__z = g_props__z['F_obs_Hb']
        Ha_obs__z = g_props__z['F_obs_Ha']
        SNHb_obs__z = Hb_obs__z / g_props__z['eF_obs_Hb']
        SNHa_obs__z = Ha_obs__z / g_props__z['eF_obs_Ha']
        balm_decr__z = Ha_obs__z/Hb_obs__z

        mask = np.isnan(balm_decr__z)
        # mask = np.bitwise_or(mask, Hb_obs__z < 1e-30)
        # mask = np.bitwise_or(mask, Ha_obs__z < 1e-20)
        mask = np.bitwise_or(mask, SNHb_obs__z < SNHb_threshold)
        mask = np.bitwise_or(mask, SNHa_obs__z < SNHa_threshold)
        balm_decr_masked__z = np.ma.masked_array(balm_decr__z, mask=mask)

        N_zones_ok = balm_decr_masked__z.count()
        N = (balm_decr_masked__z < balm_threshold).astype(int).sum()
        balm_decr_less_frac = 1.0 * N/N_zones_ok
        N_masked = np.ma.count_masked(balm_decr_masked__z)

        Hb_obs_int = g_int_props['F_obs_Hb']
        Ha_obs_int = g_int_props['F_obs_Ha']
        balm_decr_int = Ha_obs_int/Hb_obs_int
        star = ''
        if balm_decr_less_frac > balm_low_frac_threshold:
            star = '*'
        print "%s%s\t%s\t%.6f\t%.6f\t%d\t%d\t%d" % (star, g['califaID'],
                                                    g['m_type_orig'],
                                                    balm_decr_less_frac,
                                                    balm_decr_int,
                                                    N_zones_ok, N,
                                                    g['N_zone'])

        m_type.append(my_morf_seq(g['m_type_orig']))
        balm.append(balm_decr_masked__z.mean())
    h5file.close()
