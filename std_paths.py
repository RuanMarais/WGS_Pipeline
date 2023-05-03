'''
**************************************************************************************************
Standard paths for the Python scripts in this directory

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

# Virulence. The VFDB database is used for this analysis https://doi.org/10.1093/nar/gkab1107
virulence_db_path = 'virulence_db/VFDB_310323'
virulence_db = 'VFDB_all'
virulence_info = 'virulence_db/VFDB_310323/VFs.xls'

# KmerFinder path
kmerfinder_path = 'kmerfinder_master/kmerfinder/kmerfinder.py'
kmerfinder_db_path = 'kmerfinder_db/kmerfinder_db/bacteria.ATG'
kmerfinder_db_tax_path = 'kmerfinder_db/kmerfinder_db/bacteria.name'

# Path to the kaptive script and database paths
kaptive_path = 'kaptive_master/Kaptive/kaptive.py'
kaptive_db_klebsiella_k_locus_primary = 'kaptive_master/Kaptive/reference_database/Klebsiella_k_locus_primary_reference.gbk'
kaptive_db_klebsiella_k_locus_variant = 'kaptive_master/Kaptive/reference_database/Klebsiella_k_locus_variant_reference.gbk'
kaptive_db_klebsiella_o_locus_primary = 'kaptive_master/Kaptive/reference_database/Klebsiella_o_locus_primary_reference.gbk'
kaptive_db_acinetobacter_k_locus_primary = 'kaptive_master/Kaptive/reference_database/Acinetobacter_baumannii_k_locus_primary_reference.gbk'
kaptive_db_acinetobacter_OC_locus_primary = 'kaptive_master/Kaptive/reference_database/Acinetobacter_baumannii_OC_locus_primary_reference.gbk'
kaptive_db_path = {'Klebsiella_pneumoniae': [('K_primary', kaptive_db_klebsiella_k_locus_primary),
                                             ('K_variant', kaptive_db_klebsiella_k_locus_variant),
                                             ('O_primary', kaptive_db_klebsiella_o_locus_primary)],
                   'Acinetobacter_baumanii': [('K_primary', kaptive_db_acinetobacter_k_locus_primary),
                                              ('OC_primary', kaptive_db_acinetobacter_OC_locus_primary)]}