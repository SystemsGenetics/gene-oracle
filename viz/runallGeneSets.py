import os
import sys, argparse

print("Running WNTBETA...")
#WNTBETA GTEX, Genes = 42
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_42s/ --sub_dir ../logs/hallmark_wnt_beta_catenin_signaling --title WNTBETA_GTEx')


#WNTBETA TCGA, Genes = 41
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_42s_panTCGA/ --sub_dir ../logs/hallmark_wnt_beta_catenin_signaling_panTCGA --title WNTBETA_TCGA')


#Notch GTEX, Genes = 32
print("Running Notch...")
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_32s/ --sub_dir ../logs/hallmark_notch_signaling --title Notch_Signaling_GTEx')

#Notch TCGA, Genes = 31
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_31s_panTCGA/ --sub_dir ../logs/hallmark_notch_signaling_panTCGA --title Notch_Signaling_panTCGA')

print("Running MYCV2...")
#MYCV2 GTEx, Genes = 58
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_58s/ --sub_dir ../logs/hallmark_myc_targets_v2 --title MYC_Targets_V2_GTEx')

#MYCV2 TCGA, Genes = 58
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_58s_panTCGA/ --sub_dir ../logs/hallmark_hedgehog_signaling_panTCGA --title MYC_Targets_V2_panTCGA')

print("Running HedgeHog...")
#HedgeHog GTEx, Genes = 36
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_36s/ --sub_dir ../logs/hallmark_hedgehog_signaling --title HedgeHog_GTEx')

#HedgeHog TCGA, Genes = 35
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_35s_panTCGA/ --sub_dir ../logs/hallmark_hedgehog_signaling_panTCGA --title HedgeHog_panTCGA')

print("Running AS...")
#Apical Surface GTEx, Genes = 44
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_44s/ --sub_dir ../logs/hallmark_apical_surface_panTCGA --title Apical_Surface_GTEx')

#Apical Surface TCGA, Genes = 43
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_43s_panTCGA/ --sub_dir ../logs/hallmark_apical_surface --title Apical_Surface_panTCGA')

print("Running Angiogenesis...")
#Angiogenesis GTEx, Genes = 36
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_36s/ --sub_dir ../logs/hallmark_angiogenesis --title Angiogenesis_GTEx')

#Angiogenesis TCGA, Genes = 36
os.system('python gene_sets_acc.py  --rand_dir ../logs/random_36s_panTCGA/ --sub_dir ../logs/hallmark_angiogenesis_panTCGA --title Angiogenesis_panTCGA')
