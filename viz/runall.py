import os
import sys, argparse

#Angiogenesis
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_angiogenesis/" --analysis "all" --num_genes 36 --hallmark "hallmark_angiogenesis" --dataset "GTEx" ')
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_angiogenesis_panTCGA/" --analysis "all" --num_genes 36 --hallmark "hallmark_angiogenesis" --dataset "panTCGA" ')

#Apical Surface
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_apical_surface/" --analysis "all" --num_genes 44 --hallmark "hallmark_apical_surface" --dataset "GTEx" ')
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_apical_surface_panTCGA/" --analysis "all" --num_genes 43 --hallmark "hallmark_apical_surface" --dataset "panTCGA" ')

#hallmark_hedgehog_signaling
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_hedgehog_signaling/" --analysis "all" --num_genes 36 --hallmark "hallmark_hedgehog_signaling" --dataset "GTEx" ')
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_hedgehog_signaling_panTCGA/" --analysis "all" --num_genes 35 --hallmark "hallmark_hedgehog_signaling" --dataset "panTCGA" ')

#MYCV2
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_myc_targets_v2/" --analysis "all" --num_genes 58 --hallmark "hallmark_myc_targets_v2" --dataset "GTEx" ')
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_myc_targets_v2_panTCGA/" --analysis "all" --num_genes 58 --hallmark "hallmark_myc_targets_v2" --dataset "panTCGA" ')

#Notch_signaling
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_notch_signaling/" --analysis "all" --num_genes 32 --hallmark "hallmark_notch_signaling" --dataset "GTEx" ')
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_notch_signaling_panTCGA/" --analysis "all" --num_genes 31 --hallmark "hallmark_notch_signaling" --dataset "panTCGA" ')

#WNT BETA
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_wnt_beta_catenin_signaling/" --analysis "all" --num_genes 42 --hallmark "hallmark_wnt_beta_catenin_signaling" --dataset "GTEx" ')
os.system('python heatmap_gen.py --graph_type "heatmap" --directory "../logs/hallmark_wnt_beta_catenin_signaling_panTCGA/" --analysis "all" --num_genes 41 --hallmark "hallmark_wnt_beta_catenin_signaling" --dataset "panTCGA" ')
