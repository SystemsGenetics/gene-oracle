import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys, argparse
import json


'''

cand_plot.py

arguments:

--subset_name expects : HedgeHog, Mycv2, Angiogenesis, Notch



'''

#will return a array of 3 values = [og_score, nc_score, c_score]
def read_file(file):
	data = []
	with open(file, 'r') as f:
		next(f)
		for line in f:
			line = line.split()
			data.append(float(line[1]))
	return data


def read_random_file(file):
	data = {}
	with open(file,'r') as f:
		next(f)
		for line in f:
			line = line.split()
			data[str(line[0])] = line[1]
		return data	



def plot(subset_name, gtex_result, tcga_result,gtex_rand_result,tcga_rand_result):
	#plt.figure(figsize=(20,10))

	#get data for randoms
	random_pantcga = read_random_file(tcga_rand_result)
	random_gtex = read_random_file(gtex_rand_result)

	#read json with gene counts 
	with open("candidate_gene_counts.json",'r') as f:
		gene_counts = json.load(f)
	
	tcga_subset_gene_counts =  [gene_counts["pantcga"][str(subset_name)]["original"],
								gene_counts["pantcga"][str(subset_name)]["noncandidate"],
								gene_counts["pantcga"][str(subset_name)]["candidate"]]
	gtex_subset_gene_counts = [gene_counts["gtex"][str(subset_name)]["original"],
								gene_counts["gtex"][str(subset_name)]["noncandidate"],
								gene_counts["gtex"][str(subset_name)]["candidate"]]

	random_og_scores =[random_gtex[str(gtex_subset_gene_counts[0])], random_pantcga[str(tcga_subset_gene_counts[0])]]
	random_noncandidate_scores = [random_gtex[str(gtex_subset_gene_counts[1])], random_pantcga[str(tcga_subset_gene_counts[1])]]
	random_candidate_scores = [random_gtex[str(gtex_subset_gene_counts[2])], random_pantcga[str(tcga_subset_gene_counts[2])]]

	print(random_og_scores)
	print(random_noncandidate_scores)
	print(random_candidate_scores)

	# get data
	n_groups = 2 # 1 for tcga, 1 for gtex
	gtex = read_file(gtex_result)
	tcga = read_file(tcga_result)
	
	og_scores = [gtex[0], tcga[0]]
	noncandidate_scores = [gtex[1], tcga[1]]
	candidate_scores = [gtex[2], tcga[2]]
	print("og_scores" + str(og_scores[0]) + str(og_scores[1]))
	print("nc_scores" + str(noncandidate_scores[0])+ str(noncandidate_scores[1]))
	print("candidate_scores" + str(candidate_scores[0]) + str(candidate_scores[1]))

	#create plot
	fig, ax = plt.subplots()
	index = np.arange(n_groups)
	bar_width = .1
	opacity = .8

	og_bars = plt.bar(index, og_scores,bar_width,color='g',label="Originals")
	r_og_bars = plt.bar(index, random_og_scores, bar_width,color='r',label="Random")

	nc_bars = plt.bar(index+bar_width +.01, noncandidate_scores, bar_width,color='b',label="NonCandidate")
	r_nc_bars = plt.bar(index+bar_width +.01, random_noncandidate_scores, bar_width,color='r')

	c_bars = plt.bar(index+bar_width+bar_width +.02, candidate_scores, bar_width,color='#a020f0',label="Candidate")
	r_c_bars = plt.bar(index+bar_width+bar_width+ .02, random_candidate_scores, bar_width,color='r')

	plt.xlabel("Dataset")
	plt.ylabel("Accuracy")
	plt.title(subset_name)
	plt.xticks(index+bar_width, ("GTEx", "panTCGA"));
	#plt.legend()
	plt.tight_layout()
	plt.savefig(str(subset_name) + "_can_analysis.png")


if __name__ == '__main__':


	parser = argparse.ArgumentParser(description='Visualize Candidate Genes vs NonCandidate for panTCGA and GTEx datasets')
	parser.add_argument('--subset_name', help='name of subset', type=str, required=True)
	parser.add_argument('--gtex_result', help='file with gtex results for a subset', required=True)
	parser.add_argument('--tcga_result', help='file with tcga results for a subset', required=True)
	parser.add_argument('--gtex_rand_result', help='file with rand results for gtex gene counts', required=True)
	parser.add_argument('--tcga_rand_result', help='file with rand results for tcga gene counts', required=True)

	args = parser.parse_args()

	plot(args.subset_name, args.gtex_result,args.tcga_result, args.gtex_rand_result,args.tcga_rand_result)