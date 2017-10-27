#/usr/bin/python

import numpy as np 
import os

# dictionary for each sample type

tissues = {
	"Adipose-Subcutaneous" : 350,
	"Adipose-Visceral" : 227,
	"Adrenal-Gland" : 145,
	"Artery-Aorta" : 224,
	"Artery-Coronary" : 133,
	"Artery-Tibial" : 332,
	"Bladder" : 11,
	"Brain-Amygdala" : 72,
	"Brain-Anterior" : 84,
	"Brain-Caudate" : 117,
	"Brain-Cerebellar-Hemisphere" : 105,
	"Brain-Cerebellum" : 125,
	"Brain-Cortex" : 114,
	"Brain-Frontal" : 108,
	"Brain-Hippocampus" : 94,
	"Brain-Hypothalamus" : 96,
	"Brain-Nucleus" : 113,
	"Brain-Putamen" : 97,
	"Brain-Spinal" : 71,
	"Brain-Substantia" : 63,
	"Breast-Mammary" : 214,
	"Cells-EBV" : 118,
	"Cells-Transformed" : 284,
	"Cervix-Ectocervix" : 6,
	"Cervix-Endocervix" : 5,
	"Colon-Sigmoid" : 149,
	"Colon-Transverse" : 196,
	"Esophagus-Gastroesophageal" : 153,
	"Esophagus-Mucosa" : 286,
	"Esophagus-Muscularis" : 247,
	"Fallopian Tube" : 6,
	"Heart-Atrial" : 194,
	"Heart-Left" : 218,
	"Kidney-Cortex" : 32,
	"Liver" : 119,
	"Lung" : 320,
	"Minor-Salivary-Gland" : 57,
	"Muscle-Skeletal" : 430,
	"Nerve-Tibial" : 304,
	"Ovary" : 97,
	"Pancreas" : 171,
	"Pituitary" : 103,
	"Prostate" : 106,
	"Skin-Not-Sun-Exposed" : 250,
	"Skin-Sun-Exposed" : 357,
	"Small-Intestine" : 88,
	"Spleen" : 104,
	"Stomach" : 193,
	"Testis" : 172,
	"Thyroid" : 323,
	"Uterus" : 83,
	"Vagina" : 96,
	"Whole-Blood" : 393
}

files = os.listdir('../datasets/hallmark_subsets/')
cols = np.load('../datasets/cols_gtex.npy')

for f in files:
	dir = f.split('.')[0]
	if not os.path.exists('../datasets/hallmark_subsets/' + dir):
		os.makedirs('../datasets/hallmark_subsets/' + dir)

	data = np.load('../datasets/hallmark_subsets/' + f)
	data = np.copy(data[:,2:])
	data = data.astype(np.float32)

	print('generating files for ' + dir + ' ')

	j = 0

	#loop through each example and save it to a h5 file
	for key in sorted(tissues.iterkeys()):
		
		#create class directory if not already there
		if not os.path.exists('../datasets/hallmark_subsets/' + dir + '/' + key):
			os.makedirs('../datasets/hallmark_subsets/' + dir + '/' + key)

		#create h5 file for each sample
		for i in range(tissues[key]):
			data[:,j].tofile('../datasets/hallmark_subsets/' + dir + '/' + key + '/' + str('%03d' % i) + '_' + cols[j] + '.dat', sep="")
			# np.savetxt('GTEx_Data/' + key + '/' + str('%03d' % i) + '_' + cols[j] + '.txt', gtex_data[:,j], fmt='%8f')
			j = j + 1

print('\nDone!\n')

