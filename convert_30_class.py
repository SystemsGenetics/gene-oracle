#/usr/bin/python

# convert 53 class GTEx datasets to 30 class datasets
import os, shutil

subs = os.listdir('./datasets/hallmark_subsets_30')
subs.sort()

for sub in subs:
	path = './datasets/hallmark_subsets_30/' + sub

####Adipose
	os.mkdir(path + '/Adipose')
	files = os.listdir(path + '/Adipose-Subcutaneous')
	for f in files:
		shutil.move(path + '/Adipose-Subcutaneous/' + f, path + '/Adipose')
	shutil.rmtree(path + '/Adipose-Subcutaneous')

	files = os.listdir(path + '/Adipose-Visceral')
	for f in files:
		shutil.move(path + '/Adipose-Visceral/' + f, path + '/Adipose')
	shutil.rmtree(path + '/Adipose-Visceral')

####Blood Vessel
	os.mkdir(path + '/Blood-Vessel')
	files = os.listdir(path + '/Artery-Aorta')
	for f in files:
		shutil.move(path + '/Artery-Aorta/' + f, path + '/Blood-Vessel')
	shutil.rmtree(path + '/Artery-Aorta')

	files = os.listdir(path + '/Artery-Coronary')
	for f in files:
		shutil.move(path + '/Artery-Coronary/' + f, path + '/Blood-Vessel')
	shutil.rmtree(path + '/Artery-Coronary')

	files = os.listdir(path + '/Artery-Tibial')
	for f in files:
		shutil.move(path + '/Artery-Tibial/' + f, path + '/Blood-Vessel')
	shutil.rmtree(path + '/Artery-Tibial')

####Brain
	os.mkdir(path + '/Brain')
	files = os.listdir(path + '/Brain-Amygdala')
	for f in files:
		shutil.move(path + '/Brain-Amygdala/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Amygdala')

	files = os.listdir(path + '/Brain-Anterior')
	for f in files:
		shutil.move(path + '/Brain-Anterior/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Anterior')

	files = os.listdir(path + '/Brain-Caudate')
	for f in files:
		shutil.move(path + '/Brain-Caudate/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Caudate')

	files = os.listdir(path + '/Brain-Cerebellar-Hemisphere')
	for f in files:
		shutil.move(path + '/Brain-Cerebellar-Hemisphere/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Cerebellar-Hemisphere')

	files = os.listdir(path + '/Brain-Cerebellum')
	for f in files:
		shutil.move(path + '/Brain-Cerebellum/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Cerebellum')

	files = os.listdir(path + '/Brain-Cortex')
	for f in files:
		shutil.move(path + '/Brain-Cortex/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Cortex')

	files = os.listdir(path + '/Brain-Frontal')
	for f in files:
		shutil.move(path + '/Brain-Frontal/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Frontal')

	files = os.listdir(path + '/Brain-Hippocampus')
	for f in files:
		shutil.move(path + '/Brain-Hippocampus/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Hippocampus')

	files = os.listdir(path + '/Brain-Hypothalamus')
	for f in files:
		shutil.move(path + '/Brain-Hypothalamus/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Hypothalamus')

	files = os.listdir(path + '/Brain-Nucleus')
	for f in files:
		shutil.move(path + '/Brain-Nucleus/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Nucleus')

	files = os.listdir(path + '/Brain-Putamen')
	for f in files:
		shutil.move(path + '/Brain-Putamen/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Putamen')

	files = os.listdir(path + '/Brain-Spinal')
	for f in files:
		shutil.move(path + '/Brain-Spinal/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Spinal')

	files = os.listdir(path + '/Brain-Substantia')
	for f in files:
		shutil.move(path + '/Brain-Substantia/' + f, path + '/Brain')
	shutil.rmtree(path + '/Brain-Substantia')

#####Breast
	os.mkdir(path + '/Breast')
	files = os.listdir(path + '/Breast-Mammary')
	for f in files:
		shutil.move(path + '/Breast-Mammary/' + f, path + '/Breast')
	shutil.rmtree(path + '/Breast-Mammary')

#####Blood
	os.mkdir(path + '/Blood')
	files = os.listdir(path + '/Cells-EBV')
	for f in files:
		shutil.move(path + '/Cells-EBV/' + f, path + '/Blood')
	shutil.rmtree(path + '/Cells-EBV')

	files = os.listdir(path + '/Whole-Blood')
	for f in files:
		shutil.move(path + '/Whole-Blood/' + f, path + '/Blood')
	shutil.rmtree(path + '/Whole-Blood')

####Cervix
	os.mkdir(path + '/Cervix')
	files = os.listdir(path + '/Cervix-Ectocervix')
	for f in files:
		shutil.move(path + '/Cervix-Ectocervix/' + f, path + '/Cervix')
	shutil.rmtree(path + '/Cervix-Ectocervix')

	files = os.listdir(path + '/Cervix-Endocervix')
	for f in files:
		shutil.move(path + '/Cervix-Endocervix/' + f, path + '/Cervix')
	shutil.rmtree(path + '/Cervix-Endocervix')

####Colon
	os.mkdir(path + '/Colon')
	files = os.listdir(path + '/Colon-Sigmoid')
	for f in files:
		shutil.move(path + '/Colon-Sigmoid/' + f, path + '/Colon')
	shutil.rmtree(path + '/Colon-Sigmoid')

	files = os.listdir(path + '/Colon-Transverse')
	for f in files:
		shutil.move(path + '/Colon-Transverse/' + f, path + '/Colon')
	shutil.rmtree(path + '/Colon-Transverse')

####Esophagus
	os.mkdir(path + '/Esophagus')
	files = os.listdir(path + '/Esophagus-Gastroesophageal')
	for f in files:
		shutil.move(path + '/Esophagus-Gastroesophageal/' + f, path + '/Esophagus')
	shutil.rmtree(path + '/Esophagus-Gastroesophageal')

	files = os.listdir(path + '/Esophagus-Mucosa')
	for f in files:
		shutil.move(path + '/Esophagus-Mucosa/' + f, path + '/Esophagus')
	shutil.rmtree(path + '/Esophagus-Mucosa')

	files = os.listdir(path + '/Esophagus-Muscularis')
	for f in files:
		shutil.move(path + '/Esophagus-Muscularis/' + f, path + '/Esophagus')
	shutil.rmtree(path + '/Esophagus-Muscularis')

####Heart
	os.mkdir(path + '/Heart')
	files = os.listdir(path + '/Heart-Atrial')
	for f in files:
		shutil.move(path + '/Heart-Atrial/' + f, path + '/Heart')
	shutil.rmtree(path + '/Heart-Atrial')

	files = os.listdir(path + '/Heart-Left')
	for f in files:
		shutil.move(path + '/Heart-Left/' + f, path + '/Heart')
	shutil.rmtree(path + '/Heart-Left')

####Kidney
	os.mkdir(path + '/Kidney')
	files = os.listdir(path + '/Kidney-Cortex')
	for f in files:
		shutil.move(path + '/Kidney-Cortex/' + f, path + '/Kidney')
	shutil.rmtree(path + '/Kidney-Cortex')

####Muscle
	os.mkdir(path + '/Muscle')
	files = os.listdir(path + '/Muscle-Skeletal')
	for f in files:
		shutil.move(path + '/Muscle-Skeletal/' + f, path + '/Muscle')
	shutil.rmtree(path + '/Muscle-Skeletal')

####Nerve
	os.mkdir(path + '/Nerve')
	files = os.listdir(path + '/Nerve-Tibial')
	for f in files:
		shutil.move(path + '/Nerve-Tibial/' + f, path + '/Nerve')
	shutil.rmtree(path + '/Nerve-Tibial')

####Salivary
	os.mkdir(path + '/Salivary-Gland')
	files = os.listdir(path + '/Minor-Salivary-Gland')
	for f in files:
		shutil.move(path + '/Minor-Salivary-Gland/' + f, path + '/Salivary-Gland')
	shutil.rmtree(path + '/Minor-Salivary-Gland')

####Skin
	os.mkdir(path + '/Skin')
	files = os.listdir(path + '/Skin-Not-Sun-Exposed')
	for f in files:
		shutil.move(path + '/Skin-Not-Sun-Exposed/' + f, path + '/Skin')
	shutil.rmtree(path + '/Skin-Not-Sun-Exposed')

	files = os.listdir(path + '/Skin-Sun-Exposed')
	for f in files:
		shutil.move(path + '/Skin-Sun-Exposed/' + f, path + '/Skin')
	shutil.rmtree(path + '/Skin-Sun-Exposed')

	files = os.listdir(path + '/Cells-Transformed')
	for f in files:
		shutil.move(path + '/Cells-Transformed/' + f, path + '/Skin')
	shutil.rmtree(path + '/Cells-Transformed')