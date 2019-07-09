#!/usr/bin/env nextflow



/**
 * Create channels for input files.
 */
DATA_TXT_FILES = Channel.fromFilePairs("${params.input.dir}/${params.input.data_txt}", size: 1, flat: true)
DATA_NPY_FILES = Channel.fromFilePairs("${params.input.dir}/${params.input.data_npy}", size: 1, flat: true)
ROWNAME_FILES = Channel.fromFilePairs("${params.input.dir}/*_rownames.txt", size: 1, flat: true)
COLNAME_FILES = Channel.fromFilePairs("${params.input.dir}/*_colnames.txt", size: 1, flat: true)
LABEL_FILES = Channel.fromFilePairs("${params.input.dir}/${params.input.labels}", size: 1, flat: true)
GMT_FILES = Channel.fromFilePairs("${params.input.dir}/${params.input.gmt_files}", size: 1, flat: true)



/**
 * Group dataset files by dataset.
 */
LABEL_FILES.into {
	LABEL_FILES_FOR_TXT;
	LABEL_FILES_FOR_NPY
}

DATA_TXT_FILES
	.map { [it[0], [it[1]]] }
	.mix(LABEL_FILES_FOR_TXT)
	.groupTuple(size: 2)
	.map { [it[0], it[1].sort()] }
	.map { [it[0], it[1][0], it[1][1]] }
	.set { DATA_TXT_COMBINED }

DATA_NPY_FILES
	.mix(ROWNAME_FILES, COLNAME_FILES)
	.groupTuple()
	.map { [it[0], it[1].sort()] }
	.mix(LABEL_FILES_FOR_NPY)
	.groupTuple(size: 2)
	.map { [it[0], it[1].sort()] }
	.map { [it[0], it[1][1], it[1][0]] }
	.set { DATA_NPY_COMBINED }

Channel.empty()
	.mix(DATA_TXT_COMBINED, DATA_NPY_COMBINED)
	.into {
		DATASETS_FOR_PHASE1_FG;
		DATASETS_FOR_PHASE1_BG;
		DATASETS_FOR_PHASE2_EVALUATE
	}



/**
 * Send GMT files to each process that uses them.
 */
GMT_FILES.into {
	GMT_FILES_FOR_PHASE1_FG;
	GMT_FILES_FOR_PHASE1_BG;
	GMT_FILES_FOR_PHASE1_SELECT
}



/**
 * The phase1_split process splits the input GMT file into chunks.
 */
process phase1_split {
	tag "${gmt_name}"

	input:
		set val(gmt_name), file(gmt_file) from GMT_FILES_FOR_PHASE1_FG

	output:
		set val(gmt_name), file("*") into PHASE1_FG_CHUNKS mode flatten

	when:
		params.phase1.enabled == true

	script:
		"""
		split -d -n l/${params.chunks} ${gmt_file} ""
		"""
}



/**
 * Combine the datasets and gene set chunks using the cartesian product
 * in order to trigger the phase1_fg process for each input and each chunk.
 */
PHASE1_FG_CHUNKS_COMBINED = DATASETS_FOR_PHASE1_FG.combine(PHASE1_FG_CHUNKS)



/**
 * The phase1_fg process performs experiments from a single chunk of a
 * GMT file.
 */
process phase1_fg {
	tag "${dataset}/${gmt_name}/${gmt_file.name}"
	label "gpu"

	input:
		set val(dataset), file(data_files), file(labels), val(gmt_name), file(gmt_file) from PHASE1_FG_CHUNKS_COMBINED

	output:
		set val(dataset), val(gmt_name), file("*.log") into PHASE1_FG_SCORE_CHUNKS

	when:
		params.phase1.enabled == true

	script:
		"""
		taskset -c 0-1 phase1-evaluate.py \
			--dataset      ${data_files[0]} \
			--labels       ${labels} \
			--model-config ${baseDir}/example/models.json \
			--model        ${params.phase1.model} \
			--gene-sets    ${gmt_file} \
			--outfile      ${gmt_file}.log
		"""
}



/**
 * The phase1_bg process performs a single chunk of background experiments.
 */
process phase1_bg {
	tag "${dataset}/${index}"
	label "gpu"

	input:
		set val(dataset), file(data_files), file(labels) from DATASETS_FOR_PHASE1_BG
		each(index) from Channel.from( 0 .. params.chunks-1 )

	output:
		set val(dataset), file("*.log") into PHASE1_BG_SCORE_CHUNKS_RAW

	when:
		params.phase1.enabled == true

	script:
		"""
		START=${params.phase1.random_min + index}
		STOP=${params.phase1.random_max}
		STEP=${params.chunks}

		taskset -c 0-1 phase1-evaluate.py \
			--dataset      ${data_files[0]} \
			--labels       ${labels} \
			--model-config ${baseDir}/example/models.json \
			--model        ${params.phase1.model} \
			--random \
			--random-range \$START \$STOP \$STEP \
			--random-iters ${params.phase1.random_iters} \
			--outfile      \$(printf "%04d" ${index}).log
		"""
}



/**
 * Use the same background scores for each GMT file.
 */
PHASE1_BG_SCORE_CHUNKS_RAW
	.combine(GMT_FILES_FOR_PHASE1_BG)
	.map { [it[0], it[2], it[1]] }
	.set { PHASE1_BG_SCORE_CHUNKS }



/**
 * Group output chunks by dataset and GMT file so that they can be merged.
 */
Channel.empty()
	.concat(PHASE1_FG_SCORE_CHUNKS, PHASE1_BG_SCORE_CHUNKS)
	.groupTuple(by: [0, 1])
	.map { [it[0], it[1], it[2].sort {it.name}] }
	.set { PHASE1_SCORE_CHUNKS }



/**
 * The phase1_merge process takes the output chunks from phase1_evaluate
 * and merges them into a score file for each dataset/GMT pair.
 */
process phase1_merge {
	publishDir "${params.output.dir}/${dataset}", mode: "copy"
	tag "${dataset}/${gmt_name}"

	input:
		set val(dataset), val(gmt_name), file(chunks) from PHASE1_SCORE_CHUNKS

	output:
		set val(dataset), val(gmt_name), file("phase1-scores.txt") into PHASE1_SCORES

	when:
		params.phase1.enabled == true

	script:
		"""
		head -n +1 ${chunks[0]} > phase1-scores.txt

		for f in ${chunks}; do
			tail -n +2 \$f >> temp
		done

		sort -V temp >> phase1-scores.txt
		"""
}



/**
 * The phase1_select process takes a score file for a dataset / GMT and selects
 * gene sets which score significantly higher over background.
 */
process phase1_select {
	publishDir "${params.output.dir}/${dataset}", mode: "copy"
	tag "${dataset}/${gmt_name}"

	input:
		set val(dataset), val(gmt_name), file(scores) from PHASE1_SCORES
		set val(gmt_name), file(gmt_file) from GMT_FILES_FOR_PHASE1_SELECT

	output:
		set val(dataset), val(gmt_name), file("phase1-genesets.txt") into PHASE1_GENESETS
		file("phase1-select-${gmt_name}.log")

	when:
		params.phase1.enabled == true

	script:
		"""
		phase1-select.py \
			--scores    ${scores} \
			--gene-sets ${gmt_file} \
			--threshold ${params.phase1.threshold} \
			--n-sets    ${params.phase1.n_sets} \
			--outfile   phase1-genesets.txt \
			> phase1-select-${gmt_name}.log
		"""
}



/**
 * Send phase 1 gene sets to each process that uses them.
 */
PHASE1_GENESETS.into {
	GENESETS_FOR_PHASE2_EVALUATE;
	GENESETS_FOR_PHASE2_SELECT
}



/**
 * The phase2_split process splits the gene set list from phase 1 into chunks.
 */
process phase2_split {
	tag "${dataset}/${gmt_name}"

	input:
		set val(dataset), val(gmt_name), file(gmt_file) from GENESETS_FOR_PHASE2_EVALUATE

	output:
		set val(dataset), val(gmt_name), file("*") into PHASE2_EVALUATE_CHUNKS mode flatten

	when:
		params.phase2.enabled == true

	script:
		"""
		split -d -n l/${params.chunks} ${gmt_file} ""
		"""
}



/**
 * Combine the datasets and gene set chunks using the cartesian product
 * in order to trigger the phase2_evaluate process for each input and each chunk.
 */
DATASETS_FOR_PHASE2_EVALUATE
	.combine(PHASE2_EVALUATE_CHUNKS, by: 0)
	.set { PHASE2_EVALUATE_CHUNKS_COMBINED }



/**
 * The phase2_evaluate process performs subset analysis on a gene set.
 */
process phase2_evaluate {
	publishDir "${params.output.dir}/${dataset}", mode: "copy"
	tag "${dataset}/${gmt_name}/${gmt_file.name}"
	label "gpu"

	input:
		set val(dataset), file(data_files), file(labels), val(gmt_name), file(gmt_file) from PHASE2_EVALUATE_CHUNKS_COMBINED

	output:
		set val(dataset), val(gmt_name), file("*_scores_*.txt") optional true into PHASE2_SCORE_CHUNKS

	when:
		params.phase2.enabled == true

	script:
		"""
		taskset -c 0-1 phase2-evaluate.py \
			--dataset      ${data_files[0]} \
			--labels       ${labels} \
			--model-config ${baseDir}/example/models.json \
			--model        ${params.phase2.model} \
			--gene-sets    ${gmt_file} \
			--logdir       .
		"""
}



/**
 * Group phase 2 scores by dataset and geneset list.
 */
PHASE2_SCORE_CHUNKS
	.groupTuple(by: [0, 1])
	.map { [it[0], it[1], it[2].flatten()] }
	.set { PHASE2_SCORES }



/**
 * The phase2_select process takes a list of gene sets selected by phase 1, as well
 * as subset scores from phase2_evaluate, and selects candidate genes for each
 * gene set.
 */
process phase2_select {
	publishDir "${params.output.dir}/${dataset}", mode: "copy"
	tag "${dataset}/${gmt_name}"

	input:
		set val(dataset), val(gmt_name), file("phase1-genesets.txt") from GENESETS_FOR_PHASE2_SELECT
		set val(dataset), val(gmt_name), file(scores) from PHASE2_SCORES

	output:
		set val(dataset), val(gmt_name), file("phase2-genesets.txt") into PHASE2_GENESETS
		file("*.png") optional true into PHASE2_PLOTS

	when:
		params.phase2.enabled == true

	script:
		"""
		phase2-select.py \
			--gene-sets phase1-genesets.txt \
			--logdir    . \
			${params.phase2.visualize ? "--visualize" : ""} \
			--outfile   phase2-genesets.txt
		"""
}
