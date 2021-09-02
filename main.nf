#!/usr/bin/env nextflow



/**
 * Create channels for input files.
 */
EMX_FILES = Channel.fromFilePairs("${params.input_dir}/${params.emx_files}", size: 1, flat: true)
LABEL_FILES = Channel.fromFilePairs("${params.input_dir}/${params.label_files}", size: 1, flat: true)
GMT_FILES = Channel.fromFilePairs("${params.input_dir}/${params.gmt_files}", size: 1, flat: true)



/**
 * Group dataset files by name.
 */
DATASETS = EMX_FILES.join(LABEL_FILES)



/**
 * Send dataset files to each process that uses them.
 */
DATASETS
    .into {
        DATASETS_FOR_PHASE1_FG;
        DATASETS_FOR_PHASE1_BG;
        DATASETS_FOR_PHASE1_SELECT;
        DATASETS_FOR_PHASE2_EVALUATE;
        DATASETS_FOR_PHASE2_RF
    }



/**
 * Send GMT files to each process that uses them.
 */
GMT_FILES.into {
    GMT_FILES_FOR_PHASE1_FG;
    GMT_FILES_FOR_PHASE1_BG;
    GMT_FILES_FOR_PHASE1_SELECT;
    GMT_FILES_FOR_PHASE2_RF
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
        params.phase1 == true

    script:
        """
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE gmt_lines=`cat ${gmt_file} | wc -l`"

        split -d -n l/${params.chunks} ${gmt_file} ""
        """
}



/**
 * Combine the datasets and gmt chunks using the cartesian product
 * in order to trigger the phase1_fg process for each dataset and each gmt chunk.
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
        set val(dataset), file(emx_file), file(label_file), val(gmt_name), file(gmt_file) from PHASE1_FG_CHUNKS_COMBINED

    output:
        set val(dataset), val(gmt_name), file("*.log") into PHASE1_FG_SCORE_CHUNKS

    when:
        params.phase1 == true

    script:
        """
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE gmt_lines=`cat ${gmt_file} | wc -l`"
        echo "#TRACE gmt_genes=`cat ${gmt_file} | wc -w`"
        echo "#TRACE n_rows=`tail -n +1 ${emx_file} | wc -l`"
        echo "#TRACE n_cols=`head -n +1 ${emx_file} | wc -w`"
        echo "#TRACE model=${params.phase1_model}"

        phase1-evaluate.py \
            --dataset      ${emx_file} \
            --labels       ${label_file} \
            --model-config ${baseDir}/example/models.json \
            --model        ${params.phase1_model} \
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
        set val(dataset), file(emx_file), file(label_file) from DATASETS_FOR_PHASE1_BG
        each(index) from Channel.from( 0 .. params.chunks-1 )

    output:
        set val(dataset), file("*.log") into PHASE1_BG_SCORE_CHUNKS_RAW

    when:
        params.phase1 == true

    script:
        """
        echo "#TRACE random_min=${params.phase1_random_min}"
        echo "#TRACE random_max=${params.phase1_random_max}"
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE index=${index}"
        echo "#TRACE n_rows=`tail -n +1 ${emx_file} | wc -l`"
        echo "#TRACE n_cols=`head -n +1 ${emx_file} | wc -w`"
        echo "#TRACE model=${params.phase1_model}"

        START=${params.phase1_random_min + index}
        STOP=${params.phase1_random_max}
        STEP=${params.chunks}

        phase1-evaluate.py \
            --dataset      ${emx_file} \
            --labels       ${label_file} \
            --model-config ${baseDir}/example/models.json \
            --model        ${params.phase1_model} \
            --random \
            --random-range \${START} \${STOP} \${STEP} \
            --random-iters ${params.phase1_random_iters} \
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
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}"

    input:
        set val(dataset), val(gmt_name), file(chunks) from PHASE1_SCORE_CHUNKS

    output:
        set val(dataset), val(gmt_name), file("phase1-evaluate-*.log") into PHASE1_SCORES

    when:
        params.phase1 == true

    script:
        """
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE gmt_lines=`cat ${chunks} | wc -l`"

        head -n +1 ${chunks[0]} > phase1-evaluate-${gmt_name}.log

        for f in ${chunks}; do
            tail -n +2 \$f >> temp
        done

        sort -V temp >> phase1-evaluate-${gmt_name}.log
        """
}



/**
 * Combine the datasets and gmt files using the cartesian product
 * in order to trigger the phase1_select process for each dataset and each gmt file.
 */
PHASE1_SELECT_INPUTS_COMBINED = DATASETS_FOR_PHASE1_SELECT.combine(GMT_FILES_FOR_PHASE1_SELECT)



/**
 * The phase1_select process takes a score file for a dataset / GMT and selects
 * gene sets which score significantly higher over background.
 */
process phase1_select {
    publishDir "${params.output_dir}/${dataset}", mode: "copy"

    input:
        set val(dataset), file(emx_file), file(label_file), val(gmt_name), file(gmt_file) from PHASE1_SELECT_INPUTS_COMBINED
        set val(dataset), val(gmt_name), file(scores) from PHASE1_SCORES

    output:
        set val(dataset), val(gmt_name), file("phase1-genesets.txt") into PHASE1_GENESETS
        file("phase1-select-${gmt_name}.log")

    when:
        params.phase1 == true

    script:
        """
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE gmt_lines=`cat ${gmt_file} | wc -l`"
        echo "#TRACE gmt_genes=`cat ${gmt_file} | wc -w`"

        phase1-select.py \
            --dataset   ${emx_file} \
            --gene-sets ${gmt_file} \
            --scores    ${scores} \
            --threshold ${params.phase1_threshold} \
            --n-sets    ${params.phase1_n_sets} \
            > phase1-select-${gmt_name}.log
        """
}



/**
 * Send phase 1 gene sets to each process that uses them.
 */
PHASE1_GENESETS.into {
    GENESETS_FOR_PHASE2_EVALUATE;
    GENESETS_FOR_PHASE2_SELECT;
    GENESETS_FOR_PHASE2_RF
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
        params.phase2 == true

    script:
        """
        split -d -l 1 ${gmt_file} ""
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
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}/${gmt_file.name}"
    label "gpu"

    input:
        set val(dataset), file(emx_file), file(label_file), val(gmt_name), file(gmt_file) from PHASE2_EVALUATE_CHUNKS_COMBINED

    output:
        set val(dataset), val(gmt_name), file("*_scores_*.txt") optional true into PHASE2_SCORE_CHUNKS

    when:
        params.phase2 == true

    script:
        """
        phase2-evaluate.py \
            --dataset      ${emx_file} \
            --labels       ${label_file} \
            --model-config ${baseDir}/example/models.json \
            --model        ${params.phase2_model} \
            --gene-sets    ${gmt_file} \
            --n-jobs       1 \
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
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}"

    input:
        set val(dataset), val(gmt_name), file("phase1-genesets.txt") from GENESETS_FOR_PHASE2_SELECT
        set val(dataset), val(gmt_name), file(scores) from PHASE2_SCORES

    output:
        set val(dataset), val(gmt_name), file("phase2-genesets.txt") into PHASE2_GENESETS
        file("*.png") optional true into PHASE2_PLOTS

    when:
        params.phase2 == true

    script:
        """
        phase2-select.py \
            --gene-sets phase1-genesets.txt \
            --logdir    . \
            --threshold ${params.phase2_threshold} \
            ${params.phase2_visualize ? "--visualize" : ""}
        """
}



/**
 * Combine the datasets and gmt files using the cartesian product
 * in order to trigger the phase2_rf process for each dataset and each gmt file.
 */
PHASE2_RF_INPUTS_COMBINED = DATASETS_FOR_PHASE2_RF.combine(GMT_FILES_FOR_PHASE2_RF)



/**
 * The phase2_rf process takes a list of gene sets selected by phase 1 and
 * uses the feature importances of a random forest to select candidate genes
 * for each gene set.
 */
process phase2_rf {
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}"

    input:
        set val(dataset), file(emx_file), file(label_file), val(gmt_name), file(gmt_file) from PHASE2_RF_INPUTS_COMBINED
        set val(dataset), val(gmt_name), file("phase1-genesets.txt") from GENESETS_FOR_PHASE2_RF

    output:
        set val(dataset), val(gmt_name), file("phase2-rf-genesets.txt") into PHASE2_RF_GENESETS
        file("*.png") optional true into PHASE2_RF_PLOTS

    when:
        params.phase2_rf == true

    script:
        """
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE gmt_lines=`cat phase1-genesets.txt | wc -l`"
        echo "#TRACE gmt_genes=`cat phase1-genesets.txt | wc -w`"
        echo "#TRACE n_rows=`tail -n +1 ${emx_file} | wc -l`"
        echo "#TRACE n_cols=`head -n +1 ${emx_file} | wc -w`"

        phase2-rf.py \
            --dataset   ${emx_file} \
            --labels    ${label_file} \
            --gene-sets phase1-genesets.txt \
            --n-jobs    1 \
            --threshold ${params.phase2_rf_threshold} \
            ${params.phase2_rf_visualize ? "--visualize" : ""}
        """
}
