#!/usr/bin/env nextflow

nextflow.enable.dsl=2



workflow {
    // create synthetic data if specified
    if ( params.make_inputs == true ) {
        make_inputs()
        emx_files    = make_inputs.out.dataset
        labels_files = make_inputs.out.labels
        gmt_files    = make_inputs.out.gmt_file
    }

    // otherwise load input files
    else {
        emx_files    = Channel.fromFilePairs("${params.input_dir}/${params.emx_files}", size: 1, flat: true)
        labels_files = Channel.fromFilePairs("${params.input_dir}/${params.labels_files}", size: 1, flat: true)
        gmt_files    = Channel.fromFilePairs("${params.input_dir}/${params.gmt_files}", size: 1, flat: true)
    }

    // group dataset files by name
    datasets = emx_files.join(labels_files)

    // combine datasets with gmt files via cartesian product
    inputs_combined = datasets.combine(gmt_files)

    // perform phase 1 if specified
    if ( params.phase1 == true ) {
        // split gmt files into chunks
        phase1_split(gmt_files)
        gmt_chunks = phase1_split.out.chunks
            .flatMap { it[1].collect { c -> [it[0], c] } }

        // combine datasets with gmt file chunks via cartesian product
        fg_chunks = datasets.combine(gmt_chunks)

        // perform foreground evaluation
        phase1_fg(fg_chunks)
        fg_scores = phase1_fg.out.score_chunks

        // perform background evaluation
        bg_chunks = Channel.from( 0 .. params.chunks-1 )
        phase1_bg(datasets, bg_chunks)

        // combine background scores with gmt files
        bg_scores = phase1_bg.out.score_chunks
            .combine(gmt_files)
            .map { [it[0], it[2], it[1]] }

        // combine fg and bg score chunks
        // group score chunks by dataset and gmt file
        score_chunks = Channel.empty()
            .concat(fg_scores, bg_scores)
            .groupTuple(by: [0, 1])
            .map { [it[0], it[1], it[2].sort {it.name}] }

        // merge score chunks
        phase1_merge(score_chunks)
        phase1_scores = phase1_merge.out.scores

        // use scores to select gene sets
        phase1_select(inputs_combined, phase1_scores)
        gene_sets = phase1_select.out.gene_sets
    }
    else {
        gene_sets = Channel.empty()
    }

    // perform phase 2 (combinatoric search) if specified
    if ( params.phase2 == true ) {
        // split gene sets into chunks
        phase2_split(gene_sets)
        geneset_chunks = phase2_split.out.chunks
            .flatMap { it[2].collect { c -> [it[0], it[1], c] } }

        // combine datasets with gene set chunks via cartesian product
        phase2_chunks = datasets.combine(geneset_chunks, by: 0)

        // evaluate gene sets via combinatoric search
        phase2_evaluate(phase2_chunks)

        // group scores by dataset and geneset list
        phase2_scores = phase2_evaluate.out.score_chunks
            .groupTuple(by: [0, 1])
            .map { [it[0], it[1], it[2].flatten()] }

        // use scores to select gene sets
        phase2_select(gene_sets, phase2_scores)
    }

    // perform phase 2 (random forest) if specified
    if ( params.phase2_rf == true ) {
        phase2_rf(inputs_combined, gene_sets)
    }
}



/**
 * The make_inputs process generates synthetic input data
 * for an example run.
 */
process make_inputs {
    publishDir "${params.output_dir}", mode: "copy"

    output:
        tuple val("example"), path("example.emx.txt"),      emit: dataset
        tuple val("example"), path("example.labels.txt"),   emit: labels
        tuple val("example"), path("example.genesets.txt"), emit: gmt_file
        path("example.tsne.png")

    script:
        """
        make-inputs.py \
            --n-samples 1000 \
            --n-genes   100 \
            --n-classes 5 \
            --n-sets    10 \
            --visualize
        """
}



/**
 * The phase1_split process splits the input gmt file into chunks.
 */
process phase1_split {
    tag "${gmt_name}"

    input:
        tuple val(gmt_name), path(gmt_file)

    output:
        tuple val(gmt_name), path("*"), emit: chunks

    script:
        """
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE gmt_lines=`cat ${gmt_file} | wc -l`"

        split -d -n l/${params.chunks} ${gmt_file} ""
        """
}



/**
 * The phase1_fg process performs experiments from a single chunk
 * of a gmt file.
 */
process phase1_fg {
    tag "${dataset}/${gmt_name}/${gmt_file.name}"
    label "gpu"

    input:
        tuple val(dataset), path(emx_file), path(labels_file), val(gmt_name), path(gmt_file)

    output:
        tuple val(dataset), val(gmt_name), path("*.log"), emit: score_chunks

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
            --labels       ${labels_file} \
            --model-config ${baseDir}/models.json \
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
        tuple val(dataset), path(emx_file), path(labels_file)
        each index

    output:
        tuple val(dataset), path("*.log"), emit: score_chunks

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
            --labels       ${labels_file} \
            --model-config ${baseDir}/models.json \
            --model        ${params.phase1_model} \
            --random \
            --random-range \${START} \${STOP} \${STEP} \
            --random-iters ${params.phase1_random_iters} \
            --outfile      \$(printf "%04d" ${index}).log
        """
}



/**
 * The phase1_merge process takes the output chunks from phase1_evaluate
 * and merges them into a score file for each dataset/GMT pair.
 */
process phase1_merge {
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}"

    input:
        tuple val(dataset), val(gmt_name), path(chunks)

    output:
        tuple val(dataset), val(gmt_name), path("phase1-evaluate-*.log"), emit: scores

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
 * The phase1_select process takes a score file for a dataset / GMT and selects
 * gene sets which score significantly higher over background.
 */
process phase1_select {
    publishDir "${params.output_dir}/${dataset}", mode: "copy"

    input:
        tuple val(dataset), path(emx_file), path(labels_file), val(gmt_name), path(gmt_file)
        tuple val(dataset), val(gmt_name), path(scores)

    output:
        tuple val(dataset), val(gmt_name), path("phase1-genesets.txt"), emit: gene_sets
        path("phase1-select-${gmt_name}.log")

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
            --n-sets    ${params.phase1_nsets} \
            > phase1-select-${gmt_name}.log
        """
}



/**
 * The phase2_split process splits the gene set list from phase 1 into chunks.
 */
process phase2_split {
    tag "${dataset}/${gmt_name}"

    input:
        tuple val(dataset), val(gmt_name), path(gmt_file)

    output:
        tuple val(dataset), val(gmt_name), path("*"), emit: chunks

    script:
        """
        split -d -l 1 ${gmt_file} ""
        """
}



/**
 * The phase2_evaluate process performs subset analysis on a gene set.
 */
process phase2_evaluate {
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}/${gmt_file.name}"
    label "gpu"

    input:
        tuple val(dataset), path(emx_file), path(labels_file), val(gmt_name), path(gmt_file)

    output:
        tuple val(dataset), val(gmt_name), path("*_scores_*.txt"), optional: true, emit: score_chunks

    script:
        """
        phase2-evaluate.py \
            --dataset      ${emx_file} \
            --labels       ${labels_file} \
            --model-config ${baseDir}/models.json \
            --model        ${params.phase2_model} \
            --gene-sets    ${gmt_file} \
            --n-jobs       1 \
            --logdir       .
        """
}



/**
 * The phase2_select process takes a list of gene sets selected by phase 1, as well
 * as subset scores from phase2_evaluate, and selects candidate genes for each
 * gene set.
 */
process phase2_select {
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}"

    input:
        tuple val(dataset), val(gmt_name), path("phase1-genesets.txt")
        tuple val(dataset), val(gmt_name), path(scores)

    output:
        tuple val(dataset), val(gmt_name), path("phase2-genesets.txt"), emit: gene_sets
        path("*.png"), optional: true

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
 * The phase2_rf process takes a list of gene sets selected by phase 1 and
 * uses the feature importances of a random forest to select candidate genes
 * for each gene set.
 */
process phase2_rf {
    publishDir "${params.output_dir}/${dataset}", mode: "copy"
    tag "${dataset}/${gmt_name}"

    input:
        tuple val(dataset), path(emx_file), path(labels_file), val(gmt_name), path(gmt_file)
        tuple val(dataset), val(gmt_name), path("phase1-genesets.txt")

    output:
        tuple val(dataset), val(gmt_name), path("phase2-rf-genesets.txt"), emit: gene_sets
        path("*.png"), optional: true

    script:
        """
        echo "#TRACE chunks=${params.chunks}"
        echo "#TRACE gmt_lines=`cat phase1-genesets.txt | wc -l`"
        echo "#TRACE gmt_genes=`cat phase1-genesets.txt | wc -w`"
        echo "#TRACE n_rows=`tail -n +1 ${emx_file} | wc -l`"
        echo "#TRACE n_cols=`head -n +1 ${emx_file} | wc -w`"

        phase2-rf.py \
            --dataset   ${emx_file} \
            --labels    ${labels_file} \
            --gene-sets phase1-genesets.txt \
            --n-jobs    1 \
            --threshold ${params.phase2_rf_threshold} \
            ${params.phase2_rf_visualize ? "--visualize" : ""}
        """
}
