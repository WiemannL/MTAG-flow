#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ------------------
// Parameters
// ------------------
params.samples  = "samples.csv"
params.outdir   = "./results"
params.mtag_dir = "./mtag"

params.snp = "SNP"
params.a1  = "A1"
params.a2  = "A2"
params.z   = "Z"
params.n   = "N"
params.maf_min = 0.01

log.info """\
MTAG-flow - Pairwise MTAG Pipeline
===================================
GWAS CSV     : ${params.samples}
Output dir   : ${params.outdir}
MTAG code    : ${params.mtag_dir}
""".stripIndent()


// ------------------
// Load CSV
// Expected columns:
// inflammation_gwas, brain_gwas, brain_trait
// ------------------
Channel
    .fromPath(params.samples)
    .splitCsv(header:true)
    // trim spaces in header and values
    .map { row -> row.collectEntries { k,v -> [(k.trim()): v?.trim()] } }
    .map { row ->
        def brain_name = row.brain_gwas.tokenize('/').last().replace('.txt','')
        def brain_trait = row.brain_trait ?: brain_name
        tuple(
            row.inflammation_gwas,
            row.brain_gwas,
            brain_trait,
            "${params.outdir}/${brain_trait}_MTAG"
        )
    }
    .set { pairwise_inputs }


// ------------------
// MTAG process
// ------------------
process RUN_MTAG {
    tag { brain_trait }

    input:
    tuple val(infl_gwas), val(brain_gwas), val(brain_trait), val(outdir)

    output:
    path "mtag_output"

    errorStrategy 'ignore'
    publishDir "${outdir}", mode: 'copy'

    script:
    """
    mkdir -p mtag_output

    python ${params.mtag_dir}/mtag.py \
        --sumstats ${infl_gwas},${brain_gwas} \
        --out mtag_output/${brain_trait} \
        --snp_name ${params.snp} \
        --a1_name ${params.a1} \
        --a2_name ${params.a2} \
        --z_name ${params.z} \
        --n_name ${params.n} \
        --maf_min ${params.maf_min} \
        --stream_stdout \
        --force
    """
}

// ----------------------
// FILTER GWAS process
// ----------------------
process FILTER_GWS {
    tag { trait }

    input:
    path mtag_dir, emit: mtag_output_dir
    val trait

    output:
    path "*.GWS.txt"
    path "*.GWS.clump.tsv"

    publishDir "${params.outdir}/${trait}_MTAG", mode: 'copy'

    script:
    """
    # Find the formatted file (assume convention: <trait>_phenotype_formatted.txt)
    mtagfile=\$(find \${mtag_dir} -name "*_phenotype_formatted.txt" | head -n1)
    if [[ -f "\$mtagfile" ]]; then
        bin/filter_gwas_gws.py \$mtagfile --out_prefix \${trait}
    else
        echo "ERROR: No *_phenotype_formatted.txt file found in \${mtag_dir}" >&2
        exit 1
    fi
    """
}

// ----------------------
// Workflow definition
// ----------------------
workflow {
    // Run MTAG for each pair
    RUN_MTAG(pairwise_inputs)

    // Filter significant loci for each MTAG output
    // Passes "mtag_output" dir and trait to the filter process
    pairwise_inputs
        .map { infl, brain, trait, outdir -> [file("${outdir}/mtag_output"), trait] }
        .set { filter_inputs }

    FILTER_GWS(filter_inputs)
}
