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
params.maf_min = 0

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
        def raw_trait  = row.brain_trait ?: brain_name
        // If brain_trait is an absolute path, extract just the filename stem
        def brain_trait = raw_trait.tokenize('/').last().replace('.txt','')
        def infl_name  = row.inflammation_gwas.tokenize('/').last().replace('.txt','').replace('_processedGRCh38_clean','').replace('_processedGRCh38_N','').replace('_processedGRCh38','')
        tuple(
            row.inflammation_gwas,
            row.brain_gwas,
            "${infl_name}_${brain_trait}"
        )
    }
    .set { pairwise_inputs }


// ------------------
// MTAG process
// ------------------
process RUN_MTAG {
    tag { brain_trait }

    input:
    tuple val(infl_gwas), val(brain_gwas), val(brain_trait)

    output:
    path "mtag_output"

    errorStrategy 'ignore'
    publishDir "${params.outdir}/${brain_trait}_MTAG", mode: 'copy'

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
//process FILTER_GWS {
  //  tag { trait }
    //conda "/home/law22/data/conda_envs/mtag"

   // input:
    //path mtag_dir
    //val trait
    //val infl_name

    //output:
    //path "*.GWS.txt"
    //path "*.GWS.clump.tsv"

    //publishDir "${params.outdir}/${trait}_MTAG", mode: 'copy'

    //script:
    //"""
    
    //# Rename and reformat trait_1 (inflammation) and trait_2 (brain)
    //for f in ${mtag_dir}/*_trait_1.txt; do
     //   base=\$(basename \$f _trait_1.txt)
       // out="${infl_name}_MTAG.txt"
        //awk 'BEGIN{OFS="\\t"} NR==1{print "SNP","CHR","BP","A1","A2","Z","N","FRQ","BETA","SE","P"} NR>1{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$12}' \$f > \$out
    //done

    //for f in ${mtag_dir}/*_trait_2.txt; do
      //  out="${trait}_MTAG.txt"
       // awk 'BEGIN{OFS="\\t"} NR==1{print "SNP","CHR","BP","A1","A2","Z","N","FRQ","BETA","SE","P"} NR>1{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$12}' \$f > \$out
    //done

//# Filter for GWS hits on brain trait
  //  mtagfile="${trait}_MTAG.txt"
    //if [[ -f "\$mtagfile" ]]; then
      //  /home/law22/data/conda_envs/mtag/bin/python ${projectDir}/bin/filter_gwas.py \$mtagfile --out_prefix ${trait}
    //else
      //  echo "ERROR: \$mtagfile not found" >&2
       // exit 1
    //fi
   // """
//}

// ----------------------
// Workflow definition
// ----------------------
workflow {
    RUN_MTAG(pairwise_inputs)

    RUN_MTAG.out
        .merge(pairwise_inputs.map { infl, brain, trait ->
            def infl_name = infl.tokenize('/').last().replace('.txt','')
            [trait, infl_name]
        })
        .set { filter_inputs }

  //  FILTER_GWS(
    //    filter_inputs.map { dir, trait, infl -> dir },
      //  filter_inputs.map { dir, trait, infl -> trait },
        //filter_inputs.map { dir, trait, infl -> infl }
    //)
}