#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.vcfs = 'path/to/vcf-dir/*.vcf.gz'
params.vep_cachedir = "~/vep_cache"
params.parse_script  = 'parse_annotations.py'

process NormalizeVCF {
    input:
        path vcf_file
    output:
        path "*.vcf.gz", emit: normalized_vcf
    script:
        """
        prefix=\$(basename ${vcf_file} .vcf.gz)
        # Filter for PASS variants, then normalize
        bcftools view -f PASS ${vcf_file} | bcftools norm -m -any -Oz -o \$prefix.normalized.vcf.gz
        """
}

process AnnotateVariants {
    publishDir "./vcfs-annotated", mode: 'symlink'
    input:
        tuple path(vcf_file), path(vep_cache)
    output:
        path("*.annotated.vcf.gz"), emit: annotated_vcf
    script:
        """
        # Extract prefix by removing the .normalized.vcf.gz suffix
        prefix=\$(basename ${vcf_file} .normalized.vcf.gz)
        # Run VEP and pipe the output through bgzip to produce a gzipped VCF
        vep -i ${vcf_file} -o /dev/stdout --force_overwrite --stats_file \$prefix.summary.html --vcf --symbol --canonical --protein --offline --assembly GRCh38 --cache --dir ${vep_cache} | bgzip -c > \$prefix.annotated.vcf.gz
        """
}

process ParseAnnotations {
    publishDir "./annotations-parsed", mode: 'symlink'
    input:
        tuple path(annotated_vcf), file(parse_script)
    output:
        path "*.tsv", emit: summary_tsv
    script:
        """
        prefix=\$(basename ${annotated_vcf} .annotated.vcf.gz)
        tabix -p vcf ${annotated_vcf}
        python ${parse_script} ${annotated_vcf} \${prefix}.summary.tsv \${prefix}.summary.protein.tsv
        """
}

workflow {
    input_vcf      = Channel.fromPath(params.vcfs)
    parse_script_ch = Channel.value(file(params.parse_script))

    // Run normalization on the input VCFs
    normalized_ch = NormalizeVCF(input_vcf)
    
    // Pair each normalized VCF with the VEP cache directory to create a tuple:
    annotated_input = normalized_ch.map { norm_vcf -> tuple(norm_vcf, file(params.vep_cachedir)) }
    
    // Pass the tuple (normalized VCF, VEP cache) to AnnotateVariants
    annotated_ch = AnnotateVariants(annotated_input)
    
    // Combine the annotated output with the parse script and run parsing
    annotated_with_parse = annotated_ch.combine(parse_script_ch)
    ParseAnnotations(annotated_with_parse)
}
