nextflow.enable.dsl = 2

params.project_root = "/mnt/gpfs1/Users/yangjinxurong/pipeline/WGS/snp/"

def chr_list = []
if (params.genome_version == 'hg19_decoy' || params.genome_version == 'b37') {
    chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
                '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
} else {
    chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                'chrX', 'chrY', 'chrM']
}

def chrSortKey(chr) {
    def c = chr.toString().replaceFirst(/^chr/, '')
    if (c == 'X') return 23
    if (c == 'Y') return 24
    if (c == 'M' || c == 'MT') return 25
    return c.toInteger()
}

workflow {
    chromosomes_ch = Channel
        .fromList(chr_list)
        .map { chr -> tuple(chr, chrSortKey(chr)) }

    if (params.skip_header_normalization) {
        normalized_files_ch = Channel
            .fromPath("${params.reheadered_gvcf_dir}/*.g.vcf.gz", checkIfExists: true)
            .map { gvcf -> [gvcf, file("${gvcf}.tbi")] }
            .flatten()
            .collect()
    } else {
        DISCOVER_GVCFS()
        BUILD_CONTIG_TEMPLATE(DISCOVER_GVCFS.out.manifest)

        raw_gvcf_ch = DISCOVER_GVCFS.out.manifest
            .splitCsv(header: true, sep: "\t")
            .map { row -> tuple(row.sample_id, file(row.gvcf)) }

        REHEADER_GVCF(raw_gvcf_ch.combine(BUILD_CONTIG_TEMPLATE.out.contigs))

        normalized_files_ch = REHEADER_GVCF.out.normalized
            .map { sample_id, gvcf, tbi -> [gvcf, tbi] }
            .flatten()
            .collect()
    }

    joint_call_in = chromosomes_ch
        .combine(normalized_files_ch.map { [it] })
        .map { chr, sort_key, files ->
            tuple(chr, sort_key, files, file("${params.bed_dir}/${chr}.bed"))
        }

    JOINT_CALL(joint_call_in)
    FILTER_JOINT_VCF(JOINT_CALL.out.raw_vcf)
    LIFTOVER_CHR(FILTER_JOINT_VCF.out.filtered_vcf)
    ANNOTATE_CHR(LIFTOVER_CHR.out.vcf)

    merge_annotated_in = ANNOTATE_CHR.out.annotated
        .map { chr, sort_key, vcf, tbi, report ->
            ['cohort', [sort_key, vcf, tbi]]
        }
        .groupTuple(by: 0)
        .map { cohort_id, records ->
            def sorted = records.sort { a, b -> a[0] <=> b[0] }
            [cohort_id, sorted.collect { it[1] }, sorted.collect { it[2] }]
        }

    MERGE_ANNOTATED_VCF(merge_annotated_in)
    // MERGE_ANNOTATED_VCF emits the merged VCF as `annotated_vcf`.
    STANDARDIZE_ANNOTATION(MERGE_ANNOTATED_VCF.out.annotated_vcf)
}

process DISCOVER_GVCFS {
    tag "discover_gvcfs"

    output:
    path "gvcf_manifest.tsv", emit: manifest

    script:
    """
    set -euo pipefail

    ${params.python} ${baseDir}/bin/reheader.py discover-gvcfs \\
      --gvcf-dir "${params.gvcf_dir}" \\
      --out-manifest gvcf_manifest.tsv
    """
}

process BUILD_CONTIG_TEMPLATE {
    tag "build_contig_template"

    input:
    path manifest

    output:
    path "contig_template.txt", emit: contigs
    path "contig_template.tsv", emit: template_tsv
    path "input_contig_audit.tsv", emit: audit
    path "excluded_contigs.tsv", emit: excluded
    path "reference_resolution.tsv", emit: reference_resolution

    script:
    """
    set -euo pipefail

    ${params.python} ${baseDir}/bin/reheader.py build-template \\
      --manifest ${manifest} \\
      --reference-fai "${params.reference_fai}" \\
      --reference-fasta "${params.reference_fasta}" \\
      --auto-detect-reference-fai "${params.auto_detect_reference_fai}" \\
      --contig-mode "${params.contig_mode}" \\
      --out-contigs contig_template.txt \\
      --out-template-tsv contig_template.tsv \\
      --out-audit-tsv input_contig_audit.tsv \\
      --out-excluded-tsv excluded_contigs.tsv \\
      --out-reference-tsv reference_resolution.tsv
    """
}

process REHEADER_GVCF {
    tag "${sample_id}"
    publishDir "${params.reheadered_gvcf_dir}/", mode: "copy", pattern: "*.g.vcf.gz*"

    input:
    tuple val(sample_id), path(gvcf), path(template_contigs)

    output:
    tuple val(sample_id), path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"), emit: normalized
    path "${sample_id}.header_report.tsv", emit: report

    script:
    """
    set -euo pipefail
    input_gvcf=\$(readlink -f ${gvcf})
    tmp_gvcf="${sample_id}.reheader.tmp.g.vcf.gz"

    ${params.python} ${baseDir}/bin/reheader.py rewrite-header \\
      --sample-id "${sample_id}" \\
      --gvcf "\${input_gvcf}" \\
      --template-contigs ${template_contigs} \\
      --out-header ${sample_id}.header.txt \\
      --report-tsv ${sample_id}.header_report.tsv

    ${params.bcftools} reheader \\
      -h ${sample_id}.header.txt \\
      -o "\${tmp_gvcf}" \\
      "\${input_gvcf}"

    rm -f ${sample_id}.g.vcf.gz ${sample_id}.g.vcf.gz.tbi
    mv "\${tmp_gvcf}" ${sample_id}.g.vcf.gz

    ${params.tabix} -f -p vcf ${sample_id}.g.vcf.gz
    """
}

process JOINT_CALL {
    tag "GLnexus_JointCall_${chr}"
    publishDir "${params.outdir}/JointCall", mode: 'copy'

    input:
    tuple val(chr), val(sort_key), path(all_files), path(chrom_bed)

    output:
    tuple val(chr), val(sort_key), path("cohort_${chr}_joint_raw.vcf.gz"), emit: raw_vcf

    script:
    """
    rm -rf GLnexus.DB
    rm -f cohort_${chr}_joint_raw.bcf cohort_${chr}_joint_raw.vcf.gz cohort_${chr}_joint_raw.vcf.gz.tbi empty_${chr}.vcf

    emit_empty_vcf() {
        echo -e "##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" > empty_${chr}.vcf
        ${params.bgzip} -f empty_${chr}.vcf -c > cohort_${chr}_joint_raw.vcf.gz
    }

    awk '{if(\$1 !~ /^chr/) print "chr"\$0; else print \$0}' ${chrom_bed} > aligned_region.bed
    ls -1 *.g.vcf.gz | sort > glnexus_inputs.list
    FIRST_GVCF=\$(head -n 1 glnexus_inputs.list)

    if ! zgrep -q "contig=<ID=${chr}" \$FIRST_GVCF; then
        echo "WARNING: Contig ${chr} not found in input gVCF headers. Generating empty safe VCF."
        emit_empty_vcf
    else
        set +e
        ${params.singularity} run --env LC_ALL=C.UTF-8,LANG=C.UTF-8 \\
            -B `pwd`:`pwd` -B /mnt/gpfs/:/mnt/gpfs/ -B /mnt/gpfs1/:/mnt/gpfs1/ \\
            ${params.glnexus} /usr/local/bin/glnexus_cli \\
            --config DeepVariant_unfiltered \\
            -b aligned_region.bed \\
            --mem-gbytes ${params.glnexus_mem_gbytes} \\
            --threads ${task.cpus} \\
            --list glnexus_inputs.list > cohort_${chr}_joint_raw.bcf
        glnexus_status=\$?
        set -e

        if [[ \$glnexus_status -ne 0 ]]; then
            exit \$glnexus_status
        elif [[ -s cohort_${chr}_joint_raw.bcf ]]; then
            ${params.bcftools} view cohort_${chr}_joint_raw.bcf -Oz -o cohort_${chr}_joint_raw.vcf.gz
        else
            echo "WARNING: GLnexus produced no variants for ${chr}. Generating empty safe VCF."
            emit_empty_vcf
        fi
    fi

    [[ -s cohort_${chr}_joint_raw.vcf.gz ]]
    """
}

process FILTER_JOINT_VCF {
    tag "Filter_JointCall_${chr}"
    publishDir "${params.outdir}/JointCall", mode: 'copy'

    input:
    tuple val(chr), val(sort_key), path(vcf)

    output:
    tuple val(chr), val(sort_key), path("cohort_${chr}_joint_filter.vcf.gz"), path("cohort_${chr}_joint_filter.vcf.gz.tbi"), emit: filtered_vcf

    script:
    """
    ${params.bcftools} filter -i '${params.joint_filter_expr}' -Oz -o cohort_${chr}_joint_filter.vcf.gz ${vcf}
    ${params.tabix} -p vcf cohort_${chr}_joint_filter.vcf.gz
    """
}

process LIFTOVER_CHR {
    tag "Genome_Space_Interceptor_${chr}"
    publishDir "${params.outdir}/JointCall/per_chromosome_final", mode: 'copy'

    input:
    tuple val(chr), val(sort_key), path(vcf), path(tbi)

    output:
    tuple val(chr), val(sort_key), path("cohort_${chr}_joint_final.vcf.gz"), path("cohort_${chr}_joint_final.vcf.gz.tbi"), emit: vcf

    script:
    if (params.genome_version == "hg19_decoy" || params.genome_version == "b37") {
        """
        ln -s ${vcf} cohort_${chr}_joint_final.vcf.gz
        ln -s ${tbi} cohort_${chr}_joint_final.vcf.gz.tbi
        """
    } else {
        """
        echo -e "chrM\tMT" > rename_map.txt
        for i in {1..22} X Y; do echo -e "chr\$i\t\$i" >> rename_map.txt; done

        ${params.bcftools} annotate --rename-chrs rename_map.txt ${vcf} -Oz -o cohort_${chr}_joint_final.vcf.gz
        ${params.tabix} -p vcf cohort_${chr}_joint_final.vcf.gz
        """
    }
}

process ANNOTATE_CHR {
    tag "Annotate_${chr}"
    publishDir "${params.outdir}/Anno/per_chromosome", mode: 'copy'

    input:
    tuple val(chr), val(sort_key), path(vcf), path(tbi)

    output:
    tuple val(chr), val(sort_key), path("cohort_${chr}_annotated.vcf.gz"), path("cohort_${chr}_annotated.vcf.gz.tbi"), path("cohort_${chr}.vep_summary.html"), emit: annotated
    path "versions.yml", emit: versions

    script:
    // The available gnomAD v4.1 genome resources cover autosomes/X/Y only.
    // chrM is still annotated by VEP, but has no gnomAD custom source.
    def gnomadChr = chr.toString().startsWith('chr') ? chr : "chr${chr}"
    def gnomadVcf = "${params.gnomad_v4_genomes_vep_dir}/${params.gnomad_v4_genomes_filename_prefix}.${gnomadChr}.nochr.vcf.gz"
    def gnomadCustom = (chr.toString() in ['chrM', 'M', 'MT']) ? '' : "--custom file=${gnomadVcf},short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF%AF_eas"
    """
    source /mnt/gpfs/Users/yangjinxurong/software/miniconda3/bin/activate /mnt/gpfs/Users/yangjinxurong/software/miniconda3/envs/vep113_env

    export PERL5LIB="${params.loftee_perl_compat_dir}:${params.loftee_dir}:\${PERL5LIB:-}"

    vep -i ${vcf} \\
        -o cohort_${chr}.vep.vcf \\
        --format vcf --merged --cache --offline --force_overwrite \\
        --dir_cache ${params.vep_homo_113} \\
        --dir_plugins ${params.vep_plugin} \\
        --fasta ${params.vep_fasta} \\
        --fork ${task.cpus} \\
        --assembly GRCh38 --distance 500 --no_escape --show_ref_allele \\
        --check_existing --exclude_predicted --mane --symbol --mane_select \\
        --canonical --numbers --af --af_1kg -hgvsg -hgvs --vcf \\
        --stats_file cohort_${chr}.vep_summary.html \\
        --plugin MaxEntScan,${params.maxentscan_dir} \\
        --plugin SpliceAI,snv=${params.spliceai_snv},indel=${params.spliceai_indel} \\
        --plugin LoF,loftee_path:${params.loftee_dir},gerp_bigwig:${params.loftee_gerp_bigwig},human_ancestor_fa:${params.loftee_human_ancestor_fa},conservation_file:${params.loftee_conservation_file} \\
        --plugin AlphaMissense,file=${params.alphamissense_tsv} ${gnomadCustom} \\
        --custom file=${params.clinvar_vcf},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF

    ${params.snpeff_java} -jar ${params.snpeff_jar} \\
        -ud 500 GRCh38.105 \\
        cohort_${chr}.vep.vcf > cohort_${chr}.vep.snpEff.vcf

    ${params.bgzip} -f -c cohort_${chr}.vep.snpEff.vcf > cohort_${chr}_annotated.vcf.gz
    ${params.tabix} -p vcf cohort_${chr}_annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$( vep --help | grep -oP 'Versions:.*' | cut -f2 -d' ' )
        snpeff: \$(${params.snpeff_java} -jar ${params.snpeff_jar} -version | grep "SnpEff" | cut -d' ' -f2 )
    END_VERSIONS
    """
}

process MERGE_ANNOTATED_VCF {
    tag "Merge_Annotated_VCFs"
    publishDir "${params.outdir}/Anno", mode: 'copy'

    input:
    tuple val(cohort_id), path(vcf_files), path(tbi_files)

    output:
    path "cohort_annotated.vcf.gz", emit: annotated_vcf
    path "cohort_annotated.vcf.gz.tbi", emit: annotated_tbi

    script:
    """
    ${params.bcftools} concat ${vcf_files.join(' ')} -Oz -o cohort_annotated.vcf.gz
    ${params.tabix} -p vcf cohort_annotated.vcf.gz
    """
}

process STANDARDIZE_ANNOTATION {
    tag "AnnotationStandardization"
    label "annotation"
    publishDir "${params.outdir}/annotation_standardization", mode: "copy", overwrite: true

    input:
    path annotated_vcf

    output:
    path "results/variant_annotation.tsv", emit: variant_annotation
    path "results/gene_variant_map.tsv", emit: gene_variant_map
    path "results", emit: results
    path "reports", emit: reports
    path "logs", emit: logs

    script:
    """
    set -euo pipefail
    : "annotation_parser_version=gnomADg_AF_eas_csq_v1"
    ${params.python} ${params.project_root}/scripts/standardize_variant_annotation.py \\
      --annotated-vcf ${annotated_vcf} \\
      --genome-version "${params.genome_version}" \\
      --max-variants "${params.annotation_max_variants}" \\
      --rare-af-threshold "${params.rare_af_threshold}" \\
      --outdir .
    """
}

workflow.onError {
    println "Pipeline stopped: ${workflow.errorMessage}"
}
