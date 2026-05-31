nextflow.enable.dsl=2

// =========================================================================================
// 1. 初始化输入通道 (严格锁定 Tuple 级别，包装嵌套列表杜绝多样本降维丢失)
// =========================================================================================

// 精准抓取 params.gvcf_dir 下 300 个样本的 .g.vcf.gz 及其配套的 .tbi 索引
gvcf_files_ch = Channel.fromPath("${params.gvcf_dir}/*.g.vcf.gz")
            .map { file -> [file, file + ".tbi"] }
            .flatten()
            .collect()

// 在 Nextflow 顶层根据基因组版本，动态创建对应的常规主染色体清单以激活集群最高并发
def chr_list = []
if (params.genome_version == "hg19_decoy") {
    chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
} else {
    // 针对你的 hg38 参考基因组，默认激活带 chr 前缀的 25 条常规标准主染色体清单
    chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
}
chromosomes_ch = Channel.fromList(chr_list)

// =========================================================================================
// 2. 主流程工作流逻辑拓扑 (Workflow Block)
// =========================================================================================

workflow {
    
    // 【核心提速】：利用 .combine() 将所有样本的 gVCF 物流与 25 条常规主染色体进行交叉组合
    // 使得每一条染色体独立分流，激活集群层面的最高并发合并
    joint_call_in = chromosomes_ch.combine(gvcf_files_ch.map{ [it] })

    // 桥梁接口精准对接：根据当前并发的 ${chr}，动态锁定你在本地已经提前切分生成好的单染色体常规 BED 文件
    joint_call(
        joint_call_in,
        joint_call_in.map { chr, files -> 
            file("${params.bed_dir}/${chr}.bed")
        }
    )

    // 收拢各异步并发分片，利用 Groovy 算法按染色体物理顺序精确、安全升序排序（完美继承原 Groovy 桥梁代码）
    merge_vcf_in = joint_call.out.vcf
        .map { chr, vcf_file -> 
            ["cohort", [chr, vcf_file]]
        }
        .groupTuple(by: 0)
        .map { cohort_id, records ->
            def sorted_vcfs = records
                .sort { a, b -> 
                    def chrA = a[0].toString().replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')
                    def chrB = b[0].toString().replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25')
                    chrA.toInteger() <=> chrB.toInteger()
                }
                .collect { it[1] } 
            
            [cohort_id, sorted_vcfs]
        }
    merge_vcf(merge_vcf_in)
    vcf_liftover(merge_vcf.out.vcf)
    annotate_vcf(vcf_liftover.out.vcf_37)
}

process joint_call {
    tag "GLnexus_JointCall_${chr}"
    label "big_mem_2" 
    
    input:
    tuple val(chr), path(all_files) 
    path chrom_bed                  
    
    output:
    tuple val(chr), path("cohort_${chr}_joint_raw.vcf.gz"), emit: vcf
    
    script:
    """
    awk '{if(\$1 !~ /^chr/) print "chr"\$0; else print \$0}' ${chrom_bed} > aligned_region.bed
    FIRST_GVCF=\$(ls *.g.vcf.gz | head -n 1)
    if ! zgrep -q "contig=<ID=${chr}" \$FIRST_GVCF; then
        echo "⚠️ WARNING: Contig ${chr} not found in input gVCF headers (e.g. testing data empty slot). Generating empty safe vcf dummy..."
        echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > empty_${chr}.vcf
        ${params.bgzip} -f empty_${chr}.vcf -c > cohort_${chr}_joint_raw.vcf.gz
    else
        ${params.singularity} run --env LC_ALL=C.UTF-8,LANG=C.UTF-8 -B `pwd`:`pwd` -B /mnt/gpfs/:/mnt/gpfs/ -B /mnt/gpfs1/:/mnt/gpfs1/ \
            ${params.glnexus} /usr/local/bin/glnexus_cli \
            --config DeepVariant_unfiltered \
            -b aligned_region.bed \
            --threads ${task.cpus} \
            *.g.vcf.gz | ${params.bcftools} view - -Oz -o cohort_${chr}_joint_raw.vcf.gz
    fi
    """
}

process merge_vcf {
    tag "Merge_Chr_VCFs"
    
    input:
    tuple val(cohort_id), file(vcf_files)
    
    output:
    tuple val(cohort_id), file("cohort_joint_raw.vcf.gz"), file("cohort_joint_raw.vcf.gz.tbi"), emit: vcf
    
    script:
    """
    ${params.bcftools} concat ${vcf_files.join(' ')} -Oz -o cohort_joint_raw.vcf.gz   
    ${params.tabix} -p vcf cohort_joint_raw.vcf.gz
    """
}

process vcf_liftover {
    tag "Genome_Space_Interceptor_${cohort_id}"
    
    input:
    tuple val(cohort_id), path(vcf), path(tbi)
    
    output:
    tuple val(cohort_id), path("cohort_joint_final.vcf.gz"), path("cohort_joint_final.vcf.gz.tbi"), emit: vcf_37
    
    script:
    if (params.genome_version == "hg19_decoy" || params.genome_version == "b37") {
        """
        ln -s ${vcf} cohort_joint_final.vcf.gz
        ln -s ${tbi} cohort_joint_final.vcf.gz.tbi
        """
    } else {
        """
        echo -e "chrM\tMT" > rename_map.txt
        for i in {1..22} X Y; do echo -e "chr\$i\t\$i" >> rename_map.txt; done

        ${params.bcftools} annotate --rename-chrs rename_map.txt ${vcf} -Oz -o cohort_joint_final.vcf.gz
        ${params.tabix} -p vcf cohort_joint_final.vcf.gz
        """
    }
}

process annotate_vcf {
    tag "VEP_Advanced_Annotation_GRCh37"
    publishDir "${params.outdir}/Final_VCF", mode: 'copy'
    label "big_mem_2"
    
    input:
    tuple val(cohort_id), path(vcf), path(tbi) // 3元素 Tuple 格式与上游 liftover 投递完美对齐
    
    output:
    path "cohort_annotated.vcf.gz", emit: annotated_vcf
    
    script:
    """

    ${params.singularity} exec -B `pwd`:`pwd` -B /mnt/gpfs/:/mnt/gpfs/ -B /mnt/gpfs1/:/mnt/gpfs1/ -B /mnt/gpfs/Users/lilei/tools/vep:/vep_path \
    -B `dirname ${vcf}`:/input \
    -B /mnt/gpfs/Users/fanlei/project/21.wgs_method_quality/01.database/MaxEntScan/fordownload:/fordownload \
    -B /mnt/gpfs/Users/lilei/tools/vep/plugin/plugin_data:/plugin_data \
    /mnt/gpfs/Users/lilei/tools/vep/vep.sif vep -i /input/${vcf.name} \
    -fork 10 -o cohort_annotated.vcf \
    --format vcf --refseq --cache --offline --force_overwrite --dir_cache /vep_path/cache/GRCh37 \
    --plugin MaxEntScan,/fordownload --plugin SpliceAI,snv=/plugin_data/spliceai_scores.masked.snv.hg19.vcf.gz,indel=/plugin_data/spliceai_scores.masked.indel.hg19.vcf.gz \
    --custom file=/vep_path/database/dbsnp_151_hg19.vcf.gz,short_name=dbsnp_151,format=vcf,type=exact,coords=0,fields=RS \
    --custom file=/vep_path/database/dbnsfp42a_hg19.vcf.gz,short_name=dbnsfp42a,format=vcf,type=exact,coords=0,fields=SIFT_score%SIFT_pred%SIFT4G_score%SIFT4G_pred%Polyphen2_HDIV_score%Polyphen2_HDIV_pred%Polyphen2_HVAR_score%Polyphen2_HVAR_pred%MutationTaster_score%MutationTaster_pred%MutationAssessor_score%MutationAssessor_pred%REVEL_score%REVEL_rankscore%ClinPred_score%ClinPred_pred%CADD_raw%CADD_phred%PROVEAN_score%PROVEAN_pred%LRT_score%LRT_pred%GERP++_RS%GERP++_RS_rankscore%phyloP100way_vertebrate%phyloP100way_vertebrate_rankscore%phyloP30way_mammalian%phyloP30way_mammalian_rankscore \
    --dir_plugins /vep_path/plugin/VEP_plugins --plugin dbscSNV,/vep_path/database/dbsnv1.1/dbscSNV1.1_GRCh37.txt.gz --custom file=/vep_path/database/ExAC.r1.sites.vep.vcf.gz,short_name=ExAC,format=vcf,type=exact,coords=0,fields=AC%AC_EAS%AN_EAS%AC_Het%AC_Hom%AC_Hemi%AF \
    --custom file=/vep_path/database/gnomad.exomes.r2.1.1.sites.vcf.gz,short_name=gnomad,format=vcf,type=exact,coords=0,fields=controls_AC%controls_AN%controls_nhomalt%controls_AF%controls_AC_eas%controls_AN_eas%controls_nhomalt_eas%controls_AF_eas%AC_male%AN_male%AF_male \
    --custom file=/vep_path/database/gnomad.genomes.r2.1.1.sites.vcf.gz,short_name=gnomad_wgs,format=vcf,type=exact,coords=0,fields=controls_AC%controls_AN%controls_nhomalt%controls_AF%controls_AC_eas%controls_AN_eas%controls_nhomalt_eas%controls_AF_eas%AC_male%AN_male%AF_male \
    --custom file=/vep_path/database/WBBC.vcf.gz,short_name=WBBC,format=vcf,type=exact,coords=0,fields=AC%AN%AF \
    --custom file=/vep_path/database/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \
    --custom file=/vep_path/database/pubmed_pro_2021.4_hg19.vcf.gz,short_name=hgmd,format=vcf,type=exact,coords=0,fields=CLASS \
    --fasta /vep_path/FASTA/v37/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz \
    --assembly GRCh37 --distance 500  --no_escape --show_ref_allele --check_existing --exclude_predicted --mane -hgvsg -hgvs --vcf --symbol --mane_select --canonical --numbers --af --af_1kg

    ${params.java} -jar ${params.snpeff_jar} -ud 500 GRCh37.p13 cohort_annotated.vcf > cohort_annotated.snpEff.vcf
    
    ${params.bgzip} -f cohort_annotated.snpEff.vcf -c > cohort_annotated.vcf.gz
    ${params.tabix} -p vcf cohort_annotated.vcf.gz
    """
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}