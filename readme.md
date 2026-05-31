cp /mnt/gpfs/Dataset/04.project/Bioinfo/92.KW/01.analysis/260422122712_Cetus_C007-B-260421-260422122725/LC260198-S01-D03-CUDI56/LC260198-S01-D03-CUDI56.combined.g.vcf.gz .
cp /mnt/gpfs/Dataset/04.project/Bioinfo/92.KW/01.analysis/260422122712_Cetus_C007-B-260421-260422122725/LC260198-S01-D03-CUDI56/LC260198-S01-D03-CUDI56.combined.g.vcf.gz.tbi .
cp /mnt/gpfs/Dataset/04.project/Bioinfo/92.KW/01.analysis/260422122712_Cetus_C007-B-260421-260422122725/LC260196-S01-D03-CUDI54/LC260196-S01-D03-CUDI54.combined.g.vcf.gz .
cp /mnt/gpfs/Dataset/04.project/Bioinfo/92.KW/01.analysis/260422122712_Cetus_C007-B-260421-260422122725/LC260196-S01-D03-CUDI54/LC260196-S01-D03-CUDI54.combined.g.vcf.gz.tbi .
cp /mnt/gpfs/Dataset/04.project/Bioinfo/92.KW/01.analysis/260422122712_Cetus_C007-B-260421-260422122725/LC260192-S01-D01-BUDI18/LC260192-S01-D01-BUDI18.combined.g.vcf.gz .
cp /mnt/gpfs/Dataset/04.project/Bioinfo/92.KW/01.analysis/260422122712_Cetus_C007-B-260421-260422122725/LC260192-S01-D01-BUDI18/LC260192-S01-D01-BUDI18.combined.g.vcf.gz.tbi .

zcat LC260192-S01-D01-BUDI18.combined.g.vcf.gz | grep "#" | head -n 1000 | grep "<NON_REF>"
zcat LC260196-S01-D03-CUDI54.combined.g.vcf.gz | grep "#" | less -SN


###
占位符不一致： 在标准 GATK 体系中，GenotypeGVCFs 默认在矩阵中疯狂寻找带有 <NON_REF> 字符串的记录 。
DeepVariant/Sentieon 的输出习惯： 如果你的 gVCF 是用 DeepVariant 或者某些特定商业版软件呼出来的，它们遵循的是官方 VCF 4.2/4.3 规范，通常会使用 <*> 替代 <NON_REF> 作为“任意等位基因（Any alternative allele）”的占位符。
GATK 逻辑崩溃： 当 GATK 4.3.0.0 遇到了带有 <*> 的 gVCF 时，它的底层代码如果不经过特殊配置，就无法识别这个符号，从而认为你输入的是普通 VCF，随即抛出了 must contain <NON_REF> as an allele 的 User Error 报错。
### 使用deepvariant配套的glnexus


To activate this environment, use:

    mamba activate glnexus_env

Or to execute a single command in this environment, use:

    mamba run -n glnexus_env mycommand

系统版本过低不能使用conda的glnexus的包
报错：glnexus_cli: /lib64/libc.so.6: version `GLIBC_2.18' not found 
于是下载glnexus_v1.4.1.sif
```sh
singularity pull glnexus_v1.4.1.sif docker://ghcr.nju.edu.cn/dnanexus-rnd/glnexus:v1.4.1 
/mnt/gpfs1/Users/yangjinxurong/software/singularity/glnexus_v1.4.1.sif
```

[300样本 gVCF + TBI 捆绑] ──> chromosomes_ch.combine() 展开 ──> [25条常规染色体并行分流]
                                                                        │
  ┌─────────────────────────────────────────────────────────────────────┘
  ▼
process joint_call (Singularity 容器运行 GLnexus v1.4.1)
  │  └─> 接收: tuple val(chr), path(all_files) 和 path(chrom_bed)
  │  └─> 动作: -b 喂入单染色体常规非N掩码 BED，绕过 CentOS 7 的老旧 GLIBC 限制
  ▼
process merge_vcf (Mamba glnexus_env 环境)
  │  └─> 接收: tuple val(cohort_id), file(vcf_files) 
  │  └─> 动作: bcftools concat -Oz 物理安全升序排序合并，产出全基因组压缩大 VCF
  ▼
process vcf_liftover (双生态过渡桥梁)
  │  └─> 动作: 批量重命名染色体名 (chr20 -> 20)，对齐到数字命名空间
  ▼
process annotate_vcf (Singularity 容器运行 VEP + snpEff)
     └─> 动作: 动态映射 dirname，利用 100% 齐全的 GRCh37背景库进行高级临床全库注释

需要生成hg19的分染色体bed文件

由于上游软件多线程并行或随机输出等原因，这三个文件中 ##contig 染色体的出现顺序彻底被打乱了（例如 192 文件的顺序是 chr6→chr1，而 196 文件的顺序是 chr12→chr13）[3]。GLnexus 在构建数据库时要求所有输入文件的染色体顺序必须绝对一致，否则就会判定为“基因组不匹配”。

帮我执行 bash run_38.sh
1. 理解这个流程正在做的事情，确保其整个逻辑，参数，数据库和软件使用的搭配是符合生信行业流程的习惯和要求
2. 监控其任务执行状态，如果报错就调试代码。把最初的代码和数据版本先备份   
3. 查杀任务的时候确保是本任务执行的路径和命令才杀掉
4 输出正确的结果给我