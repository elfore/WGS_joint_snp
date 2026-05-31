pipedir="/mnt/gpfs1/Users/yangjinxurong/projects/WGS/joint_snp"
source /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/bin/activate /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/envs/nextflow/
nextflow -q -log ${pipedir}/log/log_38 run ${pipedir}/main_38.nf -params-file ${pipedir}/params_38.yaml -c ${pipedir}/nextflow.config -bg -with-trace -resume

# pkill -u yangjinxurong -f nextflow