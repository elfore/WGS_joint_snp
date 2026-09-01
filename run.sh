pipeDir="/mnt/gpfs1/Users/yangjinxurong/projects/WGS/snp"
source /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/bin/activate /mnt/gpfs1/Users/yangjinxurong/software/miniconda3/envs/nextflow/
nextflow -q \
    -log log/log \
    run ${pipeDir}/main.nf \
    -params-file params.yaml \
    -c ${pipeDir}/nextflow.config \
    -bg \
    -with-trace \
    -resume