sample=$1
python complexsv.py call \
    --sv_fn=/home/chenlingxi/mnt/chenlingxi/workspace/Bio_Projects/scDNA/SVCalling/svaba/$sample/${sample}.svaba.sv.vcf \
    --sample_name=$sample \
    --out_dir=$sample
