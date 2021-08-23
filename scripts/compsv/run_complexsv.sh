sample=$1
wdr=$2
python complexsv.py call \
    --sv_fn=$wdr/SVCalling/svaba/$sample/${sample}.svaba.sv.vcf \
    --sample=$sample \
    --out_dir=$wdr/SVCalling/svaba/$sample/complexsv
