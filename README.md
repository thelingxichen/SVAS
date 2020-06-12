# SVAS
SVAS( Somatic Variant Analysis Suite)
Author: Lingxi Chen

Email: chanlingxi@gmail.com

## Description
ComplexSV provides scripts for complex structure variation analsyis. Currently, it provides two functionalities:
+ [linkage_heatmap.py] Generate read shared linkages for given SV event.
+ [read_support.py] Generate supporting split reads for given SV event.

The output files can be visuliazed through https://btdraw.com/ SV modules.

## Depencency

User can install the required packages via `requirements.txt`.
```
python3
pyvcf==0.6.8
pysam==0.15.3
numpy==1.15.4
docopt==0.6.2
```

## Installation
```
git clone https://github.com/paprikachan/ComplexSV.git
```
Then, please add this directory to your `PATH`.
```
export PATH=$PWD/ComplexSV/:$PATH
```

## Linkage Heatmap

### General usage
```
$python3 linkage_heatmap.py -h
Linkage Heatmap.

Usage:
    linkage_heatmap.py call --bam_fn=IN_FILE --sv_fn=IN_FILE [--out_dir=IN_DIR] [--tenx=STR]
    linkage_heatmap.py -h | --help

Options:
    -h --help           Show this screen.
    --version           Show version.
    --bam_fn=IN_FILE    Path of tumor bam file.
    --sv_fn=IN_FILE     Path of tumor sv vcf file.
    --out_dir=OUT_DIR   Path of out directory, [default: ./]
    --tenx=BOOL         Boolean value, True if the protocal is 10x, otherwise False.
```
User can check the usability with demo case.
```
cd ComplexSV
python3 linkage_heatmap.py call --bam_fn=demo_data/tumor.bam --sv_fn=demo_data/tumor.sv.vcf --out_dir=demo_data/ --tenx=True
```

### Understanding Output

The demo output file `demo_data/10x.txt` starts with "sv" section.
```
#sv
chr11	76139879	+	chr2	63328429	-   VARTYPE=BND:TRX-tt
```
It stores the breakpoints information of SV event, ordered by:

| chrom_5p |  bkpos_5p |  strand_5p | chrom_3p |  bkpos_3p |  strand_3p | meta_info |
|---|---|---|---|---|---|---|
| chr11  | 76139879  | + | chr2  | 63328429  |- | VARTYPE=BND:TRX-tt |

where `chrom_5p`, `bkpos_5p`, `strand_5p` respectively stands for the chromosome, position, strand of 5' breakpoint, and `chrom_3p`, `bkpos_3p`, `strand_3p` respectively stands for the chromosome, position, strand of 3' breakpoint. `meta_info` stores the meta information, such as the variation type of the SV.

The second section is "heatmap" section, it starts with the comment:
```
#heatmap    linkage_type=10x barcode
```
where `linkage_type` specifies the type of read linkage used to count, now we support to show "10x barcode" and "Pair End" linkage type.

For a SV event with two breakpoints, expand the breakpoint with `1000bp`, we can generate four region pairs and their shared linkage count matrix.

Region pair 1: vertical region is `(bkpos_5p-1000bp, bkpos_5p+1000bp)` and horizontal region is `(bkpos_5p-1000bp, bkpos_5p+1000bp)`.
```
v=chr2:63327429-63329429	h=chr2:63327429-63329429	resolution=100	axis=0,0
63.0,34.0,35.0,32.0,21.0,...
34.0,82.0,56.0,37.0,36.0,...
35.0,56.0,104.0,47.0,40.0,...
32.0,37.0,47.0,84.0,40.0,...
21.0,36.0,40.0,40.0,83.0,...
...
```
Region pair 2: vertical region is `(bkpos_5p-1000bp, bkpos_5p+1000bp)` and horizontal region is `(bkpos_3p-1000bp, bkpos_3p+1000bp)`.
```
v=chr2:63327429-63329429	h=chr11:76138879-76140879	resolution=100	axis=1,0
8.0,14.0,7.0,8.0,11.0,...
9.0,17.0,13.0,14.0,13.0,...
14.0,16.0,17.0,22.0,15.0,...
9.0,14.0,11.0,13.0,13.0,...
8.0,13.0,13.0,16.0,13.0,...
...
```
Region pair 3: vertical region is `(bkpos_3p-1000bp, bkpos_3p+1000bp)` and horizontal region is `(bkpos_5p-1000bp, bkpos_5p+1000bp)`.
```
v=chr11:76138879-76140879	h=chr2:63327429-63329429	resolution=100	axis=0,1
8.0,9.0,14.0,9.0,8.0,...
14.0,17.0,16.0,14.0,13.0,...
7.0,13.0,17.0,11.0,13.0,...
8.0,14.0,22.0,13.0,16.0,...
11.0,13.0,15.0,13.0,13.0,...
...
```
Region pair 4: vertical region is `(bkpos_3p-1000bp, bkpos_3p+1000bp)` and horizontal region is `(bkpos_3p-1000bp, bkpos_3p+1000bp)`.
```
v=chr11:76138879-76140879	h=chr11:76138879-76140879	resolution=100	axis=1,1
80.0,52.0,27.0,42.0,33.0,...
52.0,140.0,56.0,61.0,50.0,...
27.0,56.0,124.0,72.0,50.0,...
42.0,61.0,72.0,152.0,57.0,...
33.0,50.0,50.0,57.0,112.0,...
...
```


## Read Support
### General usage
```
$python3 read_support.py -h
Read Support.

Usage:
    read_support.py call --bam_fn=IN_FILE --sv_fn=IN_FILE [--out_dir=OUT_DIR]
    read_support.py -h | --help

Options:
    -h --help           Show this screen.
    --version           Show version.
    --bam_fn=IN_FILE    Path of tumor bam file.
    --sv_fn=IN_FILE     Path of tumor sv vcf file.
    --out_dir=OUT_DIR   Path of out directory, [default: ./]
```
User can check the usability with demo case.
```
cd ComplexSV
python3 read_support.py call --bam_fn=demo_data/tumor.bam --sv_fn=demo_data/tumor.sv.vcf --out_dir=demo_data/
```

### Understanding Output

The demo output file `demo_data/simu.junc.reads.txt` stores supporting split reads for given SV. It starts with "sv" section.

| chrom_5p |  bkpos_5p | chrom_3p |  bkpos_3p |  inner_ins | splits_read_num |
|---|---|---|---|---|---|
| chr11  | 100000  | chr1 | 101000  | Inner-Ins:TACCGATAT  |10 | 

where `inner_ins` refers to the detected inner insertion, and `splits_read_num` refers to the number of supporting split read pairs.

Then, it list all supporting split reads, ordered by:
+ `read_query_name`, the read id;
+ `read_flags`, it consists of flag `J`, `P`, `r`, `R`, `d`, where

  |Flag|Description|
  |---|---|
  |`J`|The current read is split read|
  |`P`|The paired read of current read is split read|
  |`r`|The current read is reverse read|
  |`R`|The paired read of current read is reverse read|
  |`d`|The current read is duplicated|

+ `read_position`, the read start position aligned to reconstructed SV haplotype.
  ```
   |<-   5'segment   ->| inner_ins |<- 3'  segment  ->|
   |xxxxxxxxxxxxxxxxxxx|iiiiiiiiiii|xxxxxxxxxxxxxxxxxx|
                  :    ^
         negative : <- 0 -> positive
                  :
  split_read:     xxxxxxxxxxxxxxxxxxxxxxxxxxxx
                  ^
  read_position: -5
  ```
+ `read_seq`, the read sequence;
+ `read_qual`, the read sequence quality.


