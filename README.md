# bam2iupac
IUPAC FASTA extraction from BAM files to stdout

## install
```
#Ubuntu: sudo apt-get install libbz2-dev
#RedHat: sudo dnf install bzip2-devel
git clone --recursive https://github.com/kullrich/bam2iupac
cd bam2iupac/htslib
make
cd ..
make HTSLIB_DIR=$PWD/htslib
```

## run
```
bam2iupac --help

SYNOPSIS
  IUPAC FASTA extraction from BAM files to stdout
USAGE
  ./bam2iupac [options] --b 1.bam --n ind1 --b 2.bam --n ind2 [...]
  OR
  ./bam2iupac [options] --bamList samples.txt
OPTIONS
  --b          BAM files
  --n          Sequence IDs
  --bamList    File containing BAM paths and Sample IDs (one per line: <BAM> <SAMPLE>)
  --r          Region ('chr:start-end' or 'chr start end') coordinates are 1-based
  --regionList File containing regions (one per line: 'chr:start-end' or 'chr start end')
  --minMQ      Minimum mapping quality (default: 0)
  --minBQ      Minimum base quality (default: 0)
  --minC       Minimum coverage (default: 0)
  --maxC       Maximum coverage (default: 9999)
  --iupacRatio IUPAC ratio (default: 0.25)
  --incMQ      Include missing mapping quality value 255 (default: False)
  --incBQ      Include missing base quality value 255 (default: False)
  --help       Show this help
  --version    Print version and exit
  --debug      Debug

         NOTE: Multiple BAM files (--b/--bamList) and regions (--r/--regionList)
               are processed in the order provided.
               Multiple regions are concatenated into a single output sequence,
               enabling direct extraction of combined intervals (exons) from a GTF annotation.

EXENAME
  bam2iupac
VERSION
  0.0.3
URL
  https://github.com/kullrich/bam2iupac
```

## example to extract IUPAC FASTA

### hg19
```
./bam2iupac \
--b https://cdna.eva.mpg.de/denisova/BAM/human/DNK02.bam --n DNK02 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00521.bam --n HGDP00521 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00542.bam --n HGDP00542 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00665.bam --n HGDP00665 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00778.bam --n HGDP00778 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00927.bam --n HGDP00927 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP01029.bam --n HGDP01029 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP01284.bam --n HGDP01284 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP01307.bam --n HGDP01307 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP0456.bam --n HGDP0456 \
--b https://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam --n T \
--b https://cdna.eva.mpg.de/denisova/Den25/BAM/Denisova25.hg19.ontarget.uniq.L35MQ25.indel_realigned.MDfixed.bam --n Denisova25 \
--r 1:10000001-10001000
```

### concatenate regions
```
./bam2iupac \
--b https://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam --n T \
--r 1:10000001-10001000 \
--r 1:10001001-10002000
```

### direct distance calculation with [literal-dists](https://github.com/kullrich/literal-dists)
```
git clone https://github.com/kullrich/literal-dists
cd literal-dists
make
```

```
./bam2iupac \
--b https://cdna.eva.mpg.de/denisova/BAM/human/DNK02.bam --n DNK02 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00521.bam --n HGDP00521 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00542.bam --n HGDP00542 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00665.bam --n HGDP00665 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00778.bam --n HGDP00778 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP00927.bam --n HGDP00927 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP01029.bam --n HGDP01029 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP01284.bam --n HGDP01284 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP01307.bam --n HGDP01307 \
--b https://cdna.eva.mpg.de/denisova/BAM/human/HGDP0456.bam --n HGDP0456 \
--b https://cdna.eva.mpg.de/denisova/alignments/T_hg19_1000g.bam --n T \
--b https://cdna.eva.mpg.de/denisova/Den25/BAM/Denisova25.hg19.ontarget.uniq.L35MQ25.indel_realigned.MDfixed.bam --n Denisova25 \
--r 1:10000001-10001000 | ./literal-dist > distances.tsv
```

