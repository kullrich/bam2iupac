# bam2iupac
IUPAC FASTA extraction from BAM files to stdout

## install
```
git clone --recursive https://github.com/kullrich/bam2iupac
cd bam2iupac/htslib
make
cd ..
make HTSLIB_DIR=htslib
```

## run
```
bam2iupac --help

SYNOPSIS
  IUPAC FASTA extraction from BAM files to stdout
USAGE
  ./bam2iupac [options] --b 1.bam --n ind1 --b 2.bam --n ind2 [...]
OPTIONS
  --b		BAM files
  --n		Sequence IDs
  --r		Region ('chr:start-end' or 'chr start end')
  --minMQ	Minimum mapping quality (default: 0)
  --minBQ	Minimum base quality (default: 0)
  --minC	Minimum coverage (default: 0)
  --maxC	Maximum coverage (default: 9999)
  --iupacRatio	IUPAC ratio (default: 0.25)
  --incMQ	Include missing mapping quality value 255 (default: False)
  --incBQ	Include missing base quality value 255 (default: False)
  --help	Show this help
  --version	Print version and exit
  --debug	Debug
EXENAME
  bam2iupac
VERSION
  0.0.1
URL
  https://github.com/kullrich/bam2iupac
```

## example to extract IUPAC-fasta
```

```

