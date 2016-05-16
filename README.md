# vcfFilter - Python tool for handling VCF files.

Table of Contents
=================

  * [1 Introduction](#1-introduction)
  * [2 Requirement and Installation](#2-requirement-and-installation)
  * [3 Usage](#3-usage)
    * [3\.1 Filter by Chromosome](#31-filter-by-chromosome)
      * [3\.1\.1 Keep all autosomes (separated by "\-")](#311-keep-all-autosomes-separated-by--)
      * [3\.1\.2 Keep chromosomes (separated by ",")](#312-keep-chromosomes-separated-by-)
      * [3\.1\.3 Keep chromosomes (separated by "," and "\-")](#313-keep-chromosomes-separated-by--and--)
    * [3\.2 Filter by Region](#32-filter-by-region)
    * [3\.3 Filter by QUAL Score](#33-filter-by-qual-score)
    * [3\.4 Filter by Filter Flag](#34-filter-by-filter-flag)
      * [3\.4\.1 Keep variants with FILTER flag: "PASS"](#341-keep-variants-with-filter-flag-pass)
      * [3\.4\.2 Keep variants with FILTER flags (separated by ",")](#342-keep-variants-with-filter-flags-separated-by-)
    * [3\.5 Filter by Genotype fields](#35-filter-by-genotype-fields)
    * [3\.6 Filter Indels](#36-filter-indels)
      * [3\.6\.1 keep only sites that contain an indel](#361-keep-only-sites-that-contain-an-indel)
      * [3\.6\.2 exclude sites that contain an indel](#362-exclude-sites-that-contain-an-indel)
    * [3\.7 Filter by ID](#37-filter-by-id)
      * [3\.7\.1 IDs were seperated by ','](#371-ids-were-seperated-by-)
      * [3\.7\.2 IDs were stored in a file](#372-ids-were-stored-in-a-file)
    * [3\.8 Filter by Physical Positions](#38-filter-by-physical-positions)
      * [3\.8\.1 Physical positions were seperated by ','](#381-physical-positions-were-seperated-by-)
      * [3\.8\.2 Physical Positions were stored in a file](#382-physical-positions-were-stored-in-a-file)
  * [4 References](#4-references)

# 1 Introduction

# 2 Requirement and Installation

# 3 Usage

## 3.1 Filter by Chromosome

### 3.1.1 Keep all autosomes (separated by "-")

```
python vcfFilter.py -vcf input.vcf -chr chr1-chr22 -out output.vcf
```

### 3.1.2 Keep chromosomes (separated by ",")

```
python vcfFilter.py -vcf input.vcf -chr chr1,chr3 -out output.vcf
```

### 3.1.3 Keep chromosomes (separated by "," and "-")

```
python vcfFilter.py -vcf input.vcf -chr chr1-chr3,chr6 -out output.vcf
```

## 3.2 Filter by Region

Keep only the first 1 Mb (1-1,000,000) region on chromosome 1:  

```
python vcfFilter.py -vcf input.vcf -region chr1:1-1000000 -out output.vcf
```

## 3.3 Filter by QUAL Score

Keep variants with phred-scaled quality score no less than 30:  

```
python vcfFilter.py -vcf input.vcf -qual 30 -out output.vcf
```

## 3.4 Filter by Filter Flag

### 3.4.1 Keep variants with FILTER flag: "PASS"

```
python vcfFilter.py -vcf input.vcf -filter PASS -out output.vcf
```

### 3.4.2 Keep variants with FILTER flags (separated by ",")

```
python vcfFilter.py -vcf input.vcf -filter PASS,VQSRTrancheINDEL99.00to99.90,
VQSRTrancheINDEL99.90to100.00,VQSRTrancheSNP99.00to99.90,
VQSRTrancheSNP99.90to100.00 -out output.vcf
```

## 3.5 Filter by Genotype fields

Keep homozygous of reference alleles in sample 001 and sample 002:  

```
python vcfFilter.py -vcf input.vcf -genotype hom-ref -ind 001,002 -out output.vcf
```

Values for **-genotype**:

value        | number of zero | A==B   | description 
-------------|----------------|--------|--------------------------------
hom-ref      | 2              | yes    | keep homozygous of refernce allele, two ref, e.g., 0/0
hom-alt      | 0              | yes    | keep homozygous of alternative allele, two same alt e.g., 1/1
het          | 1              | no     | keep heterozygous, one ref and one alt e.g., 0/1, or 0/2
het-alt      | 0              | no     | keep individual has two copy of different alternative alleles, two different alt, e.g., 1/2
not-hom-ref  | 0 or 1         | yes/no | keep individual does not have two copy of reference allele, one ref or no ref, e.g., 0/1 or 1/1 or 1/2
not-two-alt  | 1 or 2         | yes/no | keep individual does not have two copy of alternative allele, at least one ref, e.g., 0/0, 0/1, 0/2
two-alt      | 0              | yes/no | keep individual has two copy of alternative allele, two same or different alt, e.g., 1/1 or 1/2
not-het      | 0 or 2         | yes/no | keep individual who is not heterozygous, hom-ref or two-alt, e.g. 0/0, 1/1, 1/2

The allele values are 0 for the reference allele, 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1 (A=0 and B=1, number of zero is 1).   

value of **-ind** is a string of individual IDs separated by ",". e.g., "-ind 1", "-ind 1,2,3"

**Note**: individual IDs must be in the header line of the input VCF file.

## 3.6 Filter Indels

If one site only has one ALT allele, "indel" here means that the length of REF and ALT allele is different. If one site has multiple ALT alleles, "indel" here means that at least one ALT allele has different length from REF allele. Site 1 and 2 will be defined as "indel" in the forllowing four sites.

```
site   REF                ALT
1      TTTTA              T,TTTTATTTA
2      A                  ATGTG,G,ATG
3      G                  C,T
4      T                  G
```

### 3.6.1 keep only sites that contain an indel

```
python vcfFilter.py -vcf input.vcf --keep-only-indels -out output.vcf
```

### 3.6.2 exclude sites that contain an indel

```
python vcfFilter.py -vcf input.vcf --remove-indels -out output.vcf
```

## 3.7 Filter by ID

### 3.7.1 IDs were seperated by ','

Multiple IDs (e.g. dbSNP rsID) can be seperated using ",". 

```
python vcfFilter.py -vcf input.vcf -ids rs1234 -out output.vcf
```

```
python vcfFilter.py -vcf input.vcf -ids rs1234,rs1235 -out output.vcf
```

### 3.7.2 IDs were stored in a file

Each row contains one ID. e.g.:

```
rs1234
rs1235
```

```
python vcfFilter.py -vcf input.vcf --ids-file ids.txt -out output.vcf
```

## 3.8 Filter by Physical Positions

### 3.8.1 Physical positions were seperated by ','

Physical positions can be seperated using ",". Each physical position includes chromosome and position that was seperated by ":".

```
python vcfFilter.py -vcf input.vcf --phy-pos chr1:1234567, -out output.vcf
```

```
python vcfFilter.py -vcf input.vcf --phy-pos chr1:1234567,chr2:9887234 -out output.vcf
```

### 3.8.2 Physical Positions were stored in a file

Each row contains two columns, the first column is chromosome and the second column is position. e.g.:

```
chr1  1234567
chr2  9887234
```

```
python vcfFilter.py -vcf input.vcf --phy-pos-file phypos.txt -out output.vcf
```

# 4 References

* [VCF (Variant Call Format) version 4.0](http://www.1000genomes.org/wiki/Analysis/vcf4.0)
* [VCF (Variant Call Format) version 4.1](http://samtools.github.io/hts-specs/VCFv4.1.pdf)
* [VCF (Variant Call Format) version 4.2](http://samtools.github.io/hts-specs/VCFv4.2.pdf)
* [VCF (Variant Call Format) version 4.3](http://samtools.github.io/hts-specs/VCFv4.3.pdf)

