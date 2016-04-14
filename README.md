# vcfFilter - Python tool for handling VCF files.

Table of Contents
=================

  * [vcfFilter \- Python tool for handling VCF files\.](#vcffilter---python-tool-for-handling-vcf-files)
  * [1 Introduction](#1-introduction)
  * [2 Requirement and Installation](#2-requirement-and-installation)
  * [3 Usage](#3-usage)
    * [3\.1 Filter by Chromosome](#31-filter-by-chromosome)
      * [3\.1\.1 Keep all autosomes (separated by "\-")](#311-keep-all-autosomes-separated-by--)
      * [3\.1\.2 Keep chromosomes (separated by ",")](#312-keep-chromosomes-separated-by-)
      * [3\.1\.3 Keep chromosomes (separated by "," and "\-")](#313-keep-chromosomes-separated-by--and--)
    * [3\.2 Filter by Position](#32-filter-by-position)
    * [3\.3 Filter by QUAL Score](#33-filter-by-qual-score)
    * [3\.4 Filter by Filter Flag](#34-filter-by-filter-flag)
      * [3\.4\.1 Keep variants with FILTER flag: "PASS"](#341-keep-variants-with-filter-flag-pass)
      * [3\.4\.2 Keep variants with FILTER flags (separated by ",")](#342-keep-variants-with-filter-flags-separated-by-)
    * [3\.5 Filter by Genotype fields](#35-filter-by-genotype-fields)
  * [4 References](#4-references)

# 1 Introduction

# 2 Requirement and Installation

# 3 Usage

## 3.1 Filter by Chromosome

### 3.1.1 Keep all autosomes (separated by "-")

```python
python vcfFilter.py -vcf input.vcf -chr chr1-chr22 -o output.vcf
```

### 3.1.2 Keep chromosomes (separated by ",")

```python
python vcfFilter.py -vcf input.vcf -chr chr1,chr3 -o output.vcf
```

### 3.1.3 Keep chromosomes (separated by "," and "-")

```python
python vcfFilter.py -vcf input.vcf -chr chr1-chr3,chr6 -o output.vcf
```

## 3.2 Filter by Position

Keep only the first 1 Mb (1-1,000,000) region on chromosome 1:  

```python
python vcfFilter.py -vcf input.vcf -pos chr1:1-1000000 -o output.vcf
```

## 3.3 Filter by QUAL Score

Keep variants with phred-scaled quality score no less than 30:  

```python
python vcfFilter.py -vcf input.vcf -qual 30 -o output.vcf
```

## 3.4 Filter by Filter Flag

### 3.4.1 Keep variants with FILTER flag: "PASS"

```python
python vcfFilter.py -vcf input.vcf -filter PASS -o output.vcf
```

### 3.4.2 Keep variants with FILTER flags (separated by ",")

```python
python vcfFilter.py -vcf input.vcf -filter PASS,VQSRTrancheINDEL99.00to99.90,VQSRTrancheINDEL99.90to100.00,VQSRTrancheSNP99.00to99.90,VQSRTrancheSNP99.90to100.00 -o output.vcf
```

## 3.5 Filter by Genotype fields

Keep homozygous of reference alleles in sample 001 and sample 002:  

```python
python vcfFilter.py -vcf input.vcf -gtp hom-ref -ind 001,002 -o output.vcf
```

Values for **-gtp**:

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

# 4 References

* [VCF (Variant Call Format) version 4.0](http://www.1000genomes.org/wiki/Analysis/vcf4.0)
* [VCF (Variant Call Format) version 4.1](http://samtools.github.io/hts-specs/VCFv4.1.pdf)
* [VCF (Variant Call Format) version 4.2](http://samtools.github.io/hts-specs/VCFv4.2.pdf)
* [VCF (Variant Call Format) version 4.3](http://samtools.github.io/hts-specs/VCFv4.3.pdf)

