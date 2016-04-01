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
