# vcfFilter - Python tool for handling VCF files.

## Table of Contents

* [vcfFilter \- Python tool for handling VCF files\.](#vcffilter---python-tool-for-handling-vcf-files)
    * [1 Filter by Chromosome](#1-filter-by-chromosome)
      * [1\.1 Keep all autosomes](#11-keep-all-autosomes)
      * [1\.2 Keep chromosome 1 and 3](#12-keep-chromosome-1-and-3)
      * [1\.3 Keep chromosome 1, 2, 3 and 6](#13-keep-chromosome-1-2-3-and-6)
    * [2 Filter by Position](#2-filter-by-position)
    * [3 Filter by QUAL Score](#3-filter-by-qual-score)


## 1 Filter by Chromosome

### 1.1 Keep all autosomes

```python
python vcfFilter.py -vcf input.vcf -chr chr1-chr22 -o output.vcf
```

### 1.2 Keep chromosome 1 and 3

```python
python vcfFilter.py -vcf input.vcf -chr chr1,chr3 -o output.vcf
```

### 1.3 Keep chromosome 1, 2, 3 and 6 

```python
python vcfFilter.py -vcf input.vcf -chr chr1-chr3,chr6 -o output.vcf
```

## 2 Filter by Position

Keep only the first 1 Mb (1-1,000,000) region on chromosome 1:  

```python
python vcfFilter.py -vcf input.vcf -pos chr1:1-1000000 -o output.vcf
```

## 3 Filter by QUAL Score

Keep variants with phred-scaled quality score no less than 30:  

```python
python vcfFilter.py -vcf input.vcf -qual 30 -o output.vcf
```
