# PyTV - Python Tools for handling VCF files.

Table of Contents
=================

  * [1 Introduction](#1-introduction)
  * [2 Requirement and Installation](#2-requirement-and-installation)
  * [3 vcfFilter](#3-vcffilter)
    * [3\.1 Variant Position Filtering](#31-variant-position-filtering)
      * [3\.1\.1 Chromosome filtering](#311-chromosome-filtering)
        * [3\.1\.1\.1 Keep all autosomes (separated by "\-")](#3111-keep-all-autosomes-separated-by--)
        * [3\.1\.1\.2 Keep chromosomes (separated by ",")](#3112-keep-chromosomes-separated-by-)
        * [3\.1\.1\.3 Keep chromosomes (separated by "," and "\-")](#3113-keep-chromosomes-separated-by--and--)
        * [3\.1\.1\.4 Exclude chromosomes](#3114-exclude-chromosomes)
      * [3\.1\.2 Region filtering](#312-region-filtering)
        * [3\.1\.2\.1 Single region](#3121-single-region)
        * [3\.1\.2\.2 Multiple regions](#3122-multiple-regions)
      * [3\.1\.3 Physical position filtering](#313-physical-position-filtering)
        * [3\.1\.3\.1 Physical positions were seperated by ','](#3131-physical-positions-were-seperated-by-)
        * [3\.1\.3\.2 Physical positions were stored in a file](#3132-physical-positions-were-stored-in-a-file)
    * [3\.2 Variant ID Filtering](#32-variant-id-filtering)
      * [3\.2\.1 IDs were seperated by ','](#321-ids-were-seperated-by-)
      * [3\.2\.2 IDs were stored in a file](#322-ids-were-stored-in-a-file)
      * [3\.2\.3 Remove variants by IDs](#323-remove-variants-by-ids)
    * [3\.3 Variant Type Filtering](#33-variant-type-filtering)
      * [3\.3\.1 Keep only sites that contain an indel](#331-keep-only-sites-that-contain-an-indel)
      * [3\.3\.2 Exclude sites that contain an indel](#332-exclude-sites-that-contain-an-indel)
    * [3\.4 FILTER Flag Filtering](#34-filter-flag-filtering)
      * [3\.4\.1 Keep variants with FILTER flag: "PASS"](#341-keep-variants-with-filter-flag-pass)
      * [3\.4\.2 Keep variants with FILTER flags (separated by ",")](#342-keep-variants-with-filter-flags-separated-by-)
    * [3\.5 INFO Field Filtering](#35-info-field-filtering)
      * [3\.5\.1 Option usage](#351-option-usage)
      * [3\.5\.2 Examples of using \-\-info option](#352-examples-of-using---info-option)
        * [3\.5\.2\.1 Remove variants located in the intergenic, downstream or upstream region](#3521-remove-variants-located-in-the-intergenic-downstream-or-upstream-region)
        * [3\.5\.2\.2 Remove variants wilth alternatice allele frequency in 1000 genomes higher than 0\.05](#3522-remove-variants-wilth-alternatice-allele-frequency-in-1000-genomes-higher-than-005)
        * [3\.5\.2\.3 Remove variants wilth alternatice allele frequency higher than 0\.05](#3523-remove-variants-wilth-alternatice-allele-frequency-higher-than-005)
        * [3\.5\.2\.4 Remove variants if the combined depth across samples is less than 200](#3524-remove-variants-if-the-combined-depth-across-samples-is-less-than-200)
        * [3\.5\.2\.5 Remove variants if the RMS mapping quality is less than 50](#3525-remove-variants-if-the-rms-mapping-quality-is-less-than-50)
    * [3\.6 QUAL Field Filtering](#36-qual-field-filtering)
    * [3\.7 Allele Filtering](#37-allele-filtering)
      * [3\.7\.1 Minimum / maximum number of alleles](#371-minimum--maximum-number-of-alleles)
      * [3\.7\.2 Biallelic / multiallelic site](#372-biallelic--multiallelic-site)
    * [3\.8 Individual Filtering](#38-individual-filtering)
    * [3\.9 Genotype Filtering](#39-genotype-filtering)
      * [3\.9\.1 Missing rate and missing count filtering](#391-missing-rate-and-missing-count-filtering)
      * [3\.9\.2 Homozygous and heterozygous filtering](#392-homozygous-and-heterozygous-filtering)
      * [3\.9\.3 Compound heterozygous filtering](#393-compound-heterozygous-filtering)
      * [3\.9\.4 Compare genotype of multiple individuals](#394-compare-genotype-of-multiple-individuals)
        * [3\.9\.4\.1 Only variants have the same genotype across specified individuals will be kept](#3941-only-variants-have-the-same-genotype-across-specified-individuals-will-be-kept)
        * [3\.9\.4\.2 Only variants have different genotype between the first individual and others will be kept](#3942-only-variants-have-different-genotype-between-the-first-individual-and-others-will-be-kept)
  * [4 References](#4-references)

# 1 Introduction

Python tools for handling VCF files.

* vcfFilter.py: Filtering of variants according to genotypes and thresholds for the 8 fixed fields in VCF file.
					
# 2 Requirement and Installation

`PyTV` uses Python 2 (Python 2.7 or higher) which is available [here](https://www.python.org/).

# 3 vcfFilter

Filtering of variants according to genotypes and thresholds for the 8 fixed fields in VCF file.

## 3.1 Variant Position Filtering

### 3.1.1 Chromosome filtering

#### 3.1.1.1 Keep all autosomes (separated by "-")

```
python vcfFilter.py --vcf input.vcf --chr chr1-chr22 --out output.vcf
```

#### 3.1.1.2 Keep chromosomes (separated by ",")

```
python vcfFilter.py --vcf input.vcf --chr chr1,chr3 --out output.vcf
```

#### 3.1.1.3 Keep chromosomes (separated by "," and "-")

```
python vcfFilter.py --vcf input.vcf --chr chr1-chr3,chr6 --out output.vcf
```

#### 3.1.1.4 Exclude chromosomes

use `--reverse` to reverse the filter. e.g., exclude chr3 and chr6

```
python vcfFilter.py --vcf input.vcf --chr chr1-chr3,chr6 --reverse --out output.vcf
```

### 3.1.2 Region filtering

#### 3.1.2.1 Single region

Keep only the first 1 Mb (1-1,000,000) region on chromosome 1:  

```
python vcfFilter.py --vcf input.vcf --region chr1:1-1000000 --out output.vcf
```

Exclude the first 1 Mb (1-1,000,000) region on chromosome 1:

```
python vcfFilter.py --vcf input.vcf --region chr1:1-1000000 --reverse -out output.vcf
```

#### 3.1.2.2 Multiple regions

```
python vcfFilter.py --vcf input.vcf --region-file regions.txt --out output.vcf
```

```
python vcfFilter.py --vcf input.vcf --region-file regions.txt --out output.vcf --reverse
```

`regions.txt` contains one or multiple rows. Each row defines a region and contains three columns, the first column is chromosome, the second column is the start of the region and the third column is the end of the region. e.g.:

```
chr1  1 1000000
chr6  33000000 40000000 
```

### 3.1.3 Physical position filtering

Includes site based on its position. `--reverse` can be used to reverse the filter.

#### 3.1.3.1 Physical positions were seperated by ','

Physical positions can be seperated using ",". Each physical position includes chromosome and position that was seperated by ":". 

```
python vcfFilter.py --vcf input.vcf --phy-pos chr1:1234567, --out output.vcf
```

```
python vcfFilter.py --vcf input.vcf --phy-pos chr1:1234567,chr2:9887234 --out output.vcf --reverse
```

#### 3.1.3.2 Physical positions were stored in a file

Each row contains two columns, the first column is chromosome and the second column is position. e.g.:

```
chr1  1234567
chr2  9887234
```

```
python vcfFilter.py --vcf input.vcf --phy-pos-file phypos.txt --out output.vcf
```

## 3.2 Variant ID Filtering

### 3.2.1 IDs were seperated by ','

Multiple IDs (e.g. dbSNP rsID) can be seperated using ",". 

```
python vcfFilter.py --vcf input.vcf --ids rs1234 --out output.vcf
```

```
python vcfFilter.py --vcf input.vcf --ids rs1234,rs1235 --out output.vcf
```

### 3.2.2 IDs were stored in a file

Each row contains one ID. e.g.:

```
rs1234
rs1235
```

```
python vcfFilter.py --vcf input.vcf --ids-file ids.txt --out output.vcf
```

### 3.2.3 Remove variants by IDs

use `--reverse` to reverse the filter.

```
python vcfFilter.py --vcf input.vcf --ids rs1234,rs1235 --reverse --out output.vcf
```

## 3.3 Variant Type Filtering

If one site only has one ALT allele, "indel" here means that the length of REF and ALT allele is different. If one site has multiple ALT alleles, "indel" here means that at least one ALT allele has different length from REF allele. Site 1 and 2 will be defined as "indel" in the forllowing four sites.

```
site   REF                ALT
1      TTTTA              T,TTTTATTTA
2      A                  ATGTG,G,ATG
3      G                  C,T
4      T                  G
```

### 3.3.1 Keep only sites that contain an indel

```
python vcfFilter.py --vcf input.vcf --keep-only-indels --out output.vcf
```

### 3.3.2 Exclude sites that contain an indel

```
python vcfFilter.py --vcf input.vcf --remove-indels --out output.vcf
```

## 3.4 FILTER Flag Filtering

### 3.4.1 Keep variants with FILTER flag: "PASS"

```
python vcfFilter.py --vcf input.vcf --filter PASS --out output.vcf
```

### 3.4.2 Keep variants with FILTER flags (separated by ",")

```
python vcfFilter.py --vcf input.vcf --filter PASS,VQSRTrancheINDEL99.00to99.90,VQSRTrancheINDEL99.90to100.00,VQSRTrancheSNP99.00to99.90,VQSRTrancheSNP99.90to100.00 --out output.vcf
```

## 3.5 INFO Field Filtering

This option filter on the presence of the flag and its value. INFO fields in VCF are encoded as a semicolon-separated series of short keys with optional values in the format:

```
<key>=<value>[,value]
```

When [annovar](http://annovar.openbioinformatics.org/en/latest/) was used to annotate the variants, multiple additional key and value pairs will be added in the INFO feild. Suppose `.` was used to reprent the missing value.   

### 3.5.1 Option usage

commands format for numeric value (only one value):

```
--info '<key>>=<value>'
--info '<key><=<value>'
--info '<key>><value>'
--info '<key><<value>'
```
commands format for numeric value or string type value (one value or multiple values seperated by ','):

```
--info '<key>=<value>[,value]'
--info '<key>!=<value>[,value]'
```

- Whick `key` can be used depends on your own vcf. 
- There are six operations `>=`, `>`, `<=`, `<`, `=`, `!=`. 
- The first four operations only accept **ONE float** value. Int value will be convert to float, but string or multiple value will cause ERROR. 
- The last two operation can accept **one or multiple float or string values**. If you want to filter by an int such as 'AC=10', you need to use `-info 'AC=10,'`.
- if there are multiple values for the key in vcf, all values must pass the cutoff.
- There is no whitespace between `key`, `operation`, and `value`.

### 3.5.2 Examples of using `--info` option

The following examples show how it works. Please be noted that annovar was used to annotate the vcf, so there are some keys may not in your vcf.

#### 3.5.2.1 Remove variants located in the intergenic, downstream or upstream region

```
python vcfFilter.py --vcf input.vcf --out output.vcf --info 'Func_refGene!=intergenic,downstream,upstream' 
```

By default, missing value will be kept, chang this by adding `--missing-value rm`

```
python vcfFilter.py --vcf input.vcf --out output.vcf --info 'Func_refGene!=intergenic,downstream,upstream' --missing-value rm
```

#### 3.5.2.2 Remove variants wilth alternatice allele frequency in 1000 genomes higher than 0.05

```
python vcfFilter.py --vcf input.vcf --info '1000g2015aug_all<=0.05' --missing-value keep --out output.vcf
```

#### 3.5.2.3 Remove variants wilth alternatice allele frequency higher than 0.05

```
python vcfFilter.py --vcf input.vcf --info 'AF<=0.05' --out output.vcf
```

Note: Currently, this only works for variants with only one alternative allele.   

#### 3.5.2.4 Remove variants if the combined depth across samples is less than 200 

```
python vcfFilter.py --vcf input.vcf --info 'DP>=200' --out output.vcf
```

#### 3.5.2.5 Remove variants if the RMS mapping quality is less than 50

```
python vcfFilter.py --vcf input.vcf --info 'MQ>=50' --out output.vcf
```

## 3.6 QUAL Field Filtering

Keep variants with phred-scaled quality score no less than 30:  

```
python vcfFilter.py --vcf input.vcf --qual 30 --out output.vcf
```

## 3.7 Allele Filtering

### 3.7.1 Minimum / maximum number of alleles

Only sites have alleles no less than (`>=`) the minimum number will be kept.

```
python vcfFilter.py --vcf input.vcf --min-alleles 2 --out output.vcf
```

Only sites have alleles no more than (`<=`) the maximum number will be kept.

```
python vcfFilter.py --vcf input.vcf --max-alleles 2 --out output.vcf
```

### 3.7.2 Biallelic / multiallelic site

A biallelic site is a specific locus in a genome that contains two observed alleles, counting the reference as one, and therefore allowing for one variant allele.

```
python vcfFilter.py --vcf input.vcf --min-alleles 2 --max-alleles 2 --out output.vcf
```

A multiallelic site is a specific locus in a genome that contains three or more observed alleles, again counting the reference as one, and therefore allowing for two or more variant alleles.

```
python vcfFilter.py --vcf input.vcf --min-alleles 3 --out output.vcf
```

## 3.8 Individual Filtering

Include or exclude certain individuals.

```
python vcfFilter.py --vcf input.vcf --keep-inds --ind 001,002 --out output.vcf
```

```
python vcfFilter.py --vcf input.vcf --remove-inds --ind 001,002 --out output.vcf
```

## 3.9 Genotype Filtering

### 3.9.1 Missing rate and missing count filtering

Exclude sites with proportion of missing data is larger than the cutoff, which is between 0 and 1, where 0 indicates no missing data allowed and 1 allows sites that are completely missing.  

```
python vcfFilter.py --vcf input.vcf --missing-rate 0.05 --out output.vcf
```

Exclude sites that the number of individuals who missing this site is larger than the cutoff. A low missing count cutoff is similar with a low missing rate cutoff (`missing_rate = missing_count / total_individuals`). 

```
python vcfFilter.py --vcf input.vcf --missing-count 10
```

### 3.9.2 Homozygous and heterozygous filtering

Keep homozygous of reference alleles in sample 001 and sample 002:  

```
python vcfFilter.py --vcf input.vcf --genotype hom-ref --ind 001,002 --missing-value keep --out output.vcf
```

Values for **--genotype**:

value        | number of zero | A==B   | description 
-------------|----------------|--------|--------------------------------
hom-ref      | 2              | yes    | keep variants that have two refernce allele, e.g., 0/0
het          | 1              | no     | keep variants that have one reference and one alternative allele, e.g., 0/1, or 0/2
hom-alt      | 0              | yes    | keep variants that have two same alternative allele, e.g., 1/1, 2/2
het-alt      | 0              | no     | keep variants that have two different alternative allele, e.g., 1/2, 1/3
two-alt      | 0              | yes/no | keep variants that have two copy of alternative allele, e.g., 1/1 or 1/2
not-hom-ref  | 0 or 1         | yes/no | keep variants that does not have two copy of reference allele, e.g., 0/1 or 1/1 or 1/2
not-het      | 0 or 2         | yes/no | keep variants that are hom-ref or two-alt, e.g., 0/0, 1/1, 1/2
not-hom-alt  | 0 or 1 or 2    | no/yes | keep variants that are not hom-alt, e.g., 0/0, 0/1, 1/2 
not-two-alt  | 1 or 2         | yes/no | keep variants that does not have two alternative alleles, e.g., 0/0, 0/1, 0/2

The allele values are 0 for the reference allele, 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1 (A=0 and B=1, number of zero is 1).   

value of `--ind` is a string of individual IDs separated by ",". e.g., `--ind 1`, `--ind 1,2,3`

**Note**: individual IDs must be in the header line of the input VCF file.

### 3.9.3 Compound heterozygous filtering

Find all heterozygous variant pairs in a gene. We need perform the gene-based annotation first, suppose we have performed the [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/user-guide/gene) `refGene` based annotation. The key in the INFO field for gene name is `Gene_refGene`. The key in the INFO field for function is `Func_refGene` and values for function are: `exonic`, `splicing`, `ncRNA_exonic`, `ncRNA_intronic`, `UTR5`, `UTR3`, `intronic`, `upstream`, `downstream`, `intergenic`.

e.g. find all variant pairs that individual 001 and 002 have two heterozygous while other individuals only have one or zero heterozygous. Only variants in exonic, splicing and ncRNA_exonic region were considered.

```
python vcfFilter.py --vcf input.vcf --comp-het --gene-key Gene_refGene --func-key Func_refGene --func-values 'exonic,splicing,ncRNA_exonic' --ind '001,002' --out output.vcf
```

**Note**          

- key and values are based on annotation, check annotated vcf to identify them.    
- include variants in which region for analysis is based on what you are looking for (functional variants? regulatory variants?).       
- precedence of function values listed above are in decreased orders       
- generally, you do not add 'intergenic' to `--func-values`      

### 3.9.4 Compare genotype of multiple individuals

#### 3.9.4.1 Only variants have the same genotype across specified individuals will be kept

e.g., Only variants have the same genotype in individual `001` and `003` is kept.

```
python vcfFilter.py --vcf input.vcf --cmp-gtp-same --ind 001,003 --missing-value keep --out output.vcf
```

#### 3.9.4.2 Only variants have different genotype between the first individual and others will be kept

e.g., Only variants have the different genotype between individual `001` and `003` will be kept.

```
python vcfFilter.py --vcf input.vcf --cmp-gtp-diff --ind 001,003 --missing-value keep --out output.vcf
```

e.g., Only variants have the different genotype between individual `001` and `003` and `004` will be kept.

```
python vcfFilter.py --vcf input.vcf --cmp-gtp-diff --ind 001,003,004 --missing-value keep --out output.vcf
```

Genotype comparisions:

```
001 cmp 003
001 cmp 004
```

**Note**: in the above example, there is no comparision between genotype of 003 and 004.

# 4 References

* [VCF (Variant Call Format) version 4.0](http://www.1000genomes.org/wiki/Analysis/vcf4.0)
* [VCF (Variant Call Format) version 4.1](http://samtools.github.io/hts-specs/VCFv4.1.pdf)
* [VCF (Variant Call Format) version 4.2](http://samtools.github.io/hts-specs/VCFv4.2.pdf)
* [VCF (Variant Call Format) version 4.3](http://samtools.github.io/hts-specs/VCFv4.3.pdf)
