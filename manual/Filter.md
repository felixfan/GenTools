# pyVcfFilter - Filtering of variants in VCF file

Table of Contents
=================

  * [pyVcfFilter](#pyvcffilter)
    * [Variant Position Filtering](#variant-position-filtering)
      * [Chromosome filtering (\-\-chr)](#chromosome-filtering---chr)
        * [Keep all autosomes (separated by "\-")](#keep-all-autosomes-separated-by--)
        * [Keep chromosomes 1 and 6 (multiple \-\-chr)](#keep-chromosomes-1-and-6-multiple---chr)
        * [Keep all autosomes and chromosome X and Y ("\-" and multiple \-\-chr)](#keep-all-autosomes-and-chromosome-x-and-y---and-multiple---chr)
        * [Exclude chromosomes](#exclude-chromosomes)
      * [Region filtering (\-\-region and \-\-region\-file)](#region-filtering---region-and---region-file)
        * [Single region](#single-region)
        * [Multiple regions](#multiple-regions)
      * [Physical position filtering (\-\-phy\-pos and \-\-phy\-pos\-file)](#physical-position-filtering---phy-pos-and---phy-pos-file)
        * [Single position](#single-position)
        * [Multiple positions were stored in a file](#multiple-positions-were-stored-in-a-file)
    * [Variant ID Filtering (\-\-id and \-\-id\-file)](#variant-id-filtering---id-and---id-file)
      * [Single ID](#single-id)
      * [Multiple IDs using multiple \-\-id](#multiple-ids-using-multiple---id)
      * [Multiple IDs using \-\-id\-file](#multiple-ids-using---id-file)
      * [Remove variants by IDs](#remove-variants-by-ids)
    * [Variant Type Filtering (\-\-keep\-only\-indels and \-\-remove\-indels)](#variant-type-filtering---keep-only-indels-and---remove-indels)
      * [Keep only sites that contain an indel](#keep-only-sites-that-contain-an-indel)
      * [Exclude sites that contain an indel](#exclude-sites-that-contain-an-indel)
    * [FILTER Flag Filtering (\-\-filter)](#filter-flag-filtering---filter)
      * [Keep variants with FILTER flag: "PASS"](#keep-variants-with-filter-flag-pass)
      * [Keep variants with multiple FILTER flags](#keep-variants-with-multiple-filter-flags)
    * [INFO Field Filtering (\-\-info)](#info-field-filtering---info)
        * [Remove variants wilth alternatice allele frequency higher than 0\.01](#remove-variants-wilth-alternatice-allele-frequency-higher-than-001)
        * [Remove variants if the combined depth across samples is less than 200](#remove-variants-if-the-combined-depth-across-samples-is-less-than-200)
        * [Remove variants if the RMS mapping quality is less than 50](#remove-variants-if-the-rms-mapping-quality-is-less-than-50)
        * [Remove variants located in the intergenic, downstream or upstream region](#remove-variants-located-in-the-intergenic-downstream-or-upstream-region)
        * [Remove variants wilth alternatice allele frequency in 1000 genomes higher than 0\.05](#remove-variants-wilth-alternatice-allele-frequency-in-1000-genomes-higher-than-005)
        * [Remove variants were predicted as 'D' by SIFT](#remove-variants-were-predicted-as-d-by-sift)
    * [QUAL Field Filtering (\-\-qual)](#qual-field-filtering---qual)
      * [Keep variants with phred\-scaled quality score no less than 30](#keep-variants-with-phred-scaled-quality-score-no-less-than-30)
    * [Allele Filtering](#allele-filtering)
      * [Minimum / maximum number of alleles](#minimum--maximum-number-of-alleles)
      * [Biallelic / multiallelic site](#biallelic--multiallelic-site)
    * [Individual Filtering](#individual-filtering)
    * [Genotype Filtering](#genotype-filtering)
      * [Missing rate and missing count filtering (\-\-missing\-rate and \-\-missing\-count)](#missing-rate-and-missing-count-filtering---missing-rate-and---missing-count)
      * [Homozygous and heterozygous filtering (\-\-genotype)](#homozygous-and-heterozygous-filtering---genotype)
        * [Keep homozygous of reference alleles in sample 'ind1' and sample 'ind2'](#keep-homozygous-of-reference-alleles-in-sample-ind1-and-sample-ind2)
        * [Keep heterozygous variants in sample 'ind1' and sample 'ind2'](#keep-heterozygous-variants-in-sample-ind1-and-sample-ind2)
      * [Compound heterozygous filtering](#compound-heterozygous-filtering)
      * [Compare genotype of multiple individuals](#compare-genotype-of-multiple-individuals)
        * [Only variants have the same genotype across specified individuals will be kept](#only-variants-have-the-same-genotype-across-specified-individuals-will-be-kept)
        * [Only variants have different genotype between the first individual and others will be kept](#only-variants-have-different-genotype-between-the-first-individual-and-others-will-be-kept)

# pyVcfFilter

Filtering of variants according to genotypes and thresholds for the 8 fixed fields in VCF file.

## Variant Position Filtering

### Chromosome filtering (`--chr`)

Filtering by **CHROM** field. `--chr` accepts one argument. You can use `--chr` multiple times to specify multiple chromosome, or you can use '-' to specify multiple chromesomes (autosomes only).  

#### Keep all autosomes (separated by "-")

```
python pyVcfFilter.py --vcf example/input.vcf --chr chr1-chr22 --out example/output.vcf
```

#### Keep chromosomes 1 and 6 (multiple `--chr`)

```
python pyVcfFilter.py --vcf example/input.vcf --chr chr1 --chr chr6 --out example/output.vcf
```

#### Keep all autosomes and chromosome X and Y ("-" and multiple `--chr`)

```
python pyVcfFilter.py --vcf example/input.vcf --chr chr1-chr22 --chr chrX --chr chrY --out example/output.vcf
```

#### Exclude chromosomes

use `--reverse` to reverse the filter. e.g., exclude chr3 and chr6

```
python pyVcfFilter.py --vcf example/input.vcf --chr chr3 --chr chr6 --reverse --out example/output.vcf
```

### Region filtering (`--region` and `--region-file`)

Filtering by **CHROM** and **POS** fields. `--region` accepts three arguments: chromosome, start, end. `--region-file` accepts one argument, which is a file name. There are multiple rows in the file, and each row has three columns: chromosome, start, end.  e.g.  

```
chr1  1000000 3000000
chr6  140000000 150000000
```

#### Single region

Keep only the 10 Mb (140Mb-150Mb) region on chromosome 6:  

```
python pyVcfFilter.py --vcf example/input.vcf --region chr6  140000000 150000000 --out example/output.vcf
```

Exclude the 10 Mb (140Mb-150Mb) region on chromosome 6: 

```
python pyVcfFilter.py --vcf example/input.vcf --region chr6  140000000 150000000 --reverse --out example/output.vcf
```

#### Multiple regions

```
python pyVcfFilter.py --vcf example/input.vcf --region-file example/regions.txt --out example/output.vcf
```

```
python pyVcfFilter.py --vcf example/input.vcf --region-file example/regions.txt --out example/output.vcf --reverse
```

### Physical position filtering (`--phy-pos` and `--phy-pos-file`)

Filtering by **CHROM** and **POS** fields. `--phy-pos`  accepts two arguments: chromosome, position. `--phy-pos-file` accepts one argument, file name of the postion information. Each row of the file contains two columns, the first column is chromosome and the second column is position. e.g.

```
chr1  7589457
chr1  12861602
```

`--reverse` can be used to reverse the filter. 

#### Single position

```
python pyVcfFilter.py --vcf example/input.vcf --phy-pos chr1 12861602 --out example/output.vcf
```

```
python pyVcfFilter.py --vcf example/input.vcf --phy-pos chr1 12861602 --reverse --out example/output.vcf
```

#### Multiple positions were stored in a file

```
python pyVcfFilter.py --vcf example/input.vcf --phy-pos-file example/phypos.txt --out example/output.vcf
```

```
python pyVcfFilter.py --vcf example/input.vcf --phy-pos-file example/phypos.txt --reverse --out example/output.vcf
```

## Variant ID Filtering (`--id` and `--id-file`)

Filtering by **ID** field. `--id` accepts one argument. You can use `--id` multiple times to specify multiple SNP ID, or you can use `--id-file` to specify a file of multiple SNP IDs. Each row of the file contains one ID. e.g.     

```
rs123
rs124
```

### Single ID

```
python pyVcfFilter.py --vcf example/input.vcf --id rs123 --out example/output.vcf
```

### Multiple IDs using multiple `--id`
```
python pyVcfFilter.py --vcf example/input.vcf --id rs123 --id rs124 --out example/output.vcf
```

### Multiple IDs using `--id-file`

```
python pyVcfFilter.py --vcf example/input.vcf --id-file example/ids.txt --out example/output.vcf
```

### Remove variants by IDs

use `--reverse` to reverse the filter.

```
python pyVcfFilter.py --vcf example/input.vcf --id rs123 --id rs124 --reverse --out example/output.vcf
```

```
python pyVcfFilter.py --vcf example/input.vcf --id-file example/ids.txt --reverse --out example/output.vcf
```

## Variant Type Filtering (`--keep-only-indels` and `--remove-indels`)

Filtering by **REF** and **ALT** fields. If one site only has one ALT allele, "indel" here means that the length of REF and ALT allele is different. If one site has multiple ALT alleles, "indel" here means that at least one ALT allele has different length from REF allele. Site 1 and 2 will be defined as "indel" in the forllowing four sites.

```
site   REF                ALT
1      TTTTA              T,TTTTATTTA
2      A                  ATGTG,G,ATG
3      G                  C,T
4      T                  G
```

### Keep only sites that contain an indel

```
python pyVcfFilter.py --vcf example/input.vcf --keep-only-indels --out example/output.vcf
```

### Exclude sites that contain an indel

```
python pyVcfFilter.py --vcf example/input.vcf --remove-indels --out example/output.vcf
```

## FILTER Flag Filtering (`--filter`)

Filtering by **FILTER** field. You can use `--filter` multiple times to specify multiple FILTER flags.

### Keep variants with FILTER flag: "PASS"

```
python pyVcfFilter.py --vcf example/input.vcf --filter PASS --out example/output.vcf
```

### Keep variants with multiple FILTER flags

```
python pyVcfFilter.py --vcf example/input.vcf --filter PASS --filter VQSRTrancheINDEL99.90to100.00 --filter VQSRTrancheSNP99.90to100.00 --out example/output.vcf
```

## INFO Field Filtering (`--info`)

Filtering by **INFO** field. This option filter on the presence of the key and its value. INFO fields in VCF are encoded as a semicolon-separated series of short keys with optional values in the format:

```
<key>=<value>[,value]
```

When [annovar](http://annovar.openbioinformatics.org/en/latest/) was used to annotate the variants, multiple additional key and value pairs will be added in the INFO feild. Suppose `.` was used to reprent the missing value.   

- Whick `key` can be used depends on your own vcf. 
- There are six operators `>=`, `>`, `<=`, `<`, `==`, `!=`. 
- The first four operators `>=`, `>`, `<=`, `<` only accept numeric value. 
- The last two operators `==`, `!=` can accept both numeric value and string value.
- if there are multiple values for the key in vcf, all values must pass the cutoff.
- multiple value can be specified using `--value` multiple times.

The following examples show how it works. **Please be noted that annovar was used to annotate the vcf, so there are some keys may not in your vcf.**

#### Remove variants wilth alternatice allele frequency higher than 0.01

```
python pyVcfFilter.py --vcf example/input.vcf --info --key 'AF' --operator '<=' --value 0.01 --out example/output.vcf
```

Note: Currently, this only works for variants with only one alternative allele.   

#### Remove variants if the combined depth across samples is less than 200 

```
python pyVcfFilter.py --vcf example/input.vcf --info --key 'DP' --operator '>=' --value 200 --out example/output.vcf
```

#### Remove variants if the RMS mapping quality is less than 50

```
python pyVcfFilter.py --vcf example/input.vcf --info --key 'MQ' --operator '>=' --value 50 --out example/output.vcf
```

#### Remove variants located in the intergenic, downstream or upstream region

```
python pyVcfFilter.py --vcf example/input.vcf --out example/output.vcf --info --key 'Func_refGene' --operator '!=' --value 'intergenic' --value 'downstream' --value 'upstream' 
```

By default, missing value will be kept, you can change this by adding `--missing-value rm`

```
python pyVcfFilter.py --vcf example/input.vcf --out example/output.vcf --info --key 'Func_refGene' --operator '!=' --value 'intergenic' --value 'downstream' --value 'upstream' --missing-value rm
```

#### Remove variants wilth alternatice allele frequency in 1000 genomes higher than 0.05

```
python pyVcfFilter.py --vcf example/input.vcf --info --key '1000g2015aug_all' --operator '<=' --value 0.05 --missing-value keep --out example/output.vcf
```

#### Remove variants were predicted as 'D' by SIFT 

```
python pyVcfFilter.py --vcf example/input.vcf --info --key 'SIFT_pred' --operator '==' --value 'D' --missing-value keep --out example/output.vcf
```

## QUAL Field Filtering (`--qual`)

Filtering by **QUAL** field.  

### Keep variants with phred-scaled quality score no less than 30

```
python pyVcfFilter.py --vcf example/input.vcf --qual 30 --out example/output.vcf
```

## Allele Filtering

Filtering by **ALT** field. 

### Minimum / maximum number of alleles

Only sites have alleles no less than (`>=`) the minimum number will be kept.

```
python pyVcfFilter.py --vcf example/input.vcf --min-alleles 2 --out example/output.vcf
```

Only sites have alleles no more than (`<=`) the maximum number will be kept.

```
python pyVcfFilter.py --vcf example/input.vcf --max-alleles 2 --out example/output.vcf
```

### Biallelic / multiallelic site

A biallelic site is a specific locus in a genome that contains two observed alleles, counting the reference as one, and therefore allowing for one variant allele.

```
python pyVcfFilter.py --vcf example/input.vcf --min-alleles 2 --max-alleles 2 --out example/output.vcf
```

A multiallelic site is a specific locus in a genome that contains three or more observed alleles, again counting the reference as one, and therefore allowing for two or more variant alleles.

```
python pyVcfFilter.py --vcf example/input.vcf --min-alleles 3 --out example/output.vcf
```

## Individual Filtering

Include or exclude certain individuals.

```
python pyVcfFilter.py --vcf example/input.vcf --keep-inds --ind 'ind1' --ind 'ind2' --out example/output.vcf
```

```
python pyVcfFilter.py --vcf example/input.vcf --remove-inds --ind 'ind1' --ind 'ind2' --out example/output.vcf
```

## Genotype Filtering

### Missing rate and missing count filtering (`--missing-rate` and `--missing-count`)

Exclude sites with proportion of missing data is larger than the cutoff, which is between 0 and 1, where 0 indicates no missing data allowed and 1 allows sites that are completely missing.  

```
python pyVcfFilter.py --vcf example/input.vcf --missing-rate 0.05 --out example/output.vcf
```

Exclude sites that the number of individuals who missing this site is larger than the cutoff. A low missing count cutoff is similar with a low missing rate cutoff (`missing_rate = missing_count / total_individuals`). 

```
python pyVcfFilter.py --vcf example/input.vcf --missing-count 10 --out example/output.vcf
```

### Homozygous and heterozygous filtering (`--genotype`)

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

**Note**: individual IDs must be in the header line of the input VCF file.

#### Keep homozygous of reference alleles in sample 'ind1' and sample 'ind2'

```
python pyVcfFilter.py --vcf example/input.vcf --genotype hom-ref --ind ind1 --ind ind2 --missing-value keep --out example/output.vcf
```

#### Keep heterozygous variants in sample 'ind1' and sample 'ind2'

```
python pyVcfFilter.py --vcf example/input.vcf --genotype het --ind ind1 --ind ind2 --missing-value rm --out example/output.vcf
```

### Compound heterozygous filtering

Find all heterozygous variant pairs in a gene. We need perform the gene-based annotation first, suppose we have performed the [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/user-guide/gene) `refGene` based annotation. The key in the INFO field for gene name is `Gene_refGene`. The key in the INFO field for function is `Func_refGene` and values for function are: `exonic`, `splicing`, `ncRNA_exonic`, `ncRNA_intronic`, `UTR5`, `UTR3`, `intronic`, `upstream`, `downstream`, `intergenic`.

e.g. find all variant pairs that individual 'ind1' and 'ind' have two heterozygous while other individuals only have one or zero heterozygous. Only variants in exonic, splicing and ncRNA_exonic region were considered.  

```
python pyVcfFilter.py --vcf example/input.vcf --comp-het --gene-key 'Gene_refGene' --func-key 'Func_refGene' --func-values 'exonic' --func-values 'splicing' --func-values 'ncRNA_exonic' --func-values 'intronic' --ind 'ind1' --ind 'ind2' --out example/output.vcf
```

**Note**          

- key and values are based on annotation, check annotated vcf to identify them.    
- include variants in which region for analysis is based on what you are looking for (functional variants? regulatory variants?).       
- precedence of function values listed above are in decreased orders       
- generally, you do not add 'intergenic' to `--func-values`      

### Compare genotype of multiple individuals

#### Only variants have the same genotype across specified individuals will be kept

e.g., Only variants have the same genotype in individual `ind1` and `ind2` is kept.

```
python pyVcfFilter.py --vcf example/input.vcf --cmp-gtp-same --ind ind1 --ind ind2 --missing-value keep --out example/output.vcf
```

#### Only variants have different genotype between the first individual and others will be kept

e.g., Only variants have the different genotype between individual `ind1` and `ind2` will be kept.

```
python pyVcfFilter.py --vcf example/input.vcf --cmp-gtp-diff --ind ind1 --ind ind2 --missing-value keep --out example/output.vcf
```

e.g., Only variants have the different genotype between individual `ind1` and `ind2` and `ind3` will be kept.

```
python pyVcfFilter.py --vcf example/input.vcf --cmp-gtp-diff --ind ind1 --ind ind2 --ind ind3 --missing-value keep --out example/output.vcf
```

Genotype comparisions:

```
ind1 cmp ind2
ind1 cmp ind3
```

**Note**: in the above example, there is no comparision between genotype of ind2 and ind3.
