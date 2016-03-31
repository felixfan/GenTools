# vcfFilter - Python tool for handling VCF files.

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
