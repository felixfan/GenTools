# -*- coding: utf-8 -*-

import argparse
import time
import sys
import operator

def cal_het(infile, outfile):
    '''
    calculates heterozygosity on a per-individual basis, only use biallelic variants on autosomes
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    het = []
    hom_alt = []
    hom_ref = []
    ind = []
    n = 0
    for r in fr:
        r = r.strip()
        if r.startswith("##"):
            pass
        elif r.startswith("#"):
            ind = r.split()[9:]
            n = len(ind)
            het = [0] * n
            hom_alt = [0] * n
            hom_ref = [0] * n
        else:
            arr = r.split()
            if (len(arr[0]) > 2 and arr[0][3:].isdigit() and int(arr[0][3:]) < 23) or (arr[0].isdigit() and int(arr[0]) < 23):
                if arr[4].find(',') == -1:
                    for i in xrange(n):
                        if arr[i+9] != ".":
                            gtp = arr[i+9][0:3]
                            if gtp.count('0') == 1:
                                het[i] += 1
                            elif gtp.count('0') == 0:
                                hom_alt[i] += 1
                            elif gtp.count('0') == 2:
                                hom_ref[i] += 1
                            else:
                                sys.exit("wrong genotype format: {}\n".format(gtp))

    fr.close()
    print('Outputting Individual Heterozygosity: Only using biallelic variants on autosomes.\n')
    fw.write('{}\t{}\t{}\t{}\t{}\n'.format('indvID', 'numHet', 'numHomRef', 'numHomAlt', 'HetRatio'))
    print('{}\t{}\t{}\t{}\t{}'.format('indvID', 'numHet', 'numHomRef', 'numHomAlt', 'HetRatio'))
    for i in xrange(n):
        hetratio = 1.0 * het[i] / hom_alt[i]
        fw.write('{}\t{}\t{}\t{}\t{:.4f}\n'.format(ind[i], het[i], hom_ref[i], hom_alt[i], hetratio))
        print('{}\t{}\t{}\t{}\t{:.4f}'.format(ind[i], het[i], hom_ref[i], hom_alt[i], hetratio))
    fw.close()
    print("\nWrite results to {}".format(outfile))

def concat_vcf(infile1, infile2, outfile):
    '''
    concatenates VCF files, input file should have same columns.
    '''
    n = 0
    m = 0
    z = 0
    fr = open(infile1)
    fw = open(outfile, 'w')
    var = []
    for r in fr:
        r = r.strip()
        arr = r.split()
        if not r.startswith("#"):
            n += 1
            var.append(arr[0]+':'+arr[1])
        fw.write('{}\n'.format(r))
    fr.close()
    fr = open(infile2)
    for r in fr:
        r = r.strip()
        arr = r.split()
        if not r.startswith("#"):
            m += 1
            tmp = arr[0]+':'+arr[1]
            if not tmp in var:
                fw.write('{}\n'.format(r))
                z += 1
    fr.close()
    fw.close()
    print("There are {} variants in {}".format(n, infile1))
    print("There are {} variants in {}, and {} variants are already in {}".format(m, infile2, m-z, infile1))
    print("Write {} variants to {}".format(n+z, outfile))

def stat_field(infile, outfile, col):
    '''
    calculate how many of each value in the field.
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    value = {}
    for r in fr:
        r = r.strip()
        if not r.startswith("#"):
            arr = r.split()
            if arr[col] in value:
                value[arr[col]] += 1
            else:
                value[arr[col]] = 1
    fr.close()
    for k,v in sorted(value.items(), key=lambda x: x[0]):
        print("{}\t{}".format(k, v))
        fw.write("{}\t{}\n".format(k, v))

def stat_info(infile, outfile, key):
    '''
    calculate how many of each value for the key in the INFO field.
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    value = {}
    for r in fr:
        r = r.strip()
        if not r.startswith("#"):
            arr = r.split()
            n += 1
            m = len(arr) - 9
            info = arr[7]
            if -1 != info.find(key):
                ar = info.split(';')
                for z in ar:
                    if z.startswith(key+'='):
                        a = z.split('=')
                        if a[1] in value:
                            value[a[1]] += 1
                        else:
                            value[a[1]] = 1
            else:
                sys.exit('{} is not in the INFO field'.format(key))
    print "There are {} individuals with {} variants in {}".format(m, n, infile)
    for k, v in sorted(value.items(), key=lambda x: x[1], reverse=True): # sort by value
        print '{}\t{}'.format(k, v)
        fw.write('{}\t{}\n'.format(k, v))
    fw.close()
    fr.close()

def cal_titv(infile, outfile):
    '''
    calculates transition transversion ratio, only use biallelic SNPs on autosomes
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    ac = []
    ag = []
    at = []
    cg = []
    ct = []
    gt = []
    ind = []
    n = 0
    for r in fr:
        r = r.strip()
        if r.startswith("##"):
            pass
        elif r.startswith("#"):
            ind = r.split()[9:]
            n = len(ind)
            ac = [0] * n
            ag = [0] * n
            at = [0] * n
            cg = [0] * n
            ct = [0] * n
            gt = [0] * n
        else:
            arr = r.split()
            if (len(arr[0]) > 2 and arr[0][3:].isdigit() and int(arr[0][3:]) < 23) or (arr[0].isdigit() and int(arr[0]) < 23):
                if len(arr[4])== 1:
                    for i in xrange(n):
                        if arr[i+9] != ".":
                            gtp = arr[i+9][0:3]
                            if gtp.count('0') == 1:
                                if arr[3] + arr[4] == 'AC' or arr[3] + arr[4] == 'CA':
                                    ac[i] += 1
                                elif arr[3] + arr[4] == 'AG' or arr[3] + arr[4] == 'GA':
                                    ag[i] += 1
                                elif arr[3] + arr[4] == 'AT' or arr[3] + arr[4] == 'TA':
                                    at[i] += 1
                                elif arr[3] + arr[4] == 'CG' or arr[3] + arr[4] == 'GC':
                                    cg[i] += 1
                                elif arr[3] + arr[4] == 'CT' or arr[3] + arr[4] == 'TC':
                                    ct[i] += 1
                                elif arr[3] + arr[4] == 'GT' or arr[3] + arr[4] == 'TG':
                                    gt[i] += 1
    fr.close()
    print('Outputting Individual transition transversion ratio: Only using biallelic SNPs on autosomes.\n')
    fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('indvID', 'AC', 'AG', 'AT', 'CG', 'CT', 'GT', "Ti", "Tv", "TiTv"))
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('indvID', 'AC', 'AG', 'AT', 'CG', 'CT', 'GT', "Ti", "Tv", "TiTv"))
    for i in xrange(n):
        ti = ag[i] + ct[i]
        tv = ac[i] + at[i] + cg[i] + gt[i]
        titv = 1.0 * ti / tv
        fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}\n'.format(ind[i], ac[i], ag[i], at[i], cg[i], ct[i], gt[i], ti, tv, titv))
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4f}'.format(ind[i], ac[i], ag[i], at[i], cg[i], ct[i], gt[i], ti, tv, titv))
    fw.close()
    print("\nWrite results to {}".format(outfile))

def convert_plink_tped(infile, outfile):
    '''
    output the genotype data in PLINK tped/tfam format, only biallelic variants will be output.
    '''
    pass

def diff_site(infile1, infile2, outfile):
    '''
    output sites that are common /unique to each file
    '''
    pass

def diff_indv(infile1, infile2, outfile):
    '''
    output individuals that are common /unique to each file
    '''
    pass

def basic_stats(infile, outfile):
    '''
    output some basic statistics: number of individuals, total number of variants, number of SNPs and indels,
    number of biallelic variants, multiple allelic variants
    '''
    pass

def run_time(starttime):
    usedtime = time.time() - starttime
    print
    print "Time used:",
    if usedtime >=60:
        ts = int(usedtime) % 60
        usedtime = int(usedtime) / 60
        tm = int(usedtime) % 60
        usedtime = int(usedtime) / 60
        th = int(usedtime) % 60
        if th > 0:
            print "%d hours"  % th,
            print "%d minutes"  % tm,
        elif tm > 0:
            print "%d minutes"  % tm,
    else:
        ts = usedtime
    print '%.2f seconds' % ts
    print "Finished at ",
    print time.strftime("%H:%M:%S %d %b %Y")

if __name__ == '__main__':
    starttime = time.time()
    desc = '''Python tool for basic statistics and manipulation of VCF file'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', action='version', version='%(prog)s 0.3.0')
    ### input
    parser.add_argument('--vcf', help='input vcf file', required=True, type=str)
    parser.add_argument('--vcf2', help='the second input vcf file', type=str)
    ### funcs
    parser.add_argument('--het', help='calculates heterozygosity on a per-individual basis', action='store_true')
    parser.add_argument('--concat', help='concatenates VCF files', action='store_true')
    parser.add_argument('--stat-chr', help='statistics based field in VCF files', action='store_true')
    parser.add_argument('--stat-filter', help='statistics based field in VCF files', action='store_true')
    parser.add_argument('--stat-info', help='statistics based on info field', action='store_true')
    parser.add_argument('--key', help='key words to be used', type=str)
    parser.add_argument('--titv', help='calculates Ti/Tv ratio on a per-individual basis', action='store_true')
    ### output
    parser.add_argument('--out', help='output file', type=str, default='output.txt')
    ### parameter
    args = vars(parser.parse_args())
    INFILE = args['vcf'] if 'vcf' in args else None
    OUTFILE = args['out']
    STATINFO = args['stat_info'] if 'stat_info' in args else False
    KEY = args['key'] if 'key' in args else None
    CONCAT = args['concat'] if 'concat' in args else False
    INFILE2 = args['vcf2'] if 'vcf2' in args else None
    HET = args['het'] if 'het' in args else False
    STATCHR = args['stat_chr'] if 'stat_chr' in args else None
    STATFILTER = args['stat_filter'] if 'stat_filter' in args else None
    TITV = args['titv'] if 'titv' in args else False
    ### log
    print "@-------------------------------------------------------------@"
    print "|     pyVcfManipulator    |      v0.3.0     |   27 Jul 2016   |"
    print "|-------------------------------------------------------------|"
    print "|  (C) 2016 Felix Yanhui Fan, GNU General Public License, v3  |"
    print "|-------------------------------------------------------------|"
    print "|    For documentation, citation & bug-report instructions:   |"
    print "|            http://felixfan.github.io/PyTV                   |"
    print "@-------------------------------------------------------------@"
    print "\n\tOptions in effect:"
    print "\t--vcf", INFILE
    if STATINFO:
        print '\t--stat-info'
        if KEY:
            print '\t--key', KEY
        else:
            sys.exit('option --key is missing')
    elif CONCAT:
        print '\t--concat'
        if INFILE2:
            print '\t--vcf2', INFILE2
        else:
            sys.exit('option --vcf2 is missing')
    elif HET:
        print '\t--het'
    elif STATCHR:
        print '\t--stat-chr'
    elif STATFILTER:
        print '\t--stat-filter'
    elif TITV:
        print '\t--titv'
    print "\t--out", OUTFILE
    print
    ### run
    if STATINFO:
        stat_info(INFILE, OUTFILE, KEY)
    elif CONCAT:
        concat_vcf(INFILE, INFILE2, OUTFILE)
    elif HET:
        cal_het(INFILE, OUTFILE)
    elif STATCHR:
        stat_field(INFILE, OUTFILE, 0)
    elif STATFILTER:
        stat_field(INFILE, OUTFILE, 6)
    elif TITV:
        cal_titv(INFILE, OUTFILE)
    else:
        print "Do nothing"
    ###time
    run_time(starttime)
