# -*- coding: utf-8 -*-

import argparse
import time
import sys

def split_str_comma(CHR):
    '''
    split string with comma
    return a list of string
    '''
    chrs = []
    if -1 == CHR.find(','):
        chrs.append(CHR)
    else:
        tmp = CHR.split(',')
        for t in tmp:
            chrs.append(t)
    return chrs

def split_str_dash(CHR):
    '''
    split string with dash (-)
    return a list of string
    '''
    chrs = []
    if -1 != CHR.find('-'):
        if -1 == CHR.find('chr') and -1 == CHR.find('CHR'):
            tmp = CHR.split('-')
            start = int(tmp[0])
            end = int(tmp[1])+1
            for t in range(start, end):
                chrs.append(str(t))
        else:
            tmp = CHR.split('-')
            start = int(tmp[0][3:])
            end = int(tmp[1][3:])+1
            for t in range(start, end):
                chrs.append(tmp[0][:3]+str(t))
    return chrs

def split_str_comma_dash(CHR):
    '''
    split string with comma (,) and/or dash (-)
    return a list of string
    '''
    chrs = []
    if -1==CHR.find(','):
        if -1 == CHR.find('-'):
           chrs.append(CHR)
        else:
            chrs.extend(split_str_dash(CHR))
    else:
        tmp = CHR.split(',')
        for it in tmp:
            if -1 == it.find('-'):
                chrs.append(it)
            else:
                chrs.extend(split_str_dash(it))
    return chrs
  
def filter_by_chr(infile, chrs,outfile):
    '''
    filter by chromosome
    ''' 
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    for r in fr:
        r = r.strip()
        if r.startswith("#"):
            fw.write("%s\n" % r)
        else:
            arr = r.split()
            n += 1
            if arr[0] in chrs:
                fw.write("%s\n" % r)
                m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
    
def filter_by_pos(infile, chr, start, end, outfile):
    '''
    filter by postion
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    for r in fr:
        r = r.strip()
        if r.startswith("#"):
            fw.write("%s\n" % r)
        else:
            arr = r.split()
            n += 1
            if arr[0] == chr and int(arr[1]) >= start and int(arr[1]) <= end:
                fw.write("%s\n" % r)
                m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
    
def filter_by_qual(infile, cutoff, outfile):
    '''
    filter by qual >= cutoff
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    for r in fr:
        r = r.strip()
        if r.startswith("#"):
            fw.write("%s\n" % r)
        else:
            arr = r.split()
            n += 1
            if float(arr[5]) >= cutoff:
                fw.write("%s\n" % r)
                m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()

def filter_by_filter(infile, flts, outfile):
    '''
    keep variants with FILTER flag: flts
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    for r in fr:
        r = r.strip()
        if r.startswith("#"):
            fw.write("%s\n" % r)
        else:
            arr = r.split()
            n += 1
            if arr[6] in flts:
                fw.write("%s\n" % r)
                m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
    
def check_allele(alleleStr, gtp, na):
    if alleleStr.startswith('./.') or alleleStr.startswith('.|.') or alleleStr == ".":  # missing vale
        if na == "keep":
            return True
        else:
            return False
    else:
        s = alleleStr.split(':')[0]
        if -1 != s.find('/') or -1 != s.find('|'): # gtp separate with '/' or '|', else return False
            if gtp == "hom-ref":
                if 2 == s.count('0'):
                    return True
                else:
                    return False
            elif gtp == "het":
                if 1 == s.count('0'):
                    return True
                else:
                    return False
            elif gtp == "hom-alt":
                if 0 == s.count('0') and s[0] == s[2]:
                    return True
                else:
                    return False
            elif gtp == "not-hom-ref":
                if s.count('0') < 2:
                    return True
                else:
                    return False
            elif gtp == "not-two-alt":
                if s.count('0') > 0:
                    return True
                else:
                    return False
            elif gtp == "het-alt":
                if 0 == s.count('0') and s[0] != s[2]:
                    return True
                else:
                    return False
            elif gtp == "two-alt":
                if 0 == s.count('0'):
                    return True
                else:
                    return False
            elif gtp == "not-het":
                if 1 != s.count('0'):
                    return True
                else:
                    return False
        else:
            return False
                    
def filter_by_genotype(infile, gtp, inds, na, outfile):
    '''
    keep variants using genotype in samples
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    for r in fr:
        r = r.strip()
        if r.startswith("##"):
            fw.write("%s\n" % r)
        elif r.startswith("#"):
            fw.write("%s\n" % r)
            myind = r.split()
        else:
            arr = r.split()
            n += 1
            flag = True
            for i in range(9,len(arr)):
                if myind[i] in inds:
                    flag = (flag and check_allele(arr[i],gtp,na))
                if flag == False:
                    break
            if flag:
                fw.write("%s\n" % r)
                m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
#######################################################
strattime = time.time()
#######################################################
parser = argparse.ArgumentParser(description='VCF Filter', prog="vcfFilter.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')
### input
parser.add_argument('-vcf', '--vcf', help='input vcf file', required=True, type=str)
### filters
parser.add_argument('-chr', '--chromosome', help='filter by chromosome', type=str)
parser.add_argument('-pos', '--position', help='filter by position', type=str)
parser.add_argument('-qual', '--qual', help='filter by qual score', type=float)
parser.add_argument('-filter', '--filter', help='filter by filter flag', type=str)
parser.add_argument('-gtp', '--genotype', help='filter by genotype', type=str,choices=["hom-ref", "hom-alt","het", "het-alt","not-hom-ref","not-two-alt","two-alt","not-het"])
###individual
parser.add_argument('-ind', '--individual', help='individual id', type=str)
###missing value
parser.add_argument('-mv', '--missingvalue', help='how to deal with missing values', default="keep", type=str, choices=["keep", "rm"])
### output
parser.add_argument('-o', '--out', help='output vcf file', default='output.vcf')
#######################################################
args = vars(parser.parse_args())
INFILE = args['vcf']
OUTFILE = args['out']
CHR = args['chromosome'] if 'chromosome' in args else None
POS = args['position'] if 'position' in args else None
QUAL = args['qual'] if 'qual' in args else None
FILTER = args['filter'] if 'filter' in args else None
GENOTYPE = args['genotype'] if 'genotype' in args else None
IND = args['individual'] if 'individual' in args else None
NA = args['missingvalue']
#######################################################
print "@-------------------------------------------------------------@"
print "|       vcfFilter     |     v1.0.0      |    14 April 2016    |"
print "|-------------------------------------------------------------|"
print "|  (C) 2015 Felix Yanhui Fan, GNU General Public License, v2  |"
print "|-------------------------------------------------------------|"
print "|    For documentation, citation & bug-report instructions:   |"
print "|          http://felixfan.github.io/vcfFilter                |"
print "@-------------------------------------------------------------@"
print "\n\tOptions in effect:"
print "\t-vcf", INFILE
if CHR:
    print "\t-chr", CHR
elif POS:
    print "\t-pos", POS
elif QUAL:
    print "\t-qual", QUAL
elif FILTER:
    print "\t-filter", FILTER
elif GENOTYPE:
    print "\t-gtp", GENOTYPE
    if IND:
        print "\t-ind", IND
    else:
        sys.exit("Error, argment -ind is missing!")
    print "\t-mv", NA
print "\t-o", OUTFILE
print
#######################################################
if CHR:
    chrs = split_str_comma_dash(CHR)
    print "keep chromosomes:", chrs
    filter_by_chr(INFILE,chrs, OUTFILE)
elif POS:
    tmp = POS.split(':')
    chr = tmp[0]
    start = int(tmp[1].split('-')[0])
    end = int(tmp[1].split('-')[1])
    print "keep region: %s, from %d to %d" % (chr, start, end)
    filter_by_pos(INFILE,chr,start,end,OUTFILE)
elif QUAL:
    cutoff = float(QUAL)
    print "keep variants with qual score no less than %f" % cutoff
    filter_by_qual(INFILE, cutoff, OUTFILE)
elif FILTER:
    flts = split_str_comma(FILTER)
    print "keep variants with FILTER flag:", flts
    filter_by_filter(INFILE, flts, OUTFILE)
elif GENOTYPE:
    print "genotype filter:"
    inds = split_str_comma(IND)
    if GENOTYPE == "hom-ref":
        print "keep homozygous of reference allele of these individuals:"
    elif GENOTYPE == "het":
        print "keep heterozygous of these individuals:"
    elif GENOTYPE == "hom-alt":
        print "keep homozygous of alternative allele of these individuals:"
    elif GENOTYPE == "not-hom-ref":
        print "keep heterozygous AND homozygous of alternative allele of these individuals:"
    elif GENOTYPE == "not-two-alt":
        print "keep heterozygous AND homozygous of reference allele of these individuals:"
    elif GENOTYPE == "het-alt":
        print "keep two different alternative alleles of these individuals:"
    print inds
    filter_by_genotype(INFILE, GENOTYPE, inds, NA, OUTFILE)
else:
    pass
###############################################################################
usedtime = time.time() - strattime
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