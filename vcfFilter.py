# -*- coding: utf-8 -*-

import argparse
import time
import sys

'''
split string with dash (-)
return a list od string
'''
def splitStrDash(CHR):
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

'''
split string with comma (,) and/or dash (-)
return a list od string
'''
def splitStrCommaDash(CHR):
    chrs = []
    if -1==CHR.find(','):
        if -1 == CHR.find('-'):
           chrs.append(CHR)
        else:
            chrs.extend(splitStrDash(CHR))
    else:
        tmp = CHR.split(',')
        for it in tmp:
            if -1 == it.find('-'):
                chrs.append(it)
            else:
                chrs.extend(splitStrDash(it))
    return chrs

'''
filter by chromosome
'''   
def chrFilterByChr(infile, chrs,outfile):
    fr = open(infile)
    fw = open(outfile, 'w')
    for r in fr:
        r = r.strip()
        if r.startswith("#"):
            fw.write("%s\n" % r)
        else:
            arr = r.split()
            if arr[0] in chrs:
                fw.write("%s\n" % r)
    fw.close()
    fr.close()
#######################################################
strattime = time.time()
#######################################################
parser = argparse.ArgumentParser(description='VCF Filter', prog="vcfFilter.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')
### input
parser.add_argument('-vcf', '--vcfFile', help='input vcf file', required=True, type=str)
### filters
parser.add_argument('-chr', '--chromosome', help='filter by chromosome', type=str)
###missing value
parser.add_argument('-mv', '--missingValue', help='how to deal with missing values', default="keep", type=str, choices=["keep", "rm"])
### output
parser.add_argument('-o', '--out', help='output vcf file', default='output.vcf')
#######################################################
args = vars(parser.parse_args())
INFILE = args['vcfFile']
OUTFILE = args['out']
CHR = args['chromosome'] if 'chromosome' in args else None
NA = args['missingValue']
#######################################################
print "@-------------------------------------------------------------@"
print "|       vcfFilter     |     v1.0.0      |     31 Mar 2016     |"
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
print "\t-mv", NA
print "\t-o", OUTFILE
print
#######################################################
if CHR:
    chrs = splitStrCommaDash(CHR)
    print "keep chromosomes:", chrs
    chrFilterByChr(INFILE,chrs, OUTFILE)
else:
    pass