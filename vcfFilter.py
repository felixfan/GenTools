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
    
def filter_by_region(infile, chr, start, end, outfile):
    '''
    filter by region
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
    
def filter_indels(infile, b, outfile):
    '''
    filter indel
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
            if -1 == arr[4].find(','): # only one alt
                if b: # remove indel
                    if len(arr[3]) == len(arr[4]):
                        fw.write("%s\n" % r)
                        m += 1
                else: # remove snp
                    if len(arr[3]) != len(arr[4]):
                        fw.write("%s\n" % r)
                        m += 1           
            else: # multiple alt
                alleles = arr[4].split(',')
                if b: # remove indel
                    flag = True
                    for allele in alleles:
                        if len(arr[3]) != len(allele):
                            flag = False
                            break
                    if flag:
                        fw.write("%s\n" % r)
                        m += 1
                else: # remove snp
                    flag = False
                    for allele in alleles:
                        if len(arr[3]) != len(allele):
                            flag = True
                            break
                    if flag:
                        fw.write("%s\n" % r)
                        m += 1        
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
def filter_by_var_ids(infile, ids, outfile):
    '''
    keep variants in ids
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
            if arr[2] in ids:
                fw.write("%s\n" % r)
                m += 1 # do NOT rm ID that already print out, since there are more variant's ID were '.'
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
def filter_by_phypos(infile, pp, outfile):
    '''
    keep variants by physical positions
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
            tmp = arr[0] + ":" + arr[1]
            if tmp in pp:
                fw.write("%s\n" % r)
                m += 1
                pp.remove(tmp) # rm variants that already print out
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
def cmp_gtp(infile, inds, na, outfile,same):
    '''
    compare genotype of multiple individuals
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    idx0 = 0
    idx = []
    for r in fr:
        r = r.strip()
        if r.startswith("##"):
            fw.write("%s\n" % r)
        elif r.startswith("#"):
            fw.write("%s\n" % r)
            myind = r.split()
            if inds[0] in myind:
                idx0 = myind.index(inds[0])
            else:
                sys.exit('%s is not in %s' % (inds[0], infile))
            for x in inds[1:]:
                if x in myind:
                    tmp = myind.index(x)
                    idx.append(tmp)
                else:
                    sys.exit('%s is not in %s' % (x, infile))
        else:
            arr = r.split()
            n += 1
            flag = True
            for i in xrange(9,len(arr)):
                if i in idx:
                    if same: # same genotype
                        if na == 'keep': # keep missing
                            if '.' != arr[idx0] and  '.' != arr[i]:
                                if arr[idx0][0:3] != arr[i][0:3]:
                                    flag = False
                                    break
                        else: # remove missing
                            if '.' == arr[idx0] or  '.' == arr[i]:
                                flag = False
                                break
                            else:
                                if arr[idx0][0:3] != arr[i][0:3]:
                                    flag = False
                                    break
                    else: # diff genotype
                        if na == 'keep': # keep missing
                            if '.' != arr[idx0] and  '.' != arr[i]:
                                if arr[idx0][0:3] == arr[i][0:3]:
                                    flag = False
                                    break
                        else: # remove missing
                            if '.' == arr[idx0] or  '.' == arr[i]:
                                flag = False
                                break
                            else:
                                if arr[idx0][0:3] == arr[i][0:3]:
                                    flag = False
                                    break
            if flag:
                fw.write("%s\n" % r)
                m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
def filter_num_alleles(infile, minc, maxc, outfile):
    '''
    filter by num of alleles
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
            n += 1
            arr = r.split()
            num = len(arr[4].split(',')) + 1
            if minc and maxc:
                if maxc >= num >= minc:
                    fw.write("%s\n" % r)
                    m += 1
                elif minc:
                    if num >= minc:
                        fw.write("%s\n" % r)
                        m += 1
            elif maxc:
                if num <= maxc:
                    fw.write("%s\n" % r)
                    m += 1
            else:
                sys.exit('please check value of --min-alleles and/or --max-alleles')
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()
def filter_by_info(infile, outfile,key,operation,values,na,stype='float'):
	'''
	filter by info
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
			ar = arr[7].split(';')
			flag = False # whether key is in info
			for z in ar:
				if z.startswith(key):
					flag = True
					a = z.split('=')
					fla = True # keep site
					if a[1] == '.': # missing value
						if na == 'keep':
							fla = True
						else:
							fla = False
					else: # non-missing value
						bs = a[1].split(',') # multiple value
						if operation == ">=":
							for b in bs:
								if float(b) < values[0]:
									fla = False
									break
						elif operation == "<=":
							for b in bs:
								if float(b) > values[0]:
									fla = False
									break
						elif operation == "<":
							for b in bs:
								if float(b) >= values[0]:
									fla = False
									break
						elif operation == ">":
							for b in bs:
								if float(b) <= values[0]:
									fla = False
									break
						elif operation == "!=":
							for b in bs:
								if 'str' == stype:
									if b in values:
										fla = False
										break
								else:
									if float(b) == values:
										fla = False
										break
						elif operation == "=":
							for b in bs:
								if 'str' == stype:
									if not b in values:
										fla = False
										break
								else:
									if float(b) != values:
										fla = False
										break
					if fla:
							fw.write("%s\n" % r)
							m += 1
				if flag: # only one key in info, do not need to check other key
					break
			if not flag: # key is not in info, wrong key?
				sys.exit('{} is not find in site {} {}'.format(key, arr[0], arr[1]))
	print "%d of %d variants were written to %s" % (m, n, outfile)
	fw.close()
	fr.close()
def check_float(s):
	try:
		return float(s)
	except ValueError:
		return False
#######################################################
strattime = time.time()
#######################################################
desc = '''Filtering of polymorphisms according to genotypes, physical positions 
					and thresholds of quality score, read depth and others.'''

parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-v', action='version', version='%(prog)s 1.0.0')
### input
parser.add_argument('-vcf', help='input vcf file', required=True, type=str)
### filters
parser.add_argument('-chr', help='filter by chromosome', type=str)
parser.add_argument('-region', help='filter by region', type=str)
parser.add_argument('-qual', help='filter by qual score', type=float)
parser.add_argument('-filter', help='filter by filter flag', type=str)
parser.add_argument('-genotype', help='filter by genotype', type=str,choices=["hom-ref", "hom-alt","het", "het-alt","not-hom-ref","not-two-alt","two-alt","not-het"])
parser.add_argument('--keep-only-indels', help='keep only indels', action='store_true')
parser.add_argument('--remove-indels', help='remove indels', action='store_true')
parser.add_argument('-ids', help='filter by ID', type=str)
parser.add_argument('--ids-file', help='filter by ID', type=str)
parser.add_argument('--phy-pos', help='filter by physical position', type=str)
parser.add_argument('--phy-pos-file', help='filter by physical position', type=str)
parser.add_argument('--cmp-gtp-same', help='compare genotype of multiple individuals', action='store_true')
parser.add_argument('--cmp-gtp-diff', help='compare genotype of multiple individuals', action='store_true')
parser.add_argument('--min-alleles', help='minimum number of alleles', type=int)
parser.add_argument('--max-alleles', help='maximum number of alleles', type=int)
parser.add_argument('-info', help='filter by info keys', type=str)
### individual
parser.add_argument('-ind', help='individual id', type=str)
### missing value
parser.add_argument('--missing-value', help='how to deal with missing values', default="keep", type=str, choices=["keep", "rm"])
### output
parser.add_argument('-out', help='output vcf file', type=str, default='output.vcf')
#######################################################
args = vars(parser.parse_args())
INFILE = args['vcf']
OUTFILE = args['out']
CHR = args['chr'] if 'chr' in args else None
REGION = args['region'] if 'region' in args else None
QUAL = args['qual'] if 'qual' in args else None
FILTER = args['filter'] if 'filter' in args else None
GENOTYPE = args['genotype'] if 'genotype' in args else None
IND = args['ind'] if 'ind' in args else None
NA = args['missing_value']
INDEL = args['keep_only_indels'] if 'keep_only_indels' in args else False
SNP = args['remove_indels'] if 'remove_indels' in args else False
IDS = args['ids'] if 'ids' in args else None
IDSFILE = args['ids_file'] if 'ids_file' in args else None
PHYPOS = args['phy_pos'] if 'phy_pos' in args else None
PHYPOSFILE = args['phy_pos_file'] if 'phy_pos_file' in args else None
CMPGTPSAME = args['cmp_gtp_same'] if 'cmp_gtp_same' in args else None
CMPGTPDIFF = args['cmp_gtp_diff'] if 'cmp_gtp_diff' in args else None
MINALLELES = args['min_alleles'] if 'min_alleles' in args else None
MAXALLELES = args['max_alleles'] if 'max_alleles' in args else None
INFO = args['info'] if 'info' in args else None
#######################################################
print "@-------------------------------------------------------------@"
print "|       vcfFilter     |     v1.0.0      |    19 May 2016      |"
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
elif REGION:
    print "\t-region", REGION
elif QUAL:
    print "\t-qual", QUAL
elif FILTER:
    print "\t-filter", FILTER
elif GENOTYPE:
    print "\t-genotype", GENOTYPE
    if IND:
        print "\t-ind", IND
    else:
        sys.exit("Error, argment -ind is missing!")
    print "\t--missing-value", NA
elif INDEL:
    print "\t--keep-only-indels"
elif SNP:
    print "\t--remove-indels" 
elif IDS:
    print "\t-ids", IDS
elif IDSFILE:
    print "\t--ids-file", IDSFILE
elif PHYPOS:
    print "\t--phy-pos", PHYPOS
elif PHYPOSFILE:
    print "\t--phy-pos-file", PHYPOSFILE
elif CMPGTPSAME:
    print "\t--cmp-gtp-same"
    if IND:
        print "\t-ind", IND
    else:
        sys.exit("Error, argment -ind is missing!")
    print "\t--missing-value", NA
elif CMPGTPDIFF:
    print "\t--cmp-gtp-diff"
    if IND:
        print "\t-ind", IND
    else:
        sys.exit("Error, argment -ind is missing!")
    print "\t--missing-value", NA
elif MINALLELES and MAXALLELES:
    print "\t--min-alleles", MINALLELES
    print "\t--max-alleles", MAXALLELES
elif MINALLELES:
    print "\t--min-alleles", MINALLELES
elif MAXALLELES:
    print "\t--max-alleles", MAXALLELES
elif INFO:
    print "\t-info", INFO
    print "\t--missing-value", NA
print "\t-out", OUTFILE
print
#######################################################
if CHR:
    chrs = split_str_comma_dash(CHR)
    print "keep chromosomes:", chrs
    filter_by_chr(INFILE,chrs, OUTFILE)
elif REGION:
    tmp = REGION.split(':')
    chr = tmp[0]
    start = int(tmp[1].split('-')[0])
    end = int(tmp[1].split('-')[1])
    print "keep region: %s, from %d to %d" % (chr, start, end)
    filter_by_region(INFILE,chr,start,end,OUTFILE)
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
elif INDEL:
    print "keep only indels"
    INDEL = not INDEL
    filter_indels(INFILE, INDEL, OUTFILE)
elif SNP:
    print "remove indels"
    filter_indels(INFILE, SNP, OUTFILE)
elif IDS:
    print "keep variants by ID"
    ids = IDS.split(',')
    filter_by_var_ids(INFILE, ids, OUTFILE)
elif IDSFILE:
    print "keep variants by ID"
    ids = []
    tf = open(IDSFILE)
    for r in tf:
        r = r.strip()
        ids.append(r)
    tf.close()
    filter_by_var_ids(INFILE, ids, OUTFILE)
elif PHYPOS:
    print "keep variants by physical position"
    pp = PHYPOS.split(',')
    filter_by_phypos(INFILE, pp, OUTFILE)
elif PHYPOSFILE:
    print "keep variants by physical position"
    pp = []
    tf = open(PHYPOSFILE)
    for r in tf:
        r = r.strip()
        arr = r.split()
        pp.append(arr[0]+":"+arr[1])
    tf.close()
    filter_by_phypos(INFILE, pp, OUTFILE)
elif CMPGTPSAME:
    print "compare genotype of multiple individuals"
    print "only variants have the same genotype across the following individuals will be kept"
    inds = split_str_comma(IND)
    if len(inds) > 1:
        print "individuals:", inds
    else:
        sys.exit('at least two individuals should be provided')
    cmp_gtp(INFILE, inds, NA, OUTFILE,True)
elif CMPGTPDIFF:
    print "compare genotype of multiple individuals"
    inds = split_str_comma(IND)
    print "only variants have the different genotype between the first individual and others will be kept"
    if len(inds) > 1:
        print "individuals:", inds[0], 'vs.', inds[1:]
    else:
        sys.exit('at least two individuals should be provided')
    cmp_gtp(INFILE, inds, NA, OUTFILE,False)
elif MINALLELES and MAXALLELES:
    print "keep sites with minimum {} alleles and maximum {} alleles".format(MINALLELES, MAXALLELES)
    filter_num_alleles(INFILE, MINALLELES, MAXALLELES, OUTFILE)
elif MINALLELES:
    print "keep sites with minimum {} alleles".format(MINALLELES)
    filter_num_alleles(INFILE, MINALLELES, MAXALLELES, OUTFILE)
elif MAXALLELES:
    print "keep sites with maximum {} alleles".format(MAXALLELES)
    filter_num_alleles(INFILE, MINALLELES, MAXALLELES, OUTFILE)
elif INFO:
    print "filter sites by INFO keys"
    if -1 != INFO.find('>='):
        tmp = INFO.split('>=')
        key = tmp[0]
        operation = '>='
        value = [float(tmp[1])]
        print "keep sites with {} {} {}".format(key, operation, value)
        filter_by_info(INFILE, OUTFILE,key,operation,value,NA)
    elif -1 != INFO.find('<='):
        tmp = INFO.split('<=')
        key = tmp[0]
        operation = '<='
        value = [float(tmp[1])]
        print "keep sites with {} {} {}".format(key, operation, value)
        filter_by_info(INFILE, OUTFILE,key,operation,value,NA)
    elif -1 != INFO.find('<'):
        tmp = INFO.split('<')
        key = tmp[0]
        operation = '<'
        value = [float(tmp[1])]
        print "keep sites with {} {} {}".format(key, operation, value)
        filter_by_info(INFILE, OUTFILE,key,operation,value,NA)
    elif -1 != INFO.find('>'):
        tmp = INFO.split('>')
        key = tmp[0]
        operation = '>'
        value = [float(tmp[1])]
        print "keep sites with {} {} {}".format(key, operation, value)
        filter_by_info(INFILE, OUTFILE,key,operation,value,NA)
    elif -1 != INFO.find('!='):
        tmp = INFO.split('!=')
        key = tmp[0]
        operation = '!='
        value = check_float(tmp[1])
        values = None
        if value:
            stype = 'float'
            values = [value]
        else:
            stype = 'str'
            values = tmp[1].split(',')
        print "keep sites with {} {} {}".format(key, operation, values)
        filter_by_info(INFILE, OUTFILE,key,operation,values,NA,stype)
    elif -1 != INFO.find('='):
        tmp = INFO.split('=')
        key = tmp[0]
        operation = '='
        value = check_float(tmp[1])
        values = None
        if value:
            stype = 'float'
            values = [value]
        else:
            stype = 'str'
            values = tmp[1].split(',')
        print "keep sites with {} {} {}".format(key, operation, values)
        filter_by_info(INFILE, OUTFILE,key,operation,values,NA,stype)
else:
    sys.exit('do nothing!')
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
