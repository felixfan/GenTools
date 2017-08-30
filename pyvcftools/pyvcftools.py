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
        if len(CHR)<3:
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

def filter_by_chr(infile, chrs,outfile, exclude=False):
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
            if not exclude:
                if arr[0] in chrs:
                    fw.write("%s\n" % r)
                    m += 1
            else:
                if not arr[0] in chrs:
                    fw.write("%s\n" % r)
                    m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()

def check_region(chrom, pos, regions):
    '''
    check whether the chrom:pos is in regions
    regions is a list, each element of it is also a list: [chr, start, end]
    '''
    flag = False
    for r in regions:
        if chrom == r[0]:
            if r[2] >= pos >= r[1]:
                flag = True
    return flag

def filter_by_region(infile, regions, outfile, exclude=False):
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
            if not exclude:
                if check_region(arr[0], int(arr[1]), regions):
                    fw.write("%s\n" % r)
                    m += 1
            else:
                if not check_region(arr[0], int(arr[1]), regions):
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
            elif gtp =="not-hom-alt":
                if 0 < s.count('0'):
                    return True
                else:
                    if s[0] != s[2]:
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
            print 'individual ids in vcf are:', myind[9:]
            for z in inds:
                if not z in myind[9:]:
                    sys.exit('Error: {} is not in the input vcf file!'.format(z))
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

def filter_by_var_ids(infile, ids, outfile, exclude=False):
    '''
    filter variants in ids
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
            if not exclude:
                if arr[2] in ids:
                    fw.write("%s\n" % r)
                    m += 1
            else:
                if not arr[2] in ids:
                    fw.write("%s\n" % r)
                    m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()

def filter_by_phypos(infile, pp, outfile, exclude=False):
    '''
    filter variants by physical positions
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
            if not exclude:
                if tmp in pp:
                    fw.write("%s\n" % r)
                    m += 1
            else:
                if not tmp in pp:
                    fw.write("%s\n" % r)
                    m += 1
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
                sys.exit('Error: %s is not in %s' % (inds[0], infile))
            for x in inds[1:]:
                if x in myind:
                    tmp = myind.index(x)
                    idx.append(tmp)
                else:
                    sys.exit('Error: %s is not in %s' % (x, infile))
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
                        elif operation == "==":
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

def filter_by_format(infile, outfile, key, operation, values, stype='float'):
    '''
    filter by format
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    idx = -1
    flag = False
    for r in fr:
        r = r.strip()
        if r.startswith("#"):
            fw.write("%s\n" % r)
        else:
            arr = r.split()
            n += 1
            if idx == -1:
                ar = arr[8].split(':')
                for z in xrange(len(ar)):
                    if ar[z] == key:
                        idx = z
                        flag = True
                if not flag:
                    sys.exit('{} is not find in format field of the input vcf file'.format(key))
            # write output
            fw.write('{}'.format(arr[0]))
            for z in arr[1:9]:
                fw.write('\t{}'.format(z))
            for z in arr[9:]:
                if z == '.':
                    fw.write('\t{}'.format('.'))
                else:
                    fla = True
                    ag = z.split(':')
                    if ag[idx]  == '.':
                        fla = False
                    else:
                        bs = ag[idx].split(',') # multiple value
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
                        fw.write("\t{}".format(z))
                    else:
                        fw.write('\t{}'.format('.'))
            fw.write('\n')
    print "%d variants were written to %s" % (n, outfile)
    fw.close()
    fr.close()

def extract_genotype(infile, pp, outfile):
    '''
    extract genotype by chr:pos
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    m = 0
    for r in fr:
        r = r.strip()
        if not r.startswith("#"):
            arr = r.split()
            tmp = "{}:{}".format(arr[0],arr[1])
            # print tmp
            if tmp in pp:
                fw.write("{}\t{}\t{}\t{}".format(arr[0], arr[1], arr[3], arr[4]))
                for a in arr[9:]:
                    if a != '.':
                        fw.write('\t{}'.format(a[:3]))
                    else:
                        fw.write('\t.')
                fw.write('\n')
                m += 1
    print "genotype of {} variants were written to {}".format(m, outfile)
    fw.close()
    fr.close()

def filter_comp_het(infile, inds, geneKey, funcKey, funcValues, outfile):
    '''
    filter by compound heterozygous
    '''
    fr = open(infile)
    fp = open(outfile + '.compHetPairs.txt', 'w')
    n = 0 # total sites
    indidx = []
    genes = {}
    gene = ''
    pp = []
    for r in fr:
        r = r.strip()
        if r.startswith("##"):
            pass
        elif r.startswith("#"):
            myind = r.split()
            print "there are {} individuals in {}\nthey are: {}".format(len(myind[9:]), infile, myind[9:])
            for x in inds:
                if x in myind:
                    idx = myind.index(x)
                    indidx.append(idx)
                else:
                    sys.exit('Error: {} is not in {}'.format(x, infile))
        else:
            arr = r.split()
            n += 1
            ar = arr[7].split(';')
            flag0 = False # whether value is in funcValues
            flag1 = False # whether func key is in info and
            flag2 = False # whether gene key is in info
            for kv in ar: # check keys in info and value in funcValues
                if kv.startswith(funcKey):
                    flag1 = True
                    a = kv.split('=')
                    if a[1] in funcValues:
                        flag0 = True
                if kv.startswith(geneKey):
                    gene = kv.split('=')[1]
                    flag2 = True
                if flag0 and flag1 and flag2:
                    break
            if flag0 and flag1 and flag2: # in the specified region of a gene
                if not gene in genes:
                    genes[gene] = []
                one = [arr[0], arr[1]] # chr & pos, record het status of other inds
                flag3 = True           # whether a var is het in all comp-het inds
                for gidx in xrange(9,len(arr)): # check ind gtp
                    if arr[gidx] != '.': # non-missing
                        gtp = arr[gidx][:3]
                        c = gtp.count('0')
                        if c == 1: # het
                            if not gidx in indidx: # non-comp-het inds
                                one.append(1)
                        else:      # non - het
                            if gidx in indidx: # comp-het inds
                                flag3 = False
                                break
                            else: # non-comp-het inds
                                one.append(0)
                    else: # missing vale
                        flag3 = False
                        break
                if flag3:
                    genes[gene].append(one)
            else:
                if not flag1:
                    sys.exit('Error: can not find "{}" in INFO field of the input vcf'.format(funcKey))
                elif not flag2:
                     sys.exit('Error: can not find "{}" in INFO field of the input vcf'.format(geneKey))
    fr.close()
    # check comp-het
    chg = [] # total genes have comp het
    chv = 0 # total number of comp het var pairs
    for k, v in genes.items():
        if len(v) > 1:
            flag5 = False # whether this gene has any comp-het pairs
            for i in xrange(len(v)):
                for j in xrange(len(v)):
                    if j > i:
                        flag4 = True # whether all non-comp-het only have zero or one het in one pair vars
                        for z in xrange(2, len(v[i])):
                            if v[i][z] + v[j][z] > 1:
                                flag4 = False
                                break
                        if flag4:
                            chrom = v[i][:2][0]
                            p1 = v[i][:2][1]
                            p2 = v[j][:2][1]
                            fp.write("{}\t{}\t{}\t{}\n".format(k, chrom, p1, p2))
                            pp.append("{}:{}".format(chrom, p1))
                            pp.append("{}:{}".format(chrom, p2))
                            chv += 1
                            flag5 = True
            if flag5:
                chg.append(k)
    pp = list(set(pp)) # all uniq var (chr:pos) in comp-het pairs
    fp.close()
    print "there are total {} variants in {}".format(n, infile)
    print "\nnumber of genes: {}".format(len(chg))
    print "those genes are: {}".format(chg)
    print "number of variants: {}".format(len(pp))
    print "number of compound heterozygous: {}".format(chv)
    print "write compound heterozygous to {}.compHetPairs.txt".format(outfile)
    print "write variants to {}".format(outfile)
    filter_by_phypos(infile, pp, outfile)
    print "write genotypes of these variants to {}.gtp.txt".format(outfile)
    extract_genotype(infile, pp, outfile+'.gtp.txt')

def filter_by_indid(infile, indids,outfile, exclude=False):
    '''
    filter by individual ID
    '''
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0 # total sites
    m = 0
    pidx = [0,1,2,3,4,5,6,7,8] # first 9 cols
    inds = []
    for r in fr:
        r = r.strip()
        if r.startswith("##"):
            fw.write("%s\n" % r)
        elif r.startswith("#"):
            myind = r.split()
            print "there are {} individuals in {}\nthey are: {}".format(len(myind[9:]), infile, myind[9:])
            if not exclude:
                for x in indids:
                    if x in myind:
                        idx = myind.index(x)
                        pidx.append(idx)
                        m += 1
                        inds.append(x)
                    else:
                        sys.exit('Error: {} is not in {}'.format(x, infile))
            else:
                for x in myind[9:]:
                    if not x in indids:
                        idx = myind.index(x)
                        pidx.append(idx)
                        m += 1
                        inds.append(x)
            if m > 0:
                print 'keep these individuals: {}'.format(inds)
            else:
                sys.exit('No individual left')
            fw.write(myind[0])
            for i in pidx[1:]:
                fw.write('\t{}'.format(myind[i]))
            fw.write('\n')
        else:
            arr = r.split()
            n += 1
            fw.write(arr[0])
            for i in pidx[1:]:
                fw.write('\t{}'.format(arr[i]))
            fw.write('\n')
    print "{} variants of {} individuals were written to {}".format(n, m, outfile)
    fw.close()
    fr.close()

def filter_missing(infile, cutoff, outfile, count=False):
    '''
    filter by missing rate / count
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
            nind = len(arr) - 9
            mind = 0
            for z in arr[9:]:
                if z == '.':
                    mind += 1
            if not count:
                if 1.0*mind/nind <= cutoff:
                    fw.write("%s\n" % r)
                    m += 1
            else:
                if mind <= cutoff:
                    fw.write("%s\n" % r)
                    m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()

def filter_parental_het(infile, inds, outfile):
    fr = open(infile)
    fw = open(outfile, 'w')
    n = 0
    m = 0
    idx = []
    myind = None
    for r in fr:
        r = r.strip()
        if r.startswith("##"):
            fw.write("%s\n" % r)
        elif r.startswith("#"):
            fw.write("%s\n" % r)
            myind = r.split()
            for x in inds:
                if x in myind:
                    tmp = myind.index(x)
                    idx.append(tmp)
                else:
                    sys.exit('Error: %s is not in %s' % (x, infile))
        else:
            arr = r.split()
            n += 1
            flag = False # no het
            for i in range(9,len(arr)):
                if myind[i] in inds:
                    flag = (flag or check_allele(arr[i],'het','rm'))
                if flag == True:
                    break
            if flag:
                fw.write("%s\n" % r)
                m += 1
    print "%d of %d variants were written to %s" % (m, n, outfile)
    fw.close()
    fr.close()

def check_float(s):
    try:
        return float(s)
    except ValueError:
        return False

def check_info_format_command(operator, value):
    '''
    help function
    '''
    values = []
    stype = 'float'
    if operator == '>=' or operator == '>':
        values = [max(map(float, value))]
    elif operator == '<=' or operator == '<':
        values = [min(map(float, value))]
    elif operator == '!=' or operator == '==':
        stype = 'str'
        values = value
    return (values, stype)

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

    desc = '''Filtering of polymorphisms according to genotypes,
    physical positions and thresholds of quality score,
    read depth and others. Basic statistics and manipulation.'''

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', action='version', version='%(prog)s 1.0.0')
    ### input
    parser.add_argument('--vcf', help='input vcf file', required=True, type=str)
    parser.add_argument('--vcf2', help='the second input vcf file', type=str)
    ### filters
    parser.add_argument('--chr', help='filter by chromosome', type=str, action='append') # multiple
    parser.add_argument('--region', help='filter by region', type=str, nargs=3) # 3 arguments
    parser.add_argument('--region-file', help='filter by region', type=str)
    parser.add_argument('--qual', help='filter by qual score', type=float)
    parser.add_argument('--filter', help='filter by filter flag', type=str, action='append')
    parser.add_argument('--genotype', help='filter by genotype', type=str,choices=["hom-ref", "hom-alt","het", "het-alt","not-hom-ref","not-two-alt","two-alt","not-het","not-hom-alt"])
    parser.add_argument('--keep-only-indels', help='keep only indels', action='store_true')
    parser.add_argument('--remove-indels', help='remove indels', action='store_true')
    parser.add_argument('--id', help='filter by ID', type=str, action='append')
    parser.add_argument('--id-file', help='filter by ID', type=str)
    parser.add_argument('--phy-pos', help='filter by physical position', type=str, nargs=2)
    parser.add_argument('--phy-pos-file', help='filter by physical position', type=str)
    parser.add_argument('--cmp-gtp-same', help='compare genotype of multiple individuals', action='store_true')
    parser.add_argument('--cmp-gtp-diff', help='compare genotype of multiple individuals', action='store_true')
    parser.add_argument('--min-alleles', help='minimum number of alleles', type=int)
    parser.add_argument('--max-alleles', help='maximum number of alleles', type=int)
    parser.add_argument('--info', help='filter by info keys', action='store_true')
    parser.add_argument('--key', help='key in info field', type=str)
    parser.add_argument('--value', help='value of key', type=str, action='append')
    parser.add_argument('--value-file', help='values of key', type=str)
    parser.add_argument('--operator', help='operator', type=str, choices=['>', '>=', '<', '<=', '==', '!='])
    parser.add_argument('--comp-het', help='filter by compound heterozygous', action='store_true')
    parser.add_argument('--gene-key', help='key for gene annotation in INFO field', type = str)
    parser.add_argument('--func-key', help='key for function annotation in INFO field', type = str)
    parser.add_argument('--func-values', help='value for function annotation in INFO field', type = str, action='append')
    parser.add_argument('--keep-inds', help='filter by individual id', action='store_true')
    parser.add_argument('--remove-inds', help='filter by individual id', action='store_true')
    parser.add_argument('--keep-inds-file', help='filter by individual id', type = str)
    parser.add_argument('--remove-inds-file', help='filter by individual id', type = str)
    parser.add_argument('--missing-rate', help='filter by missing rate', type=str)
    parser.add_argument('--missing-count', help='filter by missing count', type=str)
    parser.add_argument('--format', help='filter by format keys', type=str)
    parser.add_argument('--parental-het', help='at least one parent is heterozygous', action='store_true')
    ### statistics
    parser.add_argument('--het', help='calculates heterozygosity on a per-individual basis', action='store_true')
    parser.add_argument('--titv', help='calculates Ti/Tv ratio on a per-individual basis', action='store_true')
    parser.add_argument('--stat-chr', help='statistics based field in VCF files', action='store_true')
    parser.add_argument('--stat-filter', help='statistics based field in VCF files', action='store_true')
    parser.add_argument('--stat-info', help='statistics based on info field', action='store_true')
    ### manipulation
    parser.add_argument('--concat', help='concatenates VCF files', action='store_true')
    ### individual
    parser.add_argument('--ind', help='individual id', type=str, action='append')
    ### missing value
    parser.add_argument('--missing-value', help='how to deal with missing values', default="keep", type=str, choices=["keep", "rm"])
    ### reverse command
    parser.add_argument('--reverse', help='reverse the filter', action='store_true')
    ### output
    parser.add_argument('--out', help='output vcf file', type=str)
    ### param
    args = vars(parser.parse_args())
    INFILE = args['vcf'] if 'vcf' in args else None
    OUTFILE = args['out']
    NA = args['missing_value']
    REVERSE = args['reverse'] if 'reverse' in args else False
    CHR = args['chr'] if 'chr' in args else None
    REGION = args['region'] if 'region' in args else None
    REGIONFILE = args['region_file'] if 'region_file' in args else None
    QUAL = args['qual'] if 'qual' in args else None
    FILTER = args['filter'] if 'filter' in args else None
    GENOTYPE = args['genotype'] if 'genotype' in args else None
    IND = args['ind'] if 'ind' in args else None
    INDEL = args['keep_only_indels'] if 'keep_only_indels' in args else False
    SNP = args['remove_indels'] if 'remove_indels' in args else False
    ID = args['id'] if 'id' in args else None
    IDFILE = args['id_file'] if 'id_file' in args else None
    PHYPOS = args['phy_pos'] if 'phy_pos' in args else None
    PHYPOSFILE = args['phy_pos_file'] if 'phy_pos_file' in args else None
    CMPGTPSAME = args['cmp_gtp_same'] if 'cmp_gtp_same' in args else None
    CMPGTPDIFF = args['cmp_gtp_diff'] if 'cmp_gtp_diff' in args else None
    MINALLELES = args['min_alleles'] if 'min_alleles' in args else None
    MAXALLELES = args['max_alleles'] if 'max_alleles' in args else None
    INFO = args['info'] if 'info' in args else False
    KEY = args['key'] if 'key' in args else None
    VALUE = args['value'] if 'value' in args else None
    VALUEFILE = args['value_file'] if 'value_file' in args else None
    OPERATOR = args['operator'] if 'operator' in args else None
    COMPHET = args['comp_het'] if 'comp_het' in args else None
    GENEKEY = args['gene_key'] if 'gene_key' in args else None
    FUNCKEY = args['func_key'] if 'func_key' in args else None
    FUNCVALUES = args['func_values'] if 'func_values' in args else None
    KEEPINDS = args['keep_inds'] if 'keep_inds' in args else False
    REMOVEINDS = args['remove_inds'] if 'remove_inds' in args else False
    KEEPINDSFILE = args['keep_inds_file'] if 'keep_inds' in args else None
    REMOVEINDSFILE = args['remove_inds_file'] if 'remove_inds' in args else None
    MISSRATE= args['missing_rate'] if 'missing_rate' in args else None
    MISSCOUNT = args['missing_count'] if 'missing_count' in args else None
    FORMAT = args['format'] if 'format' in args else None
    PARENTALHET = args['parental_het'] if 'parental_het' in args else None
    STATINFO = args['stat_info'] if 'stat_info' in args else False
    CONCAT = args['concat'] if 'concat' in args else False
    INFILE2 = args['vcf2'] if 'vcf2' in args else None
    HET = args['het'] if 'het' in args else False
    STATCHR = args['stat_chr'] if 'stat_chr' in args else None
    STATFILTER = args['stat_filter'] if 'stat_filter' in args else None
    TITV = args['titv'] if 'titv' in args else False
    ###log
    print "@-------------------------------------------------------------@"
    print "|        pyVCFTools     |      v1.0.1      |    18 Jul 2017   |"
    print "|-------------------------------------------------------------|"
    print "|  (C) 2017 Felix Yanhui Fan, GNU General Public License, v2  |"
    print "|-------------------------------------------------------------|"
    print "|    For documentation, citation & bug-report instructions:   |"
    print "|            http://felixfan.github.io/GenTools               |"
    print "@-------------------------------------------------------------@"
    print "\n\tOptions in effect:"
    print "\t--vcf", INFILE
    if CHR:
        print "\t--chr", CHR
        if REVERSE:
            print "\t--reverse"
    elif REGION:
        print "\t--region", REGION
        if REVERSE:
            print "\t--reverse"
    elif REGIONFILE:
        print "\t--region-file", REGIONFILE
        if REVERSE:
            print "\t--reverse"
    elif QUAL:
        print "\t--qual", QUAL
    elif FILTER:
        print "\t--filter", FILTER
    elif GENOTYPE:
        print "\t--genotype", GENOTYPE
        if IND:
            print "\t--ind", IND
        else:
            sys.exit("Error, argment --ind is missing!")
        print "\t--missing-value", NA
    elif INDEL:
        print "\t--keep-only-indels"
    elif SNP:
        print "\t--remove-indels"
    elif ID:
        print "\t--ids", ID
        if REVERSE:
            print "\t--reverse"
    elif IDFILE:
        print "\t--ids-file", IDFILE
        if REVERSE:
            print "\t--reverse"
    elif PHYPOS:
        print "\t--phy-pos", PHYPOS
        if REVERSE:
            print "\t--reverse"
    elif PHYPOSFILE:
        print "\t--phy-pos-file", PHYPOSFILE
        if REVERSE:
            print "\t--reverse"
    elif CMPGTPSAME:
        print "\t--cmp-gtp-same"
        if IND:
            print "\t--ind", IND
        else:
            sys.exit("Error, argment --ind is missing!")
        print "\t--missing-value", NA
    elif CMPGTPDIFF:
        print "\t--cmp-gtp-diff"
        if IND:
            print "\t--ind", IND
        else:
            sys.exit("Error, argment --ind is missing!")
        print "\t--missing-value", NA
    elif MINALLELES and MAXALLELES:
        print "\t--min-alleles", MINALLELES
        print "\t--max-alleles", MAXALLELES
    elif MINALLELES:
        print "\t--min-alleles", MINALLELES
    elif MAXALLELES:
        print "\t--max-alleles", MAXALLELES
    elif INFO:
        print "\t--info"
        if KEY and OPERATOR and (VALUE or VALUEFILE):
            print "\t--key", KEY
            print "\t--operaror", OPERATOR
            if VALUEFILE:
                print "\t--value-file", VALUEFILE
            else:
                print "\t--value", VALUE
        else:
            sys.exit('--key, --operator and --value or --value-file are needed!')
        print "\t--missing-value", NA
    elif COMPHET:
        print "\t--comp-het"
        if IND:
            print "\t--ind", IND
        else:
            sys.exit("Error, argment --ind is missing!")
        if GENEKEY:
            print "\t--gene-key", GENEKEY
        else:
            sys.exit("Error, argment --gene-key is missing!")
        if FUNCKEY:
            print "\t--func-key", FUNCKEY
        else:
            sys.exit("Error, argment --func-key is missing!")
        if FUNCVALUES:
            print "\t--func-values", FUNCVALUES
        else:
            sys.exit("Error, argment --func-values is missing!")
    elif KEEPINDS:
        print "\t--keep-inds"
        if IND:
            print "\t--ind", IND
        else:
            sys.exit("Error, argment --ind is missing!")
    elif KEEPINDSFILE:
        print "\t--keep-inds-file"
    elif REMOVEINDS:
        print "\t--remove-inds"
        if IND:
            print "\t--ind", IND
        else:
            sys.exit("Error, argment --ind is missing!")
    elif REMOVEINDSFILE:
        print "\t--remove-inds-file"
    elif MISSRATE:
        print "\t--missing-rate", MISSRATE
    elif MISSCOUNT:
        print "\t--missing-count", MISSCOUNT
    elif FORMAT:
        print "\t--format", FORMAT
    elif PARENTALHET:
        print "\t--parental-het"
        if IND:
            print "\t--ind", IND
        else:
            sys.exit("Error, argment --ind is missing!")
    elif STATINFO:
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
    if CHR:
        chrs = []
        for c in CHR:
            if -1 == c.find('-'):
                chrs.append(c)
            else:
                cs = split_str_comma_dash(c)
                chrs.extend(cs)
        if REVERSE:
            print "exclude variants on chromosomes:", chrs
            filter_by_chr(INFILE,chrs, OUTFILE, True)
        else:
            print "keep variants on chromosomes:", chrs
            filter_by_chr(INFILE,chrs, OUTFILE)
    elif REGION:
        chrom = REGION[0]
        start = int(REGION[1])
        end = int(REGION[2])
        regions = [[chrom,start,end]]
        if not REVERSE:
            print "keep variants in region: from {} to {} on chromosome {}".format(start, end, chrom)
            filter_by_region(INFILE,regions,OUTFILE)
        else:
            print "exclude variants in region: from {} to {} on chromosome {}".format(start, end, chrom)
            filter_by_region(INFILE,regions,OUTFILE, True)
    elif REGIONFILE:
        regions = []
        tf = open(REGIONFILE)
        for r in tf:
            r = r.strip()
            arr = r.split()
            regions.append([arr[0], int(arr[1]), int(arr[2])])
        tf.close()
        if not REVERSE:
            print "keep variants in regions: {}".format(regions)
            filter_by_region(INFILE,regions,OUTFILE)
        else:
            print "exclude variants in regions: {}".format(regions)
            filter_by_region(INFILE,regions,OUTFILE, True)
    elif QUAL:
        cutoff = float(QUAL)
        print "keep variants with qual score no less than %f" % cutoff
        filter_by_qual(INFILE, cutoff, OUTFILE)
    elif FILTER:
        print "keep variants with FILTER flag:", FILTER
        filter_by_filter(INFILE, FILTER, OUTFILE)
    elif GENOTYPE:
        print "genotype filter:"
        if GENOTYPE == "hom-ref":
            print "keep variants that are homozygous of reference allele in these individuals:"
        elif GENOTYPE == "het":
            print "keep variants that are heterozygous in these individuals:"
        elif GENOTYPE == "hom-alt":
            print "keep variants that have two same alternative allele in these individuals:"
        elif GENOTYPE == "not-hom-ref":
            print "keep variants that are not homozygous of reference allele in these individuals:"
        elif GENOTYPE == "not-two-alt":
            print "keep variants that do not have two alternative alleles in these individuals:"
        elif GENOTYPE == "het-alt":
            print "keep variants that have two different alternative alleles in these individuals:"
        elif GENOTYPE == "two-alt":
            print "keep variants that have two alternative alleles in these individuals:"
        elif GENOTYPE == "not-hom-alt":
            print "keep variants do not have two same alternative alleles in these individuals:"
        elif GENOTYPE == "not-het":
            print "keep variants are not heterozygous in these individuals:"
        print IND
        filter_by_genotype(INFILE, GENOTYPE, IND, NA, OUTFILE)
    elif INDEL:
        print "keep only indels"
        INDEL = not INDEL
        filter_indels(INFILE, INDEL, OUTFILE)
    elif SNP:
        print "remove indels"
        filter_indels(INFILE, SNP, OUTFILE)
    elif ID or IDFILE:
        ids = []
        if ID:
            ids = ID
        else:
            tf = open(IDFILE)
            for r in tf:
                r = r.strip()
                ids.append(r)
            tf.close()
        if not REVERSE:
            print "keep variants by ID"
            filter_by_var_ids(INFILE, ids, OUTFILE)
        else:
            print "exclude variants by ID"
            filter_by_var_ids(INFILE, ids, OUTFILE, True)
    elif PHYPOS or PHYPOSFILE:
        pp = []
        if PHYPOSFILE:
            tf = open(PHYPOSFILE)
            for r in tf:
                r = r.strip()
                arr = r.split()
                pp.append(arr[0]+":"+arr[1])
            tf.close()
        else:
            pp = [PHYPOS[0] + ":" + PHYPOS[1]]
        if not REVERSE:
            print "keep variants by physical position"
            filter_by_phypos(INFILE, pp, OUTFILE)
        else:
            print "exclude variants by physical position"
            filter_by_phypos(INFILE, pp, OUTFILE, True)
    elif CMPGTPSAME:
        print "compare genotype of multiple individuals"
        print "only variants have the same genotype across the following individuals will be kept"
        if len(IND) > 1:
            print "individuals:", IND
        else:
            sys.exit('at least two individuals should be provided')
        cmp_gtp(INFILE, IND, NA, OUTFILE,True)
    elif CMPGTPDIFF:
        print "compare genotype of multiple individuals"
        print "only variants have the different genotype between the first individual and others will be kept"
        if len(IND) > 1:
            print "individuals:", IND[0], 'vs.', IND[1:]
        else:
            sys.exit('at least two individuals should be provided')
        cmp_gtp(INFILE, IND, NA, OUTFILE,False)
    elif MINALLELES and MAXALLELES:
        if MINALLELES < 2:
            sys.exit('value for --min-alleles should > 1')
        if MAXALLELES < 2:
            sys.exit('value for --max-alleles should > 1')
        if MAXALLELES < MINALLELES:
            sys.exit('value for --max-alleles should >= value for --min-alleles')
        print "keep sites with minimum {} alleles and maximum {} alleles".format(MINALLELES, MAXALLELES)
        filter_num_alleles(INFILE, MINALLELES, MAXALLELES, OUTFILE)
    elif MINALLELES:
        if MINALLELES < 2:
            sys.exit('value for --min-alleles should > 1')
        print "keep sites with minimum {} alleles".format(MINALLELES)
        filter_num_alleles(INFILE, MINALLELES, MAXALLELES, OUTFILE)
    elif MAXALLELES:
        if MAXALLELES < 2:
            sys.exit('value for --max-alleles should > 1')
        print "keep sites with maximum {} alleles".format(MAXALLELES)
        filter_num_alleles(INFILE, MINALLELES, MAXALLELES, OUTFILE)
    elif INFO:
        print "filter sites by INFO keys"
        if VALUEFILE:
            stype = 'float'
            if OPERATOR in ('==', '!='):
                stype = 'str'
            values = []
            f = open(VALUEFILE)
            for r in f:
                r = r.strip()
                values.append(r)
            f.close()
            print "keep sites with {} {} values in {}".format(KEY, OPERATOR, VALUEFILE)
            filter_by_info(INFILE, OUTFILE,KEY,OPERATOR,values,NA,stype)
        else:
            values, stype = check_info_format_command(OPERATOR, VALUE)
            print "keep sites with {} {} {}".format(KEY, OPERATOR, values)
            filter_by_info(INFILE, OUTFILE,KEY,OPERATOR,values,NA,stype)
    elif COMPHET:
        print "filter by compound hetrozygous"
        print "find compound hetrozygous in these individuals: {}".format(IND)
        filter_comp_het(INFILE, IND, GENEKEY, FUNCKEY, FUNCVALUES, OUTFILE)
    elif KEEPINDS:
        filter_by_indid(INFILE, IND,OUTFILE)
    elif REMOVEINDS:
        filter_by_indid(INFILE, IND,OUTFILE, True)
    elif KEEPINDSFILE or REMOVEINDSFILE:
        MYIND = []
        tmp = KEEPINDSFILE if KEEPINDSFILE else REMOVEINDSFILE
        f = open(tmp)
        for r in f:
            r = r.strip()
            MYIND.append(r)
        f.close()
        if KEEPINDSFILE:
            filter_by_indid(INFILE, MYIND,OUTFILE)
        else:
            filter_by_indid(INFILE, MYIND,OUTFILE, True)
    elif MISSRATE:
        MISSRATE = float(MISSRATE)
        if 1 >= MISSRATE >= 0:
            print "exclude sites with proportion of missing data is larger than {}".format(MISSRATE)
            filter_missing(INFILE, MISSRATE, OUTFILE)
        else:
            sys.exit('missing rate should be between 0 and 1')
    elif MISSCOUNT:
        MISSCOUNT = int(MISSCOUNT)
        if MISSCOUNT > 0:
            print "exclude sites with more than {} missing genotypes over all ndividuals".format(MISSCOUNT)
            filter_missing(INFILE, MISSCOUNT, OUTFILE, True)
        else:
            sys.exit('missing count should be a positive int')
    elif FORMAT:
        print "filter sites by FORMAT keys"
        key, operation, values, stype = check_info_format_command(FORMAT)
        print "keep sites with {} {} {}".format(key, operation, values)
        filter_by_format(INFILE, OUTFILE,key,operation,values,stype)
    elif PARENTALHET:
        print "at least one parent is heterozygous"
        if len(IND) != 2:
            sys.exit('two and only two individual IDs are needed')
        filter_parental_het(INFILE, IND,OUTFILE)
    elif STATINFO:
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
        sys.exit('do nothing!')
    ###run time
    run_time(starttime)
