#!/usr/bin/python

'''
Phase the PC vcf file into variants septerating at different taxonomic levels.
'''

import sys,argparse
from collections import defaultdict
from sets import Set
import numpy as np



#-----------------------------------------------------
# Step 0
# Define classes and functions
#-----------------------------------------------------


class SnpObj(object):
    """
    An object holding information on SNPs
    """
    def __init__(self):
        """set the object, defining the dictionary."""
        dict = defaultdict(list)
        self.conversion_dict = dict
    def set_dict(self, lines):
        """ Create the conversion dictionary from input lines """
        for line in lines:
            line = line.rstrip()
            split_line = line.split("\t")
            old_gene_id = split_line[0]
            new_gene_id = split_line[2]
            conv_dict = self.conversion_dict
            conv_dict[old_gene_id] = new_gene_id
            self.conversion_dict = conv_dict


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument('--inp_vcf',required=True,type=str,help='input vcf file')
conf = ap.parse_args()


with open(conf.inp_vcf) as f:
    vcf_lines = f.readlines()


#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------

# Isolates accoiated with crown rot, leather rot, apple and P. idaei
def phase_snp_line(line):
    # cr = ['414', '12420', '15_13', '15_7', '2003_3', '4032', '4040', '404', '415', '416', 'P421', 'PC13_15', '10300']
    # lr = ['11-40', '17-21']
    # md = ['62471', 'P295', 'R36_14']
    # pi = ['371', 'SCRP370', 'SCRP376']
    isolates = ['11-40', '12420', '15_13', '15_7', '17-21', '2003_3', '371', '4032', '404', '4040', '414', '415', '416', '62471', 'P295', 'P421', 'PC13_15', 'R36_14', 'SCRP370', 'SCRP376']
    taxonomies = ['lr', 'cr', 'cr', 'cr', 'lr', 'cr', 'pi', 'cr', 'cr', 'cr', 'cr', 'cr', 'cr', 'md', 'md', 'cr', 'cr', 'md', 'pi', 'pi']
    summarise_genotype = {
        '0/0' : 'ref',
        '0/1' : 'var_het',
        '1/1' : 'var_homo',
        '0/2' : 'polymophic',
        '1/2' : 'polymophic',
        '2/2' : 'polymophic'
    }
    # phase_dict = defaultdict(set)
    # summary_set = set()
    # split_line = line.split("\t")
    # for isolate, taxonomy, element in zip(isolates, taxonomies, split_line[9:]):
    #     # print "\t".join([isolate, taxonomy, element])
    #     genotype = element.split(":")[0]
    #     summary_set.add(":".join([taxonomy, summarise_genotype[genotype]]))
    #     phase_dict[taxonomy].add(summarise_genotype[genotype])
    # # print(phase_dict)
    # if 'var' in str(phase_dict['cr']) and 'var' not in str(phase_dict['md']) and 'var' not in str(phase_dict['pi']):
    #     print 'crown rot private'
    # elif 'var' not in str(phase_dict['cr']) and 'var' in str(phase_dict['md']) and 'var' not in str(phase_dict['pi']):
    #     print 'apple private'
    # elif 'var' not in str(phase_dict['cr']) and 'var' not in str(phase_dict['md']) and 'var' in str(phase_dict['pi']):
    #     print 'P. idaei private'
    # elif 'var' not in str(phase_dict['cr']) and 'var' not in str(phase_dict['md']) and 'var' not in str(phase_dict['pi']) and 'var' in str(phase_dict['lr']):
    #     print 'leather rot SNP'
    # elif 'var' not in str(phase_dict['cr']) and 'var_homo' in str(phase_dict['md']) and 'var_homo' in str(phase_dict['pi']):
    #     print 'crown rot private'
    # elif 'var_het' in str(phase_dict['cr']) and 'var' not in str(phase_dict['md']) and 'var_homo' in str(phase_dict['pi']):
    #     print 'cactorum private, apple fixed, crown rot unfixed'
    # elif 'var' in str(phase_dict['cr']) and 'var_homo' in str(phase_dict['md']) and 'var_homo' in str(phase_dict['pi']):
    #     print 'crown rot private, unfixed'
    # else:
    #     print(phase_dict)
    #     exit()
    phase_dict = defaultdict(set)
    summary_set = set()
    split_line = line.split("\t")
    for isolate, taxonomy, element in zip(isolates, taxonomies, split_line[9:]):
        # print "\t".join([isolate, taxonomy, element])
        genotype = element.split(":")[0]
        # summary_set.add(":".join([taxonomy, summarise_genotype[genotype]]))
        phase_dict[taxonomy].add(genotype)
    heterozygocity_dict = defaultdict()
    # print phase_dict.keys()
    status_list = []
    # for taxon in phase_dict.keys():
    for taxon in ['cr', 'md', 'pi']:
        if any(x in phase_dict[taxon] for x in ['0/2', '1/2', '2/2']):
            status = 'triallelic SNP'
        elif len(phase_dict[taxon]) > 1:
            status = 'var_unfixed'
        elif all(x == '0/0' for x in phase_dict[taxon]):
            # heterozygocity_dict[taxon] = 'ref'
            status = 'ref'
        elif any(x == '0/1' for x in phase_dict[taxon]):
            # heterozygocity_dict[taxon] = 'var_unfixed'
            status = 'var_unfixed'
        elif all(x == '1/1' for x in phase_dict[taxon]):
            # heterozygocity_dict[taxon] = 'var_fixed'
            status = 'var_fixed'
        else:
            # heterozygocity_dict[taxon] = 'monkeys'
            # status = 'monkeys'
            print(taxon)
            print(phase_dict[taxon])
            exit()
        status_list.append("=".join([taxon, status]))


        het_to_phase = {
        'cr=var_fixed,md=ref,pi=ref' : 'impossible - crown rot private fixed',
        'cr=var_unfixed,md=ref,pi=ref' : 'crown rot private unfixed',
        'cr=ref,md=var_fixed,pi=ref' : 'apple private fixed',
        'cr=ref,md=var_unfixed,pi=ref' : 'apple private unfixed',

        'cr=ref,md=ref,pi=var_fixed' : 'species private fixed',
        'cr=ref,md=ref,pi=var_unfixed' : 'P. idaei private unfixed',


        'cr=var_fixed,md=var_fixed,pi=ref' : 'impossible - species private fixed',
        'cr=var_unfixed,md=var_fixed,pi=ref' : 'P. cactorum private unfixed',
        'cr=var_fixed,md=var_unfixed,pi=ref' : 'impossible - P. cactorum private unfixed',
        'cr=var_unfixed,md=var_unfixed,pi=ref' : 'P. cactorum private unfixed',

        'cr=ref,md=var_fixed,pi=var_fixed' : 'crown rot private fixed',
        'cr=ref,md=var_unfixed,pi=var_fixed' : 'P. cactorum private crown rot fixed',
        'cr=ref,md=var_fixed,pi=var_unfixed' : 'ancestral variation differentially fixed',
        'cr=ref,md=var_unfixed,pi=var_unfixed' : 'ancestral variation crown rot fixed',

        'cr=var_fixed,md=ref,pi=var_fixed' : 'impossible',
        'cr=var_unfixed,md=ref,pi=var_fixed' : 'P. cactorum private apple fixed',
        'cr=var_fixed,md=ref,pi=var_unfixed' : 'impossible',
        'cr=var_unfixed,md=ref,pi=var_unfixed' : 'ancestral variation apple fixed',

        'cr=var_unfixed,md=var_fixed,pi=var_fixed' : 'crown rot private unfixed',
        'cr=var_unfixed,md=var_fixed,pi=var_unfixed' : 'ancestral variation apple fixed',
        'cr=var_unfixed,md=var_unfixed,pi=var_fixed' : 'P. cactorum private unfixed',

        'cr=var_unfixed,md=var_unfixed,pi=var_unfixed' : 'ancestral variation unfixed',

        'cr=ref,md=ref,pi=ref' : 'leather rot variant'
        }
    if "triallelic SNP" in str(status_list):
        phase = "triallelic SNP"
    else:
        # print(",".join(status_list))
        phase = het_to_phase[",".join(status_list)]
    print "\t".join([phase, line])

    # if 'var' in str(heterozygocity_dict['cr']) and 'var' not in str(heterozygocity_dict['md']) and 'var' not in str(heterozygocity_dict['pi']):
    #     print 'crown rot private'
    # elif 'var' not in str(heterozygocity_dict['cr']) and 'var' in str(heterozygocity_dict['md']) and 'var' not in str(heterozygocity_dict['pi']):
    #     print 'apple private'
    # elif 'var' not in str(heterozygocity_dict['cr']) and 'var' not in str(heterozygocity_dict['md']) and 'var' in str(heterozygocity_dict['pi']):
    #     print 'P. idaei private'
    # elif 'var' not in str(heterozygocity_dict['cr']) and 'var' not in str(heterozygocity_dict['md']) and 'var' not in str(heterozygocity_dict['pi']) and 'var' in str(heterozygocity_dict['lr']):
    #     print 'leather rot SNP'
    # elif 'var' not in str(heterozygocity_dict['cr']) and 'var_homo' in str(heterozygocity_dict['md']) and 'var_homo' in str(heterozygocity_dict['pi']):
    #     print 'crown rot private'
    # elif 'var_het' in str(heterozygocity_dict['cr']) and 'var' not in str(heterozygocity_dict['md']) and 'var_homo' in str(heterozygocity_dict['pi']):
    #     print 'cactorum private, apple fixed, crown rot unfixed'
    # elif 'var' in str(heterozygocity_dict['cr']) and 'var_homo' in str(heterozygocity_dict['md']) and 'var_homo' in str(heterozygocity_dict['pi']):
    #     print 'crown rot private, unfixed'
    # else:
    #     print(heterozygocity_dict)
        # exit()






for line in vcf_lines:
    line = line.rstrip()
    if line.startswith('#'):
        continue
    # print line
    phase_snp_line(line)
    # exit()
