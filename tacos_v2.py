#!/usr/bin/env python
#-*- coding: UTF-8 -*-
from __future__ import division
import sys
import argparse
import operator
from tabulate import tabulate
import re
import plotext as plt
import pysam
import six

'''
Run directory:
/Users/francisco/Documents/trabajo_nyu/tvag_project/tstableri_pacbio/pysplicing_tstabs
python tacos_v2.py -f pacbio_hic_CHAR29_renamed.fasta -sj tvCHAR29_rnaseqSJ.out.tab -b tvCHAR29_rnaseqAligned.sortedByCoord.out.bam -o test -5m GTWBNNH -3m DBYHWNMHDYAG
'''

program = 'pysplice.py'

parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, 
    description="\n\n\tDe novo splicing assembler for Trichomonas species\n---------------------\nTested on python 3.8.12\n\n")
requiredNamed = parser.add_argument_group('Mandatory arguments')
requiredNamed.add_argument("-f", "--input", dest='input', required=True, type=str, help='Input fasta file (Genome reference)')
requiredNamed.add_argument("-sj", "--input_sj", dest='sj_star', required=True, type=str, help='SJ.out.tab from STAR mapping')
requiredNamed.add_argument("-b", "--input_bam", dest='bam_f', required=True, type=str, help='BAM file from STAR mapping')
requiredNamed.add_argument("-o", "--output", dest='output', required=True, type=str, help='Prefix: Prefix for output files')
requiredNamed.add_argument("-5m", "--5p-motif", dest='motif5p', required=True, type=str, help='String: Splicing motif at 5p (upper case nucleotides only, no spaces)')
requiredNamed.add_argument("-3m", "--3p-motif", dest='motif3p', required=True, type=str, help='String: Splicing motif at 3p (upper case nucleotides only, no spaces)')

args = parser.parse_args()

##########################################################################################
#Load fasta and plotext summary
sequence = ""
fasta_dir = {}
fas_len = {}
lens_f = []

with open(args.input, 'r') as f:
	for line in f:
		if line.startswith('>'):
			if sequence:
				fas_len[seq_id] = len(sequence)
				fasta_dir[seq_id] = sequence
				lens_f.append(len(sequence))
				sequence = ""
			seq_id = line.split()[0].rstrip()[1:]
		else:
			sequence = sequence + line.strip()
	fasta_dir[seq_id] = sequence
	lens_f.append(len(sequence))
	fas_len[seq_id] = len(sequence)

print ('\n\t----------------- Fasta Summary -----------------\n')
cum_list=[]
y = 0 

len_sorted = sorted(lens_f, reverse=True)
for x in range(0,len(len_sorted)):
	y+=len_sorted[x]
	cum_list.append(y)

plt.plot(cum_list)
plt.scatter(cum_list)
plt.figsize(60, 20) 
plt.grid(True)
plt.title("Cumulative length of the genome reference")
plt.xlabel("Contig_number")
plt.show()

##########################################################################################
#Reverse complement sequence, function

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def reverse_complement(seq):    
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

##########################################################################################
#Extract intron sequences

forward_frag = {}
reverse_frag = {} 
undef_frag = {}
intron_coverage = {}
sj_entries = 0
intron_lens = []

with open(args.sj_star, 'r') as sj_file:
	for line in sj_file:
		line = line.rstrip('\n').rstrip('\r')
		col = line.split('\t')
		sj_entries += 1
		intron_lens.append(int(col[2]) - int(col[1]))
		intron_coverage[col[0] + ':' + str(col[1]) + '-' + str(col[2])] = str(col[6]) + ',' + str(col[7])
		if int(col[3]) == 1:
			forward_frag[col[0] + ':' + str(col[1]) + '-' + str(col[2])] = (fasta_dir[col[0]][int(col[1]) - 1 : int(col[2])]).upper()
		elif int(col[3]) == 2:
			reverse_frag[col[0] + ':' + str(col[1]) + '-' + str(col[2])] = (reverse_complement(fasta_dir[col[0]][int(col[1]) - 1 : int(col[2])])).upper()
		else:
			undef_frag['undef_plus' + col[0] + ':' + str(col[1]) + '-' + str(col[2])] = (fasta_dir[col[0]][int(col[1]) - 1 : int(col[2])]).upper()
			undef_frag['undef_mins' + col[0] + ':' + str(col[1]) + '-' + str(col[2])] = (reverse_complement(fasta_dir[col[0]][int(col[1]) - 1 : int(col[2])])).upper()

print ('\n\t----------------- SJs, Summary -----------------\n')
print('SJs, total entries: ' + str(sj_entries))
print('SJ entries in forward: ' + str(len(forward_frag.keys())))
print('SJ entries in reverse: ' + str(len(reverse_frag.keys())))
print('SJ entries, undefined strand: ' + str(len(undef_frag.keys())))

print ('\n\t----------------- SJs, intron Summary -----------------\n')
print('Minimum intron length: ' + str(min(intron_lens)))
print('Maximum intron length: ' + str(max(intron_lens)))

plt.hist(intron_lens,bins=1000)
plt.figsize(60, 20)     #######para la version mas reciente de python "plotsize" debe reemplazarse por "plot_size"
plt.grid(True)
plt.title("Intron length distribution")
plt.xlabel("Length")
plt.show()
print('\n\n')

##########################################################################################
#Search motifs at intron coordinates

#make regular expressions from nucleotide sequences
def make_degenerate_regex(motif_seq, molecule='dna'):
	if not isinstance( motif_seq, six.string_types ):
		raise TypeError( "motif_seq must be a string!" )
	if molecule == 'dna':
		degenerate_code = { "A":"A", "B":"[CGT]", "C":"C", "D":"[AGT]", "G":"G", "H":"[ACT]", "K":"[GT]", "M":"[AC]", "N":"[ACGT]", "R":"[AG]", "S":"[GC]", "T":"T", "V":"[ACG]", "W":"[AT]", "Y":"[CT]" }
	elif molecule == 'rna':
		degenerate_code = { "A":"A", "B":"[CGU]", "C":"C", "D":"[AGU]", "G":"G", "H":"[ACU]", "K":"[GU]", "M":"[AC]", "N":"[ACGU]", "R":"[AG]", "S":"[GC]", "U":"U", "V":"[ACG]", "W":"[AU]", "Y":"[CU]" }
	else:
		raise ValueError( "make_degenerate_regex requires molecule to be dna or rna" )
	regex_string = ''
	idx = 0
	while idx < len( motif_seq ):
		curr = motif_seq[idx]
		count = 1
		for next_idx in range( idx+1, len( motif_seq ) ):
			next = motif_seq[next_idx]
			if next == curr:
				count += 1
			else:
				break
		regex_string += degenerate_code[curr]
		if count > 1:
			idx = idx + count - 1
			regex_string += "{%s}" % ( count )
		idx += 1
	return regex_string

regex_1 = re.compile(make_degenerate_regex(args.motif5p))
regex_2 = re.compile(make_degenerate_regex(args.motif3p))

intron_l_pass = []
intron_cov_pass = []

out_valid_sjs = open(args.output + '_valid_SJs.txt', 'w')
out_invalid_sjs = open(args.output + '_invalid_SJs.txt', 'w')
out_valid_sjs.write('IntronID\tIntronLength\tCov(S,M)\n')
out_invalid_sjs.write('IntronID\tIntronLength\tCov(S,M)\n')

ct_frwrd = 0
valid_final_sjs = []

#filter valid SJs
for k, v in forward_frag.items():
	motif5 = str(v[0:7])
	branch_motif = v[-12:]
	match_res1 = regex_1.search(motif5)
	match_res2 = regex_2.search(branch_motif)
	if match_res1 and match_res2:
		valid_final_sjs.append(k)
		intron_l_pass.append(len(v))
		intron_cov_pass.append(int(intron_coverage[k].split(',')[0]) + int(intron_coverage[k].split(',')[1]))
		ct_frwrd += 1
		out_valid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\t' + str(match_res1) + '\t' + str(match_res2) + '\n')
	else:
		out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')
		
ct_rev = 0
for k, v in reverse_frag.items():
	motif5 = str(v[0:7])
	branch_motif = v[-12:]
	match_res1 = regex_1.search(motif5)
	match_res2 = regex_2.search(branch_motif)
	if match_res1 and match_res2:
		valid_final_sjs.append(k)
		ct_rev += 1
		intron_l_pass.append(len(v))
		intron_cov_pass.append(int(intron_coverage[k].split(',')[0]) + int(intron_coverage[k].split(',')[1]))
		out_valid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\t' + str(match_res1) + '\t' + str(match_res2) + '\n')
	else:
		out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')

ct_undef = 0
for k, v in undef_frag.items():
	motif5 = str(v[0:7])
	branch_motif = v[-12:]
	match_res1 = regex_1.search(motif5)
	match_res2 = regex_2.search(branch_motif)
	if match_res1 and match_res2:
		valid_final_sjs.append(k)
		intron_l_pass.append(len(v))
		intron_cov_pass.append(int(intron_coverage[k].split(',')[0]) + int(intron_coverage[k].split(',')[1]))
		ct_undef += 1
		if k.startswith('undef_plus'):
			k = k[10:]
		elif k.startswith('undef_mins'):
			k = k[10:]
		else:
			k = k
		out_valid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\t' + str(match_res1) + '\t' + str(match_res2) + '\n')
	else:
		try:
			out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')
		except KeyError:
			k = k[10:]
			out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')

##########################################################################################
#Summary and results

print('SJ entries in forward passing the filtering process: ' + str(ct_frwrd))
print('SJ entries in reverse passing the filtering process: ' + str(ct_rev))
print('SJ entries, unknown strand passing the filtering process: ' + str(ct_undef))

print ('\n\t----------------- Valid SJs, Summary -----------------\n')

print('Motif to search at 5p region: ' +  str(args.motif5p) + ', regex: ' + str(make_degenerate_regex(args.motif5p)))
print('Motif to search at 3p region: ' +  str(args.motif3p) + ', regex: ' + str(make_degenerate_regex(args.motif3p)))
print('\n')

print('Minimum intron length: ' + str(min(intron_l_pass)))
print('Maximum intron length: ' + str(max(intron_l_pass)))

plt.hist(intron_l_pass, bins=1000)
plt.figsize(60, 20)     #######para la version mas reciente de python "plotsize" debe reemplazarse por "plot_size"
plt.grid(True)
plt.title("Intron length distribution")
plt.xlabel("Length")
plt.show()
print('\n')
out_valid_sjs.close()

print('Valid SJs coordinates saved as: ' + args.output + '_valid_SJs.txt')
print('Invalid SJs coordinates saved as: ' + args.output + '_invalid_SJs.txt\n')

##########################################################################################
#BAM filtering

bamfile = pysam.AlignmentFile(args.bam_f, "rb")
bam_out = pysam.AlignmentFile(args.output + '_outBAM.filtered.bam', "wb", template=bamfile)

ct_xs = 0
disc_sj = 0
out_splicing3plusint = open('reads_3ormoreIntrons.txt', 'w')

for read in bamfile.fetch():
	if read.has_tag("XS"):
		ct_xs += 1
		try:
			bamsj_id = read.to_string().split()[2] + ':' + str(read.get_tag("jI")[0]) + '-' + str(read.get_tag("jI")[1])
			if len(read.get_tag("jI")) >=3:
				out_splicing3plusint.write(str(read) + '\n')
			if bamsj_id in valid_final_sjs:
				#print(read.to_string().split()[2] + '\t' + str(read.get_tag("jI")) + '\t' + str(read.get_tag("jI")[0]))
				bam_out.write(read)
			else:
				disc_sj += 1
		except IndexError:
			bamsj_id = read.to_string().split()[2] + ':' + str(read.get_tag("jI")[0])
			if bamsj_id in valid_final_sjs:
				#print(read.to_string().split()[2] + '\t' + str(read.get_tag("jI")) + '\t' + str(read.get_tag("jI")[0]))
				bam_out.write(read)
			else:
				disc_sj += 1
	else:
		bam_out.write(read)

bamfile.close()
bam_out.close()
out_splicing3plusint.close()

print('Splicing entries in the BAM file: ' + str(ct_xs))
print('Splicing entries discarded: ' + str(disc_sj))
print('Reads with 3 or more introns (if any) are saved as: reads_3ormoreIntrons.txt')
print('Filtered BAM: ' + args.output + '_outBAM.filtered.bam\n\n')
