#!/usr/bin/env python
#-*- coding: UTF-8 -*-
from __future__ import division
import sys
import argparse
import operator
import re
import pysam
import six


program = 'tacos_v3.py'

parser = argparse.ArgumentParser(prog=program, formatter_class=argparse.RawTextHelpFormatter, add_help=False,
    description="""\n\n\tThis python script is part of the pipeline: Trichomonad Assembler for COmplex Splicing. 
    
    It was developed for Trichomonas vaginalis, but it can be used to filter chimeric alignments in other species caused by non-conserved intron motifs.
    The 5' and 3' intron motifs can have a max. length of 7 and 12 nucleotides, respectively.

    ---------------------
    
    This version was tested on python 3.8.12
    
    """)
requiredNamed = parser.add_argument_group('Mandatory arguments')
requiredNamed.add_argument("-f", "--input", dest='input', required=True, type=str, help='Input fasta file (Genome reference)')
requiredNamed.add_argument("-sj", "--input_sj", dest='sj_star', required=True, type=str, help='SJ.out.tab from STAR mapping')
requiredNamed.add_argument("-b", "--input_bam", dest='bam_f', required=True, type=str, help='BAM file from STAR mapping')
optionalNamed = parser.add_argument_group('Optional arguments')
optionalNamed.add_argument("-o", "--output", dest='output', default='filtered', type=str, help='Prefix: Prefix for output files, default=filtered')
optionalNamed.add_argument("-5m", "--5p-motif", dest='motif5p', default='GTWYDNH', type=str, help='String: Splicing motif at 5p (upper case nucleotides only, no spaces, up to 7 nucleotides), default:GTWYDNH')
optionalNamed.add_argument("-3m", "--3p-motif", dest='motif3p', default='DYYWAHMHDYAG', type=str, help='String: Splicing motif at 3p (upper case nucleotides only, no spaces, up to 12 nucleotides), default:TYTAAYHWNCAG')
optionalNamed.add_argument("-v", "--verbose", dest='extended_output', choices=['y', 'yes', 'n', 'no', 'Y', 'YES', 'N', 'NO', 'Yes', 'No'], default='yes', type=str, help='Verbose mode, (choices= y, yes, n, no, Y, YES, N, NO, Yes, No), default: yes')
optionalNamed.add_argument("-h", "--help", action="help", help="Show this help message and exit")

args = parser.parse_args()


verbose_mode = args.extended_output.lower() in ['y', 'yes']
if verbose_mode:
	import plotext as plt

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

cum_list=[]
y = 0 

len_sorted = sorted(lens_f, reverse=True)
for x in range(0,len(len_sorted)):
	y+=len_sorted[x]
	cum_list.append(y)

print ('\n\n\t----------------- Fasta Summary -----------------\n')
if verbose_mode:
	plt.plot(cum_list)
	plt.scatter(cum_list)
	plt.figsize(60, 20) 
	plt.grid(True)
	plt.title("Cumulative length of the genome reference")
	plt.xlabel("Contig_number")
	plt.show()

else:
	print(f'Sequences in fasta file: {len(fasta_dir.keys())}')	
	print(f'Max. len in fasta file: {len_sorted[0]}')	
	print(f'Min. len in fasta file: {len_sorted[-1]}')	

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

if verbose_mode:
	plt.hist(intron_lens,bins=1000)
	#replace "plotsize" for "plot_size"?
	plt.figsize(60, 20)     
	plt.grid(True)
	plt.title("Intron length distribution")
	plt.xlabel("Length")
	plt.show()
	print ('\n\t----------------- SJs, Summary -----------------\n')
	print(f'SJs, total entries: {sj_entries}')
	print(f'SJ entries in forward: {len(forward_frag.keys())}' )
	print(f'SJ entries in reverse: {len(reverse_frag.keys())}' )
	print(f'SJ entries, undefined strand: {len(undef_frag.keys())}' )
	
	print ('\n\t----------------- SJs, intron Summary -----------------\n')
	print(f'Minimum intron length: {min(intron_lens)}')
	print(f'Maximum intron length: {max(intron_lens)}')
else:
	print('\n\n')
	print ('\n\t----------------- SJs, Summary -----------------\n')
	print(f'SJs, total entries: {sj_entries}')
	print(f'SJ entries in forward: {len(forward_frag.keys())}' )
	print(f'SJ entries in reverse: {len(reverse_frag.keys())}' )
	print(f'SJ entries, undefined strand: {len(undef_frag.keys())}' )
	
	print ('\n\t----------------- SJs, intron Summary -----------------\n')
	print(f'Minimum intron length: {min(intron_lens)}')
	print(f'Maximum intron length: {max(intron_lens)}')

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
		#out_valid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\t' + str(match_res1) + '\t' + str(match_res2) + '\n')
		out_valid_sjs.write('\t'.join([k, str(len(v)),str(intron_coverage[k]),str(match_res1),str(match_res2)]) + '\n')

	else:
		#out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')
		out_invalid_sjs.write('\t'.join([k,str(len(v)),str(intron_coverage[k])]) + '\n')
		
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
		#out_valid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\t' + str(match_res1) + '\t' + str(match_res2) + '\n')
		out_valid_sjs.write('\t'.join([k,str(len(v)),str(intron_coverage[k]),str(match_res1),str(match_res2)]) + '\n')
	else:
		#out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')
		out_invalid_sjs.write('\t'.join([k,str(len(v)),str(intron_coverage[k])]) + '\n')

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
		#out_valid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\t' + str(match_res1) + '\t' + str(match_res2) + '\n')
		out_valid_sjs.write('\t'.join([k,str(len(v)),str(intron_coverage[k]),str(match_res1),str(match_res2)]) + '\n')
	else:
		try:
			#out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')
			out_invalid_sjs.write('\t'.join([k,str(len(v)),str(intron_coverage[k])]) + '\n')

		except KeyError:
			k = k[10:]
			#out_invalid_sjs.write(k + '\t' + str(len(v)) + '\t' + str(intron_coverage[k]) + '\n')
			out_invalid_sjs.write('\t'.join([k,str(len(v)),str(intron_coverage[k])]) + '\n')

##########################################################################################
#Summary and results

print(f'SJ entries in forward passing the filtering process: {ct_frwrd}')
print(f'SJ entries in reverse passing the filtering process: {ct_rev}')
print(f'SJ entries, unknown strand passing the filtering process: {ct_undef}')

print ('\n\t----------------- Valid SJs, Summary -----------------\n')
print(f'Motif to search at 5p region: {args.motif5p}') #+ ', regex: ' + str(make_degenerate_regex(args.motif5p)))
print(f'Motif to search at 3p region: {args.motif3p}') #+ ', regex: ' + str(make_degenerate_regex(args.motif3p)))
print('\n')

print(f'Minimum intron length: {min(intron_l_pass)}')
print(f'Maximum intron length: {max(intron_l_pass)}')

if verbose_mode:
	plt.hist(intron_l_pass, bins=1000)
	#replace "plotsize" for "plot_size"?
	plt.figsize(60, 20)
	plt.grid(True)
	plt.title("Intron length distribution")
	plt.xlabel("Length")
	plt.show()
	print('\n')
	
out_valid_sjs.close()
out_invalid_sjs.close()

print(f'Valid SJs coordinates saved as: {args.output}_valid_SJs.txt')
print(f'Invalid SJs coordinates saved as: {args.output}_invalid_SJs.txt\n')

##########################################################################################
#BAM filtering

print('\nFiltering BAM file ... ')
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

print('\nFiltering BAM file completed ... \n\n')

bamfile.close()
bam_out.close()
out_splicing3plusint.close()

print ('\t----------------- Final summary -----------------\n')
print(f'Splicing entries in the BAM file: {ct_xs}')
print(f'Splicing entries discarded: {disc_sj}')
print(f'Reads with 3 or more introns (if any) are saved as: reads_3ormoreIntrons.txt')
print(f'Filtered BAM: {args.output}_outBAM.filtered.bam\n\n')
