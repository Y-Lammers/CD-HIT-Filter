#!/usr/bin/env python

# Author: Youri Lammers
# Contact: youri.lammers@naturalis.nl / youri.lammers@gmail.com

# Filter the cluster output from CD-hit based on
# the minimum number of sequences per cluster

# command: CD-hit_filter.py -f [cluster .fasta file] -c [.clstr cluster file] -m [min # seqs per cluster]

# the output of the CD-hit_filter.py script is a fasta file [input cluster .fasta file with the min # seqs in the file name]
# with the OTU sequences for the clusters with more than or equal number of sequences to the user specified minimum.

# import modules used by the script
import argparse, os, itertools

# set argument parser
parser = argparse.ArgumentParser(description = 'Filter the output from CD-hit based on the minimum number of read per cluster.\nThe filtered output fasta file produced has the same name as the input file with  _min_[minimum size from -c argument].fasta attachted to the name.')

parser.add_argument('-f', '--fasta', metavar='.fasta file', dest='fasta', type=str,
			help='The .fasta file containing the clusters produced by CD-hit.')
parser.add_argument('-c', '--cluster', metavar='.clstr file', dest='cluster', type=str,
			help='The .clstr file producec by CD-hit that contains the cluster information.')
parser.add_argument('-m', '--minimum', metavar='minimum size', dest='minimum', type=int,
			help='The minimum cluster size.')

args = parser.parse_args()

def read_clstr():

	# parse through the .clstr file and create a dictionary
	# with the sequences per cluster

	# open the cluster file and set the output dictionary
	cluster_file, cluster_dic = open(args.cluster), {}

	# parse through the cluster file and store the cluster name + sequences in the dictionary
	cluster_groups = (x[1] for x in itertools.groupby(cluster_file, key=lambda line: line[0] == '>'))
	for cluster in cluster_groups:
		name = cluster.next().strip()
		seqs = [seq.split('>')[1].split('...')[0] for seq in cluster_groups.next()]
		cluster_dic[name] = seqs

	# return the cluster dictionary
	return cluster_dic


def clstr_filt(cluster_dic):

	# filter the cluster in the cluster_dic based
	# on the minimum cluster size

	# grab a list of the cluster_dic keys
	# and create an empty list for the filtered sequence ids
	clusters, filtered_dic = list(cluster_dic), {}

	# parse through the keys and check the cluster size,
	# if the size is > minimum size, add the cluster to the
	# filtered dic
	for cluster in clusters:
		if len(cluster_dic[cluster]) >= args.minimum:
			for sequence in cluster_dic[cluster]:
				filtered_dic[sequence] = cluster

	# return the filtered cluster_dic
	return filtered_dic

def parse_fasta(filtered_dic):

	# parse through the cluster fasta sequences
	# create a new _filtered_[minimum size].fasta file with
	# the fasta sequences that are larger than the minimum size

	# open the fasta file and set lists for the sequences
	seq_file, seq_list, = open(args.fasta), []

	# parse through the fasta file and obtain the sequence
	seq_groups = (x[1] for x in itertools.groupby(seq_file, key=lambda line: line[0] == '>'))
	for header in seq_groups:
        	header = header.next().strip()
	        sequence = ''.join(seq_line.strip() for seq_line in seq_groups.next())
        	seq_list.append([header,sequence])

	# close the sequence file
	seq_file.close()

	# open the filter_output file
	filter_file = open('{0}_min_{1}.fasta'.format(os.path.splitext(args.fasta)[0], args.minimum), 'w')

	# loop through the sequence list
	while len(seq_list) > 0:
	
		# retrieve a fasta sequence from the list
		seq = seq_list.pop()
		
		# check if the sequence is present in the filtered_list
		if seq[0][1:] in filtered_dic:
			
			# if present, write the header / sequence to the filter_file
			filter_file.write('{0}\n{1}\n'.format(filtered_dic[seq[0][1:]], '\n'.join([seq[1][i:i+60] for i in range(0, len(seq[1]), 60)])))

	# close the filtered results file
	filter_file.close()


def main():

	# obtain a dictionary with the clusters and sequences
	cluster_dic = read_clstr()

	# get a list of sequence in larger clusters
	filtered_dic = clstr_filt(cluster_dic)
	
	# write the filtered sequence to a new file
	parse_fasta(filtered_dic)

if __name__ == '__main__':
	main()
