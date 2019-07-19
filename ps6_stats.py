#!/usr/bin/env python3
import re
import argparse

def get_args():
	parser = argparse.ArgumentParser(description="A program to introduce yourself")
	parser.add_argument("-i", "--file_name", help="Specify the file name", required=True)
	parser.add_argument("-k", "--kmer_len", help="Specify the file name", required=True)
	return parser.parse_args()
args = get_args()

k = args.kmer_len
INPUT_FILE = args.file_name
kmer_coverage_lst = []       #for DATA SHOWCASE: median coverage depth

#EXTRACT FASTA HEADER AND SEQ
INPUT_FILE_1 = open(INPUT_FILE, "r")
for line in INPUT_FILE_1 :
    line = line.strip()
    if ">" in line:
        #GET CONTIG LEN
        kmer_len = re.findall("[0-9]+_c", line)[0]     #from header extract k-mer length of each contig
        kmer_len = kmer_len.strip("_c")
        real_len = (int(kmer_len) + (int(k) - 1))
        #GET COVERAGE
        kmer_coverage = re.findall("ov_.+", line)[0]       #from header extract the k-mer coverage for the contig
        kmer_coverage = kmer_coverage.strip("ov_")
        kmer_coverage_lst.append(float(kmer_coverage))             #for DATA SHOWCASE: median coverage depth
INPUT_FILE_1.close()


#SHOWCASE DATA: 1.number of contigs, 2.total length of the genome, 3.largest contig,
#               4.average contig length, 5.mean depth of contig's coverage, 6.N50, 7.bucket distribution
#number of contigs
INPUT_FILE_2 = open(INPUT_FILE, "r")
num_contigs = 0
for line in INPUT_FILE_2:
    if ">" in line:
        num_contigs += 1
print("The number of contigs is: " + str(num_contigs))
INPUT_FILE_2.close()

#total length of the genome
INPUT_FILE_3 = open(INPUT_FILE, "r")
total_num_bp = 0
len_contig_lst = []
for line in INPUT_FILE_3:
    if ">" not in line:
        seq = line.strip()
        len_contig_lst.append(len(seq))   #will use for avg contig len
        total_num_bp += len(seq)
print("The total length of the genome assembly is: " + str(total_num_bp) + " basepairs")
INPUT_FILE_3.close()

#largest contig
len_contig_lst.sort(reverse = True)
longest_contig = len_contig_lst[0]
print("The largest contig is: "+ str(longest_contig) + " basepairs")

#average contig length
sum_contig_len = 0
for contig_len in len_contig_lst:
    sum_contig_len += contig_len
average_contig_len = (sum_contig_len/num_contigs)
print("The average contig length is: " + str(average_contig_len) + " basepairs")

#mean depth of contig's coverage
sum_kmer_coverage = 0
for kmer_coverage in kmer_coverage_lst:
    sum_kmer_coverage += kmer_coverage
mean_cov_depth = (sum_kmer_coverage/num_contigs)
print("The mean coverage depth is: " +  str(mean_cov_depth))

#N50
#find half the sum of the contig lengths
tot_sum = 0
for contig_len in len_contig_lst:
    tot_sum += contig_len
    half_tot_sum = tot_sum/2
#find the contig len associated with the half of total sum
running_sum = 0
for i in range(len(len_contig_lst)):
    running_sum += len_contig_lst[i]
    if running_sum >= half_tot_sum:
        break
print("The N50 if this genome assembly is: " + str(len_contig_lst[i]))

#BUCKET MAKING ALGORITHM
distrib_dict = {}
#make a starter ditionary with keys as hundreds (0, 100, 200, etc.) and values as zero
for i in range(0, (longest_contig), 100):
    distrib_dict[i] = 0
#finding the hundred associated with the contig len in the len contig lst
for contig_len in len_contig_lst:
    for hundred in distrib_dict.keys():
 #once the contig len is larger than the hundred add to bucket (ex: 1142 > 1100, add to 1100 bucket)
        if contig_len >= hundred:
            bucket = hundred
 #add 1 into the bucket
    if bucket in distrib_dict:
        distrib_dict[bucket] += 1
    else:
        distrib_dict[bucket] = 1

#PRINT BUCKET DISTRIBUTION
print("# Contig length" + "\t" + "Number of contigs in this category")
for bucket, freq in distrib_dict.items():
    print(str(bucket) + "\t" + str(freq))
