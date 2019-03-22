#!/usr/bin/python
#Feb. 19, 2019
#This script inspects sam file and creates summary of all reads to detect sequence variations

import sys, subprocess
from optparse import OptionParser

dict_count={}
dict_fasta={}
	
	
def main():
	#parse user input
	parser=OptionParser()
	usage = "usage: python %prog [options] arg1 arg2"
	parser.add_option("-f","--samfile",action="store",type="string",dest="sam_input",help="specify the sam file to be analyzed")
	parser.add_option("-a","--fastafile",action="store",type="string",dest="fasta_input",help="specify the fasta file used for mapping")
	parser.add_option("-s","--start",action="store",type="int",dest="start",default=10,help="position of sequences to start analysis")
	parser.add_option("-e","--end",action="store",type="int",dest="end",default=70,help="position of sequences to end analysis")
	parser.add_option("-c","--summarycutoff",action="store",type="int",dest="mapped_cutoff",default=10,help="reads cutoff for summary of mapped reads")
	parser.add_option("-u","--unmappedcutoff",action="store",type="int",dest="unmapped_cutoff",default=100,help="read cutoff for summary of unmapped reads")
	(option,args)=parser.parse_args()
	if not option.sam_input:
		parser.error ("option -f is required (sam file needed)")
	if not option.fasta_input:
		parser.error ("option -a is required (fasta file needed)")
	sam_name=option.sam_input

	#create file handles
	unmapped_reads=open(sam_name+"_unmapped_reads.txt","w")
	mapped_reads=open(sam_name+"_mapped_reads.txt","w")
	mapped_reads_summary=open(sam_name+"_mapped_reads_summary.txt","w")
	
	with open(option.sam_input) as fp:
		for line in (line for line in fp if not line.startswith('@')):
			line=line.strip('\n').split('\t')
			ref_name=line[2]
			start_pos=line[3]
			if int(start_pos) <=option.start:
				start,end=where_to_start(start_pos, option.start,option.end)
				seq=line[9][start:end]
				full_name=ref_name+":"+seq
				count_occurrence(full_name,dict_count)
		give_output(dict_count, unmapped_reads, mapped_reads,option.unmapped_cutoff)
		sorting(unmapped_reads,mapped_reads)
		run()
		read_fasta(option.fasta_input,">",option.start,option.end)
		summarize(mapped_reads,option.mapped_cutoff,mapped_reads_summary)

#Determine which position to segment each read
def where_to_start(actual_start,where_to_start,where_to_end):
	if actual_start == "0":
		begin=int(where_to_start-1)
		stop=int(where_to_end-1)
	if actual_start != "0":
		begin=abs(where_to_start-int(actual_start))
		stop=where_to_end-int(actual_start)
	return(begin,stop)

#Count how many times each read occurs
def count_occurrence(item,dictionary):
	if item not in dictionary.keys():
		dict_count[item]=1
	elif item in dictionary.keys():
		dict_count[item]=dictionary[item]+1
	return (dictionary)

#write results to output files
def give_output(dictionary,output1,output2,cutoff):
	for key,value in dictionary.items():
			if key.startswith("*") and value>cutoff:
				output1.write("%s %s\n" %(key,value))
			if not key.startswith("*"):
				output2.write("%s %s\n" %(key,value))
	output1.close()
	output2.close()

#sort sequences based on number of occurrences in descending order
def sorting(input1, input2):
	bash=open("commands.sh","w")
	opening="#!/bin/bash"
	command1="sort -k2 -n "+input1.name + ">"+input1.name+"_sorted.txt"
	string="awk -F'[: ]' '{print $1\" \"$2\" \"$3}'"
	command2="cat "+input2.name+"| "+string +" | sort -k1,1 -k3rn  > "+input2.name +"_sorted.txt"
	command3="rm "+input1.name+" "+input2.name
	bash.write("%s\n%s\n%s\n%s\n" %(opening,command1, command2,command3))
	bash.close()


#Run the bash command generated in "sorting"
def run():
	subprocess.call(['chmod', '+x', './commands.sh'])
	subprocess.call(['./commands.sh'],shell=True)

#compare each sequence to the reference
def summarize(input_file,cutoff_of_sig,output):
	file_name=input_file.name+"_sorted.txt"
	with open(file_name,"r") as unsum_result:
		for line in unsum_result:
			line=line.strip('\n').split(" ")
			if int(line[2])>cutoff_of_sig:
				ref_seq=dict_fasta[line[0]]
				mis=compare_two_seq(ref_seq,line[1])
				if mis:
					output.write("%s\t%s\t%s\t%s\n" %(line[0],line[1],line[2],mis))
				if not mis:
					output.write("%s\t%s\t%s\n" %(line[0],line[1],line[2]))
	output.close()

#compare each read to the reference sequence
def compare_two_seq (seq_a,seq_b):
	mismatches=[]
	for pos in range(0,min(len(seq_a),len(seq_b))):
		if seq_a[pos]!=seq_b[pos]:
			string=str(int(pos)+1)+":"+seq_a[pos]+">"+seq_b[pos]
			mismatches.append(string)
	return ("\t".join (mismatches))

#read fasta file into memory
def read_fasta (input_file,symbol,start,end):
	with open(input_file) as amplicons:
		for line in amplicons:
			if line.startswith(symbol):
				name=line[1:-1]
				continue
			dict_fasta[name]=line.strip()[start-1:end-1]
	return(dict_fasta)

if __name__=="__main__":
	main()

