#!/usr/bin/python
# coding=utf-8

"#Aligning, trimming and realigning"

import os

from Bio import SeqIO


#Return a list of the entries in the directory given
names_fasta = os.listdir("/your_directory") 




path="path/to/your_directory"

# Add name of folder where the first alignments are added
align_folder ="/Name_align_folder/" 
# Add name of folder where the final alignments are stored
final_folder ="/Name_final_folder/"

# Go through - each fasta file  
for name in names_fasta:
	j = name.split(".")
	
	if len(j)>1 and j[1] == 'fasta':
		print j[0]
		
		## Opening the file to get the number of sequences
		count_seq = 0
		for record in SeqIO.parse(path+'/'+name, "fasta"):
			count_seq = count_seq + 1
		
		print "number of sequences: ", count_seq

		print "First alignment"
		o_name = j[0] + "_align.fasta"
		
		# align
		command = "clustalo -i "+path+"/"+name+" -o "+path+align_folder+o_name+" --infmt=fasta --outfmt=fasta --force"		

		#print command
		os.system(command)

		print "Trimming"
		
		"# Trimming using Gblocks"
		set_b2= str(round(count_seq/2)+1)
		command_Gblocks = "Gblocks "+path+align_folder+o_name+" -t=p -b2="+set_b2+" -b3=12 -b4=2 -b5=h  -e=Gtrim "

		os.system(command_Gblocks)

		# align
		print "Second alignment"
		o_name_end_aut = j[0] + "_final.fasta"
		
		command3 = "clustalo -i "+path+align_folder+j[0]+"_align.fastaGtrim -o "+path+final_folder+o_name_end_aut+" --outfmt=fasta --force"
		
		os.system(command3)
		
print "the end"



#--------------------------------------------------------
# Gblocks
#--------------------------------------------------------
# Parameters:
# -t= Type Of Sequence (Protein, DNA, Codons) - p
# -b0= (This option does not appear in the menu) Minimum Length Of An Initial Block (Same as Minimum Length Of A Block) 
# -b1= Minimum Number Of Sequences For A Conserved Position (50% of the number of sequences + 1) - will keep (default) - KEEP can't lower
# -b2= Minimum Number Of Sequences For A Flank Position (85% of the number of sequences) - will keep (default) - bigger values -> decrease selected number og positions - set to 50%
# -b3= Maximum Number Of Contiguous Nonconserved Positions (8) - will keep (default) - big num -> increase selected number of positions
# -b4= Minimum Length Of A Block (10) - 5 (NOT default) -> bigger values -> decrease the selcted number of positions -> Set to 12
# -b5= Allowed Gap Positions (None, With Half, All) - h (NOT default) (default = None)
# -b6= Use Similarity Matrices (Yes, No) - y (default)
# -s= Selected Blocks (Yes, No) - y (default)
# -p= Results And Parameters File (Yes, Text, Short Text, No) - y (default)
# -v= (Only visible in the extended saving options) Characters Per Line In Results And Parameters File (60) - will keep (default)
# -n= (Only visible in the extended saving options) Nonconserved Blocks (Yes, No) - n (default)
# -u= (Only visible in the extended saving options) Ungapped Alignment (Yes, No) - y (NOT default)
# -k= (Only visible in the extended saving options) Mask File With The Selected Blocks (Yes, No) - n (default)
# -d= (Only visible in the extended saving options) Postscript File With The Selected Blocks (Yes, No) - n (default)
# -a= (Only available with paths files) Concatenated Blocks From Alignments In Batch (Yes, No) - n (default)
# -c= (Only available with paths files) Concatenated Input Alignments In Batch (Yes, No) - n (default)
# -w= (Only available with paths files) Concatenated Ungapped Alignments In Batch (Yes, No) - n (default)
# -e= Generic File Extension (-gb) _gb (NOT default)


