#!/usr/bin/python
# coding=utf-8



# importing modules
# to read and write csv files, sys - to take arguments
import csv, sys
import argparse
import MySQLdb as mdb
from datetime import datetime
import time
start_time =datetime.now()

"######### NOTES before Running prox_fam ##########"
# Pre-requisite tables: 
# countfams - has to be made before the script can be run 
# smurt table

# In SMG_proximity
# gff_table - named gff_ultimo is needed
 
# IN SMG_family three tables have to be set: 
# intermediate_table = 'ik_prox_fam' # intermediate table with clust+prox genes coupled to hfam
# hfam_table = 'hfam' # table containing the proteins and the protein family they belong to
# count_fam = 'countfams' # table containing the protein families, organisms and the number of homologs in each organism. 

# to run this script:
#python SMG_prox_fam_ver3_Pip.py -c TRUE -rr TRUE -de TRUE 



# connecting to the database
host = "XXX"
user = "user_name"
password = "***"
dbname = "database"
db = mdb.connect( host,user,password,dbname )
cursor = db.cursor()
# connect to DB
db = mdb.connect( host,user,password,dbname )
cursor = db.cursor()



"####### USER INPUT #######"

parser = argparse.ArgumentParser()
parser.add_argument("-clear", "-c", required=False, default='FALSE',
					 help="Clear and replace data if org_id is already found in output table.")
parser.add_argument("-plimit", "-pl", required=False, default=0,
					 help="Set the number of genes to to include before and after the cluster.")
# Input/output for SMG_prox2
parser.add_argument("-input", "-i", required=False, default='smurf',
					 help="Set the secondary metabolite Table to use as input.")
parser.add_argument("-output", "-o", required=False, default='clust_prox',   # ok
					 help="Set the Table to use as output for SMG_prox2.")
# Input/output for SMG_fam2
parser.add_argument("-input_fam", "-if", required=False, default='clust_prox', 
 					 help="Set the Table to use as input.")
parser.add_argument("-output_fam", "-of", required=False, default='prox_fam_count',
 					 help="Set the Table to use as output.")

parser.add_argument("-delete", "-de", required=False, default='TRUE',
					 help="Delete the ik_prox_xx table if it already exists.")

# add arument for choosing if new table for SMG_prox should be made or not
parser.add_argument("-rerun", "-rr", required=False, default='FALSE',
					 help="Rerun SMG_prox to create new ik_clust_prox table or use already existing.")



# User input arguments
# https://docs.python.org/2/howto/argparse.html#id1
# https://docs.python.org/2/library/argparse.html
args = parser.parse_args()
print args.clear

# Set variables
clear_flag = args.clear					# Set clear_flag variable if TRUE the organisms already in the output table is deleted and then added from the new 
prox_limit = args.plimit 				# Set the number of genes to include before and after a cluster
input_table = args.input 				# Set the input table with the SM gene clusters ()
output_table = args.output 				# Set the output table
input_table_fam = args.input_fam 		# Set the input table with the SM gene clusters
output_table_fam = args.output_fam 		# Set the output table
clear_prox_fam = args.delete 			# Set flag - delete ik_prox_fam table if TRUE for SMG_prox
rerun = args.rerun 						# Set if the SMG_prox should be rerun


"########### CHECK NEEDED TABLES ARE THERE AND UPDATED ###############"

if cursor.execute("SHOW TABLES LIKE 'countfams';"):
	# CHECK if the input orgs are there in countfams table
	cursor.execute("SELECT distinct org_id from countfams;") # get the org_ids from countfams
	countfams_org = cursor.fetchall()
	org_countfams = [x[0] for x in countfams_org]

	cursor.execute( "SELECT distinct org_id from " + input_table +";") # get the org_ids from the input data (smurf cluster)
	input_org = cursor.fetchall()
	org_input = [x[0] for x in input_org]		# save org_id from input as a list in org_input 
	
	for i in range(0,len(org_input)):
		if org_input[i] not in org_countfams:
			print "Organism not found in countfams: ", org_input[i]
			#sys.exit("ERROR: The table ik_countfams_resistance need to be updated there is an organism from the input data that is not in ik_countfams_xx.")
else:
	sys.exit("ERROR: The table ik_countfams_resistancepipeline does not exist and have to be created .")
	# If the table is not present you can use the code below to create the table using mysql
	# CREATE TABLE countfams
	# SELECT * FROM (	
	# SELECT hfam, org_id, count(hfam) as count FROM hfam
	# group by org_id, hfam) as tt;


"# Running SMG_prox_ver3 - this will find the genes in proximity of all the clusters (in case smurf predictions are not accurate."
if rerun == 'TRUE':
	# calling function to create table with all SMG and genes in the proximity
	print "SMG_prox will be rerun and a new table of " +output_table+ " will be generated."
	from SMG_prox_ver3_Pip import SMG_proximity
	SMG_proximity(clear_flag = clear_flag, prox_limit= prox_limit, input_table = input_table, output_table= output_table)
else:
	if not cursor.execute("SHOW TABLES LIKE \'" + output_table + "\';"):
		sys.exit("ERROR: The table: " +output_table+ " does not exist and have to be created set -rerun / -rr to TRUE.")
	else:
		print "The existing table of " +output_table+ " will be used in the subsequent part to find protein families."

"# Running SMG_fam_ver3 - finding the protein family for each gene and counting the number of this protein family in clusters, outside clusters, for each organism"
"# SMG_family uses an intermediate table - rember to check if this is correct"
from SMG_fam_ver3_Pip import SMG_family
SMG_family(clear_prox_fam = clear_prox_fam, input_table_fam=input_table_fam, output_table_fam= output_table_fam )

# Make sure to check the tables used in the SMG_family 
#	intermediate_table = 'ik_prox_fam_v3' 
#	hfam_table = 'hfam_SL_ID50_SC130_noScerNeuAzo'
#	count_fam = 'ik_countfams_Flavi'


# closing database
db.close()

print "THE END"
print "# INFO Total Runtime for SMG_prox_fam: ", (datetime.now()-start_time)

