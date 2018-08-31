#!/usr/bin/python
# coding=utf-8

# importing modules
# to read and write csv files, sys - to take arguments
import csv, sys
import argparse
import MySQLdb as mdb
from datetime import datetime
import time


def SMG_proximity(clear_flag, prox_limit, input_table, output_table):
	start_time =datetime.now()

	
	print "#---------------------------------------------------------------------------#"
	print "#                   Creating a table with  SM genes    "
	print "#                           and genes in proximity. "
	print "#---------------------------------------------------------------------------#"



	# connecting to the database
	host = "XXX"
	user = "user_name"
	password = "***"
	dbname = "database"
	# connect to DB
	db = mdb.connect( host,user,password,dbname )
	cursor = db.cursor()


	"#################### CHECK STARTING DATA ##########################"

	cursor.execute("SELECT org_id, count(org_id) FROM " + input_table +
		" group by org_id;")
	input_check = cursor.fetchall()

	print "######### Initial check of input data #############"
	print "The number of organisms is: ", len(input_check)
	org_input = [x[0] for x in input_check]		# save org_id from input as a list in org_input 
	print "The org_ids in the input is:", org_input

	"##### create MySQL table for results #####"
	# MySQL query to create new table (just an empty table to put the results in but with columns designed for result)- clust_prox
	# if the table is not already found
	if not cursor.execute("SHOW TABLES LIKE \'"+output_table+"\';"):
		print "# Creating table for the output:", output_table
		cursor.execute("CREATE TABLE `"+output_table+"` (\
			`u_clust_id` varchar(50) NOT NULL DEFAULT '',\
	  		`clust_backbone` varchar(50) DEFAULT '',\
	  		`clust_size` int(11) DEFAULT NULL,\
	  		`clust_type` varchar(50) DEFAULT NULL,\
	  		`org_id` int(11) NOT NULL,\
	  		`gff_seqorigin` varchar(100) NOT NULL DEFAULT '',\
	  		`ustart` int(11) DEFAULT NULL,\
	  		`uend` int(11) DEFAULT NULL,\
	  		`protein_id` varchar(100) NOT NULL DEFAULT '',\
	  		`gene_order` bigint(84) DEFAULT NULL,\
	  		`end_flag` int(11) NOT NULL,\
	  		PRIMARY KEY (`u_clust_id`,`protein_id`),\
	  		KEY `org_protein` (`org_id`,`protein_id`)\
			) ENGINE=MyISAM DEFAULT CHARSET=latin1;")
		db.commit()

	"##### CHECK DATA input and output #####"
	# is input org already found in the output table 
	cursor.execute("SELECT distinct org_id from " + output_table + " ;") # get the org_ids for org already found in the output table
	output_check =cursor.fetchall()
	org_output_start = [x[0] for x in output_check]

	# checking time for running this part - deleting
	Del_start_time =datetime.now()
	# check if the input org_ids are found in the output already
	org_del = []
	for i in range(0,len(org_input)):
		if org_input[i] in org_output_start:
			org_del.append(org_input[i])
	print "to be deleted:", org_del
	# If org already found - then they will be deleted if clear_falg is TRUE
	if org_del:
		if clear_flag == 'FALSE':
			sys.exit("ERROR: The input organisms are already found in the table, to drop these and run again set -clear -c to TRUE.")
		else:
			print "The organisms already in the table", org_del, "will be removed from table."
			query = "DELETE FROM "+ output_table +" WHERE org_id in (%s)" % ','.join(['%s'] * len(org_del))

			org_del = tuple(org_del)
			
			try:
				cursor.execute(query, org_del)
				db.commit()
			except  mdb.Error, e:
				sys.exit("# ERROR mysql %d: %s" % ( e.args[0],e.args[1] ))

	print "# INFO Runtime for deletion of org in output table: ", (datetime.now()- Del_start_time)

	"#################### Getting the proximity genes #########################"
	# This is done since smurf predictions might not be correct and could be too strict, the genenes in close prximity could also be cluster gene 
	# - so potentially we also want to investigate them for possible resistance
	
	"#####  Create look up table for geneorder limites #####"
	Lookup_start_time =datetime.now()

	# getting the min and max gene order pr cluster plus the max geneorder for the scaffold it is located on
	# Columns: org_id, clust_backbone, gff_origin, ming (min gene order), maxg (max geneorder), s_max (max geneorder for scaffold)
	query_str1 = ("SELECT tb.org_id, clust_backbone, tb.gff_seqorigin, ming, maxg, max(geneorder) as s_max FROM (select ta.org_id, clust_backbone, gff_seqorigin, min(geneorder) as ming, max(geneorder) as maxg from " + 
		input_table + " as ta join gff_ultimo on sm_protein_id=protein_id and ta.org_id=gff_ultimo.org_id " +
		"group by ta.org_id, clust_backbone) as tb " +
		"join gff_ultimo on tb.org_id=gff_ultimo.org_id AND tb.gff_seqorigin = gff_ultimo.gff_seqorigin " +
		"group by tb.org_id, clust_backbone;")
	cursor.execute(query_str1)
	cluster_border = cursor.fetchall()

	print "# INFO Runtime for creating lookup table for geneorder: ", (datetime.now()- Lookup_start_time)


	"##### For each cluster get genes in and near cluster #####"
	values_to_insert = []
	counter = 0
	total_counter = 0
	Cre_start_time =datetime.now()

	for i in range(0,len(cluster_border)):
		counter += 1
		total_counter += 1
		flag = '0'
		
		# checking if the start of cluster is close to beginning of scaffold 
		# - if it is closer than the prox_limit (how many extra genes to include)
		#  then the genes ming (start of clust+prox) will be set to one, else it will be start of clust - proxlimit
		if cluster_border[i][3] < prox_limit:
			ming = '1'
			flag = '1'
		else:
			ming = str(cluster_border[i][3] - prox_limit)

		# checking if the end of cluster is near the end of the scaffold
		if cluster_border[i][4] + prox_limit > cluster_border[i][5]:
			maxg = str(cluster_border[i][5])
			flag = '2'
		else:	
			maxg = str(cluster_border[i][4] + prox_limit)

		# uniqe cluster column with org_id, gff_seqorigin, backbone concatenated - unique for each cluster	
		newcol = [str(cluster_border[i][0]) + "_" + str(cluster_border[i][2])+ "_" + (cluster_border[i][1])]
		
		# get all the cluster genes + proximity genes plus the info needed for the output table
		# 9 columns: clust_backbone, clust_size, clust_type, org_id, gff_seqorigin, ustart, uend, protein_id, geneorder
		query_str2 =("SELECT tA.clust_backbone, tA.clust_size, tA.clust_type, tb.org_id, tb.gff_seqorigin, tb.ustart, tb.uend, tb.protein_id, tb.geneorder FROM ( " +
				"SELECT  tt.org_id, tt.gff_seqorigin, gff_ultimo.ustart, gff_ultimo.uend, gff_ultimo.protein_id, gff_ultimo.geneorder FROM (\
					select distinct gff_ultimo.gff_seqorigin, org_id from gff_ultimo where gff_seqorigin= \'" + str(cluster_border[i][2]) + "\' and org_id= \'" + str(cluster_border[i][0]) + "\') as tt " +
				"join gff_ultimo on tt.gff_seqorigin=gff_ultimo.gff_seqorigin and gff_ultimo.org_id = tt.org_id " +
				"where gff_ultimo.geneorder between \'" + ming + "\' and \'" + maxg + "\') as tb " +
				"Left JOIN " + input_table + " as tA on tb.org_id = tA.org_id AND tb.protein_id = tA.sm_protein_id;")

		cursor.execute(query_str2)
		b_genes = cursor.fetchall()
		# For each line (protein id) make a line including the unique id (newcol) columns from table and the end flag
		# and correct the format
		for item in b_genes:
			line = [newcol, list(item), list(flag)]
			line = sum(line, [])	# making the list flat - one list out of the list of lists
			values_to_insert.append(tuple(line))

		# load data in bundels 
		if counter == 10 or total_counter == len(cluster_border):
		 	print "INFO: loading to table, total count: %s" %total_counter

			stmt = "INSERT INTO "+ output_table +" (u_clust_id, clust_backbone, clust_size, clust_type, org_id, gff_seqorigin, ustart, uend, protein_id, gene_order, end_flag) \
			VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);" 

			try:
				cursor.executemany(stmt, values_to_insert)
				db.commit()
			except  mdb.Error, e:
				print values_to_insert
				sys.exit("# ERROR mysql %d: %s" % ( e.args[0],e.args[1] ))

			#sys.exit()
		 	counter = 0				# restart the counter
		 	values_to_insert = []	# empty list 

	print "# INFO Runtime for creating output table: ", (datetime.now()- Cre_start_time)




	"############# CHECK DATA ###############"

	print "CHECK"
	# getting the org_ids in the output table
	cursor.execute("SELECT distinct org_id from " + output_table + " ;") # get the org_ids for all org found in the output table
	output_check2 =cursor.fetchall()
	org_output2 = [x[0] for x in output_check2]

	# closing database
	db.close()

	
	org_found = []
	org_notfound = []
	org_already = []

	for i in range(0,len(org_input)):
		if org_input[i] in org_output2:
			org_found.append(org_input[i])
		else:
			org_notfound.append(org_input[i])

	for i in range(0,len(org_output2)):
		if org_output2[i] not in org_input:
			org_already.append(org_output2[i])

	if len(org_input) == len(org_found):
		print "Success all organisms from the input is found in the output table. These are:", org_input
	else:
		print "ERROR: Not all input organisms are found in the output. Organisms not found:", org_notfound

	if org_already != []:
		print "Organisms found already in the output data, that hasn't been updated: ", org_already


	print "Made it to the end of SMG_proximity"
	print "# INFO Total Runtime for SMG_proximity: ", (datetime.now()-start_time)


	return 



