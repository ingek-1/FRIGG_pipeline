#!/usr/bin/python
# coding=utf-8

# importing modules
# to read and write csv files, sys - to take arguments
import csv, sys
import argparse
import MySQLdb as mdb
from datetime import datetime
import time


def SMG_family (clear_prox_fam, input_table_fam, output_table_fam):

	start_time =datetime.now()

	print "#---------------------------------------------------------------------------#"
	print "#                     Adding protein families to SM genes "
	print "#                          and genes in proximity. "
	print "#---------------------------------------------------------------------------#"

	# variables:
	intermediate_table = 'prox_fam' # intermediate table with clust+prox genes coupled to hfam
	hfam_table = 'hfam'
	count_fam = 'countfams' 


	# connecting to the database
	host = "XXX"
	user = "user_name"
	password = "***"
	dbname = "database"
	# connect to DB
	db = mdb.connect( host,user,password,dbname )
	cursor = db.cursor()



	"#################### CHECK STARTING DATA ##########################"

	if cursor.execute("SHOW TABLES LIKE \'" + input_table_fam + "\';"):
		print "check the input table exists: ", input_table_fam
	else:
		sys.exit("ERROR: The table " + input_table_fam +" does not exist and have to be created to use the function: SMG_family.")
	
	# save input org - to be used at the end to check input and output match
	cursor.execute("SELECT org_id, count(org_id) FROM " + input_table_fam +	" group by org_id;")
	input_check = cursor.fetchall()
	org_input = [x[0] for x in input_check]		# save org_id from input as a list in org_input 


	"###### Check if the intermediate table exists #######"
	if cursor.execute("SHOW TABLES LIKE \'" + intermediate_table + "\';"):
		if clear_prox_fam == 'FALSE':
			sys.exit("ERROR: The intermediate table ("+ intermediate_table +") already exists, to delete this and run again set -delete/ -d to TRUE.")
		else:
			print "The intermediate table ("+ intermediate_table +") already exists and will be will be deleted."
			query = "DROP TABLE "+ intermediate_table +";"
			
			try:
				cursor.execute(query)
				#cursor.executemany(stmt, values_to_insert)
				db.commit()
			except  mdb.Error, e:
				sys.exit("# ERROR mysql %d: %s" % ( e.args[0],e.args[1] ))

	"############### Coupling cluster genes from input to hfam plus counting number  #################"
	print "Creating MySQL table: " + intermediate_table
	# generating a table with the cluster genes + prox genes and the hfam and the number of hfam per organism
	query = ""
	query = ("CREATE TABLE " + intermediate_table +
	" SELECT ta.org_id, org_name, u_clust_id, clust_backbone, clust_type, ta.protein_id, T_hfam.hfam, count as count_proteins from "+ input_table_fam + " as ta\
	join " + hfam_table +" as T_hfam ON T_hfam.org_id =ta.org_id  AND ta.protein_id = T_hfam.protein_id \
	join "+ count_fam + " as countfams on countfams.hfam = T_hfam.hfam AND ta.org_id = countfams.org_id;")
	cursor.execute(query)
	print "# INFO Runtime for creating "+ intermediate_table +" table: ", (datetime.now()- start_time)


	"############ CHECK if output table is already found - and drop it if it is found #############"
	if cursor.execute("SHOW TABLES LIKE \'" +output_table_fam +"\';"):
		print "The output table " +output_table_fam +" already exists and will be will be deleted."
		query = "DROP TABLE " +output_table_fam +";"
			
		try:
			cursor.execute(query)
			db.commit()
		except  mdb.Error, e:
			sys.exit("# ERROR mysql %d: %s" % ( e.args[0],e.args[1] ))



	# generating table prox_fam_count with added counts for # genes in cluster, # genes pr hfam pr org and # genes outside cluster, # organism with this hfam

	"########### Create output table############"
	print "Creating output table in MySQL: " + output_table_fam
	# This is counting the number of genes 
	# 
	query = ""
	query = ("CREATE TABLE " + output_table_fam + 
				" SELECT intT.org_id, org_name, intT.u_clust_id, clust_backbone, clust_type, intT.protein_id, tb.hfam, clust_S_count, clust_count, hfam_count, count_anywhere, count_org, size_hfam FROM "+ intermediate_table +" as intT \
				left join ( \
					select *, hfam_count-clust_count as count_anywhere from ( \
						select hfam, org_id, count(distinct org_id, protein_id) as clust_count, count_proteins as hfam_count from "+ intermediate_table +" \
						group by hfam, org_id ) as ta) as tb \
				ON tb.hfam =intT.hfam AND tb.org_id = intT.org_id \
				LEFT JOIN ( \
					SELECT hfam, count(distinct org_id) as count_org FROM ( \
						select * from "+ intermediate_table +" \
						group by hfam, org_id ) as tc  \
					group by hfam) as te \
				ON te.hfam = intT.hfam \
				LEFT JOIN ( \
					SELECT org_id, u_clust_id, hfam, clust_S_count FROM ( \
						SELECT *, count( hfam) as clust_S_count FROM "+ intermediate_table +" \
						group by u_clust_id, hfam) as tf) as tg \
				ON tg.hfam = intT.hfam AND tg.u_clust_id = intT.u_clust_id \
				LEFT JOIN ( \
					SELECT hfam, count(distinct org_id, protein_id) as size_hfam FROM " + hfam_table + " \
					group by hfam) as th ON intT.hfam = th.hfam;")

	cursor.execute(query)


	##### COUNTS ######
	print "\n# INFO: Background of the data \n"
	# Number of organisms
	cursor.execute("SELECT count(distinct org_id) FROM " + output_table_fam + ";")
	number_orgs = cursor.fetchall()
	print "# Number of org found in data: ", number_orgs[0][0]
	
	# Number of genes
	cursor.execute("SELECT count(distinct org_id, protein_id) FROM " + output_table_fam + " ;")
	number_genes = cursor.fetchall()
	print "# Number of genes found in data: ", number_genes[0][0]

	# Number of hfam in the data
	cursor.execute("SELECT count(distinct hfam) FROM " + output_table_fam + " ;")
	number_hfam = cursor.fetchall()
	print "# Number of hfam found in data: ", number_hfam[0][0]

	# Number of backbone genes in the data
	cursor.execute("SELECT count(distinct clust_backbone) FROM " + output_table_fam + " WHERE clust_backbone = protein_id;")
	number_bb = cursor.fetchall()
	print "# Number of backbone genes found in the data: ", number_bb[0][0]

	# Data overview - per hfam
	print"\n# INFO Dataoverview per hfam"
	# core hfams - where the hfam is found in every organism in the set
	query = ""
	query = ("SELECT count(*), group_concat(hfam)  FROM ( \
				SELECT *, count(distinct org_id) as c_org FROM ( \
					SELECT *, count( distinct org_id, protein_id) as ce FROM " + output_table_fam + "  \
					group by hfam, org_id) as ta \
				group by hfam) as tb \
			where c_org = (SELECT count(distinct org_id) from " + output_table_fam + ");")
	cursor.execute(query)
	core_hfams = cursor.fetchall()

	print "# INFO: core hfam - hfams found in every organism in the set"
	print "The number of core hfams: ", core_hfams[0][0]
	print "Core hfams: ", core_hfams[0][1]
	print "The percentage of core hfams of total hfams: ", float(core_hfams[0][0])/float(number_hfam[0][0])*100
	
	# unique hfams
	query = ""
	query = ("SELECT count(*) FROM ( \
				SELECT *, count(distinct org_id, protein_id) as count_spsp FROM " + output_table_fam + " \
				where clust_count = 1 AND hfam_count = 1 \
				group by hfam) as ta \
			where count_spsp = 1;")
	cursor.execute(query)
	N_unique_hfams = cursor.fetchall()

	print "# INFO: unique hfam - hfams found only in one organism in the set"
	print "The number of unique hfams: ", N_unique_hfams[0][0]	
	print "The percentage of unique hfams of total hfams: ", float(N_unique_hfams[0][0])/float(number_hfam[0][0])*100

	# Data overview - per org
	print"\n# INFO Dataoverview per organism"
	# hfams per org
	cursor.execute("SELECT org_id, count(distinct hfam) FROM " + output_table_fam + "  group by org_id;")
	hfam_per_org = cursor.fetchall()
	print "# INFO: hfams per organism"
	for i in range(0,len(hfam_per_org)):
		print "Organism: ", hfam_per_org[i][0], ", # Number distinct hfams per org: ", hfam_per_org[i][1]

	# paralogs hfam
	query = ""
	query = ("SELECT org_id, count(distinct hfam) FROM ( \
				SELECT * FROM " + output_table_fam + "  \
				where hfam_count > 1) as ta \
			group by org_id;")
	cursor.execute(query)
	paralogs_per_org = cursor.fetchall()
	print "# INFO: Number of paralog hfam per organism (hfam found more than once in the org)"
	for i in range(0,len(paralogs_per_org)):
		print "Organism: ", paralogs_per_org[i][0], ", # Number distinct hfams per org: ", paralogs_per_org[i][1]


	# Species specific paralogs - found in only in one org but more than once
	query = ""
	query = ("SELECT org_id, count(*)  FROM ( \
				SELECT *, count(distinct org_id, protein_id) as count_hfam_events FROM " + output_table_fam + "  \
				where hfam_count > 1 \
				group by hfam) as ta \
			where count_hfam_events = 1 \
			group by org_id;")
	cursor.execute(query)
	count_SpSp_Paralog = cursor.fetchall()

	total_SpSp_paralogs = 0
	print "# INFO: Species specific paralogs - found in only one org but more than once in that org"
	for i in range(0,len(count_SpSp_Paralog)):
		print "Organism: ", count_SpSp_Paralog[i][0], ", # Species specific paralogs: ", count_SpSp_Paralog[i][1]
		total_SpSp_paralogs = total_SpSp_paralogs + count_SpSp_Paralog[i][1]
	print "Total species specific paralogs: ", total_SpSp_paralogs

	# Species specific singletons - only found once in one org
	query = ""
	query = ("SELECT org_id, count(*) FROM ( \
				SELECT *, count(distinct org_id, protein_id) as count_spsp FROM " + output_table_fam + "  \
				where clust_count = 1 AND hfam_count = 1 \
				group by hfam) as ta \
			where count_spsp = 1 \
			group by org_id;")
	cursor.execute(query)
	count_singletons = cursor.fetchall()

	total_singletons = 0
	print "# INFO: Singletons - species specific hfam only found once in one org"
	for i in range(0,len(count_singletons)):
		print "Organism: ", count_singletons[i][0], ", # Singletons: ", count_singletons[i][1]
		total_singletons = total_singletons + count_singletons[i][1]
	print "The total number of singletons: ", total_singletons


	# closing database
	db.close()

	