#!/usr/bin/python

import os
from argparse import ArgumentParser, FileType ##for options handling
from ete3 import Tree


parser = ArgumentParser(description='This function retrieves the lastest NCBI taxonomy for all or a subset of subtrees, and adds french translations. It returns trees in Extended Newick Format (NHX)')
parser.add_argument('--updatedb', nargs='?', const='True', default='True', help='Should the NCBI taxonomy db be updated? Default is True', choices=['True','False'])
parser.add_argument('--simplify', nargs='?', const='True', default='False', help='Should the tree be simplified by removing environmental and unindentified species? Default is False', choices=['True','False'])
parser.add_argument('--taxid', nargs='+', default='1', help='List of taxids for which to return a tree, space-separated. If none is provided, the whole NCBI taxonomy tree is returned and the tree file is called "All.tre". If one or more is provided, the ouput files are named "*taxid*.tre". As a reminder, taxid for Bacteria=2, Eucaryotes=2759 and Archaea=2157.')

args = parser.parse_args()

def updateDB():
	print('Updating databases...')
	os.system("wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -N")
	os.system("tar xvzf taxdump.tar.gz -C taxo/")
	#unzip taxref NOTE: taxref has to be downloaded by hand beforehand ?
	os.system("unzip -o taxo/TAXREF_INPN_v11.zip -d taxo/")

def simplify(arbre):
	initialSize = len(arbre)
	for n in arbre.traverse():
		if (n.is_leaf()==True) and (n.rank=='no rank'):
			n.detach()
		else:
			if ('Unclassified' in n.sci_name) or ('unclassified' in n.sci_name) or ('uncultured' in n.sci_name) or ('Uncultured' in n.sci_name) or ('unidentified' in n.sci_name) or ('Unidentified' in n.sci_name) or ('environmental' in n.sci_name) or ('sp.' in n.sci_name):
				n.detach()
	print("Tree HAS BEEN simplified")
	finalSize = len(arbre)
	diffInSize = (initialSize-finalSize)
	print(str(diffInSize) + " tips have been removed (" + str(round(float(diffInSize)/float(initialSize)*100, 2)) + '%)')
	print("FINAL TREE SIZE: " + str(finalSize))
	return arbre

def getTheTrees(): 
	##DOWNLOAD taxdump and store in taxo folder
	##DOWNLOAD TAXREF BY HAND! and put it in taxo/

	class Trans:
		def __init__(self):
			self.common_name_FR = []


	print("Getting french translations...")
	TRANS = {} ##translations in french
	with open("taxo/TAXREFv11.txt") as f:  
		for line in f:
			sciname = line.split("\t")[14]
			comnameFR = line.split("\t")[19]
			if (sciname not in TRANS and line.split("\t")[19]!=''):
				TRANS[sciname] = Trans()
			if (line.split("\t")[19]!=''):
				TRANS[sciname].common_name_FR.append(comnameFR)

	#get translation of ranks
	print("Getting rank names in french...")
	RANKS = {}
	with open("taxo/ranks_FR.txt") as f:  
		for line in f:
			rank_en = line.split("\t")[0]
			rank_fr = line.split("\t")[1].rstrip() ##to remove \n
			RANKS[rank_en] = rank_fr


	class Taxid:
		def __init__(self):
			self.sci_name = ""
			self.authority = ""
			self.synonym = ""
#			self.common_name = ""
			self.common_name = []
#			self.common_name_FR = ""
			self.common_name_FR = []

	cpt = 0
	cptfr = 0
	ATTR = {} ##here we will list attribute of each species per taxid
	print("Reading NCBI taxonomy...")
	with open("taxo/names.dmp") as f:  
		for line in f:		
			taxid = line.split("|")[0].replace("\t","")
			tid_val = line.split("|")[1].replace("\t","")
			tid_type = line.split("|")[3].replace("\t","")
			if (taxid not in ATTR):
				ATTR[taxid] = Taxid()
			if (tid_type=="scientific name"):
				ATTR[taxid].sci_name = tid_val
				#and get translation in french (if any)
				if tid_val in TRANS:
					ATTR[taxid].common_name_FR = TRANS[tid_val].common_name_FR
					cptfr += 1
			if (tid_type=="authority"):
				if (ATTR[taxid].authority!=""):
					ATTR[taxid].authority = ATTR[taxid].authority + ", " + tid_val
				else:
					ATTR[taxid].authority = tid_val
			if (tid_type=="synonym"):
				if (ATTR[taxid].synonym!=""):
					ATTR[taxid].synonym = ATTR[taxid].synonym + ", " + tid_val
				else:
					ATTR[taxid].synonym = tid_val
			if (tid_type=="common name"):
				cpt +=1
				ATTR[taxid].common_name.append(tid_val)
				# if (ATTR[taxid].common_name!=""):
				# 	ATTR[taxid].common_name = ATTR[taxid].common_name + ", " + tid_val
				# else: 
				# 	ATTR[taxid].common_name = tid_val


	T = {}

	###New gettrees
	filepath = 'taxo/nodes.dmp'  
	print("Building the NCBI taxonomy tree...")
	with open(filepath) as fp:  
		first_line = fp.readline() ## remove the 1 | 1 edge
		for line in fp:
			dad = line.split("|")[1].replace("\t","")
			son = line.split("|")[0].replace("\t","")
			rank = line.split("|")[2].replace("\t","")
			if (dad not in T):
				T[dad] = Tree()
				T[dad].name = dad
#				T[dad].rank = rank
#				T[dad].rank_FR = RANKS[rank]
				T[dad].taxid = dad
				T[dad].sci_name = ATTR[dad].sci_name
				T[dad].common_name = ATTR[dad].common_name
				T[dad].synonym = ATTR[dad].synonym
				T[dad].authority = ATTR[dad].authority
				T[dad].common_name_FR = ATTR[dad].common_name_FR
			if (son not in T):
				T[son] = Tree()
				T[son].name = son
				T[son].rank = rank
				T[son].rank_FR = RANKS[rank]
				T[son].taxid = son
				T[son].sci_name = ATTR[son].sci_name
				T[son].common_name = ATTR[son].common_name
				T[son].synonym = ATTR[son].synonym
				T[son].authority = ATTR[son].authority
				T[son].common_name_FR = ATTR[son].common_name_FR
			else:
				if (hasattr(T[son], 'rank')==False):
					T[son].rank = rank
					T[son].rank_FR = RANKS[rank]
			T[dad].add_child(T[son])
	return T


if (args.updatedb=='True'):
	updateDB()

T = getTheTrees() #Get the Whole tree, for real.

if (len(args.taxid)==1):
	print("Writing the whole NCBI taxonomy tree to All.tre...")
	tout=T[args.taxid]
	tout.write(outfile = "All.tre", features = ["name", "taxid", "sci_name","common_name","rank", "authority","synonym","common_name_FR", "rank_FR"], format_root_node=True)
else:
	print("Writing trees...")
	for i in args.taxid:
		out= i + ".tre"	
		print("  " + out)
		t=T[i]
		t.write(outfile = out, features = ["name", "taxid", "sci_name","common_name","rank", "authority","synonym","common_name_FR", "rank_FR"], format_root_node=True)

