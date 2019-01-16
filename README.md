# Lifemap-VectorTiles

In this version of Lifemap, we strictly separate the retrieval of the NCBI taxonomy, and the generation of all formatted data necessary for Lifemap. 
This means that Lifemap can then be generated for any newick formatted tree, as long as the tree is in the correct format.

##get_NCBI_taxo_FR.py
This function downloads the latest NCBI taxonomy, get french translations from INPN (zip file must already be present in the taxo folder), and outputs an NHX tree file for either the whole NCBI taxonomy, or for each subtree indicated in the command line with --taxid.

For help with this funciton, type
 get_NCBI_taxo_FR.py -h

##Newick2json.py
This function takes as input a tree in NHX format (like the one generated by get_NCBI_taxo_FR.py) and outputs all json files necessary for Lifemap to be created (see after). 

For help with this funciton, type
 Newick2json.py -h