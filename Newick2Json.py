#!/usr/bin/python

#########
#WORK VERSION OF THE NEW WAY OF GENERATING LIFEMAP-RELATED FILES.
#BETTER VERSION MIGHT COME SOON (less ncbi-specific)
#########
import sys
import os
from argparse import ArgumentParser, FileType ##for options handling
import numpy as np
from ete3 import Tree
#from ete3 import NCBITaxa
import psycopg2 ##for postgresql connection
#import cPickle as pickle
from getTrees_fun import getTheTrees

parser = ArgumentParser(description='Convert a tree in NHX format (newick extended format) into a collection of json files necessary for Lifemap visualization.')
parser.add_argument('group', help='Group to look at. Can be 1,2 or 3 for Archaea, Eukaryotes and Bacteria respectively', choices=['1','2','3'])
parser.add_argument('start', help='index of the first node met in the tree', type=int)
parser.add_argument('--lang', nargs='?', const='EN', default='EN', help='Language chosen. FR for french, EN (default) for english', choices=['EN','FR'])
parser.add_argument('--updatedb', nargs='?', const='True', default='True', help='Should the NCBI taxonomy db be updated ?', choices=['True','False'])
parser.add_argument('--simplify', nargs='?', const='True', default='False', help='Should the tree be simplified by removing environmental and unindentified species?', choices=['True','False'])

args = parser.parse_args()
#print args

##update db (if requested?)
def updateDB():
	print 'Updating databases...'
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
	print "Tree HAS BEEN simplified"
	finalSize = len(arbre)
	diffInSize = (initialSize-finalSize)
	print  str(diffInSize) + " tips have been removed (" + str(round(float(diffInSize)/float(initialSize)*100, 2)) + '%)'
	print  "FINAL TREE SIZE: " + str(finalSize)
	return arbre


if (args.updatedb=='True'):
	updateDB()

##get arguments
groupnb = args.group ##will be written

T = getTheTrees()

#print sys.argv[1];
starti = args.start;
print "Downloading tree..."
if (groupnb=="1"):
	#with open('ARCHAEA.pkl', 'rb') as input:
	#t = pickle.load(input)
	t = T['2157'].detach()
	#SIMPLIFY THE TREE IF REQUESTED	
	if args.simplify=="True":
		t = simplify(t)
	#t = Tree("ARCHAEA")
	print "Archaeal tree loaded..."
	##and we save it
	t.write(outfile="ARCHAEA", features = ["name", "taxid"], format_root_node=True)
	t.x = 6.0;
	t.y = 9.660254-10.0;
	t.alpha = 30.0;
	t.ray = 10.0;
	starti = starti;
if (groupnb=="2"):
	# with open('EUKARYOTES.pkl', 'rb') as input:
	# 	t = pickle.load(input)
	t = T['2759'].detach()
	#SIMPLIFY THE TREE IF REQUESTED	
	if args.simplify=="True":
		t = simplify(t)
	print "Eukaryotic tree loaded"
	t.write(outfile="EUKARYOTES", features = ["name", "taxid"], format_root_node=True)
	t.x = -6.0;
	t.y = 9.660254-10.0;
	t.alpha = 150.0;
	t.ray = 10.0;
	starti = starti;
if (groupnb=="3"):
	# with open('BACTERIA.pkl', 'rb') as input:
	# 	t = pickle.load(input)
	t = T['2'].detach()
	#SIMPLIFY THE TREE IF REQUESTED	
	if args.simplify=="True":
		t = simplify(t)
	print "Bacterial tree loaded"
	t.write(outfile="BACTERIA", features = ["name", "taxid"], format_root_node=True)
	t.x = 0.0;
	t.y = -11.0;
	t.alpha = 270.0;
	t.ray = 10.0;
	starti = starti;

t.zoomview = np.ceil(np.log2(30/t.ray));


#specis and node ids
nbsp = len(t)
# spid = starti
# ndid = starti + nbsp
# rootnb = ndid+1
maxZoomView=0

##FUNCTIONS
def rad(deg):
    return((deg*np.pi)/180);
def halfCircle(x,y,r,start,end,nsteps):
    rs = np.linspace(start,end,num=nsteps)
    xc = x+r*np.cos(rs)
    yc = y+r*np.sin(rs)
    return(xc,yc)
def ellipse(x,y,r, alpha, nsteps):
    start=0
    end=np.pi+start
    rs = np.linspace(start,end,num=nsteps)
    a = r
    b = float(r)/6 ##Change this value to change the shape of polygons. This controls how flat is the elliptic side of the polygon. The other side is always a half cricle. 
    xs = a*np.cos(rs)
    ys = b*np.sin(rs)
    ##rotation
    xs2 = x+(xs*np.cos(alpha)-ys*np.sin(alpha))
    ys2 = y+(xs*np.sin(alpha)+ys*np.cos(alpha))
    return(xs2,ys2)
def HalfCircPlusEllips(x,y,r,alpha, start, end,nsteps):
        circ = halfCircle(x,y,r,start,end, nsteps)
        elli = ellipse(x,y,r,alpha,nsteps)
        return (np.concatenate((circ[0], elli[0])),np.concatenate((circ[1], elli[1])))
    ##write json for search

def cleanit(str): #this function cleans the strings that need to be.
    return(str.replace('"','\\"'))

def writeGeojsonNode(node):
    NodeTipJS.write("    {\n");
    NodeTipJS.write("      \"type\": \"Feature\",\n")
    NodeTipJS.write("      \"properties\": {\n")
    NodeTipJS.write("        \"taxid\":\"%s\",\n" % (node.taxid))
    NodeTipJS.write("        \"sci_name\":\"%s\",\n" % (cleanit(node.sci_name)))
    NodeTipJS.write("        \"common_name_en\":\"%s\",\n" % (cleanit(node.common_name)))
    NodeTipJS.write("        \"common_name_fr\":\"%s\",\n" % (cleanit(node.common_name_FR)))
    NodeTipJS.write("        \"authority\":\"%s\",\n" % (cleanit(node.authority)))
    NodeTipJS.write("        \"synonym\":\"%s\",\n" % (cleanit(node.synonym)))
    NodeTipJS.write("        \"rank_en\":\"%s\",\n" % (node.rank))
    NodeTipJS.write("        \"rank_fr\":\"%s\",\n" % (node.rank_FR))
    NodeTipJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    NodeTipJS.write("        \"isnode\":\"%d\",\n" % (node.is_leaf==False))
    NodeTipJS.write("        \"istip\":\"%d\",\n" % (node.is_leaf==True))
    NodeTipJS.write("        \"nbdesc\":%d\n" % (node.nbdesc))
    NodeTipJS.write("      },\n") #close properties
    NodeTipJS.write("      \"geometry\":{\n")
    NodeTipJS.write("        \"type\": \"Point\",\n")
    NodeTipJS.write("        \"coordinates\": [%.20f,%.20f]\n" % (node.x, node.y))
    NodeTipJS.write("      }\n") #close geometry
    NodeTipJS.write("    },\n") #close the whole feature

def writeGeojsonLines(node):
    BranchJS.write("    {\n");
    BranchJS.write("      \"type\": \"Feature\",\n")
    BranchJS.write("      \"properties\": {\n")
    BranchJS.write("        \"zoom\":%d\n" % (node.zoomview))
    BranchJS.write("      },\n") #close properties
    BranchJS.write("      \"geometry\":{\n")
    BranchJS.write("        \"type\": \"LineString\",\n")
    BranchJS.write("        \"coordinates\": [ [ %.20f,%.20f ],[ %.20f,%.20f ] ]\n" % (node.up.x, node.up.y, node.x, node.y))
    BranchJS.write("      }\n") #close geometry
    BranchJS.write("    },\n") #close the whole feature

def writeGeojsonPolyg(node):
    sci_name = node.sci_name
    sci_name = sci_name.replace('"','\\"')
    common_name = node.common_name
    common_name = common_name.replace('"','\\"')
    common_name_FR = node.common_name_FR
    common_name_FR = common_name_FR.replace('"','\\"')
    ##new attributes
    authority = node.authority
    authority = authority.replace('"','\\"')
    synonym = node.synonym
    synonym = synonym.replace('"','\\"')
    polyg = HalfCircPlusEllips(node.x,node.y,node.ray,rad(node.alpha) + np.pi/2, rad(node.alpha) - np.pi/2, rad(node.alpha) + np.pi/2, 30)
    polygcenter = (np.mean(polyg[0]),np.mean(polyg[1]));
    PolygonJS.write("    {\n");
    PolygonJS.write("      \"type\": \"Feature\",\n")
    PolygonJS.write("      \"properties\": {\n")
    PolygonJS.write("        \"taxid\":\"%s\",\n" % (node.taxid))
    PolygonJS.write("        \"sci_name\":\"%s\",\n" % (cleanit(node.sci_name)))
    PolygonJS.write("        \"common_name_en\":\"%s\",\n" % (cleanit(node.common_name)))
    PolygonJS.write("        \"common_name_fr\":\"%s\",\n" % (cleanit(node.common_name_FR)))
    PolygonJS.write("        \"authority\":\"%s\",\n" % (cleanit(node.authority)))
    PolygonJS.write("        \"synonym\":\"%s\",\n" % (cleanit(node.synonym)))
    PolygonJS.write("        \"rank_en\":\"%s\",\n" % (node.rank))
    PolygonJS.write("        \"rank_fr\":\"%s\",\n" % (node.rank_FR))
    PolygonJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    PolygonJS.write("        \"isnode\":\"%d\",\n" % (node.is_leaf==False))
    PolygonJS.write("        \"istip\":\"%d\",\n" % (node.is_leaf==True))
    PolygonJS.write("        \"nbdesc\":%d\n" % (node.nbdesc))
    PolygonJS.write("      },\n") #close properties
    PolygonJS.write("      \"geometry\":{\n")
    PolygonJS.write("        \"type\": \"Polygon\",\n")
    PolygonJS.write("        \"coordinates\": [ [            ")
    for i in range(0,59):
        PolygonJS.write(" [ %.20f, %.20f ],"  % (polyg[0][i], polyg[1][i]))
    PolygonJS.write(" [ %.20f, %.20f ]\n"  % (polyg[0][0], polyg[1][0])) #to close the polygon
    PolygonJS.write("             ] ]\n")
    PolygonJS.write("      }\n") #close geometry
    PolygonJS.write("    },\n") #close the whole feature
    #We write also the centre of the polygon as a point feature, in a separate file
    RankNamesJS.write("    {\n");
    CladeNameJS.write("      \"type\": \"Feature\",\n")
    CladeNameJS.write("      \"properties\": {\n")
    CladeNameJS.write("        \"taxid\":\"%s\",\n" % (node.taxid))
    CladeNameJS.write("        \"sci_name\":\"%s\",\n" % (sci_name))
    CladeNameJS.write("        \"common_name_en\":\"%s\",\n" % (common_name))
    CladeNameJS.write("        \"common_name_fr\":\"%s\",\n" % (common_name_FR))
    CladeNameJS.write("        \"authority\":\"%s\",\n" % (authority))
    CladeNameJS.write("        \"synonym\":\"%s\",\n" % (synonym))
    CladeNameJS.write("        \"rank_en\":\"%s\",\n" % (node.rank))
    CladeNameJS.write("        \"rank_fr\":\"%s\",\n" % (node.rank_FR))
    CladeNameJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    CladeNameJS.write("        \"isnode\":\"%d\",\n" % (node.is_leaf==False))
    CladeNameJS.write("        \"istip\":\"%d\",\n" % (node.is_leaf==True))
    CladeNameJS.write("        \"nbdesc\":%d\n" % (node.nbdesc))
    CladeNameJS.write("      },\n") #close properties
    CladeNameJS.write("      \"geometry\":{\n")
    CladeNameJS.write("        \"type\": \"Point\",\n")
    CladeNameJS.write("        \"coordinates\": [%.20f,%.20f]\n" % (polygcenter[0], polygcenter[1]))
    CladeNameJS.write("      }\n") #close geometry
    CladeNameJS.write("    },\n") #close the whole feature
    #We write also the rank names along a line.
    RankNamesJS.write("    {\n");
    RankNamesJS.write("      \"type\": \"Feature\",\n")
    RankNamesJS.write("      \"properties\": {\n")
    RankNamesJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    RankNamesJS.write("        \"rank_en\":\"%s\",\n" % (node.rank))
    RankNamesJS.write("        \"rank_fr\":\"%s\"\n" % (node.rank_FR))
    RankNamesJS.write("      },\n") #close properties
    RankNamesJS.write("      \"geometry\":{\n")
    RankNamesJS.write("        \"type\": \"LineString\",\n")
    RankNamesJS.write("        \"coordinates\": [ \n")
    RankNamesJS.write("            [ %.20f, %.20f ]"  % (polyg[0][i], polyg[1][i]))
    for i in range(36,45):
        RankNamesJS.write(", [ %.20f, %.20f ]"  % (polyg[0][i], polyg[1][i]))
    RankNamesJS.write("\n          ] \n")    
    RankNamesJS.write("      }\n") #close geometry
    RankNamesJS.write("    },\n") #close the whole feature



def TerminateFiles(file):
    ##remove unwanted last character(,) of json files
    consoleexex = 'head -n -1 ' + file + ' > temp.txt ; mv temp.txt '+ file;
    os.system(consoleexex);
    temp = open(file, "a");
    temp.write("    }\n  ]\n}")
    temp.close()

NodeTipJS = open("Tipsnodes.json", "w")
BranchJS = open("Branches.json", "w")
CladeNameJS = open("CladesNames.json", "w")
RankNamesJS = open("RankNames.json", "w")
PolygonJS = open("Polygons.json", "w")

NodeTipJS.write("{\n  \"type\": \"FeatureCollection\",\n  \"features\": [\n");
BranchJS.write("{\n  \"type\": \"FeatureCollection\",\n  \"features\": [\n");
CladeNameJS.write("{\n  \"type\": \"FeatureCollection\",\n  \"features\": [\n");
RankNamesJS.write("{\n  \"type\": \"FeatureCollection\",\n  \"features\": [\n");
PolygonJS.write("{\n  \"type\": \"FeatureCollection\",\n  \"features\": [\n");

print "Tree traversal..."
for n in t.traverse():
    special = 0
    n.dist=1.0
    tot = 0.0
    child = n.children
    ##NEW  -->|
    if ((len(child)==1)&(len(n)>1)):
        special=1
    if ((len(child)==1)&(len(n)==1)):
        special=2
    ## |<-- NEW
    for i in child:
        tot = tot + np.sqrt(len(i));
    nbdesc = len(n);
    ##remove special chars in names
    ####IF --LANG IS SET TO FR, WE CHGANGE HERE THE RANK AND COMMON NAMES
    # if (args.lang=='FR'):
    #     n.common_name = n.common_name_FR
    #     n.rank = n.rank_FR   
    #####OK
    ## we create a 'long' common name. the common name going to db is only the first of the list 
    ## n.common_name_long = ', '.join(n.common_name)
    n.common_name = n.common_name[0] if len(n.common_name)>0 else ""
    n.common_name = n.common_name.replace("'","''");
    n.common_name_FR = n.common_name_FR[0] if len(n.common_name_FR)>0 else ""
    n.common_name_FR = n.common_name_FR.replace("'","''");
    n.rank = n.rank.replace("'","''");
    n.rank_FR = n.rank_FR.replace("'","''");

    n.sci_name = n.sci_name.replace("'","''")
    #add parenthesis to the common name
    if n.common_name!='':
        n.common_name = "(" + n.common_name + ")"
    n.nbdesc = nbdesc;
    nbsons = len(child);
    angles = [];
    ray = n.ray;
    for i in child:
        #i.ang = 180*(len(i)/float(nbdesc))/2;
        i.ang = 180*(np.sqrt(len(i))/tot)/2; #using sqrt we decrease difference between large and small groups
        angles.append(i.ang);
        if (special==1):
            i.ray = ray-(ray*20)/100
        else:
            if (special==2):
                i.ray = ray-(ray*50)/100
            else:
                i.ray = (ray*np.sin(rad(i.ang))/np.cos(rad(i.ang)))/(1+(np.sin(rad(i.ang))/np.cos(rad(i.ang))));
        i.dist = ray - i.ray;
    ang = np.repeat(angles, 2);
    ang = np.cumsum(ang);
    ang = ang[0::2];
    ang = [i-(90-n.alpha) for i in ang];
    cpt = 0
    for i in child:
        i.alpha = ang[cpt];
        i.x = n.x + i.dist*np.cos(rad(i.alpha));
        i.y = n.y + i.dist*np.sin(rad(i.alpha));
        i.zoomview = np.ceil(np.log2(30/i.ray))
        if i.zoomview <= 0:
            i.zoomview = 0
        if maxZoomView<i.zoomview:
            maxZoomView = i.zoomview
        cpt = cpt+1;
    #we write node info
    writeGeojsonNode(n)
    writeGeojsonPolyg(n)
    if n.is_root()==False:
        writeGeojsonLines(n)

NodeTipJS.close()
BranchJS.close()
CladeNameJS.close()
RankNamesJS.close()
PolygonJS.close()
TerminateFiles("Tipsnodes.json")
TerminateFiles("Branches.json")
TerminateFiles("CladesNames.json")
TerminateFiles("RankNames.json")
TerminateFiles("Polygons.json")

print "Tree traversal... DONE (first one)"

#################################
#       WRITE JSON FILES        #
#################################
jsonfile = 'TreeFeatures'+groupnb+'.json';
json = open(jsonfile, "w");
json.write("[\n");
def writejsonNode(node):
    sci_name = node.sci_name
    sci_name = sci_name.replace('"','\\"')
    common_name = node.common_name
    common_name = common_name.replace('"','\\"')
	##new attributes
    authority = node.authority
    authority = authority.replace('"','\\"')
    synonym = node.synonym
    synonym = synonym.replace('"','\\"')
    json.write("  {\n");
    json.write("    \"taxid\":\"%s\",\n" % (node.taxid))
    json.write("    \"sci_name\":\"%s\",\n" % (sci_name))
    json.write("    \"common_name\":\"%s\",\n" % (common_name))
    #new functions: add authority and synonym
    json.write("    \"authority\":\"%s\",\n" % (authority))
    json.write("    \"synonym\":\"%s\",\n" % (synonym))
    #end
    json.write("    \"rank\":\"%s\",\n" % (node.rank))
    json.write("    \"zoom\":\"%d\",\n" % (node.zoomview+4))
    json.write("    \"nbdesc\":\"%d\",\n" % (node.nbdesc))
    json.write("    \"all\":\"%s | %s | %s | %s\",\n" % (sci_name, common_name, node.rank, node.taxid))
    json.write("    \"coordinates\": [%.20f,%.20f],\n" % (node.y, node.x))
    json.write("    \"lat\": \"%.20f\",\n" % (node.y))
    json.write("    \"lon\": \"%.20f\"\n" % (node.x))
    json.write("  },\n")




print "Tree traversal 2... "
##LAST LOOP TO write coords of polygs and JSON file
for n in t.traverse():
    #save all trees to disk
#    out="trees/" + str(n.taxid) + ".tre";
#    n.write(outfile=out, features=["taxid","sci_name","common_name","rank"]);
    ##we finish writing in the database here.
    if n.is_root()==False:
        ndid = ndid+1
        writeosmWays(n, ndid)
    if n.is_leaf()==False:
        indexes = np.linspace(ndid + 1,ndid+63,num=63)
        writeosmpolyg(n, indexes)
        ndid = ndid+63
    writejsonNode(n)
    writeGeojsonLines(n)
##after this, node.nbgenomes should be ok.
print "Tree traversal 2... DONE "
conn.commit();

##we add the way from LUCA to the root of the subtree 
ndid=ndid+1
command = "INSERT INTO lines (id, branch, zoomview, ref, way) VALUES(%d,'TRUE', '4','%s',ST_Transform(ST_GeomFromText('LINESTRING(0 -4.226497, %.20f %.20f)', 4326), 900913));" % (ndid, groupnb, t.x, t.y);
cur.execute(command);
conn.commit()



print "DONE!"
out = open("tempndid", "w")
out.write("%d" % ndid) ##we store the max id so that we start from there for next group.
print ("Max zoom view : %d" % (maxZoomView));

        
##remove unwanted last character(,) of json file
json.close()
consoleexex = 'head -n -1 ' + jsonfile + ' > temp.txt ; mv temp.txt '+ jsonfile;
os.system(consoleexex);
json = open(jsonfile, "a");
json.write("\t}\n]\n")
json.close()









def writeosmWays(node, id):
    #Create branch names
    Upsci_name = node.up.sci_name;
    Upcommon_name = node.up.common_name;
    Downsci_name = node.sci_name;
    Downcommon_name = node.common_name;
    left = Upsci_name +  " " + Upcommon_name;
    right = Downsci_name + " " + Downcommon_name;
    if (node.x >= node.up.x): #we are on the right
        wayName = "\u2190  " + left + "     -     " + right + "  \u2192"
    else: #we are on the left
        wayName = "\u2190  " + right + "     -     " + left + "  \u2192"
    command = "INSERT INTO lines (id, branch, zoomview, ref, name, way) VALUES(%d,'TRUE',%d,'%s',E'%s',ST_Transform(ST_GeomFromText('LINESTRING(%.20f %.20f, %.20f %.20f)', 4326), 900913));" % (id, node.zoomview, groupnb, wayName, node.up.x, node.up.y, node.x, node.y);
    cur.execute(command);
    ##conn.commit();
