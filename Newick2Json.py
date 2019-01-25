#!/usr/bin/python

#########
# This code allows generating the JSON files that are required by tippecanoe to build a complete Lifemap view
# It takes an input a tree in NHX format. 
# Attributes associated to the nodes and tips of the tree become the attributes of the json file. 
#########

import sys
import os
from argparse import ArgumentParser, FileType ##for options handling
import numpy as np
from ete3 import Tree
import time




def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return an open file handle

parser = ArgumentParser(description='Convert a tree in NHX format (newick extended format) into a collection of json files necessary for Lifemap visualization. These files will serve as input for "tippecanoe"')
parser.add_argument("-i", dest="filename", required=True, help="input tree file in NHX format", metavar="TreeFile", type=lambda x: is_valid_file(parser, x))
parser.add_argument('--lang', nargs='?', const='EN', default='EN', help='language chosen. FR for french, EN (default) for english', choices=['EN','FR'])

args = parser.parse_args()

sys.stdout.write("\nLoading tree... \r")
sys.stdout.flush()
t = Tree(args.filename) #read the input tree.
nbsp = len(t) ## get nb of tips 
sys.stdout.write("Loading tree... DONE [the tree has %d tips] \n" % nbsp)
sys.stdout.flush()

t.x = 6.0;
t.y = 9.660254-10.0;
t.alpha = 30.0;
t.ray = 30.0;
t.zoomview = np.ceil(np.log2(30/t.ray));
maxZoomView=0 ##

##FUNCTIONS
#getattr(t,n)

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
    NodeTipJS.write("      \"tippecanoe\": { \"maxzoom\" : %d, \"minzoom\" : %d },\n" % (32 if node.zoomview>26 else node.zoomview+5, 1 if node.zoomview < 3 else node.zoomview-2))
    NodeTipJS.write("      \"properties\": {\n")
    for f in node.features:
        NodeTipJS.write("        \"%s\":\"%s\",\n" % (f, cleanit(str(getattr(node, f)))))
    ## We add the Lifemap-specific info: zoomview, istip, isnode,
    NodeTipJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    NodeTipJS.write("        \"isnode\":\"%s\",\n" % (str(node.is_leaf()==False)))
    NodeTipJS.write("        \"istip\":\"%s\",\n" % (str(node.is_leaf()==True)))
    NodeTipJS.write("        \"nbdesc\":%d\n" % (node.nbdesc))
    NodeTipJS.write("      },\n") #close properties
    ## And finally the geometry
    NodeTipJS.write("      \"geometry\":{\n")
    NodeTipJS.write("        \"type\": \"Point\",\n")
    NodeTipJS.write("        \"coordinates\": [%.20f,%.20f]\n" % (node.x, node.y))
    NodeTipJS.write("      }\n") #close geometry
    NodeTipJS.write("    },\n") #close the whole feature

def writeGeojsonLines(node):
    BranchJS.write("    {\n");
    BranchJS.write("      \"type\": \"Feature\",\n")
##    BranchJS.write("      \"tippecanoe\": { \"maxzoom\" : %d, \"minzoom\" : %d },\n" % (32 if node.zoomview>26 else node.zoomview+5, 1 if node.zoomview < 3 else node.zoomview-2))
    BranchJS.write("      \"properties\": {\n")
    BranchJS.write("        \"zoom\":%d\n" % (node.zoomview))
    BranchJS.write("      },\n") #close properties
    BranchJS.write("      \"geometry\":{\n")
    BranchJS.write("        \"type\": \"LineString\",\n")
    BranchJS.write("        \"coordinates\": [ [ %.20f,%.20f ],[ %.20f,%.20f ] ]\n" % (node.up.x, node.up.y, node.x, node.y))
    BranchJS.write("      }\n") #close geometry
    BranchJS.write("    },\n") #close the whole feature

def writeGeojsonPolyg(node):
    polyg = HalfCircPlusEllips(node.x,node.y,node.ray,rad(node.alpha) + np.pi/2, rad(node.alpha) - np.pi/2, rad(node.alpha) + np.pi/2, 30)
    polygcenter = (np.mean(polyg[0]),np.mean(polyg[1]));
    PolygonJS.write("    {\n");
    PolygonJS.write("      \"type\": \"Feature\",\n")
    PolygonJS.write("      \"tippecanoe\": { \"maxzoom\" : %d, \"minzoom\" : %d },\n" % (32 if node.zoomview>26 else node.zoomview+5, 1 if node.zoomview < 3 else node.zoomview-2))
    PolygonJS.write("      \"properties\": {\n")
    for f in node.features:
        PolygonJS.write("        \"%s\":\"%s\",\n" % (f, cleanit(str(getattr(node, f)))))
    PolygonJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    PolygonJS.write("        \"isnode\":\"%s\",\n" % (str(node.is_leaf()==False)))
    PolygonJS.write("        \"istip\":\"%s\",\n" % (str(node.is_leaf()==True)))
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
    for f in node.features:
        CladeNameJS.write("        \"%s\":\"%s\",\n" % (f, cleanit(str(getattr(node, f)))))    
    CladeNameJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    CladeNameJS.write("        \"isnode\":\"%s\",\n" % (str(node.is_leaf()==False)))
    CladeNameJS.write("        \"istip\":\"%s\",\n" % (str(node.is_leaf()==False)))
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
    for f in node.features:
        RankNamesJS.write("        \"%s\":\"%s\",\n" % (f, cleanit(str(getattr(node, f)))))    
    RankNamesJS.write("        \"zoom\":%d,\n" % (node.zoomview))
    CladeNameJS.write("        \"nbdesc\":%d\n" % (node.nbdesc))
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



currsp = 0 ##count nb of species already visited

for n in t.traverse():
    special = 0
    n.dist=1.0
    tot = 0.0
    child = n.children
    ##NEW: we deal with singletons -->| 
    if ((len(child)==1)&(len(n)>1)):
        special=1
    if ((len(child)==1)&(len(n)==1)):
        special=2
    ## |<-- NEW
    for i in child:
        tot = tot + np.sqrt(len(i));
    nbdesc = len(n);
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
    if n.is_leaf():
        ##progress...
        currsp+=1
        sys.stdout.write("Tree traversal... " + str(np.ceil(currsp/nbsp*100)) + "% \r")
        sys.stdout.flush()
        ##

sys.stdout.write("Tree traversal... " + "DONE          " + "\n")
sys.stdout.flush()

sys.stdout.write("Closing output files... \r")
sys.stdout.flush()
sys.stdout.write("Closing output files... DONE\n")
sys.stdout.flush()

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

print("INFO: Lifemap requires at least %d zoom levels to visualize the whole tree.\n" % maxZoomView)

##we MAY need to add a branch before the root... 
