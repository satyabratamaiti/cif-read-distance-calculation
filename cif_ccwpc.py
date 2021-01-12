#!/usr/bin/env python
# coding: utf-8

import getopt
import sys
import math


class Atom:
    def __init__(self,atype,resid,resname,chain_id,coor_x,coor_y,coor_z):
        self.atype=atype
        self.resid=resid
        self.resname=resname
        self.chain_id=chain_id
        self.coor_x=coor_x
        self.coor_y=coor_y
        self.coor_z=coor_z
       
        
    def findDistance(self,another_atom):
        del_x = self.coor_x - another_atom.coor_x
        del_x = del_x * del_x
        del_y = self.coor_y - another_atom.coor_y
        del_y = del_y * del_y
        del_z = self.coor_z - another_atom.coor_z
        del_z = del_z * del_z
        val = math.sqrt(del_x + del_y + del_z)
        return val
    
  
       

        
    def details(self):
        return (self.atype + ":" + self.resid + ":" + self.resname + ":" + self.chain_id + ":" + str(self.coor_x) + ":" + str(self.coor_y) + ":" + str(self.coor_z))
        

def cif_file_operation(cfile,resid1,resid2,chain1,chain2):
    
    fp=open(cfile,'r')
    
    #attribute mapper generation
    attribute_mapper={}
    c=0
    for line in fp:
        if line.startswith("_atom_site."):
            attrb=line.strip().split(".")[-1]
            attribute_mapper.update( {attrb : c} )
            c=c+1
    
    
    atom1=None
    atom2=None
    atom3=None
    atom4=None
    
    
    
    fp=open(cfile,'r')
    for line in fp:
       
        if line.startswith("ATOM"):
            attrbs=line.split()
            atype=attrbs[attribute_mapper["label_atom_id"]].replace("\"","")
            resid=attrbs[attribute_mapper["label_seq_id"]]
            resname=attrbs[attribute_mapper["label_comp_id"]]
            chain_id=attrbs[attribute_mapper["auth_asym_id"]]
            coor_x=attrbs[attribute_mapper["Cartn_x"]]
            coor_y=attrbs[attribute_mapper["Cartn_y"]]
            coor_z=attrbs[attribute_mapper["Cartn_z"]]
            
            
            if atype=="N3" and resid==resid1 and chain_id==chain1:
                atom1=(Atom(atype,resid,resname,chain_id,float(coor_x),float(coor_y),float(coor_z)))
                
                
                
            if atype=="O2" and resid==resid2 and chain_id==chain2:
                atom2=(Atom(atype,resid,resname,chain_id,float(coor_x),float(coor_y),float(coor_z)))
                
            if atype=="O2" and resid==resid1 and chain_id==chain1:
                atom3=(Atom(atype,resid,resname,chain_id,float(coor_x),float(coor_y),float(coor_z)))
                
                
                
            if atype=="N3" and resid==resid2 and chain_id==chain2:
                atom4=(Atom(atype,resid,resname,chain_id,float(coor_x),float(coor_y),float(coor_z)))
    
    
    dist1 = atom1.findDistance(atom2)
    dist2 = atom3.findDistance(atom4)
    #print(dist1) 
    #print(dist2)
    return dist1, dist2
    #return dist2
                

#cif_file_operation("4u4u.cif","74","105","3","3")

import os.path

fp=open("BLCCWPC.csv",'r')
fw=open("distance.txt",'w')

for line_ in fp:
    line_1=line_.split(",")[-1].replace("\"","").replace("\n","")
    line_2=line_1.split("_")
    cif=line_2[0]+".cif"
    resid1=line_2[5]
    resid2=line_2[7]
    chain1=line_2[6]
    chain2=line_2[8]
    
    
    
    try :
        if os.path.exists(cif):
            dist = cif_file_operation(cif,resid1,resid2,chain1,chain2)
#            print (str(dist[0]) + "++")
#            print (str(dist[1]) + "--")
            if dist[0] <= 3.5:
                fw.write(cif + "\t" + str(dist[0]) + "\n")
            if dist[1] <= 3.5:
                fw.write(cif + "\t" + str(dist[1]) + "\n")
            
        
    except Exception :
        pass
    

fw.close()

