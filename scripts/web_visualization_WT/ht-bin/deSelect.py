#!/usr/bin/env python
import pdb,sys,os
import json
import sqlite3
import cgi,cgitb,json
from collections import defaultdict


cgitb.enable(format="text")

form=cgi.FieldStorage()
cells=form.getvalue("cells1")
cells2=form.getvalue("cells2")


print("Content-Type: text/html\n") 

#cells="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17"
#cells=""
#cells2="21,22,23,24,25,26,27,28,29,30,31,32,33,34"
#print(cells)
#print(cells)

cells=[] if cells==None else [item for item in cells.split(",") if item!=""]
cells2=[] if cells2==None else [item for item in cells2.split(",") if item!=""]

#pdb.set_trace()


if len(cells)==0 and len(cells2)==0:
	sys.exit(0)
elif len(cells)==0 and len(cells2)!=0:
	cells2=[int(item) for item in cells2]
	cells=cells2;
	ucells=[item for item in range(len(rows)) if item not in cells]
elif len(cells)!=0 and len(cells2)==0:
	cells=[int(item) for item in cells]
	ucells=[item for item in range(len(rows)) if item not in cells]
else:
	cells=[int(item) for item in cells]
	cells2=[int(item) for item in cells2]
	ucells=cells2

# use TF-DNA info to add regulators
# get TF-DNA dictionary
# TF->DNA
def getdTD(tfDNA):
	dTD = {}
	with open(tfDNA, 'r') as f:
		tfRows = f.readlines()
		tfRows = [item.strip().split() for item in tfRows]
		for row in tfRows:
			itf = row[0].upper()
			itarget = row[1].upper()
			if itf not in dTD:
				dTD[itf] = [itarget]
			else:
				dTD[itf].append(itarget)
	return dTD

def getdDT(dTD):
	gene_tf_dict = defaultdict(lambda: [])
	for key, val in dTD.items():
		for v in val:
			gene_tf_dict[v.upper()] += [key.upper()]
	return gene_tf_dict


tfDNA_file = "ht-bin/tfDNA_predicted_100.txt.update"  # may need to change to absolute path
dTD = getdTD(tfDNA_file)
dDT = getdDT(dTD)
# TF-DNA end

Z=[]
ct=0

conn=sqlite3.connect("ht-bin/ex1.db")
cur=conn.cursor()
cur.execute("SELECT * FROM extable")
# rows=cur.fetchall()  # here all genes are read in; only suitable for small datasets

for i in cur:
	ei=json.loads(i[-1])
	sei=[ei[item] for item in cells]
	uei=[ei[item] for item in ucells]
	sei=sum(sei)/len(sei)
	uei=sum(uei)/len(uei)
	ediff=sei-uei

	tfs = ','.join(dDT[i[1]])

	if ediff>0.6:
		#pdb.set_trace()
		Z.append([i[1],uei,sei,ediff,tfs])
	ct+=1
	#print(ct)

Z.sort(key=lambda x:x[3],reverse=True)
print(json.dumps(Z))
