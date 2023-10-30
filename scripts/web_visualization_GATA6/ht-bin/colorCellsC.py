#!/usr/bin/env python
import pdb,sys,os
import cgi,cgitb,json
import sqlite3


cgitb.enable()
form=cgi.FieldStorage()
geneinput=form.getvalue("cellGeneName")
geneinput=geneinput.upper()

#geneinput="NKX2-1"
# read in expression
conn=sqlite3.connect("ht-bin/ex1.db")
cur=conn.cursor()
cur.execute("SELECT * FROM extable WHERE name=?",(geneinput,))

rows=cur.fetchall()
Z=json.loads(rows[0][-1])

print("Content-Type: text/html\n")   
print(json.dumps(Z))
