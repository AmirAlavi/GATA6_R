#!/usr/bin/env python
import sqlite3,sys,pdb
from File import *
import numpy as np
import json

name = sys.argv[1]
name_db = name + '.db'

conn=sqlite3.connect(name_db)
c=conn.cursor()
sql_table="""
CREATE TABLE IF NOT EXISTS extable (
	id integer PRIMARY KEY,
	name text,
	expression text
);
"""

c.execute(sql_table)
# the original expression input (same format as SCDIFF input)
ex=np.array(TabFile(sys.argv[2]).read("\t"))
ex=ex.transpose()

for i in range(4,len(ex)):
	gi=ex[i][0].upper()  # make sure gene name is all upper-case in the database
	ei=str([float(item) for item in ex[i][1:]])
	c.execute("INSERT INTO extable (name,expression) VALUES (?,?)", [gi,ei])
conn.commit()

