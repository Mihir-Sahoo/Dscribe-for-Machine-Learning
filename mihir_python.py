# This is code to convert co-ordinates to xyz format and generate eigenvalues of Coulomb matrix
import h5py
import numpy as np 
import array
import math
import os
import pandas as pd
from io import StringIO
import csv
import ase.io
from numpy import linalg as LA
import csv
from pandas import *
from dscribe.descriptors import CoulombMatrix

with open("lf.txt","r") as data :
    emails = data.readlines()
#    print(emails)
#    x= len(emails)
#    print(x)
    batchs = int(len(emails)/24999)
#    print(batchs)
    for id,log in enumerate(emails):
#        print(id)
#        print(log)
        fileid = id/batchs
#        print(fileid)
        file=open("cl-{file}.txt".format(file=int(fileid)),'a+')
        file.write(log)

with open("lf.txt", "r") as fin:
    lines = fin.readlines()
f=open('ene.txt', 'a')
f.writelines("%s\n" % lines[i] for i in range(2,1324900,53))
f.close()
e = np.loadtxt(f"ene.txt", dtype=float)
#print(e)
#print(type(e))
#print(e.shape)
#n = m.reshape(4,3)
#print(n)
g = np.array(e)
#print(type(k))
#print(k)
dfg = pd.DataFrame(g)

for a in range (0,24999,1):
    b_file = open(f"cl-{a}.txt", "r")

    lines = b_file.readlines()
    b_file.close()

    del lines[0]  # delete first line
    del lines[0]  # delete 2nd line
    del lines[0]  # delete 3rd line
    new_file = open(f"cl-{a}.txt", "w+")

    for line in lines:
         new_file.write(line)

new_file.close()

for a in range (0,24999,1):
    m = np.loadtxt(f"cl-{a}.txt", dtype=float)
#    print(m)
#    print(type(m))
#    print(m.shape)
    s = ['C' for i in range(50)]
#    print(s)
    z = np.c_[s, m]
#    print(z)
#    print(type(z))   

    df = pd.DataFrame(z)
    df.to_csv( open(f'cl-{a}.csv', 'w'), index = False, header= False)
# df.to_csv('file.txt', index=False, header=False )
    csv.writer(open(f'cl-{a}.xyz', 'w+'), delimiter='\t').writerows(csv.reader(open(f'cl-{a}.csv')))

    src1=open(f"cl-{a}.xyz","r")
    fline=" \n"    #Prepending string
    oline=src1.readlines()
#Here, we prepend the string we want to on first line
    oline.insert(0,fline)
    src1.close()

    src1=open(f"cl-{a}.xyz","w")
    src1.writelines(oline)
    src1.close()


    src=open(f"cl-{a}.xyz","r")
    fline="50\n"    #Prepending string
    oline=src.readlines()
 #Here, we prepend the string we want to on first line
    oline.insert(0,fline)
    src.close()


#We again open the file in WRITE mode
    src=open(f"cl-{a}.xyz","w")
    src.writelines(oline)
    src.close()

#    print(src)

    atoms = ase.io.read(f'cl-{a}.xyz')

    cm = CoulombMatrix(n_atoms_max = 50, permutation="sorted_l2").create(atoms)
    cm1 = cm.reshape(50, 50)
#    print(cm1)
    df2 = pd.DataFrame(cm1)
    df2.to_csv( open(f'cl-cm-{a}.csv', 'w'), index = False, header= False)

    m1 = np.array(cm1)
    w, v = np.linalg.eig(m1)

    df3 = pd.DataFrame(w)
    df3_tr = df3.transpose()
    df3_tr.to_csv( open(f'cl-e.val{a}.csv', 'w'), index = False)

    df4 = pd.DataFrame(v)
    df4.to_csv( open(f'cl-e.vec{a}.csv', 'w'), index = False, header= False)
#    print(cm1.shape)

fout = open("e_val.csv", "a")
for line in open("cl-e.val0.csv"):
    fout.write(line)
for num in range (1, 24999, 1):
    f = open("cl-e.val"+str(num)+".csv")
    k = f.readlines()
    for line in k[1:]:
        fout.write(line)
    f.close()
fout.close()

abc = pd.read_csv("e_val.csv")
np.savetxt("foo2.csv", abc, delimiter=",", fmt ='%1.4f')
df9 = pd.read_csv("foo2.csv", header=None)
df9["50"] = "0"   # adding new column of value 0
df9["51"] = "0"   # adding new column of value 0
df9["52"] = "0"   # adding new column of value 0
df9["53"] = "0"   # adding new column of value 0
df9["54"] = "0"   # adding new column of value 0
df9.to_csv("final.csv", index= False)

aa = pd.read_csv("final.csv")
aa["Cluster-Id"] = [q for q in range(0,24999,1)]  # adding a new column named "cluster-Id" at the end with values 0-9
aa.to_csv("final-1.csv", index=False)

data = pd.read_csv("final-1.csv")


col_name = "Cluster-Id"
first_col = data.pop(col_name)


data.insert(0, col_name, first_col)

l = data

l.to_csv('a-1.csv', index = None)

s5 = pd.read_csv("a-1.csv")

headerlist = ["Cl-Id",  "eig0", "eig1", "eig2", "eig3", "eig4", "eig5", "eig6", "eig7", "eig8", "eig9", "eig10", "eig11", "eig12", "eig13", "eig14", "eig15", "eig16", "eig17", "eig18", "eig19", "eig20", "eig21", "eig22", "eig23", "eig24", "eig25", "eig26", "eig27", "eig28", "eig29", "eig130", "eig31", "eig32", "eig33", "eig34", "eig35", "eig36", "eig37", "eig38", "eig39", "eig40", "eig41", "eig42", "eig43", "eig44", "eig45", "eig46", "eig47", "eig48", "eig49", "eig50", "eig51", "eig52", "eig53", "eig54", "Energy"]

dif = pd.DataFrame(s5)
result = pd.concat([dif, dfg], axis=1)   # joining energy dataframe columnwise
result.to_csv( open(f'full.csv', 'w'), float_format ='%.4f', header = headerlist, index = False)
