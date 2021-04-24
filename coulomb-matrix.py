# Python code for creating eigenvalues Cooulomb matrix Descriptors for machine learning model
# Created by Mihir Ranjan Sahoo(mrs10@iitbbs.ac.in)
# for any queries please contact me
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

# The pos variable contains the postion array (40001, 55, 3)
# The ene variable contains the energy array (40001) 

def read_data(temp):

    f = h5py.File(f"{temp}.h5", "r")

    pos = np.swapaxes(f["positions"],0,2)
    ene = f['energies']
    return (pos,ene)
    f.close()

pos, ene = read_data(30)

for a in range (0, 40001, 10000):
 m = np.array(pos[a,:,:])
 s = ['C' for i in range(55)]
 z = np.c_[s, m]
 print(a)
 print(z)
 print(ene[a])
 
# df = pd.DataFrame(pos[a,:,:])
 df = pd.DataFrame(z)
 df.to_csv( open(f'30-{a}.csv', 'w'), index = False, header= False)
# df.to_csv('file.txt', index=False, header=False )
 csv.writer(open(f'30-{a}.xyz', 'w+'), delimiter='\t').writerows(csv.reader(open(f'30-{a}.csv')))

 src1=open(f"30-{a}.xyz","r")
 fline=" \n"    #Prepending string
 oline=src1.readlines()
 #Here, we prepend the string we want to on first line
 oline.insert(0,fline)
 src1.close()

 src1=open(f"30-{a}.xyz","w")
 src1.writelines(oline)
 src1.close()


 src=open(f"30-{a}.xyz","r")
 fline="55\n"    #Prepending string
 oline=src.readlines()
 #Here, we prepend the string we want to on first line
 oline.insert(0,fline)
 src.close()


#We again open the file in WRITE mode
 src=open(f"30-{a}.xyz","w")
 src.writelines(oline)
 src.close()

 print(src)


 atoms = ase.io.read(f'30-{a}.xyz')

 cm = CoulombMatrix(n_atoms_max = 55, permutation="sorted_l2").create(atoms)
 cm1 = cm.reshape(55, 55)
 print(cm1)
 df2 = pd.DataFrame(cm1)
 df2.to_csv( open(f'30-cm-{a}.csv', 'w'), index = False, header= False)

 m1 = np.array(cm1)
 w, v = np.linalg.eig(m1)

 df3 = pd.DataFrame(w)
 df3_tr = df3.transpose()
 df3_tr.to_csv( open(f'30-e.val{a}.csv', 'w'), index = False)

 df4 = pd.DataFrame(v)
 df4.to_csv( open(f'30-e.vec{a}.csv', 'w'), index = False, header= False)
print(cm1.shape)

fout = open("e_val.csv", "a")
for line in open("30-e.val0.csv"):
    fout.write(line)
for num in range (10000, 40001, 10000):
    f = open("30-e.val"+str(num)+".csv")
    k = f.readlines()
    for line in k[1:]:
        fout.write(line)
    f.close()
fout.close()


#ab = pd.read_csv("e_val.csv")
arr1 = [ene[a] for a in range(0, 40001, 10000)]


abc = pd.read_csv("e_val.csv")
z = np.c_[abc, arr1]
#print(z)
np.savetxt("foo2.csv", z, delimiter=",", fmt ='%1.4f')
df9 = pd.read_csv("foo2.csv", header=None)
df9.to_csv("final.csv", index= False)

aa = pd.read_csv("final.csv")
aa["temp"] = "30"
aa.to_csv("final.csv", index=False)

data = pd.read_csv("final.csv")

data.head()

col_name = "temp"
first_col = data.pop(col_name)

data.head()

data.insert(0, col_name, first_col)

l = data.head()

l.to_csv('a-1.csv', index = None)

s5 = pd.read_csv("a-1.csv")
s3 = [ "30-"+str(a) for a in range(0, 40001, 10000)]

df = pd.DataFrame(s3)
df.to_csv('myfile.csv', header = ["str-id"], index = None)
print(s3)

s4 = pd.read_csv("myfile.csv")
z2=  np.c_[s4, s5]
print(z2.dtype)

#headerlist = ["str-Id", "Temp", ["eig"+str(b) for b in range(0,55)], "Enegry"]
headerlist = ["Str-Id", "Temp", "eig0", "eig1", "eig2", "eig3", "eig4", "eig5", "eig6", "eig7", "eig8", "eig9", "eig10", "eig11", "eig12", "eig13", "eig14", "eig15", "eig16", "eig17", "eig18", "eig19", "eig20", "eig21", "eig22", "eig23", "eig24", "eig25", "eig26", "eig27", "eig28", "eig29", "eig30", "eig31", "eig32", "eig33", "eig34", "eig35", "eig36", "eig37", "eig38", "eig39", "eig40", "eig41", "eig42", "eig43", "eig44", "eig45", "eig46", "eig47", "eig48", "eig49", "eig50", "eig51", "eig52", "eig53", "eig54", "Energy"]

dif = pd.DataFrame(z2)
dif.to_csv( open(f'final-10.csv', 'w'), float_format ='%.4f', header = headerlist, index = False)
