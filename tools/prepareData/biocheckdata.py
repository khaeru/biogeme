# script to check that a biogeme estimation datafile has no error
# Authors: Ricardo Hurtubia and Michel Bierlaire
# Wed Mar 14 08:30:18 2012

from string import *
import sys
import random
import re

from math import *

file = str(sys.argv[1])
argc = len(sys.argv)
if (argc == 1):
	print("Syntax: ",sys.argv[0]," mydata.dat [rowsMerged=1]")
elif (argc == 2):
	mergedRows = 1
else:
	mergedRows = int(sys.argv[2])
	if (mergedRows < 1):
		print("The number of merged rows is incorrect: ",mergedRows)
		sys.exit() 
	else:
		print("Each ",mergedRows, "consecutive rows are merged") ;
Input1 = open(str(sys.argv[1]), 'r')

print("Check if the file ",str(sys.argv[1])," is complying with biogeme's requirements.") 
print("Reading data")
data_1 = {}
tab_data = []
len_data=0
for line in Input1:
	dataLine = str.rstrip(line)
	data_1[line] = re.compile('\s').split(dataLine)
	row = data_1[line]
#	print(row)
	tab_data.append(row)
	if (len_data == 0):
		print(len(row), "headers: ",row) 
	len_data=len_data+1
	if len_data == 500000:	
		print("500000 lines read")
	if len_data == 1000000:	
		print("1000000 lines read")
	if len_data == 1500000:	
		print("1500000 lines read")	
	if len_data == 2000000:	
		print("more than 2000000 lines read")
	
print(" ")

print(len_data, "lines")
print(len(tab_data[0]), "columns")

err=0
for i in range (len(tab_data)):
# for i in range (3):
	if i>0:
		if (mergedRows * len(tab_data[i])) >len(tab_data[0]):
			print("Length [",i,"]: ", len(tab_data[i]))
			print("Length [0]: ", len(tab_data[0]))
			print("error in line", i+1," (more columns (",mergedRows * len(tab_data[i]),") than headers (",len(tab_data[0]),"))")
			err = 1
		if (mergedRows * len(tab_data[i])) <len(tab_data[0]):
			print("error in line", i+1," (less columns (",mergedRows * len(tab_data[i]),") than headers (",len(tab_data[0]),"))")
			err = 1
		
		for j in range(len(tab_data[i])):
			
			x=tab_data[i][j]
			try:
				y=float(x)
				
			except:	
				print("error in line", i+1, " (column ",j+1," contains text: ",x,")")
				err=1

if err==0:
	print("data check finalized, no errors.")
else:
	print("The file does not comply with biogeme's requirements")
