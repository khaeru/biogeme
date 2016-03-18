# script to extract a sample from a data set
# Author: Michel Bierlaire
# Sun Nov 30 22:00:11 2014

import sys
import random


filename = str(sys.argv[1])
f = open(str(sys.argv[1]), 'r')
data_filename = 'sampled_'+filename
d = open(data_filename,'w') 

percent = 1.0/float(str(sys.argv[2]))
print("Sample ",sys.argv[2],"% of the data")
first = 1
total = 0 ;
accepted = 0 ;
for line in f:
    if first:
        first = 0
        print(line,end='',file=d) 
    else:
        total += 1 
        r = random.random() 
        if (random.random() <= percent):
            accepted += 1
            print(line,end='',file=d) 
f.close()
d.close()
print("Actual sample: ", 100.0 * accepted/total, "%") 
print("Sampled data file: ",data_filename)
