import sys
import numpy as N


file = sys.argv[1]
fh = open(file)
data = []

for line in fh:
    if(len(line.split())==4):
        vec = [float(x) for x in line.split()]
        data.append(vec)
data = N.array(data)
fh.close()
print "Correlation Coefficient"
print N.corrcoef([data[:,0],data[:,2]])[0,1]
print N.corrcoef([data[:,1],data[:,3]])[0,1]
print "Fit"
print N.polyfit(data[:,0],data[:,2],1)[0]
print N.polyfit(data[:,1],data[:,3],1)[0]
