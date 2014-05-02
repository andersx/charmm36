import sys

file = sys.argv[1]
fh=open(file,'r')

for line in fh:
    vpp = []
    s = line.split()
    vpp.append(float(s[1])*0.001)
    vpp.append(float(s[2])*4.184)
    vpp.append(float(s[3])*4.184)
    vpp.append(float(s[4])*4.184)
    vpp.append(float(s[5])*4.184)
    vpp.append(float(s[6])*0.1)
    string = "%4s " % s[0]
    for kk in vpp:
        string += " %8.3f "% kk
    print string
fh.close
