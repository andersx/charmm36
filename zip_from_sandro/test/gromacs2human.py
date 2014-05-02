import sys

j2cal = 0.239005736
f=sys.argv[1]

x = []
fh = open(f,'r')
for line in fh:
    try:
        v = [j2cal*float(p) for p in line.split()]
        x.extend(v)
    except:
        pass
fh.close()

print " Bond    Angles  Dihedr   I-Dihe   CMAP     LJ       COUL     TOT  "
print "%-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f" % (x[0],x[1],x[2], x[3],x[4], x[5]+x[7], x[6]+ x[8], x[9])


g=sys.argv[2]

y = []
fh = open(g,'r')
for line in fh:
    try:
        v = [float(p) for p in line.split()[2:]]
        y.extend(v)
    except:
        pass
fh.close()

print " Bond    Angles  Dihedr   I-Dihe   CMAP     LJ      Coul     TOT  "
print "%-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f" % (y[3],y[4]+y[5], y[6],y[7], y[8],y[11], y[12], y[0])

h=sys.argv[3]

z = []
fh = open(h,'r')
for line in fh:
    try:
        v = [float(p) for p in line.split()[1:]]
        z.extend(v)
    except:
        pass
fh.close()
ss = sum(z)
idhe = 0.0
cmap = 0.0
print " Bond    Angles  Dihedr   I-Dihe   CMAP     LJ      Coul     TOT  "
print "%-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f" % (z[4],z[2], z[3],idhe, cmap,z[1], z[0], ss)
