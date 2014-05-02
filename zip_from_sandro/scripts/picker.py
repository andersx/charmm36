import sys


ph_file = sys.argv[1]
gr_file = sys.argv[2]

ph_fh = open(ph_file)
for line in ph_fh:
    stripped = line.strip()
    split = stripped.split()
    if("Total" in split and "14" in split):
        print split[2],
    if("Total" in split and "sr" in split):
        print split[2],
ph_fh.close()

gr_fh = open(gr_file)
for line in gr_fh:
    if("@" not in line and "#" not in line and len(line.split())>0):
        print line.split()[1],
        print line.split()[2]
gr_fh.close()
