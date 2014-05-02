import sys


ph_file = sys.argv[1]
gr_file = sys.argv[2]
ch_file = sys.argv[3]

ph_fh = open(ph_file)
for line in ph_fh:
    stripped = line.strip()
    split = stripped.split()
    if("Total" in split and len(split)==2):
        v_ph = float(split[1])
ph_fh.close()

gr_fh = open(gr_file)
for line in gr_fh:
    if("@" not in line and "#" not in line and len(line.split())>0):
        v_gr = float(line.split()[1]) + float(line.split()[2])
gr_fh.close()


ch_fh = open(ch_file)
for line in ch_fh:
    if("ENER EXTERN>" in line):
        v_ch = float(line.split()[3])
ch_fh.close()

print v_ph, v_gr, v_ch*4.186
