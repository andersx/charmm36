import sys


ph_file = sys.argv[1]
ch_file = sys.argv[2]

ph_fh = open(ph_file)
for line in ph_fh:
    stripped = line.strip()
    split = stripped.split()
    if("EEF1" in split and len(split)==4):
        v_ph = float(split[2])
ph_fh.close()


ch_fh = open(ch_file)
for line in ch_fh:
    if("ENER EXTERN>" in line):
        v_ch = float(line.split()[5])
ch_fh.close()

print v_ph, v_ch
