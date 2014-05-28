from string import split

def read_charm(filename):

    f = open(filename, "r")

    e = dict()

    for line in f.readlines():

        if "ENER INTERN>" in line:
            e["bond"]       = float(split(line)[2])
            e["angle"]      = float(split(line)[3])
            e["urey"]       = float(split(line)[4])
            e["torsion"]    = float(split(line)[5])
            e["imptor"]     = float(split(line)[6])


        if "ENER CROSS>" in line:
            e["cmap"]       = float(split(line)[2])

        if "ENER EXTERN>" in line:
            e["vdw"]        = float(split(line)[2])
            e["coul"]       = float(split(line)[3])
            e["solv"]       = float(split(line)[5])


    f.close()

    return e

def read_phaistos(filename):

    f = open(filename, "r")

    e = dict()

    for line in f.readlines():

        if "bond-stretch" in line and "kcal" in line:
            e["bond"]       = round(float(split(line)[3]),5)

        if "angle-bend E" in line and "kcal" in line:
            e["angle"]       = round(float(split(line)[3]),5)
        if "urey" in line and "kcal" in line:
            e["urey"]       = round(float(split(line)[3]),5)
        if "torsion" in line and "kcal" in line:
            e["torsion"]       = round(float(split(line)[3]),5)
        if "imptor" in line and "kcal" in line:
            e["imptor"]       = round(float(split(line)[3]),5)
        if "CMAP" in line and "kcal" in line:
            e["cmap"]       = round(float(split(line)[3]),5)
        if "vdW-total" in line and "kcal" in line:
            e["vdw"]       = round(float(split(line)[3]),5)
        if "Coul-total" in line and "kcal" in line:
            e["coul"]       = round(float(split(line)[3]),5)
        if "EEF1" in line and "kcal" in line:
            e["solv"]       = round(float(split(line)[3]),5)

    f.close()

    return e

aas = "fywdenqhkrgastcvlimp"

tolerance = 0.000011


for l1 in aas:
    for l2 in aas:

        seq = "a" + l1 + l2 + "a"

        e1 = read_charm("tetra/clean_" + seq + ".charmm")
        e2 = read_phaistos("tetra/clean_" + seq + ".phaistos")

        for key in e1.keys():

            diff = e1[key] - e2[key]

            if diff > tolerance:
                print seq, key, e1[key] - e2[key]


for seq in ["apo_lfabp", "ci2", "ff_domain", "hr4660b",
            "protein_g", "rhodopsin", "top7", "wr73"]:

    e1 = read_charm("proteins/clean_" + seq + ".charmm")
    e2 = read_phaistos("proteins/clean_" + seq + ".phaistos")
    
    for key in e1.keys():
    
        diff = e1[key] - e2[key]
    
        if diff > tolerance:
            print seq, key, e1[key] - e2[key]


for seq in ["clean_AAAAAAAAAAAA",
            "clean_ACCCCCCCCCCA",
            "clean_ADDDDDDDDDDA",
            "clean_AEEEEEEEEEEA",
            "clean_AFFFFFFFFFFA",
            "clean_AGGGGGGGGGGA",
            "clean_AHHHHHHHHHHA",
            "clean_AIIIIIIIIIIA",
            "clean_AKKKKKKKKKKA",
            "clean_ALLLLLLLLLLA",
            "clean_AMMMMMMMMMMA",
            "clean_ANNNNNNNNNNA",
            "clean_APPPPPPPPPPA",
            "clean_AQQQQQQQQQQA",
            "clean_ARRRRRRRRRRA",
            "clean_ASSSSSSSSSSA",
            "clean_ATTTTTTTTTTA",
            "clean_AVVVVVVVVVVA",
            "clean_AWWWWWWWWWWA",
            "clean_AYYYYYYYYYYA",]:
    e1 = read_charm("deca/" + seq + ".charmm")
    e2 = read_phaistos("deca/" + seq + ".phaistos")
    
    for key in e1.keys():
    
        diff = e1[key] - e2[key]
    
        if diff > tolerance:
            print seq, key, e1[key] - e2[key]
