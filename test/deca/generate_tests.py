from fragbuilder import Peptide
aas = "FYWDENQHKRGASTCVLIMP"

for l1 in aas:

    # Define sequence
    seq = "A" + l1 * 10 + "A"
    print seq

    # Create peptide
    pep = Peptide(seq, nterm="charged", cterm="charged")

    # Sample (and set) backbone angles for residue 2
    for i in range(10):

        pep.sample_bb_angles(2 + i)
        print pep.get_bb_angles(2 + i)

        # Sample (and set) chi angles for residue 2
        try:
            pep.sample_chi_angles(2 + i)
            print pep.get_chi_angles(2 + i)
        except:
            pass

    # MMFF94 optimization
    pep.optimize()

    # Write files
    filename = seq + ".pdb"
    pep.write_pdb(filename, QUIET=False)


