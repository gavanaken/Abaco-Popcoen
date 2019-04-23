#! /usr/bin/env python
import os
import re
import sys
import subprocess
import math
import random
import time
import socket
import hashlib
import mdtraj
import numpy as np
import pylab
import tempfile

import popcoen

import PoPE_features
import Popcoen_Reliability_number

print(sys.version_info)



def Popcoen_all___from_structure_to_prediction(fn):

    if os.path.isfile(fn) != True:
        sys.exit("ERROR: tempfile does not exist or not readable. Filename:",fn)
    testit = open(fn)
    testitT = testit.read()
    testit.close()
    Natoms = testitT.count("ATOM") + testitT.count("HETATM")
    if Natoms < 10:
	if Natoms == 0:
	    Natoms="no"
	else:
	    Natoms="not enough"
        sys.exit("ERROR: pdb-file contains %s ATOMS. Without ATOMs no Popcoen!" %Natoms)


    try:
	centroid = mdtraj.load_pdb(fn)
    except:
	sys.exit("ERROR: pdb-file could not be loaded (using mdtraj). Please check for proper pdb format.")



    print "INFO: Start Popcoen Step1: Feature measurement."
    x, AAnames = PoPE_features.measure_all_features(centroid)



    # measure reliability number
    reliability_number_lambda = Popcoen_Reliability_number.calculate_reliability_number(x)



    # Drop unnecessary columns and convert AA-code to binary code. 
    X_is_2D = np.array( [False, True, True, True, True, True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False] )
    X_is_important = np.array( [True, False, False, True, True, True, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, True, True, True, False, False] )

    x_AAs = x[:, np.logical_and(X_is_2D, X_is_important)]
    AAs = np.arange(21)
    x = np.hstack((x[:, np.logical_and(np.logical_not(X_is_2D), X_is_important)], np.equal(x_AAs[:, :, None], AAs[None, None, :]).reshape(x_AAs.shape[0], x_AAs.shape[1] * AAs.shape[0])))



    print "INFO: Start Popcoen Step2: Neural network evaluation."
    Si_hat = popcoen.predict(x)


    # Compression error rescaling
    rs_a = 0.49252789 # old inverse rescaling gave: 1./1.69138373
    rs_b = 1.10525533 # old inverse rescaling gave: 2.18795683 / 1.69138373
    Si_predict = rs_a * Si_hat + rs_b


    S_PC = Si_predict.sum()

    print
    print "-----------------------"
    print "  Popcoen prediction:  "
    print "-----------------------"
    print
    print "Partial entropies S_i  (i.e. configurational entropy of single residue) in k_B:"
    assert len(Si_predict) == len(AAnames)
    for out in zip(AAnames, Si_predict):
	AAn, Sip = out
	print "Si-predict:", '{:>4}'.format(str(AAn.resSeq)), AAn.name, ":", '{: .3f}'.format(Sip)
    print
    print "S_PC = sum Si =", '{: .3f}'.format(S_PC)
    print
    print "Reliability-number lambda = %.3lf" % reliability_number_lambda
    print
    print "Note: The configurational entropy of a protein equals S_PC up to a chain-dependent additive constant."
    print "      See Popcoen publication for details. We plan to include the constant to Popcoen Version 2."
    print "      The constant cancels for entropy differences between distinct conformations of given chain."
    print
    print "      Popcoen entropies are reliable for lambda > ~0.5 and very reliable for lambda > 0.8."
    print
    print "Thank you for using Popcoen. For questions/recommendations/bug-reports, please contact Martin Goethe."
    print
    print "INFO: Stop. No more Popcoen."





if __name__ == "__main__":

    # expect to get pdb file content piped into
    pdbContentPipedin = sys.stdin.read()
    print "I recieved data piped into popcoenMain program."

    # save to tmp file, load via mdtraj, rm tmpfile
    t = tempfile.NamedTemporaryFile(mode='w+b', bufsize=-1, suffix='.pdb', prefix='tmpfile_mdtrajInput_', dir="./")
    t.write(pdbContentPipedin)
    t.seek(0)

    pdbfilename = t.name
    Popcoen_all___from_structure_to_prediction(pdbfilename)

