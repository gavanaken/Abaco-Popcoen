#! /usr/bin/env python
import os
import sys
import pickle
import mdtraj
import numpy as np
import math
import itertools
import random
import collections

import pylab




def measure_residue_contact_numbers(centroid, rcs = [2,4,6,8,10,12,14,16,18,20,22]):

    assert centroid.n_frames == 1

    gibweiter1, gibweiter2 =  mdtraj.compute_contacts(centroid, contacts="all", scheme="ca", periodic=False)
    # factor 10 for angstroem. [0] for only frame
    ContactMap = 10. * mdtraj.geometry.squareform(gibweiter1, gibweiter2)[0]

    CN = []
    for rc in rcs:
	# We do not normalize this. This is learned by neural network
	#Vol = 4. * np.pi/3 * rc*rc*rc
	#print Vol
	D = ContactMap < rc
	CN.append( D.sum(axis=1) )
	#print CN
    CN = np.transpose(CN)
    CN = CN.astype(float)

    return(CN, rcs)





# -----------------------------------------------------------------------------------------------





def measure_dist_CA_MassCenter(centroid):

    assert centroid.n_frames == 1

    MassCenter = mdtraj.compute_center_of_mass(centroid)[0]

    CAs    = [a.index for a in centroid.topology.atoms if a.name == 'CA']
    CA_pos = centroid.xyz[0, CAs, :]
    CA_pos = CA_pos.astype(float)
    # factor 10 for ang
    dist = [ 10. * np.linalg.norm(X - MassCenter) for X in CA_pos ]

    dist = np.array(dist)
    return(dist)



def measure_Geometry(centroid):

    assert centroid.n_frames == 1

    dist = measure_dist_CA_MassCenter(centroid)
    #print dist

    # we use all heavy atoms for r_g, not only the CA atoms
    #def filter(x):
    #    if x == "CA":
    #        return(1.0)
    #    else:
    #        return(0.0)
    #atomMASS = np.array( [filter(str(a.name)) for a in centroid.topology.atoms] )
    atomMASS = None
    r_g = 10. * ( mdtraj.compute_rg(centroid, masses=atomMASS)[0] )
    r_g = np.array( [r_g for i in range(len(dist)) ] )
    #print "rg",r_g

    Geometry = np.column_stack((r_g * r_g , dist * dist))

    return (Geometry)





# -----------------------------------------------------------------------------------------------





def initialize_StandardSASA_values_Naccess(ShowInfo = True):


    NACCESS_file_content = """STANDARD ACCESSIBILITES FOR PROBE 1.40 AND RADII vdw.radii 
ATOM S   2  ALA  107.95   0.0  69.41   0.0   0.00   0.0  69.41   0.0  38.54   0.0  71.38   0.0  36.58   0.0
ATOM S   2  CYS  134.28   0.0  96.75   0.0   0.00   0.0  96.75   0.0  37.53   0.0  97.93   0.0  36.35   0.0
ATOM S   2  ASP  140.39   0.0  48.00   0.0  54.69   0.0 102.69   0.0  37.70   0.0  49.24   0.0  91.15   0.0
ATOM S   2  GLU  172.25   0.0  59.10   0.0  75.64   0.0 134.74   0.0  37.51   0.0  60.29   0.0 111.96   0.0
ATOM S   2  PHE  199.48   0.0 164.11   0.0   0.00   0.0 164.11   0.0  35.37   0.0 165.25   0.0  34.23   0.0
ATOM S   2  GLY   80.10   0.0  32.33   0.0   0.00   0.0  32.33   0.0  47.77   0.0  37.55   0.0  42.55   0.0
ATOM S   2  HIS  182.88   0.0  96.01   0.0  51.07   0.0 147.08   0.0  35.80   0.0  97.15   0.0  85.73   0.0
ATOM S   2  ILE  175.12   0.0 137.96   0.0   0.00   0.0 137.96   0.0  37.16   0.0 139.14   0.0  35.98   0.0
ATOM S   2  LYS  200.81   0.0 115.38   0.0  47.92   0.0 163.30   0.0  37.51   0.0 116.57   0.0  84.24   0.0
ATOM S   2  LEU  178.63   0.0 141.12   0.0   0.00   0.0 141.12   0.0  37.51   0.0 142.31   0.0  36.32   0.0
ATOM S   2  MET  194.15   0.0 156.64   0.0   0.00   0.0 156.64   0.0  37.51   0.0 157.84   0.0  36.32   0.0
ATOM S   2  ASN  143.94   0.0  44.98   0.0  61.26   0.0 106.24   0.0  37.70   0.0  46.23   0.0  97.72   0.0
ATOM S   2  PRO  136.13   0.0 119.90   0.0   0.00   0.0 119.90   0.0  16.23   0.0 120.95   0.0  15.19   0.0
ATOM S   2  GLN  178.50   0.0  51.03   0.0  89.96   0.0 140.99   0.0  37.51   0.0  52.22   0.0 126.28   0.0
ATOM S   2  ARG  238.76   0.0  76.60   0.0 124.65   0.0 201.25   0.0  37.51   0.0  77.80   0.0 160.97   0.0
ATOM S   2  SER  116.50   0.0  46.89   0.0  31.22   0.0  78.11   0.0  38.40   0.0  48.55   0.0  67.95   0.0
ATOM S   2  THR  139.27   0.0  74.54   0.0  27.17   0.0 101.70   0.0  37.57   0.0  75.72   0.0  63.55   0.0
ATOM S   2  VAL  151.44   0.0 114.28   0.0   0.00   0.0 114.28   0.0  37.16   0.0 115.47   0.0  35.97   0.0
ATOM S   2  TRP  249.36   0.0 187.67   0.0  23.60   0.0 211.26   0.0  38.10   0.0 189.67   0.0  59.69   0.0
ATOM S   2  TYR  212.76   0.0 135.35   0.0  42.03   0.0 177.38   0.0  35.38   0.0 136.50   0.0  76.26   0.0"""


    # Load Standard values used in Naccess
    StandardSASA_values_Naccess = dict()

    T = NACCESS_file_content.split("\n")

    assert T[0] == "STANDARD ACCESSIBILITES FOR PROBE 1.40 AND RADII vdw.radii "
    T = T[1:]
    for l in T:
	tmp = l[:]
	while tmp.find("  ") != -1:
	    tmp = tmp.replace("  "," ")
	tmp = tmp.split()
	assert tmp[0] == "ATOM"
	RES = tmp[3]
	StandardSASA = float(tmp[4])
	StandardSASA_values_Naccess[RES] = StandardSASA
    del T,l,tmp,RES,StandardSASA

    # if an AA type is not in list, use mean value for SASA
    mean_val = round( np.mean( [StandardSASA_values_Naccess[X] for X in StandardSASA_values_Naccess] ), 3) + 0.000888888888888

    # ADD weird AAs by their brothers
    MoDEL_AAs_Substitution = []
    MoDEL_AAs_Substitution.append(["AR0","ARG"])
    MoDEL_AAs_Substitution.append(["ASH","ASP"])
    MoDEL_AAs_Substitution.append(["CYM","CYS"])
    MoDEL_AAs_Substitution.append(["CYX","CYS"])
    MoDEL_AAs_Substitution.append(["GLH","GLU"])
    MoDEL_AAs_Substitution.append(["HID","HIS"])
    MoDEL_AAs_Substitution.append(["HIE","HIS"])
    MoDEL_AAs_Substitution.append(["HIP","HIS"])
    MoDEL_AAs_Substitution.append(["LYN","LYS"])
    # For SASA, Selenomethionine has same properties as Methionine
    MoDEL_AAs_Substitution.append(["MSE","MET"])


    for aat in MoDEL_AAs_Substitution:
	StandardSASA_values_Naccess[aat[0]] = StandardSASA_values_Naccess[aat[1]]
    del MoDEL_AAs_Substitution


    # convert to defaultdict, so that unknown keys are given mean-SASA value
    StandardSASA_values_Naccess = StandardSASA_values_Naccess.items()
    StandardSASA_values_Naccess = collections.defaultdict(lambda: mean_val, StandardSASA_values_Naccess)

    if ShowInfo:
	print "# INFO: Loaded StandardSASA_values_Naccess defaultdict.  KeyError gives mean value",mean_val,"- Defined Keys are:", StandardSASA_values_Naccess.items()

    return(StandardSASA_values_Naccess)











StandardSASA_values_Naccess = "not_initialized_yet"

def measure_SASA_percentage(centroid, Show_INFOS = False):

    assert centroid.n_frames == 1
    global StandardSASA_values_Naccess

    if str(StandardSASA_values_Naccess) == "not_initialized_yet":
	StandardSASA_values_Naccess = initialize_StandardSASA_values_Naccess(ShowInfo = Show_INFOS)

    # factor 100 to come from pm^2 tp ang^2
    SASA = 100. * mdtraj.shrake_rupley(centroid, probe_radius=0.14, n_sphere_points=960, mode="residue")[0]
    total_SASA = SASA.sum()

    # standard values for normalization
    resnames = [str(a.name) for a in centroid.topology.residues]
    assert len(SASA) == len(resnames)
    StdVal = np.array( [StandardSASA_values_Naccess[X] for X in resnames] )

    SASA_nmd = SASA / StdVal

    if Show_INFOS:  # show test
	for X in zip(SASA, StdVal, SASA_nmd):
	    print "%.2lf %.2lf %.2lf" % X


    assert SASA_nmd.shape == (centroid.n_residues,)
    SASA_back = np.empty((centroid.n_residues,2))
    SASA_back[:,0] = SASA_nmd     # first col is SASA percentage of residue
    SASA_back[:,1] = total_SASA   # secnd col is always total_SASA in ang^2

    return(SASA_back)





# -----------------------------------------------------------------------------------------------





AAtype_dict = ""

def residue_surroundings(centroid, Show_INFOS = False):

    assert centroid.n_frames == 1
    global AAtype_dict

    if len(AAtype_dict) == 0:   # if dictionary not initialized yet

	Basic_AAs = ['LEU', 'ALA', 'GLY', 'VAL', 'LYS', 'GLU', 'SER', 'ASP', 'THR', 'ILE', 'ASN', 'ARG', 'PRO', 'PHE', 'GLN', 'TYR', 'HIS', 'MET', 'CYS', 'TRP']
	assert len(Basic_AAs) == 20
	Basic_AAs = dict( [[X[1], X[0]+1] for X in enumerate(Basic_AAs)] )
	AAtype_dict = Basic_AAs
	# Hence labels 1,2,3,...,19,20 are used for Basic_AAs

	# the following rare subspecies are not distinguished.
	# We treat them as if they were a specific basic AA
	# The asociated dictionary is the following, which can be extendet by other unusual AA-names.
	# These are treated than as its substitution, and not as unknown AA
	Subs_rare = {'CYM':'CYS', 'MSE':'MET', 'GLH':'GLU', 'LYN':'LYS', 'AR0':'ARG', 'ASH':'ASP' }
	for X in Subs_rare:
	    AAtype_dict[X] = AAtype_dict[Subs_rare[X]]

	# Residue names which are not in the dictionary yet, are given the value 21
	unknown = 21
	AAtype_dict = collections.defaultdict(lambda: unknown, AAtype_dict)

	if Show_INFOS:
	    print "Initialized AAtype_dict to ", AAtype_dict
	# end initialization of AAtype_dict


    Residues = [str(X.name) for X in centroid.topology.residues]
    NResidues = len(Residues)
    Residues = [AAtype_dict[X] for X in Residues]

    # empty positions are given value zero
    zero3 = [0, 0, 0]
    Residues = zero3 + Residues + zero3
    Residues = np.array(Residues)
    assert Residues.shape == (centroid.n_residues+6,)


    neigh = np.empty((centroid.n_residues, 7))
    for col in range(7):
	neigh[:,col] = Residues[col:col + centroid.n_residues]
    neigh = neigh.astype(float)

    # neigh[i,3+delta] = code of amino acid of residue i+delta where delta in {-3,-1,..,3}
    # if a residue is not known (not included in AAtype_dict) then value = 21
    # if aa does not exist (at chain-ends) then value = 0

    return (neigh)


