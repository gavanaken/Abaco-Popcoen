#! /usr/bin/env python
import os
import sys
import pickle
import mdtraj
import numpy as np
import math
import itertools
import random
import pylab


import PoPE_feature_measure
import PoPE_feature_measure2



def reinitialize_Residue_numbering(centroid):
    top_df, bonds_inimportance_here = centroid.topology.to_dataframe()
    ResSeqDict = top_df.resSeq.unique()
    ResSeqDict = dict(zip(ResSeqDict, range(len(ResSeqDict))))
    top_df.resSeq = top_df.resSeq.apply(lambda X: ResSeqDict[X])
    centroid.topology = centroid.topology.from_dataframe(top_df, bonds_inimportance_here)




def measure_all_features(centroid):


	# Check that no negative residue-Ids
	top_df, bonds_inimportance_here = centroid.topology.to_dataframe()
	if np.all(top_df.resSeq.unique() >= 0) == False:  # there are residues with negative resnumbers
	    sys.exit("ERROR: Negative residue-Ids encountered. Mdtraj version (1.8.0) cannot handle this (bug). Please rename the residue-Ids to be non-negative.")


	# drop if more than one frame
	if centroid.n_frames > 1:
	    print "WARNING: pdb-file contains more than one MODEL. Popcoen calculation for first model of file."
	    centroid = centroid[0]


	# drop atoms which dont belong to protein
	protein_atoms = centroid.top.select("protein")
	all_atoms = centroid.top.select("all")
	if len(protein_atoms) == 0:
	    sys.exit("ERROR: pdb-file does not contain protein atoms. Popcoen (=Prediction of Protein Configurational Entropy) can predict entropy only for proteins.")
	if np.all ( protein_atoms == all_atoms ) == False:
	    centroid = centroid.atom_slice(protein_atoms)
	    print "WARNING: pdb-file contains atoms not-belonging to protein. Popcoen drops all non-protein atoms."


	# obtain smallest chain ID
	top_df, bonds_inimportance_here = centroid.topology.to_dataframe()
	chainIDs = sorted(list(top_df.chainID.unique()))
	smallestchainID = chainIDs[0]
	if smallestchainID != 0:
	    print "WARNING: Untypical chainID=", smallestchainID
	del top_df, bonds_inimportance_here, chainIDs
	# drop all chains except first one
	protein_atoms = centroid.top.select("chainid == %i" % smallestchainID)
	all_atoms = centroid.top.select("all")
	if np.all ( protein_atoms == all_atoms ) == False:
	    centroid = centroid.atom_slice(protein_atoms)
	    print "WARNING: pdb-file contains multiple chains. All but first chain dropped."



	# check that no resSeq has multiple resnames. This happens 
	# because of 27th character of ATOM record (AChar=Code for insertion of residues.)
	# which cannot be handeled by mdtraj. Poor mdtraj!!!
	top_df, bonds_inimportance_here = centroid.topology.to_dataframe()
	tmpdf = top_df[["resSeq", "resName"]].copy().drop_duplicates()
	tmpdf = tmpdf.groupby("resSeq").count()
	tmpdf = tmpdf[tmpdf.resName > 1].index.values
	if len(tmpdf) > 0:
	    sys.exit("ERROR: pdb-file contains some residue-Id with multiple residue-names. This is likely because of the AChar (column 27 of ATOM record) which cannot be handeled by mdtray (version 1.8.0), and hence neither by popcoen. Please rename the residue-ids using only numbers, no AChar. Check following Residue-Ids: %s" % str(tmpdf))
	del tmpdf



	top_df, bonds_inimportance_here = centroid.topology.to_dataframe()
	# check that there is sufficent SC information. If not, Popcoen will give incorrect output as chi angles cannot be calculated.
	different_SC_atoms = set(top_df[top_df.element != "H"].name)
	for BBa in ["C", "O", "CA", "N"]:
	    try:
		different_SC_atoms.remove(BBa)
	    except KeyError:
		sys.exit("ERROR: No coordinates found for Backbone atom type %s. Information required for Popcoen calculation." %BBa)
	if len(different_SC_atoms) < 10:
	    sys.exit("ERROR: Insufficent information about side-chain atoms. Found only these types: %s." %sorted(list(different_SC_atoms)) )
	# if there are Se atoms: set Se to S atoms
	elements_contained = top_df.element.unique()
	if "Se" in elements_contained:
	    top_df.element = top_df.element.apply(lambda X: X if X != "Se" else "S")
	    centroid.topology = centroid.topology.from_dataframe(top_df, bonds_inimportance_here)
	    print "WARNING: pdf-file contained Se atoms which were set to sulfur."
	del top_df, bonds_inimportance_here



	# check that each amino acid has a CA atom
	if centroid.n_residues != len(centroid.top.select("name CA")):
	    # first we drop those residues without CA atom
	    CAs = centroid.top.select("name CA")
	    CAs = centroid.atom_slice(CAs)
	    Residues_with_CA = [X.residue.resSeq for X in CAs.topology.atoms]
	    top_df, bonds_inimportance_here = centroid.topology.to_dataframe()
	    Selected = top_df[top_df.resSeq.isin(Residues_with_CA)].index.tolist()
	    centroid = centroid.atom_slice(Selected)
	    # then we reset the resSeq values such that it is again zero..n-1 counting
	    print "WARNING: pdf-file contains residues without CA atom. Those residues were dropped and Residue-numbering was changed."
	    reinitialize_Residue_numbering(centroid)
	    if centroid.n_residues != len(centroid.top.select("name CA")):
		sys.exit("ERROR: Very strange error: Why does the number of residues and CA atoms still does not agree??? Please contact the developers and give them pdb file for reducibility.")



	# are there also hydrogens in data? If yes, we later can calculate Hbonding
	top_df, bonds_inimportance_here = centroid.topology.to_dataframe()
	Hydrogens_exist_for_Hbonding_calculation = False
	if "H" not in top_df.element.unique():
	    print "WARNING: Structure does not contain hydrogens. Hence, Hbonding cannot be calculated."
	else:
	    fraction_H = (top_df.element == "H").astype(float).sum()/centroid.n_residues
	    if fraction_H < 3.0:
		print "WARNING: There are some hydrogens but probably not all of them. Hence, Hbonding cannot be calculated."
	    else:
		Hydrogens_exist_for_Hbonding_calculation = True



	NResidues = centroid.n_residues
	AAnames = [AAname for AAname in centroid.topology.residues]



	# measure hydrogen bonding before we take out hydrogens
	if Hydrogens_exist_for_Hbonding_calculation:
	    Hbonding = PoPE_feature_measure2.Hbonding_Infos(centroid)
	    assert Hbonding.shape == (NResidues, 8)
	else:
	    sys.exit("ERROR: Insufficient Hydrogen Information found. Please preprocess your structure and add hydrogens. Popcoen developers recommend to use FoldX. Commandline:  ./foldx --command=RepairPDB --pdb=????.pdb --pdbHydrogens=true")



	if Hydrogens_exist_for_Hbonding_calculation:
	    # now we take out hydrogens so that remaining features are exactly equal for structures with and without hydrogens
	    nonhydrogen_atoms = top_df[top_df.element != "H"].index.values
	    centroid = centroid.atom_slice(nonhydrogen_atoms)
	    #print "INFO: We reduce structure to heavy atoms:",centroid
	    testitagain = sorted(list(set([str(a.name[0:1]) for a in centroid.topology.atoms])))
	    assert "H" not in testitagain
	    if testitagain not in ( ['C', 'N', 'O', 'S'] , ['C', 'N', 'O'] ):
		print "WARNING: Unusual element encountered,",testitagain



	# measure all other features except hbonding.
	# they will be stacked together differently (in more pedagogical order)
	NResfeature = np.zeros((NResidues,1), dtype=float) + float(NResidues)
	assert NResfeature.shape == (NResidues, 1)

	# This not anymore as we have GyrationTensorStuff now
	#Geometry = PoPE_feature_measure.measure_Geometry(centroid)
	#assert Geometry.shape == (NResidues, 2)

	ContactNumbers, rcs = PoPE_feature_measure.measure_residue_contact_numbers(centroid)
	assert ContactNumbers.shape == (NResidues, 11)

        SASAnmd = PoPE_feature_measure.measure_SASA_percentage(centroid)
	assert SASAnmd.shape == (NResidues, 2)

        AAneighbors = PoPE_feature_measure.residue_surroundings(centroid)
        assert AAneighbors.shape == (NResidues, 7)

	GyrationTensorStuff = PoPE_feature_measure2.measure_GyrationTensorStuff(centroid)
	assert GyrationTensorStuff.shape == (NResidues, 6)

	Angles_of_Me_and_Neighbors = PoPE_feature_measure2.measure_mean_angles_of_neighborhood(centroid)
	assert Angles_of_Me_and_Neighbors.shape == (NResidues, 49)


	# Pack them together:
	if Hydrogens_exist_for_Hbonding_calculation:
	    All_Features = np.hstack((NResfeature, AAneighbors, GyrationTensorStuff, ContactNumbers, SASAnmd, Angles_of_Me_and_Neighbors, Hbonding))
	    assert All_Features.shape == (NResidues, 1+7+6+11+2+49+8)
	else:
	    All_Features = np.hstack((NResfeature, AAneighbors, GyrationTensorStuff, ContactNumbers, SASAnmd, Angles_of_Me_and_Neighbors))
	    assert All_Features.shape == (NResidues, 1+7+6+11+2+49)



	return(All_Features, AAnames)

