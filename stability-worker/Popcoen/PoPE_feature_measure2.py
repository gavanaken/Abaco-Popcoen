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







def measure_GyrationTensorStuff(centroid, Consistency_check = False):

    assert centroid.n_frames == 1


    pos  = 10. * centroid.xyz[0,:,:].astype(float)
    mean = np.mean(pos, axis=0)

    if Consistency_check:
	MassCenter = 10. * mdtraj.compute_center_of_mass(centroid)[0]
	diff_mean_masscenter = mean - MassCenter

    pos  = pos - mean

    rg  = np.resize(0.0, (3,3))
    for X in pos:
        for a1 in range(3):
            for a2 in range(3):
                rg[a1,a2] += X[a1] * X[a2]
    rg = rg / len(pos)

    eigvals, vects = np.linalg.eigh(rg)
    vects = vects.T   # npw each row contains an eigenvector. See test below

    # the eigvals are sorted, are they? If not, then numpy has changed format. Do sort them by hand then
    assert eigvals[0] <= eigvals[1]
    assert eigvals[1] <= eigvals[2]

    # give defined direction to vectrors. abitrary but defined
    vn = []
    for v in vects:
	if v[0] < 0.0:
	    vn.append(-v)
	else:
	    vn.append(v)
    vects = np.array(vn)
    del vn


    if Consistency_check: 
	# both test should be rg^2 in angstroem^2 units
	test = eigvals.sum()
	test2 = np.power( 10.*mdtraj.compute_rg(centroid), 2.0)
	assert np.fabs((test-test2)/test2) < 1e-5

	# test eigenvector properties
	for l,v in zip(eigvals,vects):
	    assert np.linalg.norm( np.dot(rg, v)  -   l * v ) < 1e-5
	    assert np.fabs( np.linalg.norm(v) - 1.0 ) < 1e-5
	    assert v[0] >= 0.0

	# test still correctly shifted
	assert np.linalg.norm ( np.mean(pos, axis=0) ) < 1e-5

    # get positions of CA atoms
    CAs    = [a.index for a in centroid.topology.atoms if a.name == 'CA']
    CA_pos = pos[CAs, :]

    # describe distance vectors in eigenvector basis
    projections = []
    for p in CA_pos:
	# test, that vectorized version correct
	#projects = []
	#for v in vects:
	    #projects.append( np.dot(p,v) )
	#print projects
	projections.append( np.dot(vects, p) )
    projections = np.array(projections)


    # first three cols are eigvals
    BACK = np.array([eigvals for tmp in range(len(projections))])
    # then projections
    BACK = np.hstack( (BACK, projections) )



    # consistency check with old geometry features.
    if Consistency_check:
	print "INFO: Im doing consistency check with old Geometry measurement"
	import PoPE_feature_measure
	Geo = PoPE_feature_measure.measure_Geometry(centroid)
	assert len(Geo) == len(BACK)
	for t in zip(Geo,BACK):
	    assert np.fabs( (t[0][0] - t[1][0:3].sum()) / t[0][0] ) < 1e-5    # rg^2 equals sum eigvals
	    verschieb = t[1][3:] + np.dot(vects, diff_mean_masscenter)
	    assert np.fabs( t[0][1] - np.dot( verschieb, verschieb ) ) < 1e-5   # dist^2 invariant under rotation. Shift because of different midpoint definitions used inj both versions


    # the BACK field contains now six colums and NumberRes rows.
    # cols 0,1,2:  Eigvals of gyration tensor, sorted from small to large
    # cols 3,4,5:  distance vector from mean position (not MassCenter but very close to it), measured in the eigen-basis given by Gyration-tensor.
    return(BACK)











def measure_mean_angles_of_neighborhood(traj):

    assert traj.n_frames == 1

    atoms = [X for X in traj.topology.atoms]
    NRes = traj.n_residues


    sets = []
    sets.append(["psi" ,  mdtraj.compute_psi  (traj, periodic=False) ])
    sets.append(["phi" ,  mdtraj.compute_phi  (traj, periodic=False) ])
    sets.append(["chi1",  mdtraj.compute_chi1 (traj, periodic=False) ])
    sets.append(["chi2",  mdtraj.compute_chi2 (traj, periodic=False) ])
    sets.append(["chi3",  mdtraj.compute_chi3 (traj, periodic=False) ])
    sets.append(["chi4",  mdtraj.compute_chi4 (traj, periodic=False) ])
    sets.append(["omega", mdtraj.compute_omega(traj, periodic=False) ])

    AnsSets = len(sets)
    assert AnsSets == 7 # ERROR: Martin, have you changed this? Not all existing angles anymore?

    AoR = np.resize(np.nan, (NRes, AnsSets))  # mean angles of residues array


    pi = np.pi
    #No shift with this: DefinitionZeroAngle = {"psi":0.0, "phi":0.0, "chi1":0.0, "chi2":0.0, "chi3":0.0, "chi4":0.0, "omega":0.0}
    DefinitionZeroAngle = {"psi":pi/2, "phi":pi, "chi1":5.*pi/3, "chi2":pi/4, "chi3":pi/3, "chi4":pi, "omega":pi/2}




    for columnindex, s in enumerate(sets):

	TorsionName = s[0]
	atomdefinitions, timeline = s[1]

	assert timeline.shape[0] == 1  # we have only one frame here, namely centroid
	timeline = timeline[0]

	#We shift the torsion angles such that little weight at pi=-pi. Then NeuralNetwork does not have to learn this confusion
	# for omega e.g.: we shift omega by pi/2, such that trans at pi/2 and cis at -pi/2. In this way, no jumps between -pi+eps to pi-eps for trans
	timeline = timeline - DefinitionZeroAngle[TorsionName]
	where_out = (timeline < -np.pi).astype(float)
	timeline = timeline + where_out * (2. * np.pi)
	assert np.all(timeline >= -np.pi)
	assert np.all(timeline <= np.pi)


	assert len(timeline) == len(atomdefinitions)
	ZIPP = zip(atomdefinitions, timeline)

	# fill into correct row of table
        for four_AtomIndices_def_torsion, meanPositionAngle in ZIPP:
            Atom2    = atoms[four_AtomIndices_def_torsion[1]] # second atom of tuple decodes the aa this torsion belongs to
	    rowindex = Atom2.residue.index
	    AoR[rowindex, columnindex] = meanPositionAngle



    # for testing use: AoR = np.array([[i,i,i,i,i,i,i] for i in range(NRes)])
    ManyCopies = []
    for i in range(NRes):

	line = []
	for shift in range(-3,4):
	    getem = i + shift

	    if (getem < 0) or (getem > NRes-1):
		putem = np.array([np.nan for tmp in range(AnsSets)])
	    else:
		putem = AoR[getem,:]
	    line += list(putem)
	ManyCopies.append(line)

    ManyCopies = np.array(ManyCopies)
    assert ManyCopies.shape == (NRes, AnsSets * 7)

    # Now row i of ManyCopies includes the seven angles psi,phi,..,omega for residue i-3, i-2, i-1, i, i+1, i+2, i+3
    # nonexisting angles (either neighbors does not exist, e.g. ResIndex=-1, or chi's which do not exist, or phi,psi at beginning/end)
    # are given np.nan
    # Die omega angles are shifted such that trans at pi/2 and cis at -pi/2


    return (ManyCopies)









def Hbonding_Infos(centroid, Show_INFOS = False):

    assert centroid.n_frames == 1

    hbonds = mdtraj.wernet_nilsson(centroid, periodic=False)[0]
    total_n_hbonds = len(hbonds)
    assert hbonds.shape == (total_n_hbonds, 3)

    # aceptor+donor counts
    adc = np.zeros(centroid.n_residues)
    for hbond in hbonds:
	a1,a2 = centroid.topology.atom(hbond[0]), centroid.topology.atom(hbond[2])
	r1,r2 = a1.residue.index, a2.residue.index
	adc[r1] += 1
	adc[r2] += 1

    nans = [np.nan, np.nan, np.nan]
    adc = np.hstack( ( nans, adc, nans ) )
    assert adc.shape == (centroid.n_residues+6,)

    neigh = np.empty((centroid.n_residues, 8))
    neigh[:,0] = total_n_hbonds
    for col in range(7):
	neigh[:,col+1] = adc[col:col + centroid.n_residues]
    neigh = neigh.astype(float)
    assert neigh[:,4].sum() == 2.*total_n_hbonds  # factor 2 because one donor + one acceptor counted per hbond

    # info about Hbonding. 
    # col0 = total number hbonds for all rows i
    # col4 + delta = number of hbonds of residue (i+delta) with delta in {-3,-2,..,3}
    # if residue does not exist (at chain-ends) then np.nan
    # hence

    return (neigh)



