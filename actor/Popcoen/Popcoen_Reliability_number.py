#! /usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd



def convert_MD_to_csv():

    sys.exit("This function was used to produce the Geometrics_of_MoDEL_trainingdata.csv file. Do you really want to call it???")

    print "I convert MD to csv"
    MD = np.load("/home/goethe/Y___NOW_POSTDOC/Arbeiten/EntropyPrediction2015/Data_PoPE_v5/X_v5_newver-FXRd.npy")
    MD = MD[:,[0,8,9,10]]
    MD = [tuple(X) for X in MD]
    MD = list(set(MD))
    assert len(MD) == 961

    DATT = []
    for s in MD:
	NRes = s[0]
	eigs = np.array(s[1:])

	assert eigs.shape == (3,)
	Vol = 4./3.* np.pi * np.prod( np.sqrt(eigs) )

	length_ratio = sorted(list(np.sqrt(eigs)))
	length_ratio = length_ratio[2] / length_ratio[0]

	DATT.append([Vol/NRes, length_ratio, NRes])


    # here we use first time correct name for density. ints the INVERSE elipthical density
    MD = pd.DataFrame(DATT, columns=["inv_elyptical_density", "length_ratio", "NRes"])
    del DATT
    MD.to_csv("Geometrics_of_MoDEL_trainingdata.csv", sep=";", index=False)
    print "csv file saved"







def load_Geometrics_of_trainingdata(Geometrics_filename = "./Geometrics_of_MoDEL_trainingdata.csv"):
    print(os.getcwd())
    MD = pd.read_csv(Geometrics_filename, sep=";")
    return(MD)






def measure_Geometrics_from_features(featuresX):

    GG = featuresX[:,[0,8,9,10]]
    GG = [tuple(X) for X in GG]
    GG = list(set(GG))
    assert len(GG) == 1  # the rows 0,8,9,10 have the same number for all AAs of a protein
    GG = GG[0]

    NRes = GG[0]
    eigs = np.array(GG[1:])

    assert eigs.shape == (3,)
    Vol = 4./3.* np.pi * np.prod( np.sqrt(eigs) )

    length_ratio = sorted(list(np.sqrt(eigs)))
    length_ratio = length_ratio[2] / length_ratio[0]

    Geo = [ [Vol/NRes, length_ratio, NRes] ]
    # here we denote Vol/NRes correctly as INVERSE_elipthical_density, and not as in the rest of the code elyptical density
    Geo = pd.DataFrame(Geo, columns=["inv_elyptical_density", "length_ratio", "NRes"])

    return(Geo)






MDvals = ""
def count_neighbors(CASP, radius_ball = 1.784124116):

    global MDvals, R1, R2, NormMatrix
    if len(MDvals) == 0:   # not jet initialized
	#print "INFO: load MDvals data"
	MDvals = load_Geometrics_of_trainingdata()
	MDvals = MDvals.as_matrix(columns=["length_ratio", "inv_elyptical_density"])

	R1 = MDvals[:,0].std()
	R2 = MDvals[:,1].std()

	NormMatrix = np.resize([R1,R2], (len(MDvals),2))


    # grep casp data
    assert "length_ratio"      in [str(X) for X in CASP.columns]
    assert "inv_elyptical_density" in [str(X) for X in CASP.columns]
    C = CASP.as_matrix(columns=["length_ratio", "inv_elyptical_density"])


    NNdat = []   # list of Number Neighbors in radius. Streched axis
    for pt in C: # only one in this version

	selt = (MDvals - pt) / NormMatrix
	selt = np.linalg.norm(selt, axis=1) < radius_ball
	NNdat.append(selt.sum())

	if 1==0:  # testplot
	    selt = np.hstack((np.reshape(selt,(len(selt),1)) , MDvals))
	    selt = np.array([X[1:] for X in selt if X[0] == 1])

	    import pylab
	    pylab.plot(MDvals[:,0], MDvals[:,1], "x")
	    pylab.plot(selt[:,0], selt[:,1], ".")
	    pylab.plot([pt[0]], [pt[1]],"o")
	    pylab.show()
	    sys.exit(0)

    assert len(NNdat) == 1
    number_points_in_streched_circle = NNdat[0]
    return(number_points_in_streched_circle)





def calculate_reliability_number(featuresX, give_more_info_back = False):

    # first calculate the two geometric properties "length_ratio" and "inverse_elyptical_density"
    Geo = measure_Geometrics_from_features(featuresX)

    # then count how many proteins of MoDEL data are in a (streched) sphere around with radius 1.7841, so that area = 10
    neigs = count_neighbors(Geo)

    # calculate reliability number from this. The fitting was done with PLOT_Reliability_number_fitting.py
    # the fitted half-life so close to 100 that we simply use 100
    kappa = 100./np.log(2)
    relia = 1. - np.exp(-neigs/kappa)

    if give_more_info_back:
	return([relia, Geo, neigs])
    else:
	return(relia)

