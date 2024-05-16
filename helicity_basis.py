#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# LIBRARIES
# managing files interactively
import sys
# parser for Les Houches Events file
import internal.lhe_parser as lhe_parser
# statistics and vectors
import math
import statistics as stat
import numpy as np
print('Libraries properly imported')

# FUNCTIONS
# function that returns angle between vectors (in radians)
def angle(v, w): return np.arccos(v.dot(w)/(np.linalg.norm(v)*np.linalg.norm(w)))
def cosine(v, w): return v.dot(w)/(np.linalg.norm(v)*np.linalg.norm(w))

# MAIN
def main():

    #############################################################
    # open and load file of events
    #############################################################

    # open and load file of events
    eventfile = sys.argv[1]
    E_thres = 500.0 # energy threshold in GeV
    F=lhe_parser.EventFile(eventfile)
    # initialize vectors for storage
    top_momenta = []
    antilep_momenta = []
    antitop_momenta = []
    lep_momenta = []
    # vectors for boosting and constructing helicity basis
    t_antit_sum = [] # 
    p = [] # beam direction
    # top rest frame helicity basis (orthonormal)
    k = [] # top direction
    r = [] # on (p,k) plane
    n = [] # orthogonal to (p,k) plane
    """# anti-top rest frame helicity basis
    antik = [] # anti-top direction
    antir = [] 
    antin = []"""
    # initialize entanglement matrices
    C_kk = 0
    C_rr = 0
    C_nn = 0
    err_kk = 0
    err_rr = 0
    err_nn = 0
    # initialize x_ij
    x_kk = []
    x_rr = []
    x_nn = []

    #############################################################
    # parse the file event by event
    #############################################################
    
    print(f"Processing the events...")
    for iev, event in enumerate(F):
        # find particles in the event file ('event' type is a list)
        bool_top = False
        bool_antitop = False
        bool_antilep = False
        bool_lep = False
        for part in event:
            ###############    top & anti-top    ################
            # top
            if part.pid == 6: 
                bool_top = True
                top_momenta.append(lhe_parser.FourMomentum(part))
            # antitop
            if part.pid == -6:
                bool_antitop = True
                antitop_momenta.append(lhe_parser.FourMomentum(part))
            ###############    anti-lepton & lepton    ################
            # anti-lepton
            if part.pid==-11 or part.pid==-13 or part.pid==-15:
                bool_antilep = True
                antilep_momenta.append(lhe_parser.FourMomentum(part))
            # lepton
            if part.pid==11 or part.pid==13 or part.pid==15:
                bool_lep = True
                lep_momenta.append(lhe_parser.FourMomentum(part))
        # check whether an event lacks one of the particles above
        if bool_top==False or bool_antitop==False or bool_antilep==False or bool_lep==False:
            print(f"Particle missing in event number {iev}!")
        
        #############################################################
        # perform boosts and build helicity bases
        # note: boosts do not commute!
        #############################################################

        # BOOST TO C.O.M.F. OF tt~

        # c.o.m. momentum
        t_antit_sum.append(top_momenta[iev]+antitop_momenta[iev])

        # perform Lorentz boosts on all momenta
        top_momenta[iev] = top_momenta[iev].boost_to_restframe(t_antit_sum[iev])
        antitop_momenta[iev] = antitop_momenta[iev].boost_to_restframe(t_antit_sum[iev])
        antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(t_antit_sum[iev])
        lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(t_antit_sum[iev])

        # DEFINE HELICITY BASES

        # beam direction
        if top_momenta[iev].pz*t_antit_sum[iev].pz>0:
            p.append(np.array([0,0,1]))
        else:
            p.append(np.array([0,0,-1]))

        # top direction
        k.append(np.array([top_momenta[iev].px, top_momenta[iev].py, top_momenta[iev].pz]))
        k[iev] = k[iev]/np.linalg.norm(k[iev]) # normalization
        # on the (k,p) plane
        theta = angle(p[iev],k[iev])
        r.append((p[iev]-k[iev]*np.cos(theta))/np.sin(theta))
        # orthogonal to k and r
        n.append(np.cross(k[iev],r[iev]))
        """
        # anti-top direction
        antik.append(np.array([antitop_momenta[iev].px, antitop_momenta[iev].py, antitop_momenta[iev].pz]))
        antik[iev] = antik[iev]/np.linalg.norm(antik[iev]) # normalization
        # on the (antik,p) plane
        antitheta = angle(p[iev],antik[iev])
        antir.append((p[iev]-antik[iev]*np.cos(antitheta))/np.sin(antitheta))
        # orthogonal to antik and antir
        antin.append(np.cross(antik[iev],antir[iev]))
        """
        # BOOST FROM C.O.M.F. TO REST FRAME OF SINGLE t / SINGLE t~

        antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(top_momenta[iev])
        lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(antitop_momenta[iev])

        #############################################################
        # prepare to calculate the C matrix
        #############################################################
    
        # get the spatial component to project on the helicity basis
        antilep_3momenta = np.array([antilep_momenta[iev].px, antilep_momenta[iev].py, antilep_momenta[iev].pz])
        lep_3momenta = np.array([lep_momenta[iev].px, lep_momenta[iev].py, lep_momenta[iev].pz])
        # calculate the x_ii
        if t_antit_sum[iev]*t_antit_sum[iev]>E_thres**2:
            continue
        x_kk.append(-9.*cosine(antilep_3momenta, k[iev])*cosine(lep_3momenta, k[iev]))
        x_rr.append(-9.*cosine(antilep_3momenta, r[iev])*cosine(lep_3momenta, r[iev]))
        x_nn.append(-9.*cosine(antilep_3momenta, n[iev])*cosine(lep_3momenta, n[iev]))

    # end of for loop
    print(f"{len(x_kk)} events found below E<{E_thres} threshold")

    # determine matrix elements
    print(f"From the kinematic calculations follows that:")
    C_kk = np.mean(x_kk)
    err_kk = stat.stdev(x_kk)/math.sqrt(len(x_kk))
    print(f"C_kk = {C_kk} +- {err_kk}")
    C_rr = np.mean(x_rr)
    err_rr = stat.stdev(x_rr)/math.sqrt(len(x_rr))
    print(f"C_rr = {C_rr} +- {err_rr}")
    C_nn = np.mean(x_nn)
    err_nn = stat.stdev(x_nn)/math.sqrt(len(x_nn))
    print(f"C_nn = {C_nn} +- {err_nn}")

    # find entanglement
    print(C_kk+C_rr) # at threshold, this quantity should be negative
    C = abs(C_kk+C_rr)-C_nn
    err_C = math.sqrt(err_kk**2+err_nn**2+err_rr**2)
    print(f"C = {C} +- {err_C}")

if __name__ == "__main__":
    main()
