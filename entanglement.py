#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# LIBRARIES
## managing files interactively
from argparse import ArgumentParser
from pathlib import Path
import sys
## parser for Les Houches Events file
import internal.lhe_parser as lhe_parser
## statistics and vectors
import math
import numpy as np # type: ignore
import statistics as stat
## check
print('\nLibraries properly imported\n')

# FUNCTIONS
## function that returns angle between vectors (in radians)
def angle(v, w): return np.arccos(v.dot(w)/(np.linalg.norm(v)*np.linalg.norm(w)))
## function that returns the cosine of the angle between vectors
def cosine(v, w): return v.dot(w)/(np.linalg.norm(v)*np.linalg.norm(w))

# MAIN
def main():

    #############################################################
    # open and load file of events interactively
    ############################################################# 

    ## give the user a description of code functionalities
    cl_parser = ArgumentParser(description="Studying the C matrices for detecting entanglement and BIs violation")
    cl_parser.add_argument("filename", help="(name of the file to load)", nargs=1)
    args = cl_parser.parse_args()

    #############################################################
    # iterate through run directories and prompt user
    ############################################################# 

    ## resolve the absolute path to the script
    file_path = Path(__file__)
    print(f"{file_path=}\n") 
    ## construct the path to the events directory (parent goes up the dir tree)
    events_path = file_path.parent.parent / "Events" 
    print(f"{events_path=}\n")
    ## check if the events directory exists
    if not events_path.exists(): 
        print("Events directory does not exist. Exiting...")
        sys.exit()
    ## print out the directories matching the glob pattern and propose to analyze
    print(f"Available runs:\n") 
    for cur_run in events_path.glob("run_*"):
        print(f"-\t{cur_run}") 
        event_file = cur_run / args.filename[0]
        event_file.resolve() 
        if not event_file.exists():
            print(f"   I am skipping this directory because it has no files")
            continue
        else:
            print(f"   {event_file=}")
            ## ask the user if the file is correct
            print(f"I am going to use the file named {event_file}. Continue?")
            while True:
                user_input = input("(yes/no): ")
                if user_input.lower() in ["yes", "y"]:
                    print("Continuing...")
                    break
                elif user_input.lower() in ["no", "n"]:
                    print("Exiting...")
                    sys.exit()
                else:
                    print("Invalid input. Please enter yes/no.")
    
        ## initialize vectors for storage
        top_momenta = []
        antilep_momenta = []
        antitop_momenta = []
        lep_momenta = []
        ## vectors for boosting and constructing helicity basis
        t_antit_sum = [] # total momentum of tt~ pair
        p = [] # beam direction
        ## helicity basis (orthonormal)
        k = [] # top direction
        r = [] # on (p,k) plane
        n = [] # orthogonal to (p,k) plane
        ## initialize entanglement matrix elements
        C_kk = 0
        C_rr = 0
        C_nn = 0
        err_kk = 0
        err_rr = 0
        err_nn = 0
        ## initialize x_ij (only diagonal terms i=j matter)
        x_kk = []
        x_rr = []
        x_nn = []

        #############################################################
        # parse the file event by event
        #############################################################
        F=lhe_parser.EventFile(event_file)
        for iev, event in enumerate(F):
            ## find particles in the event file ('event' type is a list)
            bool_top = False
            bool_antitop = False
            bool_antilep = False
            bool_lep = False
            for part in event:
                ################    top & anti-top    ################
                # top
                if part.pid == 6: 
                    bool_top = True
                    top_momenta.append(lhe_parser.FourMomentum(part))
                # antitop
                if part.pid == -6:
                    bool_antitop = True
                    antitop_momenta.append(lhe_parser.FourMomentum(part))
                ################    anti-lepton & lepton    ################
                # anti-lepton
                if part.pid==-11 or part.pid==-13 or part.pid==-15:
                    bool_antilep = True
                    antilep_momenta.append(lhe_parser.FourMomentum(part))
                # lepton
                if part.pid==11 or part.pid==13 or part.pid==15:
                    bool_lep = True
                    lep_momenta.append(lhe_parser.FourMomentum(part))
                ## check whether an event lacks one of the particles above
            if bool_top==False or bool_antitop==False or bool_antilep==False or bool_lep==False:
                print(f"Particle missing in event number {iev}!")
            
            #############################################################
            # perform boosts and build helicity bases
            # note: boosts do not commute!
            #############################################################

            # BOOST TO C.O.M.F. OF tt~

            ## c.o.m. momentum
            t_antit_sum.append(top_momenta[iev]+antitop_momenta[iev])

            ## perform Lorentz boosts on all momenta
            top_momenta[iev] = top_momenta[iev].boost_to_restframe(t_antit_sum[iev])
            antitop_momenta[iev] = antitop_momenta[iev].boost_to_restframe(t_antit_sum[iev])
            antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(t_antit_sum[iev])
            lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(t_antit_sum[iev])

            # DEFINE HELICITY BASES

            ## beam direction
            if top_momenta[iev].pz*t_antit_sum[iev].pz>0:
                p.append(np.array([0,0,1]))
            else:
                p.append(np.array([0,0,-1]))

            ## top direction
            k.append(np.array([top_momenta[iev].px, top_momenta[iev].py, top_momenta[iev].pz]))
            k[iev] = k[iev]/np.linalg.norm(k[iev]) # normalization
            ## on the (k,p) plane
            theta = angle(p[iev],k[iev])
            r.append((p[iev]-k[iev]*np.cos(theta))/np.sin(theta))
            ## orthogonal to k and r
            n.append(np.cross(k[iev],r[iev]))

            # BOOST FROM C.O.M.F. TO REST FRAME OF SINGLE t / SINGLE t~

            antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(top_momenta[iev])
            lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(antitop_momenta[iev])

        #############################################################
        # calculate the C matrix
        #############################################################
        
            ## get the spatial component to project on the helicity basis
            antilep_3momenta = np.array([antilep_momenta[iev].px, antilep_momenta[iev].py, antilep_momenta[iev].pz])
            lep_3momenta = np.array([lep_momenta[iev].px, lep_momenta[iev].py, lep_momenta[iev].pz])
            ## calculate the x_ii
            if t_antit_sum[iev]*t_antit_sum[iev]>E_thres**2:
                continue
            x_kk.append(-9.*cosine(antilep_3momenta, k[iev])*cosine(lep_3momenta, k[iev]))
            x_rr.append(-9.*cosine(antilep_3momenta, r[iev])*cosine(lep_3momenta, r[iev]))
            x_nn.append(-9.*cosine(antilep_3momenta, n[iev])*cosine(lep_3momenta, n[iev]))
        
        print(f"{len(x_kk)} events found below E<{E_thres} threshold")

        C_kk = np.mean(x_kk)
        err_kk = stat.stdev(x_kk)/math.sqrt(len(x_kk))
        print(f"{C_kk} +- {err_kk}")
        C_rr = np.mean(x_rr)
        err_rr = stat.stdev(x_rr)/math.sqrt(len(x_rr))
        print(f"{C_rr} +- {err_rr}")
        C_nn = np.mean(x_nn)
        err_nn = stat.stdev(x_nn)/math.sqrt(len(x_nn))
        print(f"{C_nn} +- {err_nn}")

        print(C_kk+C_rr)

        print(f"C = {abs(C_kk+C_rr)-C_nn}")

if __name__ == "__main__":
    main()