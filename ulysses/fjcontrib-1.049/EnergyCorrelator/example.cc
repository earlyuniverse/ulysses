// Example showing usage of energy correlator classes.
//
// Compile it with "make example" and run it with
//
//   ./example < ../data/single-event.dat
//
// Copyright (c) 2013-2016
// Andrew Larkoski, Lina Necib, Gavin Salam, and Jesse Thaler
//
// $Id: example.cc 1097 2018-01-05 00:04:20Z linoush $
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <sstream>
#include "EnergyCorrelator.hh" // In external code, this should be fastjet/contrib/EnergyCorrelator.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);
void analyze(const vector<PseudoJet> & input_particles);

//----------------------------------------------------------------------
int main(){

    //----------------------------------------------------------
    // read in input particles
    vector<PseudoJet> event;
    read_event(event);
    cout << "# read an event with " << event.size() << " particles" << endl;

    //----------------------------------------------------------
    // illustrate how this EnergyCorrelator contrib works

    analyze(event);

    return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){
    string line;
    while (getline(cin, line)) {
        istringstream linestream(line);
        // take substrings to avoid problems when there are extra "pollution"
        // characters (e.g. line-feed).
        if (line.substr(0,4) == "#END") {return;}
        if (line.substr(0,1) == "#") {continue;}
        double px,py,pz,E;
        linestream >> px >> py >> pz >> E;
        PseudoJet particle(px,py,pz,E);

        // push event onto back of full_event vector
        event.push_back(particle);
    }
}

////////
//
//  Main Routine for Analysis
//
///////

void analyze(const vector<PseudoJet> & input_particles) {

    /////// EnergyCorrelator /////////////////////////////

    // Initial clustering with anti-kt algorithm
    JetAlgorithm algorithm = antikt_algorithm;
    double jet_rad = 1.00; // jet radius for anti-kt algorithm
    JetDefinition jetDef = JetDefinition(algorithm,jet_rad,E_scheme,Best);
    ClusterSequence clust_seq(input_particles,jetDef);
    vector<PseudoJet> antikt_jets  = sorted_by_pt(clust_seq.inclusive_jets());

    for (int j = 0; j < 2; j++) { // Two hardest jets per event
        if (antikt_jets[j].perp() > 200) {

            PseudoJet myJet = antikt_jets[j];

            // various values of beta
            vector<double> betalist;
            betalist.push_back(0.1);
            betalist.push_back(0.2);
            betalist.push_back(0.5);
            betalist.push_back(1.0);
            betalist.push_back(1.5);
            betalist.push_back(2.0);

            // various values of alpha
            vector<double> alphalist;
            alphalist.push_back(0.1);
            alphalist.push_back(0.2);
            alphalist.push_back(0.5);
            alphalist.push_back(1.0);


            // checking the two energy/angle modes
            vector<EnergyCorrelator::Measure> measurelist;
            measurelist.push_back(EnergyCorrelator::pt_R);
            measurelist.push_back(EnergyCorrelator::E_theta);
            //measurelist.push_back(EnergyCorrelator::E_inv);

            vector<string> modename;
            modename.push_back("pt_R");
            modename.push_back("E_theta");
            //modename.push_back("E_inv");

            for (unsigned int M = 0; M < measurelist.size(); M++) {

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelator:  ECF(N,beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s %14s %14s %14s %15s\n","beta", "N=1 (GeV)", "N=2 (GeV^2)", "N=3 (GeV^3)", "N=4 (GeV^4)", "N=5 (GeV^5)");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelator ECF0(0,beta,measurelist[M]);
                    EnergyCorrelator ECF1(1,beta,measurelist[M]);
                    EnergyCorrelator ECF2(2,beta,measurelist[M]);
                    EnergyCorrelator ECF3(3,beta,measurelist[M]);
                    EnergyCorrelator ECF4(4,beta,measurelist[M]);
                    EnergyCorrelator ECF5(5,beta,measurelist[M]);

                    printf("%7.3f %14.2f %14.2f %14.2f %14.2f %15.2f \n",beta,ECF1(myJet),ECF2(myJet),ECF3(myJet),ECF4(myJet),ECF5(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorRatio:  r_N^(beta) = ECF(N+1,beta)/ECF(N,beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s %14s %14s %14s %15s \n","beta", "N=0 (GeV)", "N=1 (GeV)", "N=2 (GeV)", "N=3 (GeV)","N=4 (GeV)");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorRatio r0(0,beta,measurelist[M]);
                    EnergyCorrelatorRatio r1(1,beta,measurelist[M]);
                    EnergyCorrelatorRatio r2(2,beta,measurelist[M]);
                    EnergyCorrelatorRatio r3(3,beta,measurelist[M]);
                    EnergyCorrelatorRatio r4(4,beta,measurelist[M]);

                    printf("%7.3f %14.4f %14.4f %14.4f %14.4f %15.4f \n",beta,r0(myJet),r1(myJet),r2(myJet),r3(myJet),r4(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorDoubleRatio:  C_N^(beta) = r_N^(beta)/r_{N-1}^(beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3", "N=4");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorDoubleRatio C1(1,beta,measurelist[M]);
                    EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M]);
                    EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M]);
                    EnergyCorrelatorDoubleRatio C4(4,beta,measurelist[M]);

                    printf("%7.3f %14.6f %14.6f %14.6f %14.6f \n",beta,C1(myJet),C2(myJet),C3(myJet),C4(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorC1:  C_1^(beta) = ECF(2,beta)/ECF(1,beta)^2 with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta","C1 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorC1 c1(beta,measurelist[M]);

                    printf("%7.3f %14.6f \n",beta,c1(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorC2:  C_2^(beta) = ECF(3,beta)*ECF(1,beta)/ECF(2,beta)^2 with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta","C2 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorC2 c2(beta,measurelist[M]);

                    printf("%7.3f %14.6f \n",beta,c2(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;


                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorD2:  D_2^(beta) = ECF(3,beta)*ECF(1,beta)^3/ECF(2,beta)^3 with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta","D2 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorD2 d2(beta,measurelist[M]);

                    printf("%7.3f %14.6f \n",beta,d2(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorGeneralizedD2:  D_2^(alpha, beta) = ECFN(3,alpha)/ECFN(2,beta)^(3*alpha/beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %18s %18s %18s %18s\n","beta","alpha = 0.100","alpha = 0.200","alpha = 0.500","alpha = 1.000");

                for (unsigned int B = 1; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    printf("%7.3f ", beta);
                    for (unsigned int A = 0; A < alphalist.size(); A++) {
                        double alpha = alphalist[A];

                        EnergyCorrelatorGeneralizedD2 d2(alpha, beta, measurelist[M]);

                        printf("%18.4f ", d2(myJet));
                    }
                    printf("\n");
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;



                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorGeneralized (angles = N Choose 2):  ECFN(N, beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %7s %14s %14s %14s\n","beta", "N=1", "N=2", "N=3", "N=4");


                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorGeneralized ECF1(-1,1, beta,  measurelist[M]);
                    EnergyCorrelatorGeneralized ECF2(-1,2, beta, measurelist[M]);
                    EnergyCorrelatorGeneralized ECF3(-1,3, beta, measurelist[M]);
                    EnergyCorrelatorGeneralized ECF4(-1,4, beta, measurelist[M]);
                    //EnergyCorrelatorGeneralized ECF5(-1, 5, beta, measurelist[M]);

                    printf("%7.3f %7.2f %14.10f %14.10f %14.10f \n", beta, ECF1(myJet), ECF2(myJet), ECF3(myJet),
                           ECF4(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" <<
                endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorGeneralized:  ECFG(angles, N, beta=1) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %7s %14s %14s %14s\n","angles", "N=1", "N=2", "N=3", "N=4");

                double beta = 1.0;
                for (unsigned int A = 1; A < 2; A++) {
                    double angle = A;

                    EnergyCorrelatorGeneralized ECF1(angle, 1, beta, measurelist[M]);
                    EnergyCorrelatorGeneralized ECF2(angle, 2, beta, measurelist[M]);
                    EnergyCorrelatorGeneralized ECF3(angle, 3, beta, measurelist[M]);
                    EnergyCorrelatorGeneralized ECF4(angle, 4, beta, measurelist[M], EnergyCorrelator::slow);

                    printf("%7.0f %7.2f %14.10f %14.10f %14.10f \n", angle, ECF1(myJet), ECF2(myJet), ECF3(myJet),
                           ECF4(myJet));

                }

                for (unsigned int A = 2; A < 4; A++) {
                    double angle = A;

                    EnergyCorrelatorGeneralized ECF3(angle, 3, beta, measurelist[M]);
                    EnergyCorrelatorGeneralized ECF4(angle, 4, beta, measurelist[M]);

                    printf("%7.0f %7s %14s %14.10f %14.10f \n", angle, " " , " " ,ECF3(myJet), ECF4(myJet));
                }

                for (unsigned int A = 4; A < 7; A++) {
                    double angle = A;

                    EnergyCorrelatorGeneralized ECF4(angle, 4, beta, measurelist[M]);
                    printf("%7.0f %7s %14s %14s %14.10f \n", angle, " ", " ", " ", ECF4(myJet) );
                }
                cout << "-------------------------------------------------------------------------------------" <<
                endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorNseries:  N_i(beta) = ECFG(i+1, 2, beta)/ECFG(i, 1, beta)^2 with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorNseries N1(1,beta,measurelist[M]);
                    EnergyCorrelatorNseries N2(2,beta,measurelist[M]);
                    EnergyCorrelatorNseries N3(3,beta,measurelist[M]);

                    printf("%7.3f %14.6f %14.6f %14.6f \n",beta,N1(myJet),N2(myJet),N3(myJet));

                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;


                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorN2:  N2(beta) = ECFG(3, 2, beta)/ECFG(2, 1, beta)^2 with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta", "N2 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorN2 N2(beta,measurelist[M]);

                    printf("%7.3f %14.6f \n",beta,N2(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;


                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorN3:  N3(beta) = ECFG(4, 2, beta)/ECFG(3, 1, beta)^2 with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta", "N3 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorN3 N3(beta,measurelist[M]);

                    printf("%7.3f %14.6f \n",beta,N3(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorMseries:  M_i(beta) = ECFG(i+1, 1, beta)/ECFN(i, 1, beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s %14s %14s \n","beta", "N=1", "N=2", "N=3");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorMseries M1(1,beta,measurelist[M]);
                    EnergyCorrelatorMseries M2(2,beta,measurelist[M]);
                    EnergyCorrelatorMseries M3(3,beta,measurelist[M]);


                    printf("%7.3f %14.6f %14.6f %14.6f \n",beta,M1(myJet),M2(myJet),M3(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;


                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorM2:  M2(beta) = ECFG(3, 1, beta)/ECFG(3, 1, beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta", "M2 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorM2 M2(beta,measurelist[M]);

                    printf("%7.3f %14.6f \n",beta,M2(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorCseries:  C_i(beta) = ECFN(i-1, beta)*ECFN(i+1, beta)/ECFN(i, beta)^2 with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %20s %20s %20s \n","beta", "N=1", "N=2", "N=3");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorCseries C1(1,beta,measurelist[M]);
                    EnergyCorrelatorCseries C2(2,beta,measurelist[M]);
                    EnergyCorrelatorCseries C3(3,beta,measurelist[M]);


                    printf("%7.3f %20.10f %20.10f %20.10f \n",beta,C1(myJet),C2(myJet),C3(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;


                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorUseries:  U_i(beta) = ECFG(i+1, 1, beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %20s %20s %20s \n","beta", "N=1", "N=2", "N=3");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorUseries U1(1,beta,measurelist[M]);
                    EnergyCorrelatorUseries U2(2,beta,measurelist[M]);
                    EnergyCorrelatorUseries U3(3,beta,measurelist[M]);


                    printf("%7.3f %20.10f %20.10f %20.10f \n",beta,U1(myJet),U2(myJet),U3(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorU1:  U1(beta) = ECFG(2, 1, beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta", "U1 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorU1 U1(beta,measurelist[M]);

                    printf("%7.3f %14.10f \n",beta,U1(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorU2:  U2(beta) = ECFG(3, 1, beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta", "U2 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorU2 U2(beta,measurelist[M]);

                    printf("%7.3f %14.10f \n",beta,U2(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;

                cout << "-------------------------------------------------------------------------------------" << endl;
                cout << "EnergyCorrelatorU3:  U3(beta) = ECFG(4, 1, beta) with " << modename[M] << endl;
                cout << "-------------------------------------------------------------------------------------" << endl;
                printf("%7s %14s \n","beta", "U3 obs");

                for (unsigned int B = 0; B < betalist.size(); B++) {
                    double beta = betalist[B];

                    EnergyCorrelatorU3 U3(beta,measurelist[M]);

                    printf("%7.3f %14.10f \n",beta,U3(myJet));
                }
                cout << "-------------------------------------------------------------------------------------" << endl << endl;


                // timing tests for the developers
                double do_timing_test = false;
                if (do_timing_test) {

                    cout << "jet with pt = " << myJet.pt() << " and " << myJet.constituents().size() << " constituents" << endl;

                    clock_t clock_begin, clock_end;
                    double num_iter;
                    double beta = 0.5;

                    cout << setprecision(6);

                    // test C1
                    num_iter = 20000;
                    clock_begin = clock();
                    EnergyCorrelatorDoubleRatio C1s(1,beta,measurelist[M],EnergyCorrelator::slow);
                    EnergyCorrelatorDoubleRatio C1f(1,beta,measurelist[M],EnergyCorrelator::storage_array);
                    cout << "timing " << C1s.description() << endl;
                    cout << "timing " << C1f.description() << endl;
                    for (int t = 0; t < num_iter; t++) {
                        C1s(myJet);
                    }
                    clock_end = clock();
                    cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C1"<< endl;

                    num_iter = 20000;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        C1f(myJet);
                    }
                    clock_end = clock();
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C1"<< endl;


                    // test C2
                    num_iter = 1000;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M],EnergyCorrelator::slow);
                        C2(myJet);
                    }
                    clock_end = clock();
                    cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C2"<< endl;

                    num_iter = 10000;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorDoubleRatio C2(2,beta,measurelist[M],EnergyCorrelator::storage_array);
                        C2(myJet);
                    }
                    clock_end = clock();
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C2"<< endl;

                    // test C3
                    num_iter = 100;
                    clock_begin = clock();

                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M],EnergyCorrelator::slow);
                        C3(myJet);
                    }
                    clock_end = clock();
                    cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C3"<< endl;

                    num_iter = 3000;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorDoubleRatio C3(3,beta,measurelist[M],EnergyCorrelator::storage_array);
                        C3(myJet);
                    }
                    clock_end = clock();
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C3"<< endl;

                    // test C4
                    num_iter = 10;
                    clock_begin = clock();

                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorDoubleRatio C4(4,beta,measurelist[M],EnergyCorrelator::slow);
                        C4(myJet);
                    }
                    clock_end = clock();
                    cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C4"<< endl;

                    num_iter = 300;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorDoubleRatio C4(4,beta,measurelist[M],EnergyCorrelator::storage_array);
                        C4(myJet);
                    }
                    clock_end = clock();
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per C4"<< endl;

                    // test N2
                    num_iter = 10;
                    clock_begin = clock();

                    num_iter = 300;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorN2 N2(beta,measurelist[M],EnergyCorrelator::storage_array);
                        N2(myJet);
                    }
                    clock_end = clock();
                    EnergyCorrelatorN2 N2test(beta,measurelist[M],EnergyCorrelator::storage_array);
                    cout << "Beta is: "<< beta << endl;
                    cout << "Result of N2: "<< N2test(myJet) << endl;
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per N2"<< endl;


                    num_iter = 300;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorN3 N3(beta,measurelist[M],EnergyCorrelator::storage_array);
                        N3(myJet);
                    }
                    clock_end = clock();
                    EnergyCorrelatorN3 N3test(beta,measurelist[M],EnergyCorrelator::storage_array);
                    cout << "Beta is: "<< beta << endl;
                    cout << "Result of N3: "<< N3test(myJet) << endl;
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per N3"<< endl;




                    num_iter = 300;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorGeneralized ECF1(2,3, beta,  measurelist[M]);
                        ECF1(myJet);
                    }
                    clock_end = clock();
                    EnergyCorrelatorGeneralized ECF1test(2,3, beta,  measurelist[M]);
                    cout << "Beta is: "<< beta << endl;
                    cout << "Result of 2e3: "<< ECF1test(myJet) << endl;
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per 2e3"<< endl;


                    num_iter = 300;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorGeneralized ECF3(2,4, beta,  measurelist[M]);
                        ECF3(myJet);
                    }
                    clock_end = clock();
                    EnergyCorrelatorGeneralized ECF2test(2,4, beta,  measurelist[M]);
                    cout << "Beta is: "<< beta << endl;
                    cout << "Result of 2e4: "<< ECF2test(myJet) << endl;
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per 2e4"<< endl;


//                    num_iter = 300;
//                    clock_begin = clock();
//                    for (int t = 0; t < num_iter; t++) {
//                        EnergyCorrelatorGeneralized ECF5(2,5, beta,  measurelist[M]);
//                        ECF5(myJet);
//                    }
//                    clock_end = clock();
//                    EnergyCorrelatorGeneralized ECF5test(2,5, beta,  measurelist[M]);
//                    cout << "Beta is: "<< beta << endl;
//                    cout << "Result of 2e5: "<< ECF5test(myJet) << endl;
//                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per 2e5"<< endl;
//

                    // test M2
                    num_iter = 10;
                    clock_begin = clock();

                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorM2 M2(beta,measurelist[M],EnergyCorrelator::slow);
                        M2(myJet);
                    }
                    clock_end = clock();
                    cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per M2"<< endl;

                    num_iter = 300;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorM2 M2(beta,measurelist[M],EnergyCorrelator::storage_array);
                        M2(myJet);
                    }
                    clock_end = clock();
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per M2"<< endl;

                    // test M3
                    num_iter = 10;
                    clock_begin = clock();

                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorMseries M3(3,beta,measurelist[M],EnergyCorrelator::slow);
                        M3(myJet);
                    }
                    clock_end = clock();
                    cout << "Slow method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per M3"<< endl;

                    num_iter = 300;
                    clock_begin = clock();
                    for (int t = 0; t < num_iter; t++) {
                        EnergyCorrelatorMseries M3(3,beta,measurelist[M],EnergyCorrelator::storage_array);
                        M3(myJet);
                    }
                    clock_end = clock();
                    cout << "Storage array method: " << (clock_end-clock_begin)/double(CLOCKS_PER_SEC*num_iter)*1000 << " ms per M3"<< endl;


                }
            }
        }
    }
}



