/****************************************************************
*****************************************************************
LSN_Exercises_08
*****************************************************************
*****************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "./prng/random.h"
#include "stats.h"

using namespace std;

// psi_trial (up to a normalization constant)
double trial_minimizer(double x, double m, double s) {
    return (exp(-(0.5 * pow((m - x), 2)) / pow(s, 2)) + exp(-(0.5 * pow((m + x), 2)) / pow(s, 2)));
};
// square modulus
double square_modulus(double x, double m, double s) {
    return pow(trial_minimizer(x, m, s), 2);
};
// second derivative
double sec_der(double x, double m, double s) {
    return ( pow(s,-4.)*( (m*m-2.*m*x-s*s+x*x)*exp(-(pow(x-m,2.))/(2*s*s)) + (m*m+2.*m*x-s*s+x*x)*exp(-(pow(x+m,2.))/(2*s*s)) ) );
};
// kinetic term
double kinetic(double x, double m, double s, double hbar, double mass) {
    return pow(trial_minimizer(x,m,s),-1.)*((-0.5)*pow(hbar,2)/mass)*sec_der(x, m, s);
};
// potential (well)
double potential(double x) {
    return pow(x,4.)-2.5*pow(x,2.);
};
// total energy 
double energy(double x, double m, double s, double hbar, double mass) {
    return (kinetic(x,m,s,hbar,mass)+potential(x));
};
// set units and mass
double hbar = 1.;
double mass = 1.;

int main(int argc, char *argv[]) {

    /****************************************************************************
    setup
    *****************************************************************************/

    // a generator is constructed
    Random rnd;
    int seed[4];
    int p1, p2;
    // read primes
    ifstream Primes("./prng/Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
    }
    Primes.close();
    // take in seeds
    ifstream input("./prng/seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
    }

    /****************************************************************************
     * Part 1
    Variational Monte Carlo code 
    M(RT)^2 algorithm 
    sample the square modulus of trial wave function $|\Psi_T^{\sigma,\mu}(x)|^2$ 
    using uniform transition probability $T(x_{new}|x_{old})$ for new candidates
    *****************************************************************************/

    // set trial parameters
    double m = 1.; // mean
    double s = 0.5;// standard deviation
    // sampling
    double q = 0.; // metropolis starting point
    vector<double> Q; // store the Markov chain
    int refused_moves = 0;
    int tot_n_of_proposals = 0;
    int wanted_moves = pow(10,4);
    double width = 2.5;

    cout << "\n=============== Part 1 =================== " << endl;
    cout << "============= do histogram ===============\n " << endl;

    for (int j = 0; j < wanted_moves; j=j) {
        double q_new_proposed = (q - 0.5 * width + width * rnd.Rannyu()); // propose a move anyway
        tot_n_of_proposals++;
        if (rnd.Rannyu() < min(1.0, square_modulus(q_new_proposed, m, s) / square_modulus(q, m, s))) {   
            q = q_new_proposed;
            Q.push_back(q);
            j++;
        } else {
            refused_moves++;
        }
    }
    Print( Q , "extracted.dat" );
    cout << "Total number of accepted moves: " << tot_n_of_proposals-refused_moves << endl; // should be = wanted_moves ...
    double acceptance = static_cast<double>(tot_n_of_proposals-refused_moves) / static_cast<double>(tot_n_of_proposals);
    cout << "Metropolis acceptance = " << acceptance * 100.0 << "%" << endl;

    /****************************************************************************
     * Part 2
    by using data blocking, the code computes the expectation value for the Hamiltonian
    *****************************************************************************/

    cout << "\n=============== Part 2 =================== " << endl;
    cout << "===== find energy for fixed mu,sigma =====\n " << endl;

    int number_of_blks = 200; // i.e. estimate the integral for number_of_blks times
    vector<double> integrals; // store the results
    for (int g=0; g<number_of_blks; g++) {
        // STORE A CHAIN of x on which to sample the integral
        q = rnd.Rannyu(m-s,m+s); // metropolis starting point (in the middle of a hill)
        Q.clear(); // store the Markov chain
        int throws_per_integration_block = 500;
        width = 2.5;
        for (int j = 0; j < throws_per_integration_block; j=j) {
            double q_new_proposed = (q - 0.5 * width + width * rnd.Rannyu()); // propose a move anyway
            if (rnd.Rannyu() < min(1.0, square_modulus(q_new_proposed, m, s) / square_modulus(q, m, s))) {   
                q = q_new_proposed;
                Q.push_back(q);
                j++;
            }
        }
        // COMPUTE THE INTEGRAL (mean)
        double integral = 0;
        for (int i = 0; i < throws_per_integration_block; i++) {
            integral = static_cast<double>(i)/static_cast<double>(i+1)*integral + energy(Q[i],m,s,hbar,mass)/static_cast<double>(i+1); 
        }
        integrals.push_back(integral); // save to vector
        if ((g+1)%20==0) {
            cout << "Block "<< g+1 << endl;
        }
    }
    // print integrals, pure as they are from the blocks
    cout << "Saving block results to file... " << endl;
    Print(integrals,"integrals.dat"); 

    // DATA BLOCKING
    cout << "Data blocking... " << endl;
    // let us build a progression of increasing blocks
    vector<double> ave_ene_integr;
    vector<double> sigma_ene_integr;
    for (int i=2; i<=number_of_blks; i++) {
        vector<double> truncated_ene_integr(integrals.begin(), integrals.begin()+i); // the vector is sliced up
        ave_ene_integr.push_back(Average(truncated_ene_integr));
        sigma_ene_integr.push_back(sqrt(Variance(truncated_ene_integr)/static_cast<double>(i-1)));
    }
    // save to file
    cout << "Saving progressive averages and errors to file... " << endl;
    Print( ave_ene_integr , "ave_ene_integr.dat" );
    Print( sigma_ene_integr , "sigma_ene_integr.dat" );

    /****************************************************************************
    * Part 3
    simulated annealing
    *****************************************************************************/

    cout << "\n=============== Part 3 ===================" << endl; 
    cout <<   "========== (SA) param. optimiz. ==========\n" << endl; 

    // starting parameters
    double starting_m = 1.0;
    double m_width = 0.2; // candidates extraction
    double starting_s = 0.3;
    double s_width = 0.2; // candidates extraction
    // store sequence of parameters (to be plotted in 2D plane later)
    vector<double> m_sequence;
    vector<double> s_sequence;

    // simulated annealing
    double temperature = 5.;
    double beta = 1./temperature; // inverse temperature
    double final_temp = 0.5*pow(10.,-1); // desired (low) temperature
    double cooling_factor = 0.97; // must multiply temp by c_fact to lower temp
    int necessary_cooling_steps = 1 + log(final_temp/temperature)/log(cooling_factor);
    int inside_steps = 500; 
    int SA_length = inside_steps*necessary_cooling_steps;

    cout << "Simulated annealing requires "<< necessary_cooling_steps <<" cooling steps"<< endl;

    // energy (for SA)
    double old_energy = 0.;
    double new_energy = 0.;
    vector<double> config_ene_sequence; // to be plotted later
    // establish the number of throws for estimating energy
    int tot_throws = 300*500;   // i.e. estimate the integral for tot_throws times
                                // NOTE: the first part of the exercise ensures that 
                                //       with such a number of throws the uncertainty will be less than 2%
    
    // ======================================================
    //        OLD ENERGY CALCULATION    =====================
    // ======================================================

    // STORE A CHAIN of x on which to sample the integral
    q = rnd.Rannyu(starting_m-starting_s,starting_m+starting_s); // metropolis starting point (in the middle of a hill)
    Q.clear(); // store the Markov chain
    width = 2.6;
    for (int j = 0; j < tot_throws; j=j) {
        double q_new_proposed = (q - 0.5 * width + width * rnd.Rannyu()); // propose a move anyway
        if (rnd.Rannyu() < min(1.0, square_modulus(q_new_proposed, starting_m, starting_s) / square_modulus(q, starting_m, starting_s))) {   
            q = q_new_proposed;
            Q.push_back(q);
            j++;
        }
    }
    // COMPUTE THE INTEGRAL (mean), i.e. THE ENERGY OF THE CONFIGURATION
    for (int i = 0; i < tot_throws; i++) {
        old_energy = static_cast<double>(i)/static_cast<double>(i+1)*old_energy + energy(Q[i],starting_m,starting_s,hbar,mass)/static_cast<double>(i+1); 
    }

    // ======================================================
    // END OF OLD ENERGY CALCULATION    =====================
    // ======================================================

    // ======================================================
    // ======================================================
    // build the chain of parameters: metropolis    
    // ======================================================
    // ======================================================
    double TOTPROPOSALS = 0;
    for (int l = 0; l < SA_length; l=l) {
        if ( (l+1)%inside_steps == 0 ) {
            beta = beta / cooling_factor;
            cout << "SA step: " << (l+1)/inside_steps << endl;
        }
        TOTPROPOSALS++;
        // PROPOSALS // propose a move anyway
        double proposed_m = (starting_m - 0.5 * m_width + m_width * rnd.Rannyu()); 
        double proposed_s = 0;
        do {
            proposed_s = (starting_s - 0.5 * s_width + s_width * rnd.Rannyu());
        } while (proposed_s==0); // sigma can never be zero!

        // ======================================================
        // COST FUNCTION COMPUTATION ============================
        // ======================================================

        // STORE A CHAIN of x on which to sample the integral
        q = rnd.Rannyu(proposed_m-proposed_s,proposed_m+proposed_s); // metropolis starting point (in the middle of a hill)
        Q.clear(); // store the Markov chain
        width = 2.6;
        for (int j = 0; j < tot_throws; j=j) {
            double q_new_proposed = (q - 0.5 * width + width * rnd.Rannyu()); // propose a move anyway
            if (rnd.Rannyu() < min(1.0, square_modulus(q_new_proposed, proposed_m, proposed_s) / square_modulus(q, proposed_m, proposed_s)) ) {   
                q = q_new_proposed;
                Q.push_back(q);
                j++;
            }
        }
        // COMPUTE THE INTEGRAL (mean), i.e. THE ENERGY OF THE CONFIGURATION
        double integral = 0;
        for (int i = 0; i < tot_throws; i++) {
            integral = static_cast<double>(i)/static_cast<double>(i+1)*integral + energy(Q[i],proposed_m,proposed_s,hbar,mass)/static_cast<double>(i+1); 
        }
        new_energy = integral;

        // ======================================================
        // END OF COST FUNCTION COMPUTATION =====================
        // ======================================================

        // decide whether to accept or not the new mu and sigma
        if ( rnd.Rannyu() < min(1.0, exp( beta*(-new_energy+old_energy) ) ) ) {
            //cout << "exp = "<< exp( beta*(-new_energy+old_energy) ) << endl;   
            starting_m = proposed_m;
            starting_s = proposed_s;
            old_energy = new_energy;
            if ( (l+1)%inside_steps == 0 ) {
                s_sequence.push_back(starting_s);
                m_sequence.push_back(starting_m);
                config_ene_sequence.push_back(old_energy);
            }
            l++;
        }
    }
    cout << "Acceptance = " <<static_cast<double>(SA_length)/TOTPROPOSALS << endl;
    cout << "Best mu = " << m_sequence.back() << endl;
    cout << "Best sigma = " << s_sequence.back() << endl;    

    // print the parameters obtained
    Print(m_sequence,"m_sequence.dat");
    Print(s_sequence,"s_sequence.dat");
    Print(config_ene_sequence,"e_sequence.dat");
    
    /****************************************************************************
    * Part 4
    use parameters obtained
    *****************************************************************************/
    
    cout << "\n=============== Part 4 ===================" << endl; 
    cout <<   "=========  use  of  parameters  ==========\n" << endl; 

    // use the parameters found
    m = m_sequence.back();
    s = s_sequence.back();
    
    number_of_blks = 200; // i.e. estimate the integral for number_of_blks times
    integrals.clear(); // store the results
    for (int g=0; g<number_of_blks; g++) {
        // STORE A CHAIN of x on which to sample the integral
        q = rnd.Rannyu(m-s,m+s); // metropolis starting point (in the middle of a hill)
        Q.clear(); // store the Markov chain
        int throws_per_integration_block = 500;
        width = 2.5;
        for (int j = 0; j < throws_per_integration_block; j=j) {
            double q_new_proposed = (q - 0.5 * width + width * rnd.Rannyu()); // propose a move anyway
            if (rnd.Rannyu() < min(1.0, square_modulus(q_new_proposed, m, s) / square_modulus(q, m, s))) {   
                q = q_new_proposed;
                Q.push_back(q);
                j++;
            }
        }
        // COMPUTE THE INTEGRAL (mean)
        double integral = 0;
        for (int i = 0; i < throws_per_integration_block; i++) {
            integral = static_cast<double>(i)/static_cast<double>(i+1)*integral + energy(Q[i],m,s,hbar,mass)/static_cast<double>(i+1); 
        }
        integrals.push_back(integral);
        if ((g+1)%20==0) {
            cout << "Block "<< g+1 << endl;
        }
    }
    // print integrals, pure as they are from the blocks
    cout << "Saving block results to file... " << endl;
    Print(integrals,"final_integrals.dat");  // DIFFERENT NAME!
    // DATA BLOCKING
    cout << "Data blocking... " << endl;
    // let us build a progression of increasing blocks
    ave_ene_integr.clear();
    sigma_ene_integr.clear();
    for (int i=2; i<=number_of_blks; i++) {
        vector<double> truncated_ene_integr(integrals.begin(), integrals.begin()+i); // the vector is sliced up
        ave_ene_integr.push_back(Average(truncated_ene_integr));
        sigma_ene_integr.push_back(sqrt(Variance(truncated_ene_integr)/static_cast<double>(i-1)));
    }
    // save to file
    cout << "Saving progressive averages and errors to file... " << endl;
    Print( ave_ene_integr , "ave_final_integr.dat" );
    Print( sigma_ene_integr , "sigma_final_integr.dat" );
    
    /****************************************************************************
    closing
    *****************************************************************************/
    cout << "\n================= END ====================" << endl; 
    rnd.SaveSeed();
    return 0;
};

/****************************************************************************
*****************************************************************************/