    /****************************************************************************
    * Part 3
    simulated annealing
    *****************************************************************************/
   
    cout << "\n=============== Part 3 ===================" << endl; 
    cout <<   "========== (SA) param. optimiz. ==========\n" << endl; 

    // ======================================================
    //  PARAMS AND VAR INITIALIZATION   =====================
    // ======================================================
    
    // starting parameters
    double starting_m = 1.;
    double m_width = 0.25; // candidates extraction
    double starting_s = 0.3;
    double s_width = 0.25; // candidates extraction
    // simulated annealing params
    double temperature = 5.;
    double beta = 1./temperature; // inverse temperature
    double final_temp = pow(10.,-2); // desired (low) temperature
    double cooling_factor = 0.99; // must multiply temp by c_fact to lower temp
    int necessary_cooling_steps = 1 + log(final_temp/temperature)/log(cooling_factor);
    int inside_steps = 100; 
    int SA_length = inside_steps*necessary_cooling_steps;
    int TOTPROPOSALS = 0;

    cout << "Simulated annealing requires "<< necessary_cooling_steps <<" cooling steps"<< endl;

    // store sequence of parameters (to be plotted in 2D plane later)
    vector<double> m_sequence;
    vector<double> s_sequence;
    // energy (for SA)
    double old_energy = 0.;
    double new_energy = 0.;
    // establish the number of throws for estimating energy
    int energy_block_length = 500;
    int energy_blocks = 300;
    int tot_throws = energy_block_length*energy_blocks;   // i.e. estimate the integral with tot_throws points
                                // NOTE: the first part of the exercise ensures that 
                                //       with such a number of throws the uncertainty will be less than 2%
    
    // ======================================================
    //        OLD ENERGY CALCULATION    =====================
    // ======================================================

    // STORE A CHAIN of x on which to sample the integral
    q = rnd.Rannyu(starting_m-starting_s,starting_m+starting_s); // metropolis starting point (in the middle of a hill)
    width = 7.0;
    for (int j = 0; j < tot_throws; j=j) {
        double q_new_proposed = (q - 0.5 * width + width * rnd.Rannyu()); // propose a move anyway
        if (rnd.Rannyu() < min(1.0, square_modulus(q_new_proposed, starting_m, starting_s) / square_modulus(q, starting_m, starting_s))) {   
            q = q_new_proposed;
            old_energy = static_cast<double>(j)/static_cast<double>(j+1)*old_energy + energy(q,starting_m,starting_s,hbar,mass)/static_cast<double>(j+1); 
            j++;
        }
    }

    // ======================================================
    //    BUILD MU AND SIGMA CHAIN: M(RT)^2    ==============    
    // ======================================================

    // print on file while running
    ofstream file1;
    ofstream file2;
    ofstream file3;
    ofstream file4;
    file1.open("live_m.dat" , ios::out);
    file2.open("live_s.dat", ios::out);
    file3.open("live_energy.dat", ios::out);
    file4.open("live_energy_error.dat", ios::out);

    for (int l = 0; l < SA_length; l=l) {

        // ======================================================
        // COOLING                   ============================
        // ======================================================

        if ( (l+1)%inside_steps == 0 ) { // cools down ONLY every inside_steps steps
            beta = beta/cooling_factor;// to allow for thermalization / burn-in period
            cout << "================ SA cooling step: " << (l+1)/inside_steps << endl;
        }
        cout << "Inside step: " << (l+1)%inside_steps << endl;

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

        vector<double> integral_array_blk; // auxiliary variable for each energy block (dimension = energy_blocks)

        for (int h = 0; h < energy_blocks ; h++ ) {
            // STORE A CHAIN of x on which to sample the integral
            q = rnd.Rannyu(proposed_m-proposed_s,proposed_m+proposed_s); // metropolis starting point (in the middle of a hill)
            width = 7.0;
            double integral = 0;
            for (int j = 0; j < energy_block_length; j=j) {
                double q_new_proposed = (q - 0.5 * width + width * rnd.Rannyu()); // propose a move anyway
                if (rnd.Rannyu() < min(1.0, square_modulus(q_new_proposed, proposed_m, proposed_s) / square_modulus(q, proposed_m, proposed_s)) ) {   
                    q = q_new_proposed;
                    if (rnd.Rannyu()<0.5) {
                        q = -q;
                    }
                    integral = static_cast<double>(j)/static_cast<double>(j+1)*integral + energy(q,proposed_m,proposed_s,hbar,mass)/static_cast<double>(j+1); 
                    j++;
                }
            }
            integral_array_blk.push_back(integral); // save the result of a single block
        }

        new_energy = Average(integral_array_blk); // once every block is done, take the average

        // ======================================================
        //   START ROAMING IN PARAMS SPACE  =====================
        // ======================================================

        // decide whether to accept or not the new mu and sigma
        if ( rnd.Rannyu() < min(1.0, exp(beta*(-new_energy+old_energy))) ) {
            //cout << "exp = "<< exp( beta*(-new_energy+old_energy) ) << endl;   
            starting_m = proposed_m; // update
            starting_s = proposed_s;
            old_energy = new_energy;
            if ( (l+1)%inside_steps == 0 ) {
                m_sequence.push_back(starting_m); // update
                s_sequence.push_back(starting_s); 
                double ene_error_istantaneous = sqrt(Variance(integral_array_blk)/(energy_block_length-1));
                file1 << starting_m << endl;
                file2 << starting_s << endl;
                file3 << old_energy << endl;
                file4 << ene_error_istantaneous << endl;
            }
            l++; // increasing j ONLY when move is accepted gives control over the number of moves!
        }
    }

    // ======================================================
    //            =====================
    // ======================================================

    cout << "Acceptance = " <<static_cast<double>(SA_length)/TOTPROPOSALS << endl;
    cout << "Best mu = " << m_sequence.back() << endl;
    cout << "Best sigma = " << s_sequence.back() << endl;  

    // close live update
    file1.close();
    file2.close();
    file3.close();
    file4.close();  
    // print the parameters obtained
    Print(m_sequence,"m_sequence.dat");
    Print(s_sequence,"s_sequence.dat");
