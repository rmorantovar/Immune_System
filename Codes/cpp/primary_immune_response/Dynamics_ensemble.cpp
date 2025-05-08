//
//  Dynamics_ensemble.cpp
//  
//
//  Created by Roberto Moran Tovar on 13.03.21.
//
//Template to run a stochastic/deterministic simulation of the antigen and bcells dynamics.

#include "../library/Immuno_functions.hpp"

#include <stdio.h>
#include <chrono>

//----------------------------------------------------------------------------------
using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[]) //argv has 1:L , 2:N , 3:T , 4:T0 , 5:alpha , 6:beta , 7:gamma , 8:N_ensemble , 9:linear , 10:type
{
    string Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Ensemble/";
    string Text_files_path2 = "/Users/robertomorantovar/Documents/Research/Evolution_Immune_System/Codes/";
    cout<<">Running simulation of the Bcells-Antigen dynamics ..."<< endl;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(r, time(NULL));
    clock_t t1,t2;
    t1=clock();
    //-----------------------------------------------------------------------------
    //Parameters:
    std::string alpha_s (argv[5]);
    double alpha = stod(alpha_s);
    std::string beta_s (argv[6]);
    double beta = stod(beta_s);;
    std::string gamma_s (argv[7]);
    double gamma = stod(gamma_s);;
    int L  = atoi(argv[1]); //length of the sequence
    int L_alphabet (20); //length of the alphabet
    long long int N = atoi(argv[2]); // number of bcells
    int T = atoi(argv[3]); //number of days for the simulation
    int T0 = atoi(argv[4]); //number of days for the simulation
    double dT = 0.5; //time step
    long long int NT = (T-T0)/dT; //number of steps
    long long int N_ensemble = atoi(argv[8]);
    long long A_0 = exp(alpha*T0);
    int linear = atoi(argv[9]);
    std::string type (argv[10]);
    

    //------------Energy Matrix------------------------------------------------------
    vector < vector < double > > MJ;
    MJ.resize(L_alphabet);
    for (int k= 0; k<L_alphabet; k++)
    {
        (MJ[k]).resize(L_alphabet);
    };

    ifstream file(Text_files_path2+"Input_files/MJ2.txt");

    //------------ Alphabet ----------------------------------------------------------
    //Array with the Alphabet
    vector < string > Alphabet;
    Alphabet.resize(L_alphabet);
    ifstream file2(Text_files_path2+"Input_files/Alphabet.txt");
    cout << "The Alphabet is :";
    for (int k = 0; k < L_alphabet; k++) {

        file2 >> Alphabet[k];
        cout << Alphabet[k] ;
    
    }
    cout << "\n";
    for (unsigned int i = 0; i < L_alphabet; i++) {
        for (unsigned int j = 0; j < L_alphabet; j++) {
            file >> MJ[i][j];
        }
    }
    //------------- Antigen -------------------------------------------------------------
    //Array with the antigen
    //string Antigen_aa;
    //cout << "Insert the Aminoacid sequence of the antigen:\n";
    //getline(cin, Antigen_aa);
    std::string Antigen_aa (argv[11]);
    
    vector < int > Antigen;
    Antigen.resize(L);
    aa_to_positions(L, L_alphabet, Alphabet, Antigen, Antigen_aa);
    
    //---------Generating Bcells ---------------------------------------------------------
    //Array with Bcells
    vector < bcell > Bcells;
    Bcells.resize(N);
    
    //Array with time series of the antigen
    vector < long double > Time_series_Antigen;
    Time_series_Antigen.resize(NT);
    
    //---------Activated linages ---------------------------------------------------------
    //Array for time series of the number of active bcell linages per time
    vector <double> N_active_linages;
    N_active_linages.resize(NT);
    //Array for total final number of active bcell linages
    vector <int> N_final_active_linages;
    //N_final_active_linages.resize(N_ensemble);
    
    //Output files
    ofstream fout (Text_files_path+"energies_ensemble_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_Linear-"+std::to_string(linear)+"_"+type+".txt"); // Energies
    
    ofstream fout_bcells (Text_files_path+"bcells_ensemble_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear)+"_"+type+".txt"); // B cells final clone size
    
    ofstream fout_N_final_active (Text_files_path+"N_final_active_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear)+"_"+type+".txt"); // B cells final clone size
    
    //print in file the time series of the average of the number of activated bcell linages.
    ofstream fout_N_active_linages (Text_files_path+"N_active_linages_ensemble_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear)+"_"+type+".txt");

    // ------------ Run ensemble of trajectories ------------
    cout << "Running ensemble of trajectories ..." << endl;
    for(int i_ensemble = 0 ; i_ensemble<N_ensemble ; i_ensemble++){
        
        //Generate bcells
        generate_Bcells(N, L, L_alphabet, Bcells);
        
        // Choose the antigen-specific bcells
        vector < bcell* > Naive;
        int n_naive = 0;
        choose_naive_Bcells(N, L, L_alphabet, MJ, Antigen, Bcells, Naive, n_naive, type, r);
        
        //initialize time series arrays
        if (linear==0) {
            Time_series_Antigen[0] = A_0;
        } else {
            Time_series_Antigen[0] = 1e3;
        }
        
        //Matrix with the time series of the antigen-specific Bcells
        /*vector<vector < long double > > Time_series_Bcells;
        Time_series_Bcells.resize(n_naive);
        for(int n= 0; n<n_naive; n++)
        {
            Time_series_Bcells[n].resize(NT);
            Time_series_Bcells[n][0] = Naive[n]->cs;
        };
        */
        // Run ODE
        ODE_ensemble(linear, alpha, beta, gamma, NT, dT, n_naive, Naive, Time_series_Antigen, N_active_linages, N_final_active_linages);
        
        for (int n = 0 ; n<n_naive ; n++){
            //print in file the energies and the activation state of the antigen-specific bcells.
            fout << Naive[n]->e << "\t" << Naive[n]->active << endl;
            //Print the final clone-size of bcells
            if(Naive[n]->active==1){
                fout_bcells << Naive[n]->cs << endl;
            }
        }
            
    }
        
    for (int t= 0; t<NT; t++)
    {
        fout_N_active_linages << N_active_linages[t]/N_ensemble << "\t";
    };
    
    for (int n=0; n<N_ensemble; n++){
        fout_N_final_active << N_final_active_linages[n] << "\t";
    }
    
    
    fout.close();
    fout_bcells.close();
    fout_N_active_linages.close();
    fout_N_final_active.close();
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
    return 0;
}
