//
//  Dynamics.cpp
//  
//
//  Created by Roberto Moran Tovar on 12.03.21.
//
//Template to run a stochastic/deterministic simulation of the antigen and bcells dynamics.

#include "../library/Immuno_functions.hpp"

#include <stdio.h>

//----------------------------------------------------------------------------------
using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[]) //argv has 1:L 2:N , 3:T , 4:T0 , 5:alpha , 6:beta , 7:gamma , 8:linear ; 9:energy_model
{
    string Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/Single_trajectory/";
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
    int T0 = atoi(argv[4]); //initial number of days for the simulation
    double dT = 0.005; //time step
    long long int NT = (T-T0)/dT; //number of steps
    long long A_0 = exp(alpha*T0);
    int linear = atoi(argv[8]);
    string energy_model (argv[9]);

    //------------Energy Matrix------------------------------------------------------
    vector < vector < double > > MJ;
    MJ.resize(L_alphabet);
    for (int k= 0; k<L_alphabet; k++)
    {
        (MJ[k]).resize(L_alphabet);
    };

    ifstream file("../Input_files/MJ2.txt");
    for (unsigned int i = 0; i < L_alphabet; i++) {
        for (unsigned int j = 0; j < L_alphabet; j++) {
            file >> MJ[i][j];
        }
    }

    //------------ Alphabet ----------------------------------------------------------
    //Array with the Alphabet
    vector < string > Alphabet;
    Alphabet.resize(L_alphabet);
    ifstream file2("../Input_files/Alphabet.txt");
    cout << "The Alphabet is :";
    for (int k = 0; k < L_alphabet; k++) {

        file2 >> Alphabet[k];
        cout << Alphabet[k] ;
    
    }
    cout << "\n";
    
    //------------- Antigen -------------------------------------------------------------
    //Array with the antigen
    string Antigen_aa;
    cout << "Insert the Aminoacid sequence of the antigen:\n";
    getline(cin, Antigen_aa);
    
    vector < int > Antigen;
    Antigen.resize(L);
    aa_to_positions(L, L_alphabet, Alphabet, Antigen, Antigen_aa);
    
    //---------Generating Bcells ---------------------------------------------------------
    //Array with Bcells
    vector < bcell > Bcells;
    Bcells.resize(N);
    //generate_Bcells(N, L, L_alphabet, Bcells);
    generate_Bcells_with_e(N, L, L_alphabet, Bcells, MJ, Antigen, energy_model, r);

    //---------Choosing antigen-specific Bcells ---------------------------------------------------------
    //Array with Naive-specific Bcells
    vector < bcell* > Naive;
    int n_naive = 0;
    choose_naive_Bcells2(N, L, L_alphabet, MJ, Antigen, Bcells, Naive, n_naive, energy_model, r);
    
    //Matrix with the time series of the antigen-specific Bcells
    vector<vector < long double > > Time_series_Bcells;
    Time_series_Bcells.resize(n_naive);
    for(int n= 0; n<n_naive; n++)
    {
        Time_series_Bcells[n].resize(NT);
        Time_series_Bcells[n][0] = Naive[n]->cs;
    };
    
    
    //Array with time series of the antigen
    vector < long double > Time_series_Antigen;
    Time_series_Antigen.resize(NT);
    if (linear==0) {
        Time_series_Antigen[0] = A_0;
    } else {
        Time_series_Antigen[0] = 1;// I used 10E2 before 
    }
    
    
    //Array for time series of the number of active bcell linages
    vector <int> N_active_linages;
    N_active_linages.resize(NT);
    
    
    cout << "Reactive B cells:" << n_naive << endl;
    
    //Output files
    
 
    ofstream fout (Text_files_path+"energies_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_Linear-"+std::to_string(linear)+"_"+energy_model+".txt");
    ofstream fout_antigen (Text_files_path+"antigen_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear)+"_"+energy_model+".txt");
    ofstream fout_bcells (Text_files_path+"bcells_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear)+"_"+energy_model+".txt");
    ofstream fout_N_active_linages (Text_files_path+"N_active_linages_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear)+"_"+energy_model+".txt");

    
    cout << "e_bar:" << mean_energy(L, L_alphabet, MJ, Antigen) << endl;
    // Run ODE
    ODE(linear, alpha, beta, gamma, NT, dT, n_naive, Naive, Time_series_Bcells, Time_series_Antigen, N_active_linages);
    
    for (int n= 0; n<n_naive; n++)
    {
        fout << Naive[n]->e << "\t" << Naive[n]->active << "\t" << Naive[n]->GC << endl;
    };
    
    //Print time series of antigen and bcells
    for(int t=0 ; t<NT; t++){
        fout_antigen << Time_series_Antigen[t] << endl;
        fout_N_active_linages << N_active_linages[t] << "\t";
        for (int n = 0 ; n<n_naive ; n++){
            fout_bcells << Time_series_Bcells[n][t] << "\t";
        }
        fout_bcells << endl;
    }
        
    
    fout.close();
    fout_antigen.close();
    fout_bcells.close();
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;

    return 0;
}
