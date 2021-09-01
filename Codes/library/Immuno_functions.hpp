//
//  Immuno_functions.hpp
//  
//  Created by Roberto Moran Tovar on 27.02.21.
//

#ifndef Immuno_functions_h
#define Immuno_functions_h


#endif /* Immuno_functions_h */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>




//Library for random number generators
#include "./random.cpp"
//There are two functions extracted from the library
//double randX(min,max): a random number between min and max
//int randIX(min,max): an integer random number between min and max

const long double N_A = 6.02214076E23;

using namespace std;

// ---------------- CLASSES ----------------

class bcell {
public:
    vector < int > seq; //vector with the positions of the aa sequence in the alphabet
    bcell();
    bcell(int const & L, int const & L_alphabet, vector< int > & seq);
    double e; //energy with respect to the current epitope.
    double cs;
    bool plasma;
    bool GC;
    bool active;
};
bcell::bcell(){

}
bcell::bcell(int const & L, int const & L_alphabet, vector< int > & seq){
    this->seq = seq;
    cs = 1.0;
    plasma = 0;
    GC = 0;
    active = 0;
    
}

// ---------------- FUNCTION ---------------

//Function to calculate the energy: Implement the Energy Matrix
double Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen, string type, gsl_rng *r)
{

    double E = 0.0;
    if (type == "MJ") {
        for(int i=0; i<L ; i++){
            E = E + MJ[Antigen[i]][sequence[i]];
        }
    }else if(type == "Random"){
        E = -56.0;
        for(int i=0; i<L ; i++){
            E = E + gsl_ran_gaussian(r, 1.17);
        }
    }
    return E;

};

//Function to calculate the energy difference due to a mutation
inline double delt( int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen, int const & pos, int const & aa)
{
    double deltaE (0.);
    deltaE = MJ[Antigen[pos]][aa] - MJ[Antigen[pos]][sequence[pos]];
    return deltaE;
};

//Function to change from aminoacids to positions
void aa_to_positions( int const & L, int const & L_alphabet, vector< string > & Alphabet,  vector< int > & sequence_pos,  string  sequence_aa)
{
    for(int i=0; i<L ;i++){
        
        for(int j=0; j<L_alphabet ;j++){
            if(sequence_aa[i] == Alphabet[j][0]){
                sequence_pos[i] = j;
            }
        }
    }
};

//Function to calculate complementary sequence
void find_complementary(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector<int> & complementary_sequence)
{
    for(int i=0; i<L ; i++){
        vector < double > v;
        v.resize(L);
        v = MJ[sequence[i]];
        int index = std::min_element(v.begin(), v.end()) - v.begin();
        complementary_sequence[i] = index;
    }
};

// Function to return energies < E_max

void high_energy_filter(){
}

// Function to generate the initial amount of sequences
void generate_Bcells(int N, int L, int L_alphabet, vector<bcell> & Bcells){
    
    //---------Array with the current Sequence-------------------------------------------
    vector < int > Sequence;
    Sequence.resize(L);
    
    for(int n =0 ; n<N ; n++){
        
        //Initiating Sequence with random sequence------------------------------------
        for (int k= 0; k<L; k++)
        {
            Sequence[k] = randIX(0,L_alphabet-1);
        };
        
        // Create a bcell and add it to the vector
        bcell bcell_i(L, L_alphabet, Sequence);
        Bcells[n]  = bcell_i;
    }
}

// Function that selects the antigen-specific naive Bcells from all the sequences
void choose_naive_Bcells(int N, int L, int L_alphabet, vector< vector<double> > const & MJ, vector< int > const & Antigen, vector<bcell> & Bcells, vector<bcell*> & Naive, int & n_naive, string type, gsl_rng *r){
    
    vector <int> MS;
    MS.resize(L);
    find_complementary(L, L_alphabet, MJ, Antigen, MS);
    double e_MS = Energy(L, L_alphabet, MJ, MS, Antigen, "MJ", r);
    double e;
    for(int n = 0 ; n<N ; n++){
        e = Energy(L, L_alphabet, MJ, Bcells[n].seq, Antigen, type, r);
        if(e<e_MS+28){
            Bcells[n].e = e;
            Naive.push_back( &Bcells[n]);
            n_naive++;
        }
    }
}

//Function that mutates any sequences in one random position
void mutate_sequence(int L, int L_alphabet, vector< int > & sequence){
    int pos = randIX(0,L-1);
    int aa = randIX(0,L_alphabet-1);
    if(sequence[pos] != aa){
        sequence[pos]=aa;
    }
}

// Function to run de set of differential equations
void ODE(int linear, double const eta, double const nu, double const gamma, long long NT, double dT, int n_naive, vector<bcell*> & Naive, vector<vector < long double > > & Time_series_Bcells, vector < long double > & Time_series_Antigen, vector < int > & N_active_linages){
    double p_GC = 0.15;
    double f = 0;
    double N_active_bcells = 0;
    int n_active_linages = 0;
    if (linear==0) {
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (eta*Time_series_Antigen[t-1] - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (nu*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active);
                if(Naive[n]->active ==0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e));
                    if(f>0.5){
                        Naive[n]->active = 1;
                        n_active_linages++;
                        double r_GC = randX(0,1);
                        // Decide if the activated linage will have plasma or GC fate.
                        if(r_GC<p_GC){
                            Naive[n]->GC = 1;
                        }
                        else{
                            Naive[n]->plasma = 1;
                        }
                    }
                }else{
                    N_active_bcells = N_active_bcells + Time_series_Bcells[n][t];
                }
            }
            N_active_linages[t] = N_active_linages[t] + n_active_linages;
        }
    } else {
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (eta*2000 - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (nu*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active);
                if(Naive[n]->active ==0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e));
                    if(f>0.5){
                        Naive[n]->active = 1;
                        n_active_linages++;
                        double r_GC = randX(0,1);
                        // Decide if the activated linage will have plasma or GC fate.
                        if(r_GC<p_GC){
                            Naive[n]->GC = 1;
                        }
                        else{
                            Naive[n]->plasma = 1;
                        }
                    }
                }else{
                    N_active_bcells = N_active_bcells + Time_series_Bcells[n][t];
                }
            }
            N_active_linages[t] = N_active_linages[t] + n_active_linages;
        }
    }
    
}

// Function to run de set of differential equations
void ODE_ensemble(int linear, double const eta, double const nu, double const gamma, long long NT, double dT, int n_naive, vector<bcell*> & Naive, vector < long double > & Time_series_Antigen, vector < double > & N_active_linages, vector <double> & N_final_active_linages){
    double f = 0;
    double N_active_bcells = 0;
    int n_active_linages_t = 0;
    int n_active_linages = 0;
    if (linear==0) {
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (eta*Time_series_Antigen[t-1] - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            n_active_linages_t = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Naive[n]->cs = Naive[n]->cs + (nu*Naive[n]->cs*dT*(Naive[n]->active));
                // This function, contrary to the one for the single dynamics, does not use the Time_series_Bcells array. It uses the variable cs of the Bcell Class.
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (nu*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active); // this uses the time_series arrays
                if(Naive[n]->active == 0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e));
                    if(f>0.5){
                        Naive[n]->active = 1;
                        n_active_linages_t++;
                        n_active_linages++;
                    }
                }else{
                    N_active_bcells = N_active_bcells + Naive[n]->cs;
                }
            }
            N_active_linages[t] = N_active_linages[t] + n_active_linages_t;
        }
        N_final_active_linages.push_back(n_active_linages);
    } else {
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (eta*2000 - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            n_active_linages_t = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Naive[n]->cs = Naive[n]->cs + (nu*Naive[n]->cs*dT*(Naive[n]->active));
                // This function, contrary to the one for the single dynamics, does not use the Time_series_Bcells array. It uses the variable cs of the Bcell Class.
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (nu*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active); // this uses the time_series arrays
                if(Naive[n]->active == 0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e));
                    if(f>0.5){
                        Naive[n]->active = 1;
                        n_active_linages_t++;
                        n_active_linages++;
                    }
                }else{
                    N_active_bcells = N_active_bcells + Naive[n]->cs;
                }
            }
            N_active_linages[t] = N_active_linages[t] + n_active_linages_t;
        }
        N_final_active_linages.push_back(n_active_linages);
    }
    
}


