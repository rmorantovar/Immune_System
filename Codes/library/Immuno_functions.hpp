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
#include <valarray>
#include <numeric>
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
    bcell(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > & seq, vector< int > const & Antigen, string energy_model, gsl_rng *r);
    double e; //energy with respect to the current epitope.
    double cs; //clone size
    bool plasma; // plasma cell?
    //bool GC; // germinal center cell?
    //bool engaged; //engaged with an antigen?
    bool active; // activated clone?
    double activation_time; // time of activation

    // Functions
    double Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen, string energy_model, gsl_rng *r);

};
bcell::bcell(){
}
bcell::bcell(int const & L, int const & L_alphabet, vector< int > & seq)
{
    this->seq = seq;
    this->cs = 1.0;
    this->plasma = 1;
    //this->GC = 0;
    //this->engaged = 0;
    this->active = 0;
    this->activation_time = -1;
}
bcell::bcell(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > & seq, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{
    this->seq = seq;
    this->cs = 1.0;
    this->plasma = 1;
    //this->GC = 0;
    //this->engaged = 0;
    this->active = 0;
    this->activation_time = -1;
    this->e = this->Energy(L, L_alphabet, MJ, seq, Antigen, energy_model, r);
}
double bcell::Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{

    double E = 0.0;
    if (energy_model == "MJ") {
        for(int i=0; i<L ; i++){
            E = E + MJ[Antigen[i]][sequence[i]];
        }
    }else if(energy_model == "Random"){
        E = -56.0;
        for(int i=0; i<L ; i++){
            E = E + gsl_ran_gaussian(r, 1.17);
        }
    }
    return E;
};
// ---------------- FUNCTION ---------------
//Function to calculate the energy: Implement the Energy Matrix
double Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & MJ, vector< int > const & sequence, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{

    double E = 0.0;
    if (energy_model == "MJ") {
        for(int i=0; i<L ; i++){
            E = E + MJ[Antigen[i]][sequence[i]];
        }
    }else if(energy_model == "Random"){
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
//Function that mutates any sequences in one random position
void mutate_sequence(int L, int L_alphabet, vector< int > & sequence)
{
    int pos = randIX(0,L-1);
    int aa = randIX(0,L_alphabet-1);
    if(sequence[pos] != aa){
        sequence[pos]=aa;
    }
}
//Function to calculate the mean energy of the associated PWM
double mean_energy(int L, int L_alphabet, vector< vector<double> > const & MJ, vector< int > const & Antigen)
{

    double mean_E = 0.0;

    for(int i=0; i<L ; i++){
        for (int j = 0; j < L_alphabet; ++j)
        {
            mean_E = mean_E + MJ[Antigen[i]][j];
        }
    }
    mean_E=mean_E/(L_alphabet);

    return mean_E;
}
//Function to calculate the minimum energy of the associated PWM
double MS_energy(int L, int L_alphabet, vector< vector<double> > const & MJ, vector< int > const & Antigen)
{

    double min_E = 0.0;

    for(int i=0; i<L ; i++){

        vector < double > v;
        v.resize(L);
        v = MJ[Antigen[i]];
        int index = std::min_element(v.begin(), v.end()) - v.begin();

        min_E = min_E + MJ[Antigen[i]][index];
    }

    return min_E;
}
// Function to return energies < E_max
void high_energy_filter(){

}
// Function to generate the initial amount of sequences
void generate_Bcells(int N, int L, int L_alphabet, vector<bcell> & Bcells)
{
    
    //---------Array with the current Sequence-------------------------------------------
    vector < int > sequence;
    sequence.resize(L);
    
    for(int n =0 ; n<N ; n++){
        
        //Initiating Sequence with random sequence------------------------------------
        for (int k= 0; k<L; k++)
        {
            sequence[k] = randIX(0,L_alphabet-1);
        };
        
        // Create a bcell and add it to the vector
        bcell bcell_i(L, L_alphabet, sequence);
        Bcells[n]  = bcell_i;
    }
}
// Function to generate the initial amount of sequences with energy
void generate_Bcells_with_e(int N, int L, int L_alphabet, vector<bcell> & Bcells, vector< vector<double> > const & MJ, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{
    
    //---------Array with the current Sequence-------------------------------------------
    vector < int > Sequence;
    Sequence.resize(L);
    double e_MS = MS_energy(L, L_alphabet, MJ, Antigen);
    for(int n =0 ; n<N ; n++){
        
        //Initiating Sequence with random sequence------------------------------------
        for (int k= 0; k<L; k++)
        {
            Sequence[k] = randIX(0,L_alphabet-1);
        };
        
        // Create a bcell and add it to the vector
        bcell bcell_i(L, L_alphabet, MJ, Sequence, Antigen, energy_model, r);
        bcell_i.e = bcell_i.e - e_MS - 27.63;
        Bcells[n]  = bcell_i;
    }
}
// Function that selects the antigen-specific naive Bcells from all the sequences
void choose_naive_Bcells(int N, int L, int L_alphabet, vector< vector<double> > const & MJ, vector< int > const & Antigen, vector<bcell> & Bcells, vector<bcell*> & Naive, int & n_naive, string energy_model, gsl_rng *r)
{
    
    vector <int> MS;
    MS.resize(L);
    find_complementary(L, L_alphabet, MJ, Antigen, MS);
    double e_MS = Energy(L, L_alphabet, MJ, MS, Antigen, "MJ", r);
    cout << "e_MS:" << e_MS <<endl;
    double e;
    for(int n = 0 ; n<N ; n++){
        e = Energy(L, L_alphabet, MJ, Bcells[n].seq, Antigen, energy_model, r);
        if(e<e_MS+35){ //Modulate this parameter properly with \epsilon_m
            Bcells[n].e = e;
            Naive.push_back( &Bcells[n]);
            n_naive++;
        }
    }
}
// Function that selects the antigen-specific naive Bcells from all the sequences when energies are already calculated
void choose_naive_Bcells2(int N, int L, int L_alphabet, vector< vector<double> > const & MJ, vector< int > const & Antigen, vector<bcell> & Bcells, vector<bcell*> & Naive, int & n_naive, string energy_model, gsl_rng *r)
{
    
    vector <int> MS;
    MS.resize(L);
    find_complementary(L, L_alphabet, MJ, Antigen, MS);
    double e_MS = Energy(L, L_alphabet, MJ, MS, Antigen, energy_model, r);
    //cout << "e_MS:" << e_MS <<endl;
    double min_e  = Bcells[0].e;
    double e;
    for(int n = 1 ; n<N ; n++){
        e = Bcells[n].e;
        if(e<min_e){ 
            min_e = e;
        }
    }
    for(int n = 0 ; n<N ; n++){
        e = Bcells[n].e;
        if(e<(min_e+15)){
            Naive.push_back( &Bcells[n]);
            n_naive++;
        }
    }
}
// Recognition probability 
double p_R_dT(double dT, double k_on, double rho_A, double gamma, double theta, double e)
{
    double r_0 = (rho_A * k_on) * (pow(gamma, theta)/(pow(gamma, theta) + pow(k_on*exp(e), theta)) ); //This rate has units of s^-1
    double p_R_dT = 1-exp(-r_0*dT*60*60*24); // dT has units of days from the main simulation. Here we correct the units.
    return p_R_dT;
}
// Function to run the set of differential equations
void ODE(int linear, double const alpha, double const beta, double const gamma, long long NT, double dT, int n_naive, vector<bcell*> & Naive, vector<vector < long double > > & Time_series_Bcells, vector < long double > & Time_series_Antigen, vector < double > & N_active_linages)
{
    double p_GC = 0.1;
    double f = 0;
    double N_active_bcells = 0;
    int n_active_linages = 0;
    if (linear==0) {
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*Time_series_Antigen[t-1] - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active);
                if(Naive[n]->active ==0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e));
                    if(f>0.5){
                        Naive[n]->active = 1;
                        n_active_linages++;
                        double r_GC = randX(0,1);
                        // Decide if the activated linage will have plasma or GC fate.
                        if(r_GC<p_GC){
                            Naive[n]->plasma = 0;
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
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*1 - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active);
                if(Naive[n]->active ==0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e));
                    if(f>0.5){
                        Naive[n]->active = 1;
                        n_active_linages++;
                        double r_GC = randX(0,1);
                        // Decide if the activated linage will have plasma or GC fate.
                        if(r_GC<p_GC){
                            Naive[n]->plasma = 0;
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
// Function to run an ensemble of trajectories
void ODE_ensemble(int linear, double const alpha, double const beta, double const gamma, long long NT, double dT, int n_naive, vector<bcell*> & Naive, vector < long double > & Time_series_Antigen, vector < double > & N_active_linages, vector <int> & N_final_active_linages)
{
    double f = 0;
    double N_active_bcells = 0;
    int n_active_linages_t = 0;
    int n_active_linages = 0;
    if (linear==0){
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*Time_series_Antigen[t-1] - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            n_active_linages_t = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Naive[n]->cs = Naive[n]->cs + (beta*Naive[n]->cs*dT*(Naive[n]->active));
                // This function, contrary to the one for the single dynamics, does not use the Time_series_Bcells array. It uses the variable cs of the Bcell Class.
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active); // this uses the time_series arrays
                if(Naive[n]->active == 0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e)); //Modulate this parameter promerly with \epsilon_MS
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
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*2000 - gamma*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            n_active_linages_t = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Naive[n]->cs = Naive[n]->cs + (beta*Naive[n]->cs*dT*(Naive[n]->active));
                // This function, contrary to the one for the single dynamics, does not use the Time_series_Bcells array. It uses the variable cs of the Bcell Class.
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active); // this uses the time_series arrays
                if(Naive[n]->active == 0){
                    f = (Time_series_Antigen[t]/N_A)/((Time_series_Antigen[t]/N_A) + exp(25+Naive[n]->e)); //Modulate this parameter promerly with \epsilon_MS
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
// Function to run an ensemble of trajectories of EF response
void EF_dynamics(int linear, double const alpha, double const beta, double const gamma, double const theta, long long NT, double dT, int n_naive, vector<bcell*> & Naive, long double Antigen, vector < double > & N_active_linages)
{
    double p_GC = 0.1;
    double r;
    //double N_active_bcells = 0; // To use if we want to kill antigen from B cells
    int n_active_linages = 0;
    double k_on = 1e6;
    if (linear==0) {
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            //Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*Time_series_Antigen[t-1] - 0*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            Antigen = exp(alpha*t*dT);
            if(Antigen<1){
                Antigen = 0;
            }
            //N_active_bcells = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active);
                if(Naive[n]->active ==0){
                    r = randX(0,1);
                    //Modulate this parameter promerly with \epsilon_MS
                    if(r < ( p_R_dT(dT, k_on, Antigen/N_A, gamma, theta, Naive[n]->e) )){
                        Naive[n]->active = 1;
                        Naive[n]->activation_time = t*dT;
                        n_active_linages++;
                        double r_GC = randX(0,1);
                        // Decide if the activated linage will have plasma or GC fate.
                        if(r_GC<p_GC){
                            Naive[n]->plasma = 0;
                        }
                    }
                }else{
                    Naive[n]->cs = exp(beta*(t*dT - Naive[n]->cs));
                    //N_active_bcells = N_active_bcells + Time_series_Bcells[n][t];
                }
            }
            N_active_linages[t] = N_active_linages[t] + n_active_linages;
        }
    } else {
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            //Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*1E13  - 0*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            Antigen = 1+alpha*1e13*t*dT;
            if(Antigen<1){
                Antigen = 0;
            }
            //N_active_bcells = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active);
                if(Naive[n]->active ==0){
                    r = ((double) rand() / (RAND_MAX));
                    //Modulate this parameter promerly with \epsilon_MS
                    if(r < ( p_R_dT(dT, k_on, Antigen/N_A, gamma, theta, Naive[n]->e) )){
                        Naive[n]->active = 1;
                        Naive[n]->activation_time = t*dT;
                        n_active_linages++;
                        double r_GC = randX(0,1);
                        // Decide if the activated linage will have plasma or GC fate.
                        if(r_GC<p_GC){
                            Naive[n]->plasma = 0;
                        }
                    }
                }else{
                    Naive[n]->cs = exp(beta*(t*dT - Naive[n]->cs));
                    //N_active_bcells = N_active_bcells + Time_series_Bcells[n][t];
                }
            }
            N_active_linages[t] = N_active_linages[t] + n_active_linages;
        }
    }   
}
// Function to run an ensemble of trajectories of EF response
void EF_dynamics_ensemble(int linear, double const alpha, double const beta, double const gamma, double const theta, long long NT, double dT, int n_naive, vector<bcell*> & Naive, vector < long double > & Time_series_Antigen, vector < double > & N_active_linages, vector <int> & N_final_active_linages)
{
    double r;
    double N_active_bcells = 0;
    int n_active_linages_t = 0;
    int n_active_linages = 0;
    double k_on = 1e6;
    if (linear==0){
        for(int t = 1; t< NT ; t++){ // for loop of time
            //Update the antigen
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*Time_series_Antigen[t-1] - 0*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            n_active_linages_t = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Naive[n]->cs = Naive[n]->cs + (beta*Naive[n]->cs*dT*(Naive[n]->active));
                // This function, contrary to the one for the single dynamics, does not use the Time_series_Bcells array. It uses the variable cs of the Bcell Class.
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active); // this uses the time_series arrays
                if(Naive[n]->active == 0){
                    r = ((double) rand() / (RAND_MAX));
                    //cout << f << " "; //Modulate this parameter promerly with \epsilon_MS
                    if(r < ( (k_on*Time_series_Antigen[t]*pow(gamma, theta)*(60*60*24))/(N_A*alpha * (pow(gamma, theta) + k_on*pow(exp(Naive[n]->e), theta)) ) ) ){
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
            Time_series_Antigen[t] = Time_series_Antigen[t-1] + (alpha*1E8 - 0*Time_series_Antigen[t-1]*N_active_bcells)*dT;
            if(Time_series_Antigen[t]<1){
                Time_series_Antigen[t] = 0;
            }
            N_active_bcells = 0;
            n_active_linages_t = 0;
            //Update Bcells
            for(int n = 0 ; n<n_naive ; n++){
                Naive[n]->cs = Naive[n]->cs + (beta*Naive[n]->cs*dT*(Naive[n]->active));
                // This function, contrary to the one for the single dynamics, does not use the Time_series_Bcells array. It uses the variable cs of the Bcell Class.
                //Time_series_Bcells[n][t] = Time_series_Bcells[n][t-1] + (beta*Time_series_Bcells[n][t-1])*dT*(Naive[n]->active); // this uses the time_series arrays
                if(Naive[n]->active == 0){
                    r = ((double) rand() / (RAND_MAX));
                    if(r < ( (k_on*Time_series_Antigen[t]*pow(gamma, theta)*(60*60*24))/(N_A*alpha * (pow(gamma, theta) + k_on*pow(exp(Naive[n]->e), theta)) ) ) ){
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
void EF_response(int linear, double const alpha, double const beta, double gamma, double const theta, double N_c, double To, double Tf, long long int NT, double dT, int n_naive, vector<bcell*> & Naive)
{
    double k_on = 1e6*24*3600; // M*days^-1
    //double k_off = 1e-8*24*3600 days^-1;
    double k_off;
    gamma = gamma*24; // days^-1
    double alpha2 = alpha;
    int z;
    double r;
    double r_GC;
    double p_GC = 0.05;

    long long int n_time = NT; //number of steps

    // Time array 
    valarray<double> time(n_time);
    time[0] = To;
    for (int t = 1; t < n_time; t++)
    {   
        time[t] = time[t-1] + dT;
    }
    // arrays for probability and cumulative probability for the time of recognition
    valarray<long double> u_on(n_time);
    valarray<long double> prob_recognition(n_time);
    valarray<double> cum_prob_recognition(n_time);
    
    for(int n = 0 ; n<n_naive ; n++){
        //getting k_off from the energy
        k_off = k_on*exp(Naive[n]->e);
        u_on = (exp(alpha2*time)/N_A)*k_on*N_c;
        double p_act = 1/(1+pow((k_off/gamma),theta));
        prob_recognition = (u_on*p_act) * exp(-((u_on/alpha2)*p_act)) * dT;
        partial_sum(begin(prob_recognition), end(prob_recognition), begin(cum_prob_recognition));
        r = randX(0,1);
        z = 1;
        while((cum_prob_recognition[z] < r) & (z < n_time)){
            z++;
        }
        if(z<n_time){
            Naive[n]->activation_time = time[z];
            Naive[n]->active = 1;
        }
        r_GC = randX(0,1);
        // Decide if the activated linage will have plasma or GC fate.
        if(r_GC<p_GC){
            Naive[n]->plasma = 0;
        }
    }

    
    //double N_active_bcells = 0; // To use if we want to kill antigen from B cells
    int n_active_linages = 0;
    
    if (linear==0) {
        
    } else {

    }   
}

