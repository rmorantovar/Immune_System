//
//  functions.hpp
//  
//  Created by Roberto Moran Tovar on 27.02.21.
//

#ifndef functions_h
#define functions_h
#endif /* functions_h */

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
    bcell(int const & L, int const & L_alphabet, vector< vector<double> > const & E_matrix, vector< int > & seq, vector< int > const & Antigen, string energy_model, gsl_rng *r);
    double e; //energy with respect to the current epitope.
    double cs; //clone size
    bool plasma; // plasma cell?
    //bool GC; // germinal center cell?
    //bool engaged; //engaged with an antigen?
    bool active; // activated clone?
    double activation_time; // time of activation

    // Functions
    double Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & sequence, vector< int > const & Antigen, string energy_model, gsl_rng *r);

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
bcell::bcell(int const & L, int const & L_alphabet, vector< vector<double> > const & E_matrix, vector< int > & seq, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{
    this->seq = seq;
    this->cs = 1.0;
    this->plasma = 1;
    //this->GC = 0;
    //this->engaged = 0;
    this->active = 0;
    this->activation_time = -1;
    this->e = this->Energy(L, L_alphabet, E_matrix, seq, Antigen, energy_model, r);
}
double bcell::Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & sequence, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{

    double E = 0.0;
    if (energy_model == "MJ2") {
        for(int i=0; i<L ; i++){
            E = E + E_matrix[Antigen[i]][sequence[i]];
        }
    }if (energy_model == "TCRen") {
        for(int i=0; i<L ; i++){
            E = E + E_matrix[Antigen[i]][sequence[i]];
        }
    }if (energy_model == "TCRen-CM") {
        for(int i=0; i<L ; i++){
            E = E + E_matrix[Antigen[i]][sequence[i]];
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
// Transpose 2D vector
vector<vector<double> > transpose(vector<vector<double> > const &b)
{
    if (b.size() == 0)
        return b;

    vector<vector<double> > trans_vec(b[0].size(), vector<double>());

    for (int i = 0; i < b.size(); i++)
    {
        for (int j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }
    return trans_vec;
}
//Function to calculate the energy: Implement the Energy Matrix
double Energy(int const & L, int const & L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & sequence, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{

    double E = 0.0;
    if (energy_model == "MJ2") {
        for(int i=0; i<L ; i++){
            E = E + E_matrix[Antigen[i]][sequence[i]];
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
inline double delt( int const & L, int const & L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & sequence, vector< int > const & Antigen, int const & pos, int const & aa)
{
    double deltaE (0.);
    deltaE = E_matrix[Antigen[pos]][aa] - E_matrix[Antigen[pos]][sequence[pos]];
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
void find_complementary(int const & L, int const & L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & sequence, vector<int> & complementary_sequence)
{
    for(int i=0; i<L ; i++){
        vector < double > v;
        v.resize(L);
        v = E_matrix[sequence[i]];
        int index = std::min_element(v.begin(), v.end()) - v.begin();
        complementary_sequence[i] = index;
    }
};
//Function to calculate the mean energy of the associated PWM
double mean_energy(int L, int L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & Antigen)
{

    double mean_E = 0.0;

    for(int i=0; i<L ; i++){
        for (int j = 0; j < L_alphabet; ++j)
        {
            mean_E = mean_E + E_matrix[Antigen[i]][j];
        }
    }
    mean_E=mean_E/(L_alphabet);

    return mean_E;
}
//Function to calculate the minimum energy of the associated PWM
double MS_energy(int L, int L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & Antigen)
{

    double min_E = 0.0;
    for(int i=0; i<L ; i++){

        vector < double > v;
        v.resize(L);
        v = E_matrix[Antigen[i]];
        int index = std::min_element(v.begin(), v.end()) - v.begin();

        min_E = min_E + E_matrix[Antigen[i]][index];
    }

    return min_E;
}
// Function to generate the initial amount of sequences with energy
void generate_Bcells_with_e(int N, int L, int L_alphabet, vector<bcell> & Bcells, vector< vector<double> > const & E_matrix, vector< int > const & Antigen, string energy_model, gsl_rng *r)
{
    //---------Array with the current Sequence-------------------------------------------
    vector < int > Sequence;
    Sequence.resize(L);
    double e_MS = MS_energy(L, L_alphabet, E_matrix, Antigen);
    for(int n =0 ; n<N ; n++){
        
        //Initiating Sequence with random sequence------------------------------------
        for (int k= 0; k<L; k++)
        {
            Sequence[k] = randIX(0,L_alphabet-1);
        };
        
        // Create a bcell and add it to the vector
        bcell bcell_i(L, L_alphabet, E_matrix, Sequence, Antigen, energy_model, r);
        bcell_i.e = bcell_i.e - e_MS - 25.0;
        Bcells[n]  = bcell_i;
    }
}
// Function that selects the antigen-specific naive Bcells from all the sequences when energies are already calculated
void choose_naive_Bcells2(int N, int L, int L_alphabet, vector< vector<double> > const & E_matrix, vector< int > const & Antigen, vector<bcell> & Bcells, vector<bcell*> & Naive, int & n_naive, string energy_model, gsl_rng *r)
{
    
    //vector <int> MS;
    //MS.resize(L);
    //find_complementary(L, L_alphabet, E_matrix, Antigen, MS);
    //double e_MS = Energy(L, L_alphabet, E_matrix, MS, Antigen, energy_model, r);
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
        if(e<(min_e+10)){ //Use 10. Use 7 for elite
            Naive.push_back( &Bcells[n]);
            n_naive++;
        }
    }
}
void EF_response(int linear, double const alpha, double const beta, double gamma, double const theta, long double const N_c, double To, double Tf, long long int NT, double dT, int n_naive, vector<bcell*> & Naive)
{
    long double k_on = 1e6*24*3600; // M*days^-1
    //double k_off = 1e-8*24*3600 days^-1;
    double k_off;
    gamma = gamma*24; // days^-1
    long double alpha2 = alpha;
    int z;
    double r;
    double r_GC;
    double p_GC = 0.05;

    long long int n_time = NT; //number of steps

    // Time array 
    valarray<long double> time(n_time);
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
        double p_act = 1/pow((1+(k_off/gamma)),theta);
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

            r_GC = randX(0,1);
            // Decide if the activated linage will have plasma or GC fate.
            if(r_GC<p_GC){
                Naive[n]->plasma = 0;
            }
        }
        
    }

    //double N_active_bcells = 0; // To use if we want to kill antigen from B cells
    int n_active_linages = 0;
    
    if (linear==0) {
        
    } else {

    }   
}

