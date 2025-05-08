//
//  generate_random_bcells.cpp
//  
//
//  Created by Roberto Moran Tovar on 15.03.21.
//
//

#include "../lib/Immuno_functions.hpp"

#include <stdio.h>

//----------------------------------------------------------------------------------
using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[]) //argv has 1:L , 2:N , 3:N_ensemble
{
    string Text_files_path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Random/";
    cout<<">Running simulation of the Bcells-Antigen dynamics ..."<< endl;
    clock_t t1,t2;
    t1=clock();
    //-----------------------------------------------------------------------------
    //Parameters:
   
    int L  = atoi(argv[1]); //length of the sequence
    int L_alphabet (20); //length of the alphabet
    long long int N = atoi(argv[2]); // number of bcells
    long long int N_ensemble = atoi(argv[3]);

    //------------Energy Matrix------------------------------------------------------
    vector < vector < double > > MJ;
    MJ.resize(L_alphabet);
    for (int k= 0; k<L_alphabet; k++)
    {
        (MJ[k]).resize(L_alphabet);
    };

    ifstream file("../Input_files/MJ2.txt");

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
    for (unsigned int i = 0; i < L_alphabet; i++) {
        for (unsigned int j = 0; j < L_alphabet; j++) {
            file >> MJ[i][j];
        }
    }
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
    
    //Output files
    ofstream fout (Text_files_path+"energies_random_L-"+std::to_string(L)+"_N-"+ std::to_string(N)+"_Antigen-"+Antigen_aa+".txt"); // Energies
    
    double e;
    // ------------ Run ensemble of trajectories ------------
    cout << "Generating random B cells ..." << endl;
    for(int i_ensemble = 0 ; i_ensemble<N_ensemble ; i_ensemble++){
        
        //Generate bcells
        generate_Bcells(N, L, L_alphabet, Bcells);
        
        for (int n = 0 ; n<N ; n++){
            e = Energy(L, L_alphabet, MJ, Bcells[n].seq, Antigen);
            //print in file the energies
            fout << e << endl;
            //Print the final clone-size of bcells
        }
    }
    
    fout.close();
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
    return 0;
}

