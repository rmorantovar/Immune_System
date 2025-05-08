//
//  antigenic_0.cpp
//  
//
//  Created by Roberto Moran Tovar on 29.03.21.
//

#include "../lib/Immuno_functions.hpp"

#include <stdio.h>

//----------------------------------------------------------------------------------
using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[]) //argv has 1:L
{
    string Text_files_path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Antigenic/";
    cout<<">Running simulation of the Bcells-Antigen dynamics ..."<< endl;
    clock_t t1,t2;
    t1=clock();
    //-----------------------------------------------------------------------------
    //Parameters:
   
    int L  = atoi(argv[1]); //length of the sequence
    int L_alphabet (20); //length of the alphabet
    int N = 2E6;

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
        
    //------------- Antigens -------------------------------------------------------------
    //Arrays with the antigens
    //string Antigen_aa;
    //cout << "Insert the Aminoacid sequence of the antigen:\n";
    //getline(cin, Antigen_aa);
    
    vector < int > Antigen1;
    Antigen1.resize(L);
    vector < int > Antigen2;
    Antigen2.resize(L);
    
    //---------Array with the Memory Bcell--------------------------------------------
    vector < int > Memory1;
    Memory1.resize(L);
    vector < int > Memory2;
    Memory2.resize(L);
    
    //aa_to_positions(L, L_alphabet, Alphabet, Antigen1, Antigen_aa);
    //aa_to_positions(L, L_alphabet, Alphabet, Antigen2, Antigen_aa);
    
    double e11;
    double e12;
    double e22;
    double e21;
    for (int M = 1; M<6; M++) {
        
        cout << "Running Hamming distance " << M << endl;
        //Output files
        ofstream fout (Text_files_path+"energies_hamming_distance_"+to_string(M)+"_.txt");
        
        for (int n = 0; n<N ; n++) {
            
            //Initializing Antigen 1 and 2
            for (int k= 0; k<L; k++){
                Antigen1[k] = randIX(0,L_alphabet-1);
                Antigen2[k] = Antigen1[k];
            };
            
            for (int m = 0; m<M; m++) {
                mutate_sequence(L, L_alphabet, Antigen2);
            }
            
            find_complementary(L, L_alphabet, MJ, Antigen1, Memory1);
            find_complementary(L, L_alphabet, MJ, Antigen2, Memory2);
            
            e11 = Energy(L, L_alphabet, MJ, Memory1, Antigen1);
            e12 = Energy(L, L_alphabet, MJ, Memory2, Antigen1);
            e22 = Energy(L, L_alphabet, MJ, Memory2, Antigen2);
            e21 = Energy(L, L_alphabet, MJ, Memory1, Antigen2);
            
            fout << e11 << "\t" << e12 << "\t" << e22 << "\t" << e21 << endl;
        }
    }
    
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
    return 0;
}

