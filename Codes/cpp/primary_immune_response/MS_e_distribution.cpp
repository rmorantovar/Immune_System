//
//  MS_e_distribution.cpp
//  
//  Template to find the distribution of MSs
//  Created by Roberto Moran Tovar on 02.03.21.
//

#include "../lib/Immuno_functions.hpp"

//----------------------------------------------------------------------------------
using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    string Text_files_path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/";
    cout<<">Generating master sequences ...\n"<< endl;
    clock_t t1,t2;
    t1=clock();
    //-----------------------------------------------------------------------------
    //Parameters: (they are fine as they are)
    int L  = atoi(argv[1]); //length of the sequence
    int L_alphabet (20);
    int nT (1); //Number of temperature points
    double T1 (.1) ; double T2 (2);
    int N0;
    sscanf(argv[2], "%d", &N0);
    std::cout << "N=" << N0 << std::endl;
    
    //------------ Energy Matrix ------------------------------------------------------
    //MJ Matrix
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
    //cout << "The Alphabet is :";

    for (int k = 0; k < L_alphabet; k++) {

        file2 >> Alphabet[k];
    
    }

    //------------- Initiating Antigen ------------------------------------------------
    //Array with the antigen
    vector < int > Antigen;
    Antigen.resize(L);
    for (int k= 0; k<L; k++){
        Antigen[k] = randIX(0,L_alphabet-1);
    };

    //---------------Initiating the Master sequence----------------------
    //Array with the current sequence
    vector < int > master_sequence;
    master_sequence.resize(L);

    double e;
    
    //Output file
    ofstream fout (Text_files_path+"MS_energies_L-"+std::to_string(L)+"_N-"+ std::to_string(N0)+".txt");
    
    for(int i = 0; i<N0 ; i++){

        for (int k= 0; k<L; k++){
            Antigen[k] = randIX(0,L_alphabet-1);
        };

        find_complementary(L, L_alphabet, MJ, Antigen, master_sequence);
        e = Energy(L, L_alphabet, MJ, master_sequence, Antigen);
        fout << e << std::endl;
    }
    fout.close();
    
    //------------------------------------------------------------------------------
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
    return 0;
}
