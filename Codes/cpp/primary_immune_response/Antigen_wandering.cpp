//
//  Antigen_wandering.cpp
//  
//
//  Created by Roberto Moran Tovar on 04.05.21.
//
//Template to run a random walk in the antigen space

#include "../lib/Immuno_functions.hpp"

#include <stdio.h>


int main(int argc, char* argv[]) //argv has 1:L ; 2:N_ensemble ;
{
    string Text_files_path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/Antigen_wandering/";
    cout<<">Running simulation of Antigen random walk ..."<< endl;
    clock_t t1,t2;
    t1=clock();
    //-----------------------------------------------------------------------------
    //Parameters:
    
    int L  = atoi(argv[1]); //length of the sequence
    int L_alphabet (20); //length of the alphabet
    long long int N_ensemble = atoi(argv[2]); // Ensemble size
    long long int T = 1000; //Random walk time
    
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
    
    //---------Array with the Memory Bcell and MS--------------------------------------------
    vector < int > Memory;
    Memory.resize(L);
    vector < int > Memory1;
    Memory1.resize(L);
    vector < int > Memory2;
    Memory2.resize(L);
    vector < int > Memory3;
    Memory3.resize(L);
    vector < int > Memory4;
    Memory4.resize(L);
    vector < int > Memory5;
    Memory5.resize(L);
    
    vector < int > MS;
    MS.resize(L);
    
    double E_memory = 0;
    double E_memory1 = 0;
    double E_memory2 = 0;
    double E_memory3 = 0;
    double E_memory4 = 0;
    double E_memory5 = 0;
    double E_MS = 0;
    
    ofstream fout (Text_files_path+"energies_Mem0-MS_L-"+std::to_string(L)+"_N-"+ std::to_string(N_ensemble)+"_antigen-"+Antigen_aa+".txt");
    ofstream fout1 (Text_files_path+"energies_Mem1-MS_L-"+std::to_string(L)+"_N-"+ std::to_string(N_ensemble)+"_antigen-"+Antigen_aa+".txt");
    ofstream fout2 (Text_files_path+"energies_Mem2-MS_L-"+std::to_string(L)+"_N-"+ std::to_string(N_ensemble)+"_antigen-"+Antigen_aa+".txt");
    ofstream fout3 (Text_files_path+"energies_Mem3-MS_L-"+std::to_string(L)+"_N-"+ std::to_string(N_ensemble)+"_antigen-"+Antigen_aa+".txt");
    ofstream fout4 (Text_files_path+"energies_Mem4-MS_L-"+std::to_string(L)+"_N-"+ std::to_string(N_ensemble)+"_antigen-"+Antigen_aa+".txt");
    ofstream fout5 (Text_files_path+"energies_Mem5-MS_L-"+std::to_string(L)+"_N-"+ std::to_string(N_ensemble)+"_antigen-"+Antigen_aa+".txt");
    ofstream fout_MS (Text_files_path+"energies_MS_L-"+std::to_string(L)+"_N-"+ std::to_string(N_ensemble)+"_antigen-"+Antigen_aa+".txt");
    
    //------- Start the random walk ---------
    cout<< ">Running the Antigen random walk…"<< endl;
    
    for (int n=0; n<N_ensemble; n++) {
        //------- Initializing arrays -------
        //for (int k= 0; k<L; k++)
        //{
        //    Antigen[k] = randIX(0,L_alphabet-1);
        //};
        
        aa_to_positions(L, L_alphabet, Alphabet, Antigen, Antigen_aa);
        
        find_complementary(L, L_alphabet, MJ, Antigen, Memory);
        find_complementary(L, L_alphabet, MJ, Antigen, Memory1);
        find_complementary(L, L_alphabet, MJ, Antigen, Memory2);
        find_complementary(L, L_alphabet, MJ, Antigen, Memory3);
        find_complementary(L, L_alphabet, MJ, Antigen, Memory4);
        find_complementary(L, L_alphabet, MJ, Antigen, Memory5);
        find_complementary(L, L_alphabet, MJ, Antigen, MS);
        
        // create imperfectly maturated memory
        mutate_sequence(L, L_alphabet, Memory1);
        
        mutate_sequence(L, L_alphabet, Memory2);
        mutate_sequence(L, L_alphabet, Memory2);
        
        for (int i = 0; i<10; i++) {
            mutate_sequence(L, L_alphabet, Memory3);
        }
        for (int i = 0; i<20; i++) {
            mutate_sequence(L, L_alphabet, Memory4);
        }
        for (int i = 0; i<100; i++) {
            mutate_sequence(L, L_alphabet, Memory5);
        }
        //---------------------------------------
        
        E_memory = Energy(L, L_alphabet, MJ, Memory, Antigen);
        E_memory1 = Energy(L, L_alphabet, MJ, Memory1, Antigen);
        E_memory2 = Energy(L, L_alphabet, MJ, Memory2, Antigen);
        E_memory3 = Energy(L, L_alphabet, MJ, Memory3, Antigen);
        E_memory4 = Energy(L, L_alphabet, MJ, Memory4, Antigen);
        E_memory5 = Energy(L, L_alphabet, MJ, Memory5, Antigen);
        E_MS = Energy(L, L_alphabet, MJ, MS, Antigen);
        
        
        fout << E_memory - E_MS << "\t";
        fout1 << E_memory1 - E_MS << "\t";
        fout2 << E_memory2 - E_MS << "\t";
        fout3 << E_memory3 - E_MS << "\t";
        fout4 << E_memory4 - E_MS << "\t";
        fout5 << E_memory5 - E_MS << "\t";
        fout_MS << E_MS << "\t";
        
        for (int t = 0; t<T; t++) {
            mutate_sequence(L, L_alphabet, Antigen);
            find_complementary(L, L_alphabet, MJ, Antigen, MS);
            E_memory = Energy(L, L_alphabet, MJ, Memory, Antigen);
            E_memory1 = Energy(L, L_alphabet, MJ, Memory1, Antigen);
            E_memory2 = Energy(L, L_alphabet, MJ, Memory2, Antigen);
            E_memory3 = Energy(L, L_alphabet, MJ, Memory3, Antigen);
            E_memory4 = Energy(L, L_alphabet, MJ, Memory4, Antigen);
            E_memory5 = Energy(L, L_alphabet, MJ, Memory5, Antigen);
            E_MS = Energy(L, L_alphabet, MJ, MS, Antigen);
            fout << E_memory - E_MS << "\t";
            fout1 << E_memory1 - E_MS << "\t";
            fout2 << E_memory2 - E_MS << "\t";
            fout3 << E_memory3 - E_MS << "\t";
            fout4 << E_memory4 - E_MS << "\t";
            fout5 << E_memory5 - E_MS << "\t";
            fout_MS << E_MS << "\t";
        }
        
        fout << endl;
        fout1 << endl;
        fout2 << endl;
        fout3 << endl;
        fout4 << endl;
        fout5 << endl;
        fout_MS << endl;
    }
    cout<< ">Random walk completed…"<< endl;

    fout.close();
    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
    fout5.close();
    
    fout_MS.close();
    
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;

    return 0;
}
