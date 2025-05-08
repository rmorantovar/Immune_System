//Template to find bad binder.
//Input: (all parameters are already set internally!)

#include "../lib/Immuno_functions.hpp"
 
//----------------------------------------------------------------------------------
using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{	
	string Text_files_path = "../../../../../Dropbox/Research/Evolution_Immune_System/Text_files/";
	cout<<">Finding a bad antigen ...\n"<< endl;
	clock_t t1,t2;
    t1=clock();
	//-----------------------------------------------------------------------------
	//Parameters: (they are fine as they are)
    int L  = atoi(argv[1]); //length of the sequence
	int L_alphabet (20);
	int nT (1); //Number of temperature points
	double T1 (.1) ; double T2 (2);
	long long int n0 (0*L), d0 (10*L); //Number of steps: initial prelude, distance between sampling points
	int N0[1] = {1E8};

	
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
	    //cout << Alphabet[k] ;
	
	}
	//cout << "\n";

	//------------- Initiating Antigen ------------------------------------------------
	//Array with the antigen
	vector < int > Antigen;
	Antigen.resize(L);
    for (int k= 0; k<L; k++){
        Antigen[k] = randIX(0,L_alphabet-1);
    };
	vector < int > Antigen_i;
	Antigen_i.resize(L);

	//---------------Initiating the Master sequence----------------------
	//Array with the current sequence
	vector < int > master_sequence;
	master_sequence.resize(L);
	
	double e = -80;
	double e_new;
	for(int i = 0; i<1E5 ; i++){

		for (int k= 0; k<L; k++){
			Antigen_i[k] = randIX(0,L_alphabet-1);
		};

		find_complementary(L, L_alphabet, MJ, Antigen_i, master_sequence);
		e_new = Energy(L, L_alphabet, MJ, master_sequence, Antigen_i);

		if(e_new>e){
			e = e_new;
			Antigen = Antigen_i;
		}
	}
    find_complementary(L, L_alphabet, MJ, Antigen, master_sequence);
    
	cout << "Antigen:";
	for (int k= 0; k<L; k++){
		cout << Alphabet[Antigen[k]];
	};
	cout << "\n";
	cout << "Master Sequence:";
	for (int k= 0; k<L; k++){
		cout << Alphabet[master_sequence[k]];
	};
	cout << "\n";
	cout << "Binding energy:"<< e << "\n";
	
	//------------------------------------------------------------------------------
	cout<< ">Simulation completed…"<< endl;
	t2= clock();
	cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
	return 0;
}
