//
//  EF_dynamics.cpp
//  
//
//  Created by Roberto Moran Tovar on 25.02.22.
//
//Template to run a stochastic simulation of the EF response.

#include "../library/functions.hpp"
#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem>

using namespace std;

/* Flag set by ‘--verbose’. */
static int linear_flag=0;
static int ensemble_flag=0;

//----------------------------------------------------------------------------------
//using namespace std;
namespace fs = std::filesystem;

// define a struct that represents a row of data
struct DataRow {
	double energy;
	int active;
	int plasma;
	double act_time;
	int ensemble_id;
    //std::string sequence;
    //std::vector<int> sequence;
};
//----------------------------------------------------------------------------------
int main(int argc, char* argv[]) //argv 
{
	string Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/";
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(r, time(NULL));
    clock_t t1,t2;
    t1=clock();
    int barWidth = 70;
	//-----------------------------------------------------------------------------
    //Parameters:
    double lambda_A;
    double lambda_B;
    double k_pr;
    double theta;
    int L_seq; //length of the sequence
    long double N_c;
    int L_alphabet (20); //length of the alphabet
    long long int N_bcs; // number of bcells
    int Tf; //number of days for the simulation
    int To;; //initial number of days for the simulation
    long long int NT;
    double dT = 0.05; //time step
    if(ensemble_flag==0){
    	dT = 0.05;
    }
    long long int N_ens = 1;
    long long A0;
    std::string energy_model;
    std::string Antigen_aa;
    //-----------------------------------------------------------------------------
    //Read flags and inputs
    int c;
	while (1)
	{
	  static struct option long_options[] =
	    {
	      /* These options set a flag. */
	      {"linear", no_argument,  &linear_flag, 1},
	      {"ensemble", no_argument,  &ensemble_flag, 1},
	      /* These options don’t set a flag.
	         We distinguish them by their indices. */
	      {"lambda_A", required_argument, 0, 'a'},
	      {"lambda_B",  required_argument, 0, 'b'},
	      {"k_pr",required_argument, 0, 'k'},
	      {"proof_reading",required_argument, 0, 'q'},
	      {"To",    required_argument, 0, 't'},
	      {"Tf",    required_argument, 0, 'T'},
	      {"energy_model",    required_argument, 0, 'E'},
	      {"N_c",    required_argument, 0, 'C'},
	      {"N_bcs",    required_argument, 0, 'B'},
	      {"Antigen_seq", required_argument, 0, 's'},
	      {"N_ens", required_argument, 0, 'N'},
	      {0, 0, 0, 0}
	    };
	  /* getopt_long stores the option index here. */
	  int option_index = 0;

	  c = getopt_long (argc, argv, "a:b:k:q:t:T:E:C:B:s:N:",
	                   long_options, &option_index);
	  /* Detect the end of the options. */
	  if (c == -1)
	    break;

	  switch (c)
	    {
	    case 0:
	      /* If this option set a flag, do nothing else now. */
	      if (long_options[option_index].flag != 0)
	      	break;
	      printf ("option %s", long_options[option_index].name);
	      if (optarg)
	        printf (" with arg %s", optarg);
	      	printf ("\n");
	      break;

	    case 'a':
	    	lambda_A = atof(optarg);
	      	break;

	    case 'b':
	    	lambda_B = atof(optarg);
	    	break;

	    case 'k':
	    	k_pr = atof(optarg);
	    	break;

	    case 'q':
	    	theta = atof(optarg);
	    	break;

	    case 't':
	    	To = atof(optarg);
			break;

	    case 'T':
	    	Tf = atof(optarg);
			break;
	    
	    case 'E':
	    	energy_model = optarg;
			break;

	    case 'C':
	    	N_c = stold(optarg);
			break;

	    case 'B':
	    	N_bcs = atoi(optarg);
			break;

		case 's':
	    	Antigen_aa = optarg;
			break;

		case 'N':
	    	N_ens = atoi(optarg);
			//printf ("option -N with value `%s'\n", optarg);
			break;

	    case '?':
			/* getopt_long already printed an error message. */
			break;

	    default:
			abort ();
	    }
	}
	//Report if antigen growth is linear
	if (linear_flag)
	puts ("Antigen grows linearly");
	
	//Report if an ensemble is performed
	if (ensemble_flag)
	puts ("Performing an ensemble of trajectories");

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
	  printf ("non-option ARGV-elements: ");
	  while (optind < argc)
	    printf ("%s ", argv[optind++]);
	  putchar ('\n');
	}
	NT = (Tf-To)/dT; //number of steps
	L_seq = Antigen_aa.length();
	A0 = exp(lambda_A*To);

	//------------Energy Matrix------------------------------------------------------
    vector < vector < double > > E_matrix_t;
    E_matrix_t.resize(L_alphabet);
    for (int k= 0; k<L_alphabet; k++)
    {
        (E_matrix_t[k]).resize(L_alphabet);
    };
    
    ifstream file("../Input_files/"+energy_model+".txt");
    for (unsigned int i = 0; i < L_alphabet; i++) {
        for (unsigned int j = 0; j < L_alphabet; j++) {
            file >> E_matrix_t[i][j];
        }
    }
    file.close();
    vector < vector < double > > E_matrix = transpose(E_matrix_t);
    //------------ Alphabet ----------------------------------------------------------
    //Array with the Alphabet
    vector < string > Alphabet;
    Alphabet.resize(L_alphabet);
    ifstream file2("../Input_files/Alphabet_"+energy_model+".txt");
    cout << "The Alphabet is :";
    for (int k = 0; k < L_alphabet; k++) {

        file2 >> Alphabet[k];
        cout << Alphabet[k] ;
    
    }
    cout << "\n";
    
    for(int k = 0 ; k<L_alphabet; k++){
    	cout << E_matrix[0][k] << endl; 
    }




}

