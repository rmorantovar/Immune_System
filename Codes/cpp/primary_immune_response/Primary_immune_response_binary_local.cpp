//
//  Primary_immune_response_binary_local.cpp
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
    file2.close();
	//------------- Antigen ----------------------------------------------------------------
    vector < int > Antigen;
    Antigen.resize(L_seq);
    aa_to_positions(L_seq, L_alphabet, Alphabet, Antigen, Antigen_aa);
    //variable with antigen size
    long double Antigen_t; // To be removed
    //---------Generating Bcells ---------------------------------------------------------
    //Array with Bcells
    vector < bcell > Bcells;
    Bcells.resize(N_bcs);
    //---------Activated linages ---------------------------------------------------------
    //Array for time series of the average number of active bcell linages per time
    //vector <double> m_bar; // To be removed 
    //m_bar.resize(NT);
    //Array for total final number of active bcell linages
    //vector <int> N_final_active_linages;
    //N_final_active_linages.resize(N_ens);

    //------------------------------------------------------------------------------------
    //-------Files-----
    //Output files
    if(ensemble_flag){ // ENSEMBLE OF TRAJECTORIES
    	// create a vector to store the data rows
    	std::vector<DataRow> data;
    	// ------------ Run ensemble of trajectories ------------
	    cout << "Running ensemble of trajectories ..." << endl;
	    for(int i_ens = 0 ; i_ens<N_ens ; i_ens++){
	    	// ----------------------------PRINTING PROGRESS BAR-----------------------------
	    	float progress = i_ens/(double)N_ens;
	    	std::cout << "[";
	    	int pos = barWidth * progress;
	    	for (int i = 0; i < barWidth; ++i) {
	        	if (i < pos) std::cout << "=";
	        	else if (i == pos) std::cout << ">";
	        	else std::cout << " ";
	    	}
	    	std::cout << "] " << int(progress * 100.0) << " %\r";
	    	std::cout.flush();			
	        //--------------------------------------------------------------------------------
	        //-----------------Generate bcells-----------------
	        generate_Bcells_with_e(N_bcs, L_seq, L_alphabet, Bcells, E_matrix, Antigen, energy_model, r);
	        
	        //-----------------Choose the antigen-specific bcells-----------------
	        vector < bcell* > Naive;
	        int n_naive = 0;
	        choose_naive_Bcells2(N_bcs, L_seq, L_alphabet, E_matrix, Antigen, Bcells, Naive, n_naive, energy_model, r);

	        // Run EF dynamics
	        EF_response(linear_flag, lambda_A, lambda_B, k_pr, theta, N_c, To, Tf, NT, dT, n_naive, Naive);

	        for (int n= 0; n<n_naive; n++)
		    {
		    	DataRow row;
		        row.energy = Naive[n]->e;
				row.active = Naive[n]->active;
				row.plasma = Naive[n]->plasma;
				row.act_time = Naive[n]->activation_time;
				row.ensemble_id = i_ens;
			    //row.sequence = Naive[n]->seq;
			    data.push_back(row);
		    }; 
	    };

	    string parameters_path = "L-"+std::to_string(L_seq)+"_Nbc-"+ std::to_string(N_bcs)+"_Antigen-"+Antigen_aa+"_lambda_A-"+std::to_string(lambda_A)+"_lambda_B-"+std::to_string(lambda_B)+"_k_pr-"+std::to_string(k_pr)+"_theta-"+std::to_string(theta)+"_Nc-"+std::to_string(N_c)+"_Linear-"+std::to_string(linear_flag)+"_N_ens-"+ std::to_string(N_ens)+"_"+energy_model;
    	fs::create_directories(Text_files_path+"Ensemble/"+parameters_path);
    	std::ofstream fout(Text_files_path+"Ensemble/"+parameters_path+"/energies_ensemble.bin", std::ios::binary);
	    for (const auto& row : data) {
	    	// write the float and int to the file
	    	fout.write(reinterpret_cast<const char*>(&row.energy), sizeof(double));
	    	//fout.write((char*)&row.energy, sizeof(double));
	    	fout.write(reinterpret_cast<const char*>(&row.active), sizeof(int));
	    	//fout.write((char*)&row.active, sizeof(int));
	    	fout.write(reinterpret_cast<const char*>(&row.plasma), sizeof(int));
	    	//fout.write((char*)&row.plasma, sizeof(int));
	    	fout.write(reinterpret_cast<const char*>(&row.act_time), sizeof(double));
	    	//fout.write((char*)&row.act_time, sizeof(double));
	    	fout.write(reinterpret_cast<const char*>(&row.ensemble_id), sizeof(int));
    	}
    	fout.close();

	    std::cout << std::endl;


    }else{ // SINGLE TRAJECTORY
    	string parameters_path = "L-"+std::to_string(L_seq)+"_Nbc-"+ std::to_string(N_bcs)+"_Antigen-"+Antigen_aa+"_lambda_A-"+std::to_string(lambda_A)+"_lambda_B-"+std::to_string(lambda_B)+"_k_pr-"+std::to_string(k_pr)+"_theta-"+std::to_string(theta)+"_Nc-"+std::to_string(N_c)+"_Linear-"+std::to_string(linear_flag)+"_N_ens-"+ std::to_string(N_ens)+"_"+energy_model;
    	fs::create_directories(Text_files_path+"Trajectories/"+parameters_path);
    	cout<<">Running simulation of the EF dynamics ..."<< endl;
    	ofstream fout (Text_files_path+"Trajectories/"+parameters_path+"/energies.txt"); // Energies, activation, fate, sequence and activation time

    	//Generate bcells
        generate_Bcells_with_e(N_bcs, L_seq, L_alphabet, Bcells, E_matrix, Antigen, energy_model, r);
        
        // Choose the antigen-specific bcells
        vector < bcell* > Naive;
        int n_naive = 0;
        choose_naive_Bcells2(N_bcs, L_seq, L_alphabet, E_matrix, Antigen, Bcells, Naive, n_naive, energy_model, r);
        
    	cout << "e_bar: " << mean_energy(L_seq, L_alphabet, E_matrix, Antigen) << endl;
    	cout << "e_min: " << MS_energy(L_seq, L_alphabet, E_matrix, Antigen) << endl;

	    // Run EF dynamics
	    EF_response(linear_flag, lambda_A, lambda_B, k_pr, theta, N_c, To, Tf, NT, dT, n_naive, Naive);
	    
	    for (int n= 0; n<n_naive; n++)
	    {
	        fout << Naive[n]->e << "\t" << Naive[n]->active << "\t" << Naive[n]->plasma << "\t" << Naive[n]->activation_time << "\t";
	        for (int i=0; i<L_seq; i++){
	        	fout  << Alphabet[Naive[n]->seq[i]];
	        }
	        fout << endl;
	    };
	    
	    fout.close();

    }
    
    //------------------------------------------------------------------------------------
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;

    return 0;
}