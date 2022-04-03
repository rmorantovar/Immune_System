//
//  EF_dynamics.cpp
//  
//
//  Created by Roberto Moran Tovar on 25.02.22.
//
//Template to run a stochastic simulation of the EF response.

#include "../library/Immuno_functions.hpp"

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
    double alpha;
    double beta;
    double gamma;
    int q;
    int L_seq; //length of the sequence
    int L;
    int L_alphabet (20); //length of the alphabet
    long long int N_bcs; // number of bcells
    int Tf; //number of days for the simulation
    int To;; //initial number of days for the simulation
    long long int NT;
    double dT = 0.5; //time step
    if(ensemble_flag==0){
    	dT = 0.001;
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
	      {"alpha", required_argument, 0, 'a'},
	      {"beta",  required_argument, 0, 'b'},
	      {"gamma",required_argument, 0, 'g'},
	      {"proof_reading",required_argument, 0, 'q'},
	      {"To",    required_argument, 0, 't'},
	      {"Tf",    required_argument, 0, 'T'},
	      {"energy_model",    required_argument, 0, 'E'},
	      {"L",    required_argument, 0, 'L'},
	      {"N_bcs",    required_argument, 0, 'B'},
	      {"Antigen_seq", required_argument, 0, 's'},
	      {"N_ens", required_argument, 0, 'N'},
	      {0, 0, 0, 0}
	    };
	  /* getopt_long stores the option index here. */
	  int option_index = 0;

	  c = getopt_long (argc, argv, "a:b:g:q:t:T:E:L:B:s:N:",
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
	    	alpha = atof(optarg);
	      	break;

	    case 'b':
	    	beta = atof(optarg);
	    	break;

	    case 'g':
	    	gamma = atof(optarg);
	    	break;

	    case 'q':
	    	q = atof(optarg);
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

	    case 'L':
	    	L = atoi(optarg);
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
	A0 = exp(alpha*To);

	//------------Energy Matrix------------------------------------------------------
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
    file.close();
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
    file2.close();
	//------------- Antigen ----------------------------------------------------------------
    vector < int > Antigen;
    Antigen.resize(L_seq);
    aa_to_positions(L_seq, L_alphabet, Alphabet, Antigen, Antigen_aa);
    //variable with antigen size
    long double Antigen_t;
    //---------Generating Bcells ---------------------------------------------------------
    //Array with Bcells
    vector < bcell > Bcells;
    Bcells.resize(N_bcs);
    //---------Activated linages ---------------------------------------------------------
    //Array for time series of the average number of active bcell linages per time
    vector <double> m_bar;
    m_bar.resize(NT);
    //Array for total final number of active bcell linages
    vector <int> N_final_active_linages;
    N_final_active_linages.resize(N_ens);

    //------------------------------------------------------------------------------------
    //-------Files-----
    //Output files
    if(ensemble_flag){ // ENSEMBLE OF TRAJECTORIES
    	string parameters_path = "L-"+std::to_string(L_seq)+"_Nbc-"+ std::to_string(N_bcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model;
    	fs::create_directories(Text_files_path+"Ensemble/"+parameters_path);
    	ofstream fout_energies (Text_files_path+"Ensemble/"+parameters_path+"/energies_ensemble.txt"); // Energies
    	ofstream fout_bcells (Text_files_path+"Ensemble/"+parameters_path+"/bcells_ensemble.txt"); // B cells final clone size
    	ofstream fout_N_final_active (Text_files_path+"Ensemble/"+parameters_path+"/N_final_active.txt"); // B cells final clone siz
    	ofstream fout_m_bar (Text_files_path+"Ensemble/"+parameters_path+"/m_bar.txt"); // time series of the average of the number of activated bcell linages.
    	// ------------ Run ensemble of trajectories ------------
	    cout << "Running ensemble of trajectories ..." << endl;
	    for(int i_ens = 0 ; i_ens<N_ens ; i_ens++){
	    	//-------------------------------------------------------
	    	// printing progress bar
	    	float progress = i_ens/N_ens;
	    	std::cout << "[";
	    	int pos = barWidth * progress;
	    	for (int i = 0; i < barWidth; ++i) {
	        	if (i < pos) std::cout << "=";
	        	else if (i == pos) std::cout << ">";
	        	else std::cout << " ";
	    	}
	    	std::cout << "] " << int(progress * 100.0) << " %\r";
	    	std::cout.flush();			
	        //-------------------------------------------------------
	        //Generate bcells
	        generate_Bcells_with_e(N_bcs, L_seq, L_alphabet, Bcells, MJ, Antigen, energy_model, r);
	        
	        // Choose the antigen-specific bcells
	        vector < bcell* > Naive;
	        int n_naive = 0;
	        choose_naive_Bcells2(N_bcs, L_seq, L_alphabet, MJ, Antigen, Bcells, Naive, n_naive, energy_model, r);

	        //initialize time series arrays
	        //if (linear_flag==0) {
	        //    //Time_series_Antigen[0] = A0;
	        //} else {
	        //    //Time_series_Antigen[0] = 1;
	        //};

	        // Run EF dynamics
	        //EF_dynamics_ensemble(linear_flag, alpha, beta, gamma, q, NT, dT, n_naive, Naive, Time_series_Antigen, m_bar, N_final_active_linages);
	        
	        for (int n = 0 ; n<n_naive ; n++){
	            //print in file the energies and the activation state of the antigen-specific bcells.
	            fout_energies << Naive[n]->e << "\t" << Naive[n]->active << endl;
	            //Print the final clone-size of bcells
	            if(Naive[n]->active==1){
	                fout_bcells << Naive[n]->cs << endl;
	            };
	        };
	    };
	    std::cout << std::endl;
	    for (int t= 0; t<NT; t++)
	    {
	        fout_m_bar << m_bar[t]/N_ens << "\t";
	    };
	    
	    for (int i=0; i<N_ens; i++){
	        fout_N_final_active << N_final_active_linages[i] << "\t";
	    }

	    fout_energies.close();
	    fout_bcells.close();
	    fout_m_bar.close();
	    fout_N_final_active.close();

    }else{ // SINGLE TRAJECTORY
    	string parameters_path = "L-"+std::to_string(L_seq)+"_Nbc-"+ std::to_string(N_bcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model;
    	fs::create_directories(Text_files_path+"Trajectories/"+parameters_path);
    	cout<<">Running simulation of the EF dynamics ..."<< endl;
    	ofstream fout (Text_files_path+"Trajectories/"+parameters_path+"/energies.txt"); // Energies, activation, fate, sequence and activation time
    	//ofstream fout_antigen (Text_files_path+"Trajectories/"+parameters_path+"/antigen.txt"); // Antigen trajectory
    	//ofstream fout_bcells (Text_files_path+"Trajectories/"+parameters_path+"/bcells.txt"); // B cell time series
    	ofstream fout_m_bar (Text_files_path+"Trajectories/"+parameters_path+"/m_bar.txt");

    	//Generate bcells
        generate_Bcells_with_e(N_bcs, L_seq, L_alphabet, Bcells, MJ, Antigen, energy_model, r);
        
        // Choose the antigen-specific bcells
        vector < bcell* > Naive;
        int n_naive = 0;
        choose_naive_Bcells2(N_bcs, L_seq, L_alphabet, MJ, Antigen, Bcells, Naive, n_naive, energy_model, r);
        //Matrix with the time series of the antigen-specific Bcells
	    //vector<vector < long double > > Time_series_Bcells;
	    //Time_series_Bcells.resize(n_naive);
	    //for(int n= 0; n<n_naive; n++)
	    //{
	    //    Time_series_Bcells[n].resize(NT);
	    //    Time_series_Bcells[n][0] = Naive[n]->cs;
	    //}

        //initialize time series arrays
        //if (linear_flag==0) {
        //    //Time_series_Antigen[0] = A0;
        //} else {
            //Time_series_Antigen[0] = 1;
        //};
    	
    	cout << "e_bar: " << mean_energy(L_seq, L_alphabet, MJ, Antigen) << endl;
    	cout << "e_min: " << MS_energy(L_seq, L_alphabet, MJ, Antigen) << endl;

	    // Run EF dynamics
	    EF_dynamics(linear_flag, alpha, beta, gamma, q, NT, dT, n_naive, Naive, Antigen_t, m_bar);
	    
	    for (int n= 0; n<n_naive; n++)
	    {
	        fout << Naive[n]->e << "\t" << Naive[n]->active << "\t" << Naive[n]->plasma << "\t" << Naive[n]->activation_time << "\t";
	        for (int i=0; i<L_seq; i++){
	        	fout  << Alphabet[Naive[n]->seq[i]];
	        }
	        fout << endl;
	    };
	    
	    //Print time series of antigen and bcells
	    //for (int n = 0 ; n<n_naive ; n++){
	    //	for(int t=0 ; t<NT; t++){
            	//fout_bcells << Time_series_Bcells[n][t] << " ";
        //    }
        //    fout_bcells << endl;
        //}
        for(int t=0 ; t<NT; t++){
        	//fout_antigen << Time_series_Antigen[t] << " ";
	        fout_m_bar << m_bar[t] << " ";
        }
        //fout_antigen << endl;
        fout_m_bar << endl;
	    
	      
	    fout.close();
	    //fout_antigen.close();
	    //fout_bcells.close();
	    fout_m_bar.close();
    }
    
    //------------------------------------------------------------------------------------
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;

    return 0;
}