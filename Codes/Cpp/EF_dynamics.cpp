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

using namespace std;

/* Flag set by ‘--verbose’. */
static int linear_flag=0;
static int ensemble_flag=0;

//----------------------------------------------------------------------------------
//using namespace std;

//----------------------------------------------------------------------------------
int main(int argc, char* argv[]) //argv 
{
	string Text_files_path = "/Users/robertomorantovar/Dropbox/Research/Evolution_Immune_System/Text_files/Dynamics/";
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(r, time(NULL));
    clock_t t1,t2;
    t1=clock();
	//-----------------------------------------------------------------------------
    //Parameters:
    double alpha;
    double beta;
    double gamma;
    int L_seq; //length of the sequence
    int L;
    int L_alphabet (20); //length of the alphabet
    long long int N_bbcs; // number of bcells
    int Tf; //number of days for the simulation
    int To;; //initial number of days for the simulation
    long long int NT;
    double dT = 0.5; //time step
    if(ensemble_flag==0){
    	dT = 0.005;
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
	      {"gamma",required_argument, 0, 'c'},
	      {"To",    required_argument, 0, 't'},
	      {"Tf",    required_argument, 0, 'T'},
	      {"energy_model",    required_argument, 0, 'E'},
	      {"L",    required_argument, 0, 'L'},
	      {"N_bbcs",    required_argument, 0, 'B'},
	      {"Antigen_seq", required_argument, 0, 's'},
	      {"N_ens", required_argument, 0, 'N'},
	      {0, 0, 0, 0}
	    };
	  /* getopt_long stores the option index here. */
	  int option_index = 0;

	  c = getopt_long (argc, argv, "a:b:c:t:T:E:L:B:s:N:",
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

	    case 'c':
	    	gamma = atof(optarg);
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
	    	N_bbcs = atoi(optarg);
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

    //---------Generating Bcells ---------------------------------------------------------
    //Array with Bcells
    vector < bcell > Bcells;
    Bcells.resize(N_bbcs);

    //Array with time series of the antigen
    vector < long double > Time_series_Antigen;
    Time_series_Antigen.resize(NT);
    
    //---------Activated linages ---------------------------------------------------------
    //Array for time series of the number of active bcell linages per time
    vector <double> N_active_linages;
    N_active_linages.resize(NT);
    //Array for total final number of active bcell linages
    vector <int> N_final_active_linages;
    //N_final_active_linages.resize(N_ens);

    //------------------------------------------------------------------------------------
    //-------Files-----
    //Output files
    if(ensemble_flag){ // ENSEMBLE OF TRAJECTORIES
    	ofstream fout (Text_files_path+"Ensemble/energies_ensemble_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt"); // Energies
    	ofstream fout_bcells (Text_files_path+"Ensemble/bcells_ensemble_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt"); // B cells final clone size
    	ofstream fout_N_final_active (Text_files_path+"Ensemble/N_final_active_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt"); // B cells final clone siz
    	//print in file the time series of the average of the number of activated bcell linages.
    	ofstream fout_N_active_linages (Text_files_path+"Ensemble/N_active_linages_ensemble_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt");
    	// ------------ Run ensemble of trajectories ------------
	    cout << "Running ensemble of trajectories ..." << endl;
	    for(int i_ens = 0 ; i_ens<N_ens ; i_ens++){
	        
	        //Generate bcells
	        generate_Bcells_with_e(N_bbcs, L_seq, L_alphabet, Bcells, MJ, Antigen, energy_model, r);
	        
	        // Choose the antigen-specific bcells
	        vector < bcell* > Naive;
	        int n_naive = 0;
	        choose_naive_Bcells2(N_bbcs, L_seq, L_alphabet, MJ, Antigen, Bcells, Naive, n_naive, energy_model, r);

	        //initialize time series arrays
	        if (linear_flag==0) {
	            Time_series_Antigen[0] = A0;
	        } else {
	            Time_series_Antigen[0] = 1;
	        };

	        // Run EF dynamics
	        EF_dynamics_ensemble(linear_flag, alpha, beta, gamma, NT, dT, n_naive, Naive, Time_series_Antigen, N_active_linages, N_final_active_linages);
	        
	        for (int n = 0 ; n<n_naive ; n++){
	            //print in file the energies and the activation state of the antigen-specific bcells.
	            fout << Naive[n]->e << "\t" << Naive[n]->active << endl;
	            //Print the final clone-size of bcells
	            if(Naive[n]->active==1){
	                fout_bcells << Naive[n]->cs << endl;
	            };
	        };
	    };

	    for (int t= 0; t<NT; t++)
	    {
	        fout_N_active_linages << N_active_linages[t]/N_ens << "\t";
	    };
	    
	    for (int i=0; i<N_ens; i++){
	        fout_N_final_active << N_final_active_linages[i] << "\t";
	    }

	    fout.close();
	    fout_bcells.close();
	    fout_N_active_linages.close();
	    fout_N_final_active.close();

    }else{ // SINGLE TRAJECTORY
    	cout<<">Running simulation of the EF dynamics ..."<< endl;
    	ofstream fout (Text_files_path+"Trajectories/energies_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt");
    	ofstream fout_antigen (Text_files_path+"Trajectories/antigen_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt");
    	ofstream fout_bcells (Text_files_path+"Trajectories/bcells_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt");
    	ofstream fout_N_active_linages (Text_files_path+"Trajectories/N_active_linages_L-"+std::to_string(L_seq)+"_N-"+ std::to_string(N_bbcs)+"_Antigen-"+Antigen_aa+"_alpha-"+std::to_string(alpha)+"_beta-"+std::to_string(beta)+"_gamma-"+std::to_string(gamma)+"_Linear-"+std::to_string(linear_flag)+"_"+energy_model+".txt");

    	//Generate bcells
        generate_Bcells_with_e(N_bbcs, L_seq, L_alphabet, Bcells, MJ, Antigen, energy_model, r);
        
        // Choose the antigen-specific bcells
        vector < bcell* > Naive;
        int n_naive = 0;
        choose_naive_Bcells2(N_bbcs, L_seq, L_alphabet, MJ, Antigen, Bcells, Naive, n_naive, energy_model, r);
        //Matrix with the time series of the antigen-specific Bcells
	    vector<vector < long double > > Time_series_Bcells;
	    Time_series_Bcells.resize(n_naive);
	    for(int n= 0; n<n_naive; n++)
	    {
	        Time_series_Bcells[n].resize(NT);
	        Time_series_Bcells[n][0] = Naive[n]->cs;
	    }
        //initialize time series arrays
        if (linear_flag==0) {
            Time_series_Antigen[0] = A0;
        } else {
            Time_series_Antigen[0] = 1;
        };
    	
    	cout << "e_bar:" << mean_energy(L, L_alphabet, MJ, Antigen) << endl;
	    // Run EF dynamics
	    EF_dynamics(linear_flag, alpha, beta, gamma, NT, dT, n_naive, Naive, Time_series_Bcells, Time_series_Antigen, N_active_linages);
	    
	    for (int n= 0; n<n_naive; n++)
	    {
	        fout << Naive[n]->e << "\t" << Naive[n]->active << "\t" << Naive[n]->GC << endl;
	    };
	    
	    //Print time series of antigen and bcells
	    for(int t=0 ; t<NT; t++){
	        fout_antigen << Time_series_Antigen[t] << endl;
	        fout_N_active_linages << N_active_linages[t] << "\t";
	        for (int n = 0 ; n<n_naive ; n++){
	            fout_bcells << Time_series_Bcells[n][t] << "\t";
	        }
	        fout_bcells << endl;
	    }
	    fout.close();
	    fout_antigen.close();
	    fout_bcells.close();
	    fout_N_active_linages.close();
    }
    
    //------------------------------------------------------------------------------------
    cout<< ">Simulation completed…"<< endl;
    t2= clock();
    cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;

    return 0;
}