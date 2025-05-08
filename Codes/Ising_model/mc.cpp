//Template to run a Monte Carlo simulation for the Ising model on a 2D square lattice
//Input: (all parameters are already set internally!)
//Output: output.txt (columns: 1 temperature , 2 energy per spin, 3 heat capacity, 4 magnetisation per spin, 5 susceptibility  )


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;
//Library for random number generators
#include "./lib/random.cpp"
//There are two functions extracted from the library
//double randX(min,max): a random number between min and max
//int randIX(min,max): an integer random number between min and max  

//----------------------------------------------------------------------------------

//Function to calculate the energy: (Try not to sum over all N^2 possible entries of the coupling matrix but rather sum only over the bonds.)
double energy(double const & h, double const & J, int const & La, int const & Lb, vector< vector<int> > const & Pool)
{
	double E (0.);
    //FOR YOU TO FILL IN:
    /*
	for(int i=0;i<La-1;i++){
		for(int j=0;j<Lb-1;j++){
			double e=Pool[i][j]*Pool[i+1][j]*J + Pool[i][j]*Pool[i][j+1]*J + Pool[i][j]*h;
			E+=e;
		};
	};
	*/
	for(int i=0;i<La-1;i++){
		double e=Pool[i][0]*Pool[i+1][0]*J + Pool[i][0]*h;
		//double e=Pool[i][Lb-1]*Pool[i+1][Lb-1]*J + Pool[i][Lb-1]*Pool[i][0]*J + Pool[i][Lb-1]*h;
		E=E+e;
	};
	/*
	for(int j=0;j<Lb-1;j++){
		double e=Pool[La-1][j]*Pool[La-1][j+1]*J + Pool[La-1][j]*Pool[0][j]*J + Pool[La-1][j]*h;
		E=E+e;
	};
	*/
	E=E + Pool[La-1][0]*Pool[0][0]*J + Pool[La-1][0]*h;
	//E=E+Pool[Lb-1][La-1]*Pool[0][La-1]*J + Pool[Lb-1][La-1]*Pool[Lb-1][0]*J + Pool[Lb-1][La-1]*h;
    //
	return -E;
};

//Function to calculate the magnetisation:
double magnetisation(double const & h, double const & J, int const & La, int const & Lb, vector< vector<int> > const & Pool)
{
	double M (0.);
    //FOR YOU TO FILL IN:
	for(int i=0;i<La;i++){
		for(int j=0;j<Lb;j++){
			M+=Pool[i][j];
		};
	};
    //
	return M;
};


//Function to calculate the energy difference due to a spin flip
//Do NOT use energy(state 1)-energy(state 2) as this is computatioLally costly: the energy difference depends on the state of the four neighbor of the spins
inline double delt(double const & h, double const & J, int const & La, int const & Lb, vector< vector<int> > const & Pool, int const & ka, int const & kb)
{
	double deltaE (0.);
    //FOR YOU TO FILL IN

    if(ka!=La-1 && kb!=Lb-1 && ka!=0 && kb!=0){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[ka+1][kb] + Pool[ka][kb]*Pool[ka-1][kb] + Pool[ka][kb]*Pool[ka][kb+1] + Pool[ka][kb]*Pool[ka][kb-1]);
  	}
  	if(ka==La-1 && kb!=Lb-1 && kb!=0){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[0][kb] + Pool[ka][kb]*Pool[ka-1][kb] + Pool[ka][kb]*Pool[ka][kb+1] + Pool[ka][kb]*Pool[ka][kb-1]);
  	}
  	if(ka!=La-1 && kb==Lb-1 && ka!=0){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[ka+1][kb] + Pool[ka][kb]*Pool[ka-1][kb] + Pool[ka][kb]*Pool[ka][0] + Pool[ka][kb]*Pool[ka][kb-1]);
  	}
  	if(kb!=Lb-1 && ka==0 && kb!=0){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[ka+1][kb] + Pool[ka][kb]*Pool[La-1][kb] + Pool[ka][kb]*Pool[ka][kb+1] + Pool[ka][kb]*Pool[ka][kb-1]);
  	}
  	if(ka!=La-1 && ka!=0 && kb==0){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[ka+1][kb] + Pool[ka][kb]*Pool[ka-1][kb] + Pool[ka][kb]*Pool[ka][kb+1] + Pool[ka][kb]*Pool[ka][Lb-1]);
  	}
  	if(ka==La-1 && kb==Lb-1){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[0][kb] + Pool[ka][kb]*Pool[ka-1][kb] + Pool[ka][kb]*Pool[ka][0] + Pool[ka][kb]*Pool[ka][kb-1]);
  	}
  	if(ka==La-1 && kb==0){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[0][kb] + Pool[ka][kb]*Pool[ka-1][kb] + Pool[ka][kb]*Pool[ka][kb+1] + Pool[ka][kb]*Pool[ka][Lb-1]);
  	}
  	if(ka==0 && kb==Lb-1){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[ka+1][kb] + Pool[ka][kb]*Pool[La-1][kb] + Pool[ka][kb]*Pool[ka][0] + Pool[ka][kb]*Pool[ka][kb-1]);
  	}
  	if(ka==0 && kb==0){
  		deltaE = -2*J*(Pool[ka][kb]*Pool[ka+1][kb] + Pool[ka][kb]*Pool[La-1][kb] + Pool[ka][kb]*Pool[ka][kb+1] + Pool[ka][kb]*Pool[ka][Lb-1]);
  	}

  	deltaE = 2*h*(Pool[ka][kb]);

    //
	return deltaE;
};

//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	cout<<">Running Monte Carlo simulation of the Ising model on a 2D square lattice..."<< endl;
	clock_t t1,t2;
    t1=clock();
	//-----------------------------------------------------------------------------
	//Parameters: (they are fine as they are)
	int La (80), Lb (1);//Number of spins in each direction:
	int N (La*Lb); //Total number of spins
	double h (1.0), J(0.0); //External field and coupling
	double T1 (.1), T2(1); //Temperature min and max
	int nT (1); //Number of temperature points
	long long int n0 (0*N), N0 (100000*N), d0 (1); //Number of steps: initial prelude, total, distance between sampling points
	//Array of spins
	vector < vector<int> > Pool;
	Pool.resize(La);
	for (int ka= 0; ka<La; ka++)
	{
		(Pool[ka]).resize(Lb);
	};
	//Initiating spins at random
	for (int ka= 0; ka<La; ka++)
	{
		for (int kb= 0; kb<Lb; kb++)
		{
			Pool[ka][kb]= 1 ;
		};
	};
	//------------------------------------------------------------------------------
	ofstream fout ("output.txt");

	for (int kT= 0; kT<nT; kT++)
	{

		//double T (T2-(T2-T1)*double(kT)/double(nT-1));
		double T = T2;
		cout<< ">T= "<< T<< endl;
		
        
		double E, M; //Energy, magnetisation
		E= energy(h,J,La,Lb,Pool);
		fout<< E<< "\t"<<0.0 <<endl;
		M= magnetisation(h,J,La,Lb,Pool);

		double avrE (0.), avrEE (0.), avrM (0.), avrMM (0.); //Averages of E, E^2, M, M^2
        
        //Starting the trajectory:
		int countData (0); //Number of data point sampled in the trajectory
		for (long long int k= 0; k < N0; k++)
		{
			//FOR YOU TO FILL-IN:
			//Pick up a spin and calculate the energy difference if it were flipped
  			int a = randIX(0,La-1);
  			int b = randIX(0,Lb-1);
  			
  			double deltaE = delt(h,J,La,Lb,Pool,a,b);
			//Decide whether to actually flip the spin or not: (Metropolis' algorithm)
			if (deltaE<0){
				Pool[a][b]=Pool[a][b]*(-1);
			}
			else{
				double rand = randX(0,1);
				if(rand < exp((-1*deltaE)/T)){
					Pool[a][b]=Pool[a][b]*(-1);
				}
			};
			//
			
			//Calculate the observables starting after n0 steps: sample data points every d0 steps (for Eq.(2))
			if (k>=n0)
			{
				if ( (k%d0) == 0)
				{
					//Increase the number of data sum
					countData++;
					E= energy(h,J,La,Lb,Pool);
					fout<< E<< "\t"<<deltaE <<endl;
					M= magnetisation(h,J,La,Lb,Pool);
					//Update the sums in Eq.(2) for energy, heat capacity(EE suffices), magnetisation and susceptibility(MM suffices)
					//FOR YOU TO FILL-IN
					avrE+=E;
					avrM+=M;
					avrEE+=(E*E);
					avrMM+=(M*M);
					//
					
				};
			};
		};
        //Divide the sums by the number of data points to get the averages
		avrE/=double(countData); avrEE/=double(countData); avrM/=double(countData); avrMM/=double(countData);
		
		//FOR YOU TO FILL IN:
		double heatcapacity= (avrEE - avrE*avrE)/(T*T) ;
        double susceptibility= (avrMM - avrM*avrM)/(T);
		//
		
		//Write to an output file:
		//fout<< T<< "\t"<< avrE/double(N)<< "\t"<< heatcapacity/double(N)<< "\t"<< avrM/double(N)<< "\t"<< susceptibility/double(N) <<endl;
	};
	fout.close();
	//------------------------------------------------------------------------------
	cout<< ">Simulation completed…"<< endl;
	t2= clock();
	cout<< "(Running time: "<< double(t2-t1)/CLOCKS_PER_SEC <<" seconds.)"<< endl;
	return 0;
}
