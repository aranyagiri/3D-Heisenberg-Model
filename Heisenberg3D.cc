/* Advanced Computational Coursework
 * Project -  Magnetisation of 3D Lattice ( Heisenberg Model ) 
 * Monte Carlo Simulation using single cluster Wolff Algorithm
 * Code written in C++ and Compiler-Software is ROOT-CERN
 * By Aranya Giri, University of Houston.
 * Output is a root file "SizeoutTemp.root" for each lattice size and temperature.
 * Inside the root file : A tree with three branches 1. Magnetisation/Spin 2. Energy 3. Cluster Size (after each step).
*/

//Necessary Header Files
#include<climits>

//ROOT-CERN Header Files
#include<TMath.h>
#include<TCanvas.h>
#include<TH1.h>
#include<TFile.h>
#include<TRandom2.h>
#include <TVector3.h>
#include <TMath.h>

//Temperature Range - Choose suitable temperature (NOte : Critical Temperature for 3D Heisenberg is 1.442929 K as per ref.)
const int Temp_Count = 25;
double Temp_Range[Temp_Count]	={1.20, 1.25, 1.30, 1.32, 1.34, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.52, 1.54, 1.56, 1.60, 1.65};

// Definition of parameters for simulation
const int    	Dim   		 	= 	18; 					//Lattice Dimension
const double    H     		 	= 	0;					//Magnetic Field
const int 	total_step	 	= 	1000000;				//1 step = 1 cluster flipping.
const int 	matrix_space 	 	= 	3; 					//3-D Matrix Lattice Structure
const int 	spin_space   		= 	3; 					//Spin has 3-D freedom
const double 	J			= 	1;					//Interaction term (>0 for Ferromagnetism)

//Definition of dummy parameters
int 		nodes 		 	=       TMath::Power(Dim, matrix_space);	//Total no. of Nodes		

//Random numbers with seed in the argument.
TRandom2* 	rnd_create_theta 	= 	new TRandom2(0);			//Random spin creation at each Site
TRandom2* 	rnd_create_phi 	 	= 	new TRandom2(0);			//Random spin creation at each site
TRandom2* 	rnd_Sitepick 	 	= 	new TRandom2(0); 			//Random site pick to form cluster around it
TRandom2* 	rnd_Vecpick 	 	= 	new TRandom2(0);			//Random unit vector creation for cluster formation
TRandom2* 	rnd_link 	 	= 	new TRandom2(0);			//Random linking probability

//3-D Array
TVector3 	Lattice[Dim][Dim][Dim];							//Each Site holds spins state having 3-D freedom		
int 		DummyLattice[Dim][Dim][Dim];						//For Cluster Formation

//Data Storage
TFile 		*file;									//To store Tree, Histograms for each temp
TTree 		*tree;									//To store the mag. energy & cluster size data in the file

//Function Definition 
void 		InitializeTheLattice();							//To initilise the Lattice sites with spins having infinite degrees of freedom in 3D
void 		PrintDummyLattice();							//Printing the dummy lattice for verification purposes
void 		Cluster_Forming(int site[], double normal[], double Temp);		//Cluster formation using wolff algorithm
double 		ProbabilityLink(TVector3 Su, TVector3 Sj, TVector3 u, double Temp);	//Probablity for neighbours to join the cluster 
void 		Cluster_Flipping(double normal[]);					//Flipping the spin state w.r.t picked normal vector
int 		Cluster_Size();								//To calculate the cluster size at each step
double 		Magnetisation();							//To calculate the magnetisation of whole lattice after each cluster flipping (i.e., step)
double  	Energy();								//Same, to calculate the energy after each step

//_________________________ Functions ________________________________

void InitializeTheLattice(){

	double roll_theta;
	double roll_phi;

	for (int i = 0; i < Dim; i++){
		for(int j=0; j<Dim; j++){
			for(int k=0; k<Dim; k++){
				
				roll_theta = rnd_create_theta->Integer(INT_MAX);
				roll_theta = (roll_theta/INT_MAX)*TMath::Pi() ;

				roll_phi = rnd_create_phi->Integer(INT_MAX);
				roll_phi = (roll_phi/INT_MAX)*2*TMath::Pi() ;

				Lattice[i][j][k].SetX( TMath::Cos(roll_phi)*TMath::Sin(roll_theta) );
				Lattice[i][j][k].SetY( TMath::Sin(roll_phi)*TMath::Sin(roll_theta) );
				Lattice[i][j][k].SetZ( TMath::Cos(roll_theta) );

			}
		}
	}

}//Function-Close

//_________________________________________________________________________________________________________________________________
void PrintDummyLattice(){

	for (int i = 0; i< Dim; i++){
		for (int j = 0; j< Dim; j++){
			for (int k = 0; k< Dim; k++){

				cout<<DummyLattice[i][j][k];
			}
			cout<<endl;
		}
		cout<<endl;
	}

}//Function-Close

//__________________________________________________________________________________________________________
void Cluster_Forming(int site[], double normal[], double Temp){

	/* 	site[3] array contains picked(reference) site location.
	 *  	normal[3] array contains the components of Spin choosen randomly.
	 *  	Su[3] is the spic state at site[3].
	 *  	u[3] is the copy of normal[3]
	 *  	delta_s is like signam function +1 if >0 or -1 if <0
	*/
	
	TVector3 Su, u;

	Su.SetX( Lattice[ site[0] ][ site[1] ][ site[2] ].X() );
	Su.SetY( Lattice[ site[0] ][ site[1] ][ site[2] ].Y() );
	Su.SetZ( Lattice[ site[0] ][ site[1] ][ site[2] ].Z() );

	u.SetX( normal[0] );
	u.SetY( normal[1] );
	u.SetZ( normal[2] );

	double dot_su = Su.X()*u.X() + Su.Y()*u.Y() + Su.Z()*u.Z(); //Dot product of Su & u.
	double delta_s = TMath::Sign(1, dot_su );

	int value = DummyLattice[ site[0] ][ site[1] ][ site[2] ];

	/* If	value = 0 : not part of cluster.
	 * 	value = 1 : part of cluster.
	 * 	value = 2 : part of cluster with neighbours being checked for cluster formation.
	 */

	if( (value == 0 ) || (value == 1) ){

		//Assign 2 in the site as a marking that it's neighbour will be checked.
		DummyLattice[ site[0] ][ site[1] ][ site[2] ] = 2;

		int neibor[ (2*matrix_space) ]; // 2 times the matrix space is # neighbours.

		//0 (along x)
		neibor[0] = 	(int) ((site[0] + 1) % Dim);
		neibor[1] = 	(int) ((site[0] - 1 + Dim) % Dim);

		//1 (along y)
		neibor[2] = 	(int) ((site[1] + 1) % Dim);
		neibor[3] = 	(int) ((site[1] - 1 + Dim) % Dim);

		//2 (along z)
		neibor[4] = 	(int) ((site[2] + 1) % Dim);
		neibor[5] = 	(int) ((site[2] - 1 + Dim) % Dim);

		//Loop over all the neighbours of selected site.
		for(int n =0; n< (2*matrix_space) ; n++){

			TVector3 Sj;		//Neighbouring site under consideration
			double delta_j;		//Signam of choosen neighbour site with u unit vector
			int x,y,z;		//Copy of location of the considered neighbouring site

			double roll_link	= 	rnd_link->Integer(INT_MAX);
			roll_link		/= 	INT_MAX;

			if(n<=1){
				x = neibor[n];
				y = site[1];
				z = site[2];
			}
			if( (n>1) && (n<=3) ){	
				x = site[0];
				y = neibor[n];
				z = site[2];
			}
			if(n>3){
				x = site[0];
				y = site[1];
				z = neibor[n];
			}

			if( DummyLattice[x][y][z] == 0 ){ //To check weather the neighbour site is already member of cluster

				Sj.SetX( Lattice[x][y][z].X() );
				Sj.SetY( Lattice[x][y][z].Y() );
				Sj.SetZ( Lattice[x][y][z].Z() );

				double dot_sj = Sj.X()*u.X() + Sj.Y()*u.Y() + Sj.Z()*u.Z(); //Dot product of Sj & u
				delta_j = TMath::Sign(1, dot_sj );

				if( delta_s == delta_j ){ //**IMP condition for the neighbouring sites**

					double link_prob = ProbabilityLink(Su, Sj, u, Temp); //Find the linking probabilty to the cluster.
					if(roll_link < link_prob){ DummyLattice[x][y][z] = 1; }

				}
			}

		}//for-n close

	}//(value == 0 ) && (value == 1)

	//Recursive function calling to grow the cluster around picked site Su.
	for(int i=0; i<Dim; i++){
		for(int j=0; j<Dim; j++){
			for(int k=0; k<Dim; k++){	

				if(DummyLattice[i][j][k] == 1){ //Pass the neighbour as new Su.

					site[0] = i;
					site[1] = j;
					site[2] = k;
					Cluster_Forming( site, normal, Temp);
				}
			}
		}
	}

	return ;

}//Function-Close

//_________________________________________________________________________________________________________________________________
double ProbabilityLink(TVector3 Su, TVector3 Sj, TVector3 u, double Temp){

	double dot_su = Su.X()*u.X() + Su.Y()*u.Y() + Su.Z()*u.Z(); //Dot product of Su & u
	double dot_sj = Sj.X()*u.X() + Sj.Y()*u.Y() + Sj.Z()*u.Z(); //Dot product of Sj & u

	double pvalue = 1 - TMath::Exp( (-2*dot_su*dot_sj)/Temp ) ;
	return pvalue;

}//Function-Close

//_________________________________________________________________________________________________________________________________
void Cluster_Flipping(double normal[]){

	TVector3 u;

	u.SetX( normal[0] );
	u.SetY( normal[1] );
	u.SetZ( normal[2] );

	for(int i = 0; i<Dim; i++){
		for(int j = 0; j<Dim; j++){
			for(int k = 0; k<Dim; k++){
				if( DummyLattice[i][j][k] == 2){
					double dot_lattice = Lattice[i][j][k].X()*u.X() + Lattice[i][j][k].Y()*u.Y() + Lattice[i][j][k].Z()*u.Z();

					Lattice[i][j][k].SetX( Lattice[i][j][k].X() - 2*( dot_lattice )*( u.X() ) );
					Lattice[i][j][k].SetY( Lattice[i][j][k].Y() - 2*( dot_lattice )*( u.Y() ) );
					Lattice[i][j][k].SetZ( Lattice[i][j][k].Z() - 2*( dot_lattice )*( u.Z() ) );
				}
			}
		}
	}

}//Function-Close

//________________________________________________________________________________________________________________________________
int Cluster_Size(){

	int count = 0;

	for(int i = 0; i<Dim; i++){
		for(int j = 0; j<Dim; j++){
			for(int k = 0; k<Dim; k++){
				if( DummyLattice[i][j][k] == 2){  count ++; }
			}
		}
	}

	return count;
}//Function-Close

//_________________________________________________________________________________________________________________________________
double Magnetisation(){

	double m[3] = {0.0};
	double mvalue = 0.0;

	for(int i = 0; i<Dim; i++){
		for(int j = 0; j<Dim; j++){
			for(int k = 0; k<Dim; k++){
				m[0] += Lattice[i][j][k].X();	
				m[1] += Lattice[i][j][k].Y();
				m[2] += Lattice[i][j][k].Z();
			}
		}
	}
	mvalue = TMath::Sqrt( ( (m[0]*m[0]) + (m[1]*m[1]) + (m[2]*m[2]) ) )/nodes;

	return mvalue;
}//Function-Close

//________________________________________________________________________________________________________________________________
double Energy(){

	double evalue = 0.0;

	for(int i = 0; i<Dim; i++){
		for(int j = 0; j<Dim; j++){
			for(int k = 0; k<Dim; k++){

				int neibor[ (2*matrix_space) ]; // 2 times the matrix space is # neighbours.

				//i (along x)
				neibor[0] = 	(int) ((i + 1) % Dim);
				neibor[1] = 	(int) ((i - 1 + Dim) % Dim);

				//j (along y)
				neibor[2] = 	(int) ((j + 1) % Dim);
				neibor[3] = 	(int) ((j - 1 + Dim) % Dim);

				//k (along z)
				neibor[4] = 	(int) ((k + 1) % Dim);
				neibor[5] = 	(int) ((k - 1 + Dim) % Dim);

				//Loop over all the neighbours of selected site.
				for(int n =0; n< (2*matrix_space) ; n++){

					int x,y,z;		//Copy of location of the considered neighbouring site

					if(n<=1){
						x = neibor[n];
						y = j;
						z = k;
					}
					if( (n>1) && (n<=3) ){	
						x = i;
						y = neibor[n];
						z = k;
					}
					if(n>3){
						x = i;
						y = j;
						z = neibor[n];
					}

					//Dot product of spins with the nearest neighbours
					double dot_lattice = Lattice[i][j][k].X()*Lattice[x][y][z].X() + Lattice[i][j][k].Y()*Lattice[x][y][z].Y() + Lattice[i][j][k].Z()*Lattice[x][y][z].Z();
					evalue += dot_lattice;

				}//for-n-loop

			}//for-k-loop
		}//for-k-loop
	}//for-i-loop

	return (-J*evalue);

}//Function-Close

//_____________________Driver Code________________________________________________________________________________________________
void Heisenberg3D(){

	for(int temp_itr = 0; temp_itr<Temp_Count; temp_itr++){ //for different temperatures.

		//Variable for different calculation
		double 		mag 		= 	0.0;
		double 		energy		= 	0.0;
		int 		csize		=	0;

		//Save the files of all temperatures in folders with Dim number
		file = new TFile(Form("Size%d/out%f.root",Dim,Temp_Range[temp_itr]), "RECREATE");
		
		tree = new TTree("tree","tree");	
		tree->Branch( "mag", 	&mag, 	"mag/D");
		tree->Branch( "energy",	&energy,"energy/D");
		tree->Branch( "csize", 	&csize, "csize/I");	

		InitializeTheLattice();			

		for(int step=0; step < total_step; step++ ){

			//initialise dummy lattice
			for(int i=0; i< Dim; i++){
				for(int j=0; j<Dim; j++){
					for(int k=0; k<Dim; k++){
						DummyLattice[i][j][k] = 0;
					}
				}
			}

			//picking up a random spin site 's' 
			//and picking up a random direction vector 'u'

			double roll[matrix_space];
			double vec;
			double NormSq = 0.0;
			int site_pick[matrix_space];		//Contain indices of picked site
			double vec_pick[spin_space];		//Contain direction of picked spin

			for(int i=0; i<matrix_space; i++){
				roll[i] 	= 	rnd_Sitepick->Integer(INT_MAX);
				site_pick[i]	= 	(int) ( Dim*(roll[i]/INT_MAX) );
			}

			for(int axis = 0; axis< spin_space; axis++ ){

				vec		= 	rnd_Vecpick->Integer(INT_MAX);
				vec_pick[axis]	=	vec/INT_MAX;
				NormSq		+=	vec_pick[axis]*vec_pick[axis];

			}

			//Make u a unit vector
			for(int axis=0; axis< spin_space; axis++){
				vec_pick[axis]	/= 	TMath::Sqrt( NormSq );
			}

			Cluster_Forming(site_pick, vec_pick, Temp_Range[temp_itr] );
			Cluster_Flipping(vec_pick);
			
			if( (step % 50000) == 0) { //notification
				cout<<"Cluster Steps Covered : "<<step<<" Magnetisation : "<<mag<<endl;
			}

			if(step > (total_step/10) ){

				csize 	= 	Cluster_Size();
				mag 	= 	Magnetisation();
				energy 	= 	Energy();

				tree->Fill();
			}

		} //for step-loop close

		//hMagnet->Draw();

		tree->Write();
		file->Close();

		cout<<"End for Temp. : "<<Temp_Range[temp_itr]<<endl;
	}//for-temp close

	cout<<"******************** End of Simulation ***************"<<endl;
}//Driver Code Close
