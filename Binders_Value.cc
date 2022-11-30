//To calculate Binders Cumulant for a Lattice Size 

//Header Files

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>

const int TotalFile = 21;
int Dim = 8;

double 	Temp_Range[TotalFile] = {1.30, 1.32, 1.34, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.52, 1.54, 1.56};

TFile *file_in;
TTree *tree;
TFile *file_out;

double mag; //per spin is in the data

double M4[TotalFile] = {0.0};
double M2[TotalFile] = {0.0};
double M1[TotalFile] = {0.0};

double M4_Error[TotalFile] = {0.0};
double M2_Error[TotalFile] = {0.0};
double M1_Error[TotalFile] = {0.0};

double U_B[TotalFile] = {0.0};
double U_B_Error[TotalFile] = {0.0};

double Temp_Error[TotalFile] = {0.0};

TGraphErrors *gr_m2;
TGraphErrors *gr_m4;
TGraphErrors *gr_ub;

//Functions_________________________________________________-

double UB_Calculation( double m4, double m2){

	double u_b = 0.0;
	u_b = 1 - ( m4/(3*m2*m2) );
	return u_b;
}
//1st-Moment
double M1_Calculation(vector <double> quant, long int totalentry, int N){

	double m1 = 0.0;

	for(int i =0; i< totalentry; i++){
		m1 += TMath::Power(quant[i],N);
	}

	m1 /= totalentry;

	return m1;

}//Function Close

//Driver Function
void Binders_Value(){

	for(int jentry=0; jentry<TotalFile; jentry++){

		vector <double> Rmag;

		file_in = new TFile( Form("Size%d/out%f.root", Dim, Temp_Range[jentry]), "READ");

		long int Entries = 0;

		tree = (TTree*)file_in->Get("tree");
		Entries = tree->GetEntries();

		tree->SetBranchAddress("mag", &mag);

		for(int entry=0; entry<Entries; entry++){
		
			tree->GetEntry(entry);
			Rmag.push_back(mag);

		}//entry for loop close

		M1[jentry] = M1_Calculation(Rmag, Entries,1);
		M2[jentry] = M1_Calculation(Rmag, Entries,2);
		M4[jentry] = M1_Calculation(Rmag, Entries,4);
		
		U_B[jentry] = UB_Calculation( M4[jentry], M2[jentry] );

	}//File for loop Close

	file_out = new TFile( Form("Binders_Value%d.root", Dim), "RECREATE" );

	gr_m2 = new TGraphErrors(TotalFile, Temp_Range, M2, Temp_Error, M2_Error);
	gr_m2->SetName("gr_m2");
	gr_m2->SetTitle("m2 Vs Temp ");
   	gr_m2->SetMarkerColor(4);
  	gr_m2->SetMarkerStyle(21);

	gr_m4 = new TGraphErrors(TotalFile, Temp_Range, M4, Temp_Error, M4_Error);
	gr_m4->SetName("gr_m4");
	gr_m4->SetTitle("m4 Vs Temp ");
   	gr_m4->SetMarkerColor(4);
  	gr_m4->SetMarkerStyle(21);

	gr_ub = new TGraphErrors(TotalFile, Temp_Range, U_B, Temp_Error, U_B_Error);
	gr_ub->SetName("gr_ub");
	gr_ub->SetTitle("u_b Vs Temp ");
   	gr_ub->SetMarkerColor(4);
  	gr_ub->SetMarkerStyle(21);

	gr_m2->Write();
	gr_m4->Write();
	gr_ub->Write();

}//Driver Function Close
