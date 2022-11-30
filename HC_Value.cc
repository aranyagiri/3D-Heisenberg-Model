//Finding the alpha exponent

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>

const int TotalFile = 6;

double Temp = 1.45;

int Dim_Range[TotalFile] = {8, 10, 12, 14, 16, 18};
const int matrix_size = 3;

TFile *file_in;
TTree *tree;
TFile *file_out;

double energy;

double C1_energy = 0.0;
double C2_energy = 0.0;

double heat_cap[TotalFile] = {0.0}; //specific heat capacity per spin
double heat_cap_Log[TotalFile] = {0.0};

double C1_energy_Error = 0.0;
double C2_energy_Error = 0.0;

double heat_cap_Error[TotalFile] = {0.0};

double Dim_Log[TotalFile] = {0.0};
double Dim_Error[TotalFile] = {0.0};

//Graph for mean Energy, Sigma of E and Heat Capacity
TGraphErrors *gr_mE, *gr_sdE, *gr_hc;

double C1_Calculation(vector <double> quant, long int totalentry){

	double C1 = 0.0;

	for(int i =0; i< totalentry; i++){
		C1 += quant[i];
	}

	C1 /= totalentry;
	return C1;

}//Function Close

double C2_Calculation(vector <double> quant, long int totalentry, double mean){

	double C2 = 0.0;

	for(int i =0; i< totalentry; i++){
		C2 += TMath::Power( (quant[i] - mean), 2);
	}
	
	C2 /= totalentry;
	return C2;

}//Function Close

double HC_Calculation(double E2, double N){

	double ChiperSpin = 0.0;

	ChiperSpin = E2*(1/(Temp*Temp));
	ChiperSpin /= N;

	return ChiperSpin;
}//Function Close

double Fit_Function(double *x, double *par){

	double fit_value = 0;

	//par[0] = m
	//par[1] = c
	
	fit_value = par[0]*x[0] + par[1];

	return fit_value;

}

//Driver Function
void HC_Value(){
	for(int jentry=0; jentry<TotalFile; jentry++){

		vector <double> Renergy;

		file_in = new TFile( Form("Size%d/out%f.root", Dim_Range[jentry], Temp), "READ");

		long int Entries = 0;

		tree = (TTree*)file_in->Get("tree");
		Entries = tree->GetEntries();

		tree->SetBranchAddress("energy", &energy);

		for(int entry=0; entry<Entries; entry++){

			tree->GetEntry(entry);
			Renergy.push_back(energy);

		}//entry for loop close

		C1_energy = C1_Calculation(Renergy, Entries);
		C2_energy = C2_Calculation(Renergy, Entries, C1_energy);

		double nodes = TMath::Power(Dim_Range[jentry], matrix_size);
		heat_cap[jentry] = HC_Calculation(C2_energy, nodes);

		heat_cap_Log[jentry] = TMath::Log( heat_cap[jentry] );
		Dim_Log[jentry] = TMath::Log(Dim_Range[jentry]);

	}//File for loop Close

	file_out = new TFile("alpha_expo.root", "RECREATE" );

	gr_hc = new TGraphErrors(TotalFile, Dim_Log, heat_cap_Log, Dim_Error, heat_cap_Error);
	gr_hc->SetName("gr_hc");
	gr_hc->SetTitle("Heat Capacity Vs Log(L) ");
	gr_hc->SetMarkerColor(4);
	gr_hc->SetMarkerStyle(21);

	TF1 *func = new TF1("fit", Fit_Function, 2, 0, 3);
	func->SetParNames("slope","y-intercept");
	func->SetParameters(1 ,1);

	//plot the fitted graph
	TCanvas *c1 = new TCanvas("c1", "c1",72,64,1275,917);
	gStyle->SetOptFit(111);
	gStyle->SetOptStat(0);
	c1->Range(1.094951,-186.2599,8.797976,1220.16);
	c1->SetFillColor(0);
	c1->SetBorderMode(0);
	c1->SetBorderSize(2);
	c1->SetTopMargin(0.06846239);
	c1->SetBottomMargin(0.1324355);
	c1->SetFrameBorderMode(0);
	c1->SetFrameLineStyle(2);
	c1->SetFrameBorderMode(0);

	gr_hc->Fit("fit","R");
	gr_hc->Draw(); //Just to visualize

	gr_hc->Write();


}//Driver Function Close
