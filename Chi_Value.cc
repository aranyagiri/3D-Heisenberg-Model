//Finding the alpha exponent

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>

const int TotalFile = 6;

double Temp = 1.44;

int Dim_Range[TotalFile] = {8, 10, 12, 14, 16, 18};
const int matrix_size = 3;

TFile *file_in;
TTree *tree;
TFile *file_out;

double mag;


double Chi[TotalFile] = {0.0}; //specific heat capacity per spin
double Chi_Log[TotalFile] = {0.0};
double Chi_Error[TotalFile] = {0.0};

double Dim_Log[TotalFile] = {0.0};
double Dim_Error[TotalFile] = {0.0};

//Graph for mean Magnetisation, Sigma of Magnetisation and Chi
TGraphErrors *gr_mM, *gr_sdM, *gr_Chi;

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
	}//
	
	C2 /= totalentry;
	return C2;

}//Function Close

double Chi_Calculation(double M2, double N){

	double ChiperSpin = 0.0;

	ChiperSpin = M2*(N/(Temp));

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
void Chi_Value(){
	for(int jentry=0; jentry<TotalFile; jentry++){

		vector <double> Rmag;

		double C1_mag = 0.0;
		double C2_mag = 0.0;
		
		file_in = new TFile( Form("Size%d/out%f.root", Dim_Range[jentry], Temp), "READ");

		long int Entries = 0;

		tree = (TTree*)file_in->Get("tree");
		Entries = tree->GetEntries();

		tree->SetBranchAddress("mag", &mag);

		for(int entry=0; entry<Entries; entry++){

			tree->GetEntry(entry);
			Rmag.push_back(mag);

		}//entry for loop close

		C1_mag = C1_Calculation(Rmag, Entries);
		C2_mag = C2_Calculation(Rmag, Entries, C1_mag );

		double nodes = TMath::Power(Dim_Range[jentry], matrix_size);
		Chi[jentry] = Chi_Calculation(C2_mag, nodes);

		Chi_Log[jentry] = TMath::Log( Chi[jentry] );
		Dim_Log[jentry] = TMath::Log(Dim_Range[jentry]);

	}//File for loop Close

	file_out = new TFile(Form("gamma_expo%f.root",Temp), "RECREATE" );

	gr_Chi = new TGraphErrors(TotalFile, Dim_Log, Chi_Log, Dim_Error, Chi_Error);
	gr_Chi->SetName("gr_Chi");
	gr_Chi->SetTitle("Magnetic Susceptiblity Vs Log(L) ");
	gr_Chi->SetMarkerColor(kRed); 
	gr_Chi->SetMarkerSize(2);
	gr_Chi->SetMarkerStyle(20); 				// 20 : Solid Circle
	gr_Chi->GetXaxis()->SetTitle(" Log(L) ");
	gr_Chi->GetXaxis()->CenterTitle();
	gr_Chi->GetYaxis()->SetTitle("Magnetic Susceptiblity");
	gr_Chi->GetYaxis()->CenterTitle();


	TF1 *func = new TF1("fit", Fit_Function, 2, 0, 3);
	func->SetParNames("gamma/nu","Scaling Factor");
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

	gr_Chi->Fit("fit","R");
	gr_Chi->Draw(); //Just to visualize

	gr_Chi->Write();


}//Driver Function Close
