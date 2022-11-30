//V-Value Calculation
//Header Files

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>

const int TotalFile = 6;
const int pro_count = 4; //m, m^2, m^3, m^4

//double 	Temp_Range[TotalFile] = {1.30, 1.32, 1.34, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.52, 1.54, 1.56};

double Temp = 1.44; //1.45 & 1.44
int Dim_Range[TotalFile] = {8, 10, 12, 14, 16, 18};

TFile *file_in;
TTree *tree;
TFile *file_out;

double mag, energy; //per spin is in the data

//Cumulants 1 For M : Magnetisation/Spin & E : Energy

double CN_mag[pro_count][TotalFile] = {0.0}; //C1, C2, C3, C4 for mag
double C1_energy[TotalFile] = {0.0}; //C1 for Energy
double C1_pro[pro_count][TotalFile] = {0.0}; //C1 for mag^{1,2,3,4}*Energy

double V1[TotalFile] = {0.0};
double V1_Error[TotalFile] = {0.0};

double Dim_Log[TotalFile] = {0.0};
double Dim_Error[TotalFile] = {0.0};

TGraphErrors *gr_V1;

//Functions_________________________________________________-
//Nth-Moment
/*
double MN_Calculation(vector <double> quant, long int totalentry, double mean, int N){

	double mN = 0.0;

	for(int i =0; i< totalentry; i++){
		mN += TMath::Power( (quant[i] - mean), N);
	}
	
	mN /= totalentry;
	return mN;	

}//Function Close
*/

//1st-Moment
double M1_Calculation(vector <double> quant, long int totalentry, int N){

	double m1 = 0.0;

	for(int i =0; i< totalentry; i++){
		m1 += TMath::Power(quant[i], N);
	}

	m1 /= totalentry;
	return m1;

}//Function Close

double V1_Calculation( double m3, double m4, double E, double m3E, double m4E){

	double v1 = 0.0;

	double m3_derivative = ( m3E - m3*E)*(-1/(Temp*Temp) );
	double m4_derivative = ( m4E - m4*E)*(-1/(Temp*Temp) );

	double m3_log = TMath::Log(m3_derivative);
	double m4_log = TMath::Log(m4_derivative);

	v1 = 4*m3_log - 3*m4_log;

	return v1;
}//Function Close

double Fit_Function(double *x, double *par ){

	double fit_value = 0;

	//par[0] = m
	//par[1] = c
	
	fit_value = par[0]*x[0] + par[1];

	return fit_value;

}

//Driver Function
void V1_Value(){

	for(int jentry=0; jentry<TotalFile; jentry++){

		vector <double> Rmag, Renergy;
		vector <double> Rpro1, Rpro2, Rpro3, Rpro4;

		file_in = new TFile( Form("Size%d/out%f.root", Dim_Range[jentry], Temp), "READ");

		long int Entries = 0;

		tree = (TTree*)file_in->Get("tree");
		Entries = tree->GetEntries();

		tree->SetBranchAddress("mag", &mag);
		tree->SetBranchAddress("energy", &energy);

		for(int entry=0; entry<Entries; entry++){

			tree->GetEntry(entry);
			Rmag.push_back(mag);
			Renergy.push_back(energy);

			Rpro1.push_back( ( ( TMath::Power(mag, 1) )*energy) );
			Rpro2.push_back( ( ( TMath::Power(mag, 2) )*energy) );
			Rpro3.push_back( ( ( TMath::Power(mag, 3) )*energy) );
			Rpro4.push_back( ( ( TMath::Power(mag, 4) )*energy) );

		}//entry for loop close

		C1_energy[jentry] = M1_Calculation(Renergy, Entries, 1);
		
		for(int i=0; i<pro_count; i++){
			CN_mag[i][jentry] = M1_Calculation(Rmag, Entries, i+1);
		}
		
		C1_pro[0][jentry] = M1_Calculation(Rpro1, Entries, 1);
		C1_pro[1][jentry] = M1_Calculation(Rpro2, Entries, 1);
		C1_pro[2][jentry] = M1_Calculation(Rpro3, Entries, 1);
		C1_pro[3][jentry] = M1_Calculation(Rpro4, Entries, 1);

		V1[jentry] = V1_Calculation( CN_mag[2][jentry], CN_mag[3][jentry], C1_energy[jentry], C1_pro[2][jentry], C1_pro[3][jentry]);

		Dim_Log[jentry] = TMath::Log(Dim_Range[jentry]);
	}//File for loop Close

	file_out = new TFile(Form("V1_nu_expo%f.root",Temp), "RECREATE" );

	gr_V1 = new TGraphErrors(TotalFile, Dim_Log, V1, Dim_Error, V1_Error);
	gr_V1->SetName("gr_V1");
	gr_V1->SetTitle("V_1 Vs Log(L) ");
        gr_V1->SetMarkerColor(kRed); 
	gr_V1->SetMarkerSize(2);
	gr_V1->SetMarkerStyle(20); 				// 20 : Solid Circle
	gr_V1->GetXaxis()->SetTitle(" Log(L) ");
	gr_V1->GetXaxis()->CenterTitle();
	gr_V1->GetYaxis()->SetTitle(" V_1 ");
	gr_V1->GetYaxis()->CenterTitle();

	TF1 *func = new TF1("fit", Fit_Function, 2, 0, 3);

	func->SetParNames("1/nu","V_1(tL^(1/nu))");
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

	gr_V1->Fit("fit","R");
	gr_V1->Draw(); //Just to visualize

	gr_V1->Write();

}
