/*

   Code to print data stored in the out%f.root files.
   Reference to code ~/3rdSem/practice.cc
   Written by Aranya Giri (2077644)

*/

//Necessary Header Files

#include<TFile.h>
#include<TTree.h>
#include<TMath.h>
#include<TH1.h>

TFile 	*file_in;	//To read the input file.
TTree	*tree;	//To read the tree inside the input file
TFile 	*file_out;	//To save the information in output file.

const int 	TotalFile = 25;

double 	Temp_Range[TotalFile] = {1.20, 1.25, 1.30, 1.32, 1.34, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50, 1.52, 1.54, 1.56, 1.60, 1.65};
int 	Dim = 18;
int 	matrix_size = 3;

int 	nodes = TMath::Power(Dim, matrix_size);

//Available file variables.
double 	mag; 	//Magnetization
int 	csize; 	//Cluster Size
double 	energy;	//Energy

double 	Error_Temp[TotalFile] = {0.0};

double 	AvgMag[TotalFile] = {0.0};	//Magnetisation Mean
double 	Error_AvgMag[TotalFile] = {0.0};	//Mag Mean Error

double 	AvgEnergy[TotalFile] = {0.0};	//Magnetisation Mean
double 	Error_AvgEnergy[TotalFile] = {0.0};	//Mag Mean Error

double AvgCSize[TotalFile] = {0.0};
double Error_AvgCSize[TotalFile] = {0.0};

double 	CorrTimeM[TotalFile] = {0.0};	//Correlation time
double 	CorrTimeE[TotalFile] = {0.0};	//Correlation time


TH1D *hmag[TotalFile];
TH1I *hcsize[TotalFile];
TH1D *henergy[TotalFile];

//Ploting

TGraphErrors *gr_mag;
TGraphErrors *gr_csize;
TGraphErrors *gr_energy;
TGraph *gr_tauM;
TGraph *gr_tauE;

//Function to find the Correlation at each time (sweep)

double Correlation_Calculator(long int time, long int timeMax, vector <double> quant){

	double chi_t = 0.0;

        double first_term = 0.0;

	for(long int tdash = 0; tdash< (timeMax - time); tdash++){
		first_term += quant[tdash]*quant[tdash + time];
	}

	first_term /= (timeMax - time);

	double second_term = 0.0;

	for(long int tdash = 0; tdash< (timeMax - time); tdash++){
		second_term += quant[tdash];
	}

	second_term /= (timeMax - time);

	double third_term = 0.0;

	for(long int tdash = 0; tdash< (timeMax - time); tdash++){
		third_term += quant[tdash + time];
	}

	third_term /= (timeMax - time);

	chi_t = first_term - (second_term*third_term) ;

	return chi_t;

}//Function Close

//Function to Calculate the Error
double Error_Calculation(double sigma, double tau, long int N){

	double error = 0.0;
	double delta_time = 1;

	if(N>0) {
		error = sigma*( TMath::Sqrt( (1 + (2*tau/delta_time) )/(N-1) ) );
		return error;
	}
	else return 0;
}

//_______________________Driver Code___________________________

void Tau_Calc(){

	for(int jentry=0; jentry < TotalFile; jentry++){
		
		double 	CorrTimeStep = 0.0;	//Correlation time
		double 	Corr_T0Step = 0.0;	//Initial Correalation at t= 0;

		vector <double> Rmag;
		vector <int>	Rcsize;
		vector <double> Renergy;

		hmag[jentry] 	=  new TH1D( Form("hmag%f",Temp_Range[jentry]), "Magnetisation",200,0,1);
		hcsize[jentry]	=  new TH1I( Form("hcsize%f", Temp_Range[jentry]), "Cluster Size", (nodes/2), 0, nodes);
		henergy[jentry]	=  new TH1D( Form("henergy%f",Temp_Range[jentry]), "Energy/Spin",200,0,1);

		file_in = new TFile( Form("Size%d/out%f.root", Dim, Temp_Range[jentry]), "READ" );

		long int Entries = 0;
		long int timeMax = 0;

		tree = (TTree*)file_in->Get("tree");
		Entries = tree->GetEntries();

		tree->SetBranchAddress("mag", &mag);
		tree->SetBranchAddress("csize", &csize);
		tree->SetBranchAddress("energy", &energy);

		vector <double> CorrM, CorrE;

		for(int entry = 0; entry < Entries ; entry++){
			
			tree->GetEntry(entry);
			Rmag.push_back(mag);
			Rcsize.push_back(csize);
			Renergy.push_back(energy);

			hmag[jentry]->Fill(mag);
			hcsize[jentry]->Fill(csize);
			henergy[jentry]->Fill(energy);

			AvgMag[jentry] += mag;
			AvgCSize[jentry] +=csize;
			AvgEnergy[jentry] += (energy/nodes) ;

		}

		//Calculating coorelation time_steps with magnetisation data
		timeMax = Entries;

		for(long int time = 0; time < timeMax; time++){

			double Chi_calc = Correlation_Calculator(time, timeMax, Rmag);
			if(Chi_calc > 0 ) CorrM.push_back( Chi_calc);
			else break;
		}

		Corr_T0Step = CorrM[0];

		for(long int time = 0; time < (CorrM.size()-1); time++){
			
			CorrM[time] = CorrM[time]/Corr_T0Step;
			CorrM[time+1] = CorrM[time+1]/Corr_T0Step;

			CorrTimeStep += 0.5*(CorrM[time]+CorrM[time+1]);
		}

		//Calculating Correlation time
		AvgMag[jentry] /= Entries;
		AvgCSize[jentry] /= Entries;

		CorrTimeM[jentry] = (CorrTimeStep*AvgCSize[jentry])/nodes ;
		
		//Error Calculation
		double StdDev_AvgMag 	= 0.0;
		double StdDev_AvgCSize 	= 0.0;
		
		for(int entry =0; entry<Entries; entry++){
 		
			StdDev_AvgMag += TMath::Power( (Rmag[entry] - AvgMag[jentry]), 2);
			StdDev_AvgCSize += TMath::Power( (Rcsize[entry] - AvgCSize[jentry] ) , 2);

		}

		StdDev_AvgMag 	= TMath::Sqrt( (StdDev_AvgMag/Entries) );
		StdDev_AvgCSize = TMath::Sqrt( (StdDev_AvgCSize/Entries) );

		Error_AvgMag[jentry] = Error_Calculation( StdDev_AvgMag, CorrTimeM[jentry], Entries);
		Error_AvgCSize[jentry] = Error_Calculation( StdDev_AvgCSize, CorrTimeM[jentry], Entries);
		
		//------------------------------
		//calculating coorelation time step with energy data
	/*
		for(long int time = 0; time < timeMax; time++){

			double Chi_calc = Correlation_Calculator(time, timeMax, Renergy);
			if(Chi_calc > 0 ) CorrE.push_back( Chi_calc);
			else break;
		}

		Corr_T0Step = CorrE[0];
		CorrTimeStep = 0.0;

		for(long int time = 0; time < (CorrE.size()-1); time++){
			
			CorrE[time] = CorrE[time]/Corr_T0Step;
			CorrE[time+1] = CorrE[time+1]/Corr_T0Step;

			CorrTimeStep += 0.5*(CorrE[time]+CorrE[time+1]);
		}
*/
		//Calculating Correlation time
		AvgEnergy[jentry] /= Entries;
/*
		CorrTimeE[jentry] = (CorrTimeStep*AvgCSize[jentry])/nodes ;
		
		//Error Calculation
		double StdDev_AvgEnergy	= 0.0;
		
		for(int entry =0; entry<Entries; entry++){
 		
			StdDev_AvgEnergy += TMath::Power( (Renergy[entry] - AvgEnergy[jentry]), 2);
		}

		StdDev_AvgEnergy = TMath::Sqrt( (StdDev_AvgEnergy/Entries) );

		Error_AvgEnergy[jentry] = Error_Calculation( StdDev_AvgEnergy, CorrTimeE[jentry], Entries);
*/
		file_in->Close();

	} //for loop for file entry close

	file_out = new TFile( Form("TauM_TauE_Mag_Energy_CSize_Calc%d.root", Dim), "RECREATE" );

	for(int i=0; i<TotalFile; i++){
		
		hmag[i]->Write();
		hcsize[i]->Write();
		henergy[i]->Write();
	
	}

	gr_tauM = new TGraph(TotalFile, Temp_Range, CorrTimeM);
	gr_tauM->SetName("gr_tauM");
	gr_tauM->SetTitle("Tau with Mag data Vs Temp ");
   	gr_tauM->SetMarkerColor(4);
  	gr_tauM->SetMarkerStyle(21);

	gr_tauE = new TGraph(TotalFile, Temp_Range, CorrTimeE);
	gr_tauE->SetName("gr_tauE");
	gr_tauE->SetTitle("Tau with Energy data Vs Temp ");
   	gr_tauE->SetMarkerColor(4);
  	gr_tauE->SetMarkerStyle(21);

	gr_mag = new TGraphErrors(TotalFile, Temp_Range, AvgMag, Error_Temp, Error_AvgMag);
	gr_mag->SetName("gr_mag");
	gr_mag->SetTitle("Mag./Spin Avg Vs Temp ");
   	gr_mag->SetMarkerColor(4);
  	gr_mag->SetMarkerStyle(21);

	gr_energy = new TGraphErrors(TotalFile, Temp_Range, AvgEnergy, Error_Temp, Error_AvgEnergy);
	gr_energy->SetName("gr_energy");
	gr_energy->SetTitle("Energy/Spin Avg Vs Temp ");
   	gr_energy->SetMarkerColor(4);
  	gr_energy->SetMarkerStyle(21);

	gr_csize = new TGraphErrors(TotalFile, Temp_Range, AvgCSize, Error_Temp, Error_AvgCSize);
	gr_csize->SetName("gr_csize");
	gr_csize->SetTitle("Cluster Size Avg Vs Temp ");
   	gr_csize->SetMarkerColor(4);
  	gr_csize->SetMarkerStyle(21);

	gr_tauM->Write();
	gr_tauE->Write();
	gr_mag->Write();
	gr_energy->Write();
	gr_csize->Write();

  	//gr->Draw("ALP");

	cout<<"----------------------End of Program--------------------"<<endl;

} //Driver Code Close
