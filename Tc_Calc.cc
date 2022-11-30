//Tc Calculation Using Binders value

TFile *file_in;

const int TotalFile = 6;

TGraphErrors *gr[TotalFile];

int Dim_Range[TotalFile] = {8, 10, 12, 14, 16, 18};

TCanvas *c =  new TCanvas("c", "c", 800, 600);

void Tc_Calc(){

	for(int jentry=0; jentry<TotalFile; jentry++){

		file_in = new TFile( Form("Binders_Value%d.root", Dim_Range[jentry]) , "READ");		
		gr[jentry] = (TGraphErrors*)file_in->Get("gr_ub");

	}

	int i = 0;
	//Matrix8
	gr[i]->SetTitle(" Binders Cumulant Vs Temperature");
	gr[i]->SetMarkerColor(kCyan+2); 
	gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerStyle(31); 				// 31 : Astrik
	gr[i]->GetXaxis()->SetTitle(" Temperature (K) ");
	gr[i]->GetXaxis()->CenterTitle();
	gr[i]->GetYaxis()->SetTitle(" Binders Cumulant ");
	gr[i]->GetYaxis()->CenterTitle();

	//Matrix10
	i = 1;
	gr[i]->SetTitle(" Binders Cumulant Vs Temperature");
	gr[i]->SetMarkerColor(kRed); 
	gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerStyle(20); 				// 20 : Solid Circle
	gr[i]->GetXaxis()->SetTitle(" Temperature (K) ");
	gr[i]->GetXaxis()->CenterTitle();
	gr[i]->GetYaxis()->SetTitle(" Binders Cumulant ");
	gr[i]->GetYaxis()->CenterTitle();

	i = 2; // Matrix12
	gr[i]->SetTitle(" Binders Cumulant Vs Temperature");
	gr[i]->SetMarkerColor(kBlue); 
	gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerStyle(21); 				// 21 : Solid Square
	gr[i]->GetXaxis()->SetTitle(" Temperature (K)");
	gr[i]->GetXaxis()->CenterTitle();
	gr[i]->GetYaxis()->SetTitle(" Binders Cumulant ");
	gr[i]->GetYaxis()->CenterTitle();

	i = 3; // Matrix14
	gr[i]->SetTitle(" Binders Cumulant Vs Temperature");
	gr[i]->SetMarkerColor(kMagenta); 
	gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerStyle(29); 				// 29 : Solid Star
	gr[i]->GetXaxis()->SetTitle(" Temperature (K) ");
	gr[i]->GetXaxis()->CenterTitle();
	gr[i]->GetYaxis()->SetTitle(" Binders Cumulant ");
	gr[i]->GetYaxis()->CenterTitle();

	i = 4; // Matrix16
	gr[i]->SetTitle(" Binders Cumulant Vs Temperature");
	gr[i]->SetMarkerColor(kGreen+2); 
	gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerStyle(34); 				// 34 : Solid Cross
	gr[i]->GetXaxis()->SetTitle(" Temperature (K) ");
	gr[i]->GetXaxis()->CenterTitle();
	gr[i]->GetYaxis()->SetTitle(" Binders Cumulant ");
	gr[i]->GetYaxis()->CenterTitle();

	i = 5; // Matrix18
	gr[i]->SetTitle(" Binders Cumulant Vs Temperature");
	gr[i]->SetMarkerColor(kOrange+1); 
	gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerStyle(34); 				// 22 : Solid Up-Traingle
	gr[i]->GetXaxis()->SetTitle(" Temperature (K) ");
	gr[i]->GetXaxis()->CenterTitle();
	gr[i]->GetYaxis()->SetTitle(" Binders Cumulant ");
	gr[i]->GetYaxis()->CenterTitle();

	/* A for Axis (Always include to only the main graph to show complete range )
	 * 
	 */
	gr[1]->Draw("APLE");
	gr[0]->Draw("LPE");
	gr[2]->Draw("LPE");
	gr[3]->Draw("LPE");
	gr[4]->Draw("LPE");
	gr[5]->Draw("LPE");


	TLegend *leg = new TLegend(0.7302937,0.7132196,0.8956723,0.8965885,NULL,"brNDC");
   	leg->SetBorderSize(1);
   	leg->SetTextFont(22);
   	leg->SetLineColor(1);
   	leg->SetLineStyle(1);
   	leg->SetLineWidth(1);
   	leg->SetFillColor(0);
   	leg->SetFillStyle(1001);

   	leg->SetHeader("3-D Heisenberg System - Tc","C"); // option "C" allows to center the header
	leg->AddEntry(gr[5],"Lattice Dim. 18 ","lfp");
	leg->AddEntry(gr[4],"Lattice Dim. 16 ","lfp");
   	leg->AddEntry(gr[3],"Lattice Dim. 14 ","lfp");
	leg->AddEntry(gr[2],"Lattice Dim. 12 ","lfp");
	leg->AddEntry(gr[1],"Lattice Dim. 10 ","lfp");
	leg->AddEntry(gr[0],"Lattice Dim. 08 ","lfp");
	leg->Draw();

}
