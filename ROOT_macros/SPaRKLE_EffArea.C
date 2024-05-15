#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH1DLog.h"

using namespace std;

int SPaRKLE_EffArea()
{
    cout << "SPaRKLE_EffArea" << endl;
    TString fname = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST/gamma_19042024_0_um.root";
    TFile *f = new TFile(fname, "READ");
    cout << "File opened" << endl;
    TTree *Edep = (TTree*)f->Get("Edep");
    cout << "Tree opened" << endl;

    Edep -> Draw("RandEnergy", "", "goff");

    double RandEnergy = 0.;
    Edep -> SetBranchAddress("RandEnergy", &RandEnergy);

    double Emin_gamma = 999999.;
    double Emax_gamma = -999999.;

    for(int i=0; i < Edep->GetEntries(); i++)
    {
        Edep -> GetEntry(i);
        if(RandEnergy < Emin_gamma) Emin_gamma = RandEnergy;
        if(RandEnergy > Emax_gamma) Emax_gamma = RandEnergy;
    }



    TH1DLog *hEdep = new TH1DLog();
    hEdep -> SetXAxis(Emin_gamma, Emax_gamma, 200);
    hEdep -> SetTitle("Energy Deposition");
    hEdep -> SetXTitle("Energy [MeV]");
    hEdep -> SetYTitle("Counts");
    hEdep -> GenerateHistogram();


    hEdep -> SetName("hNoSelection");
    hEdep -> GenerateHistogram();
    TH1D *hNoSelection = (TH1D*) hEdep -> GetHistogram();
    hEdep -> SetName("hSelection_bottom");
    hEdep -> GenerateHistogram();
    TH1D *hSelection_bottom = (TH1D*) hEdep -> GetHistogram();

    hEdep -> SetName("hSelection_noBottom");
    hEdep -> GenerateHistogram();
    TH1D *hSelection_noBottom = (TH1D*) hEdep -> GetHistogram();

    hNoSelection -> SetLineColor(kBlack);
    hSelection_bottom -> SetLineColor(kRed);
    hSelection_noBottom -> SetLineColor(kBlue);

    hNoSelection -> SetMarkerColor(kBlack);
    hSelection_bottom -> SetMarkerColor(kRed);
    hSelection_noBottom -> SetMarkerColor(kBlue);



    double Threshold_Calo = 0.04;
    double Threshold_Plastic = 0.04;

    Edep -> Draw("RandEnergy>>hNoSelection", "", "goff");
    


    TString Condition1 = "(Ed_DrilledVeto <" + std::to_string(Threshold_Plastic) + ")";
    TString Condition2 = "(Ed_BottomVeto <" + std::to_string(Threshold_Plastic) + ")";
    TString Condition3 = "(Ed_Calo0 >" + std::to_string(Threshold_Calo) + ")";

    cout << Condition1 << endl;
    cout << Condition2 << endl;
    cout << Condition3 << endl;

    TString Selection_noBottom = Condition1 + " && " + Condition3;
    TString Selection_Bottom = Condition1 + " && " + Condition2 + " && " + Condition3;


    Edep -> Draw("RandEnergy>>hSelection_noBottom", Selection_noBottom, "goff");
    Edep -> Draw("RandEnergy>>hSelection_bottom", Selection_Bottom, "goff");

    
    double Areagen = 20.*20.;

    for(int i=0; i < hNoSelection->GetNbinsX(); i++)
    {
        double binContent = hNoSelection -> GetBinContent(i);
        double binError = hNoSelection -> GetBinError(i);
        
        double binContent_bottom = hSelection_bottom -> GetBinContent(i);
        double binError_bottom = hSelection_bottom -> GetBinError(i);
        
        double binContent_noBottom = hSelection_noBottom -> GetBinContent(i);
        double binError_noBottom = hSelection_noBottom -> GetBinError(i);

        double effArea_bottom = binContent_bottom / binContent * Areagen;
        double effArea_noBottom = binContent_noBottom / binContent * Areagen;

        double effArea_bottom_error = effArea_bottom * TMath::Sqrt(TMath::Power(binError_bottom/binContent_bottom, 2) + TMath::Power(binError/binContent, 2));
        double effArea_noBottom_error = effArea_noBottom * TMath::Sqrt(TMath::Power(binError_noBottom/binContent_noBottom, 2) + TMath::Power(binError/binContent, 2));

        hSelection_bottom -> SetBinContent(i, effArea_bottom);
        hSelection_bottom -> SetBinError(i, effArea_bottom_error);

        hSelection_noBottom -> SetBinContent(i, effArea_noBottom);
        hSelection_noBottom -> SetBinError(i, effArea_noBottom_error);

    }


    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    //hNoSelection -> Draw("");
    hSelection_bottom -> Draw("E1 same");
    hSelection_noBottom -> Draw("E1 same");


    std::string outname = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST/EffArea.csv";

    ofstream outfile;
    outfile.open(outname);

    outfile << "#emin_bin, emax_bin, effArea_bottom, effArea_bottom_error, effArea_noBottom, effArea_noBottom_error" << endl;

    for(int i=0; i < hNoSelection->GetNbinsX(); i++)
    {
        double emin = hNoSelection -> GetBinLowEdge(i);
        double emax = hNoSelection -> GetBinLowEdge(i+1);
        double effArea_bottom = hSelection_bottom -> GetBinContent(i);
        double effArea_bottom_error = hSelection_bottom -> GetBinError(i);
        double effArea_noBottom = hSelection_noBottom -> GetBinContent(i);
        double effArea_noBottom_error = hSelection_noBottom -> GetBinError(i);

        // Remove NaN
        if(effArea_bottom != effArea_bottom) effArea_bottom = 0.;
        if(effArea_bottom_error != effArea_bottom_error) effArea_bottom_error = 0.;
        if(effArea_noBottom != effArea_noBottom) effArea_noBottom = 0.;
        if(effArea_noBottom_error != effArea_noBottom_error) effArea_noBottom_error = 0.;


        outfile << emin << ", " << emax << ", " << effArea_bottom << ", " << effArea_bottom_error << ", " << effArea_noBottom << ", " << effArea_noBottom_error << endl;
    }

    outfile.close();

    hSelection_bottom -> SetStats(0);
    hSelection_noBottom -> SetStats(0);

    hSelection_bottom -> SetTitle("Efficiency Area");
    hSelection_bottom -> SetXTitle("Energy [MeV]");
    hSelection_bottom -> SetYTitle("Effective Area [cm^2]");
    

    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg -> AddEntry(hSelection_bottom, "Bottom Veto VETO ON", "lpe");
    leg -> AddEntry(hSelection_noBottom, "Bottom Veto VETO OFF", "lpe");

    leg -> Draw();




    gPad -> SetLogx();
    gPad -> SetLogy();
    gPad -> SetGrid();

    










    return 0;

    


}