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
#include "TChain.h"


#include "TH1DLog.h"
#include "TH2DLog.h"

using namespace std;


void Normalize_TH2D(TH2D *h)
{
    for(int i = 0; i < h->GetNbinsX(); i++)
    {
        double Counts = 0;
        for(int j = 0; j < h->GetNbinsY(); j++)
        {
            Counts += h->GetBinContent(i+1, j+1);

        }

        if(Counts == 0) continue;
        for(int j = 0; j < h->GetNbinsY(); j++)
        {
            double c = h->GetBinContent(i+1, j+1);
            double s = h->GetBinError(i+1, j+1);

            h->SetBinContent(i+1, j+1, c/Counts);
            h->SetBinError(i+1, j+1, s/Counts);
        }
    }

    return;
}








int Photons()
{
    TString filename = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST_Acc/gamma_2.root";
    TFile *f = new TFile(filename);
    TTree *Edep = (TTree*)f->Get("Edep");


    // Set the resolution of the detector
    double GAGG_Res = 0.1; // %
    double Plastic_Res = 0.3; // %
    double SD_Res = 6e-3; //MeV


    double Threshold_GAGG = 0.01; // MeV
    double Threshold_Plastic = 0.05; // MeV
    double Threshold_SD = 0.04; // MeV

    Edep->SetAlias("gCalo_A1",Form("(Calo_A1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_A1*%g)", GAGG_Res));
    Edep->SetAlias("gCalo_A2", Form("(Calo_A2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_A2*%g)", GAGG_Res));
    Edep->SetAlias("gCalo_B1", Form("(Calo_B1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_B1*%g)", GAGG_Res));
    Edep->SetAlias("gCalo_B2", Form("(Calo_B2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_B2*%g)", GAGG_Res));

    Edep->SetAlias("gSD1", Form("(SD1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g)", SD_Res));
    Edep->SetAlias("gSD2", Form("(SD2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g)", SD_Res));

    Edep->SetAlias("gVT", Form("(VT + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VT)",Plastic_Res));
    Edep->SetAlias("gVB", Form("(VB + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VB)", Plastic_Res));
    Edep->SetAlias("gVL0", Form("(VL0 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL0)", Plastic_Res));
    Edep->SetAlias("gVL1", Form("(VL1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL1)", Plastic_Res));
    Edep->SetAlias("gVL2", Form("(VL2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL2)", Plastic_Res));
    Edep->SetAlias("gVL3", Form("(VL3 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL3)", Plastic_Res));




    TString GoodEvents_gPID_1 = Form("(gCalo_A1 > %g) && (gSD1 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_PID_1 = Form("(Calo_A1 > %g) && (SD1 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);

    TString GoodEvents_gPID_2 = Form("(gCalo_A2 > %g) && (gSD1 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_PID_2 = Form("(Calo_A2 > %g) && (SD1 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);

    TString GoodEvents_photons_gPID_1 = Form("(gCalo_A1 > %g) && (gSD1 < %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_photons_PID_1 = Form("(Calo_A1 > %g) && (SD1 < %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_photons_gPID_2 = Form("(gCalo_A2 > %g) && (gSD1 < %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_photons_PID_2 = Form("(Calo_A2 > %g) && (SD1 < %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);

    TString GoodEvents_photons_gPID_1_NoVB  = Form("(gCalo_A1 > %g) && (gSD1 < %g) && (gVT < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_photons_PID_1_NoVB   = Form("(Calo_A1 > %g) && (SD1 < %g) && (VT < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_photons_gPID_2_NoVB  = Form("(gCalo_A2 > %g) && (gSD1 < %g) && (gVT < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_photons_PID_2_NoVB   = Form("(Calo_A2 > %g) && (SD1 < %g) && (VT < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);


    TString GoodEvents_gGBM_1 = Form("(gCalo_B1 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_GBM_1 = Form("(Calo_B1 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_gGBM_2 = Form("(gCalo_B2 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_GBM_2 = Form("(Calo_B2 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);

    TString GoodEvents_gGBM_1_NoVB = Form("(gCalo_B1 > %g) && (gVT < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_GBM_1_NoVB = Form("(Calo_B1 > %g) && (VT < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG,  Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_gGBM_2_NoVB = Form("(gCalo_B2 > %g) && (gVT < %g)  && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_GBM_2_NoVB = Form("(Calo_B2 > %g) && (VT < %g)  && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);




    // X axis
    double Egen_0 = 10e-3;
    double Egen_1 = 100.;
    int Nx = 200;


    double E_det_0 = 10e-3;
    double E_det_1 = 500.;
    int Ny = 200;

    TH2DLog *h_migr = new TH2DLog();
    h_migr -> SetXAxis(Egen_0, Egen_1, Nx);
    h_migr -> SetYAxis(E_det_0, E_det_1, Ny);
    h_migr -> SetTitle("");


    // Calo_A1
    h_migr -> SetName("h_A1");
    h_migr -> GenerateHistogram();
    TH2D *h_A1 = (TH2D*)h_migr->GetHistogram();

    // Calo_A2
    h_migr -> SetName("h_A2");
    h_migr -> GenerateHistogram();
    TH2D *h_A2 = (TH2D*)h_migr->GetHistogram();

    // Calo_B1
    h_migr -> SetName("h_B1");
    h_migr -> GenerateHistogram();
    TH2D *h_B1 = (TH2D*)h_migr->GetHistogram();

    // Calo_B2
    h_migr -> SetName("h_B2");
    h_migr -> GenerateHistogram();
    TH2D *h_B2 = (TH2D*)h_migr->GetHistogram();


    cout << "Filling Migration Matrices" << endl;
    Edep -> Draw("Calo_A1:Egen>>h_A1", GoodEvents_photons_PID_1_NoVB.Data(), "goff");
    
    cout << "Filling Migration Matrices" << endl;
    Edep -> Draw("Calo_A2:Egen>>h_A2", GoodEvents_photons_PID_2_NoVB.Data(), "goff");
    
    cout << "Filling Migration Matrices" << endl;
    Edep -> Draw("Calo_B1:Egen>>h_B1", GoodEvents_GBM_1_NoVB.Data(), "goff");
    
    cout << "Filling Migration Matrices" << endl;
    Edep -> Draw("Calo_B2:Egen>>h_B2", GoodEvents_GBM_2_NoVB.Data(), "goff");



    cout << "Normalizing Migration Matrices" << endl;
    Normalize_TH2D(h_A1);
    Normalize_TH2D(h_A2);
    Normalize_TH2D(h_B1);
    Normalize_TH2D(h_B2);


    h_A1 -> SetStats(0);
    h_A2 -> SetStats(0);
    h_B1 -> SetStats(0);
    h_B2 -> SetStats(0);

    h_A1 -> SetTitle("");
    h_A1 -> SetXTitle("E_{gen} [MeV]");
    h_A1 -> SetYTitle("E_{det} Calo_A1 [MeV]");

    h_A2 -> SetTitle("");
    h_A2 -> SetXTitle("E_{gen} [MeV]");
    h_A2 -> SetYTitle("E_{det} Calo_A2 [MeV]");

    h_B1 -> SetTitle("");
    h_B1 -> SetXTitle("E_{gen} [MeV]");
    h_B1 -> SetYTitle("E_{det} Calo_B1 [MeV]");

    h_B2 -> SetTitle("");
    h_B2 -> SetXTitle("E_{gen} [MeV]");
    h_B2 -> SetYTitle("E_{det} Calo_B2 [MeV]");

    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
    c2 -> Divide(2, 2);
    c2 -> cd(1);
    gPad -> SetLogx();
    gPad -> SetLogy();
    h_A1 -> Draw("colz");
    c2 -> cd(2);
    h_B1 -> Draw("colz");
    gPad -> SetLogx();
    gPad -> SetLogy();
    c2 -> cd(3);
    h_B2 -> Draw("colz");
    gPad -> SetLogx();
    gPad -> SetLogy();
    c2 -> cd(4);
    h_A2 -> Draw("colz");
    gPad -> SetLogx();
    gPad -> SetLogy();

    h_A1 -> SaveAs("A1.root");
    h_A2 -> SaveAs("A2.root");
    h_B1 -> SaveAs("B1.root");
    h_B2 -> SaveAs("B2.root");











    double E0 = 10e-3; // MeV
    double E1 = 500.; // MeV
    int N = 100;
    double Lx = 12.;
    double Ly = 12.;

    double Area = Lx * Ly;


    // NO FILTERING
    TH1DLog *h_gamma = new TH1DLog();
    h_gamma -> SetXAxis(E0, E1, N);
    h_gamma -> SetName("h_gamma_0");
    h_gamma -> SetTitle("");
    h_gamma -> SetXTitle("Egen");
    h_gamma -> SetYTitle("Effective Area [cm^{2}]");

    h_gamma -> GenerateHistogram();
    TH1D* h_gamma_0 = (TH1D*)h_gamma->GetHistogram();
    h_gamma_0 -> SetStats(0);


    // FILTERING CONSIDERING THE TOTAL EVENT CONFINEMENT
    h_gamma = new TH1DLog();
    h_gamma -> SetXAxis(E0, E1, N);
    h_gamma -> SetName("h_gamma_1");
    h_gamma -> SetTitle("");
    h_gamma -> SetXTitle("Egen");
    h_gamma -> SetYTitle("Effective Area [cm^{2}]");

    h_gamma -> GenerateHistogram();
    TH1D* h_gamma_1 = (TH1D*)h_gamma->GetHistogram();
    h_gamma_1 -> SetStats(0);

    // FILTERING NOT CONSIDERING THE BOTTOM VETO CONFINEMENT
    h_gamma = new TH1DLog();
    h_gamma -> SetXAxis(E0, E1, N);
    h_gamma -> SetName("h_gamma_2");
    h_gamma -> SetTitle("");
    h_gamma -> SetXTitle("Egen");
    h_gamma -> SetYTitle("Effective Area [cm^{2}]");

    h_gamma -> GenerateHistogram();
    TH1D* h_gamma_2 = (TH1D*)h_gamma->GetHistogram();
    h_gamma_2 -> SetStats(0);


    // FILL THE HISTOGRAMS
    TString Condition = Form("(%s) || (%s) || (%s) || (%s)", GoodEvents_GBM_1.Data(), GoodEvents_GBM_2.Data(), GoodEvents_photons_PID_1.Data(), GoodEvents_photons_PID_2.Data());
    cout << "Condition: " << Condition.Data() << endl;

    cout << "Filling histograms..." << endl;
    Edep -> Draw("Egen>>h_gamma_0", "", "goff");
    cout << "Filling histograms..." << endl;
    Edep -> Draw("Egen>>h_gamma_1", Condition.Data(), "goff");

    Condition = Form("(%s) || (%s) || (%s) || (%s)", GoodEvents_GBM_1_NoVB.Data(), GoodEvents_GBM_2_NoVB.Data(), GoodEvents_photons_PID_1_NoVB.Data(), GoodEvents_photons_PID_2_NoVB.Data());
    cout << "Condition: " << Condition.Data() << endl;

    cout << "Filling histograms..." << endl;
    Edep -> Draw("Egen>>h_gamma_2", Condition.Data(), "goff");


    for(int i = 0; i < h_gamma_0->GetNbinsX(); i++)
    {
        double N_gen = h_gamma_0->GetBinContent(i+1);
        double N_det_1 = h_gamma_1->GetBinContent(i+1);
        double N_det_2 = h_gamma_2->GetBinContent(i+1);

        double sigma_N_gen = sqrt(N_gen);
        double sigma_N_det_1 = sqrt(N_det_1);
        double sigma_N_det_2 = sqrt(N_det_2);

        double eff_area_1 = N_det_1 / N_gen * Area;
        double eff_area_2 = N_det_2 / N_gen * Area;

        double sigma_eff_area_1 = eff_area_1 * sqrt(pow(sigma_N_det_1/N_det_1, 2) + pow(sigma_N_gen/N_gen, 2));
        double sigma_eff_area_2 = eff_area_2 * sqrt(pow(sigma_N_det_2/N_det_2, 2) + pow(sigma_N_gen/N_gen, 2));


        h_gamma_1 -> SetBinContent(i+1, eff_area_1);
        h_gamma_1 -> SetBinError(i+1, sigma_eff_area_1);

        h_gamma_2 -> SetBinContent(i+1, eff_area_2);
        h_gamma_2 -> SetBinError(i+1, sigma_eff_area_2);
    }

    TCanvas *c = new TCanvas("c", "c", 800, 800);

    h_gamma_1 -> SetLineColor(kRed);
    h_gamma_1 -> Draw();
    h_gamma_2 -> SetLineColor(kBlue);
    h_gamma_2 -> Draw("same");

    TLegend *l = new TLegend(0.1, 0.7, 0.4, 0.9);
    l -> AddEntry(h_gamma_1, "w confinement", "l");
    l -> AddEntry(h_gamma_2, "w/o confinement", "l");
    l -> Draw();

    gPad -> SetGrid();
    gPad -> SetLogx();
    gPad -> SetLogy();



    return 0;
}