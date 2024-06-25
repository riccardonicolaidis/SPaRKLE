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

using namespace std;


int Charged_Particles()
{

    // Filenames
    // Insert the paths of the files to be analyzed
    vector<TString> filenames;
    filenames.push_back("/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST_Acc/e.root");
    filenames.push_back("/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST_Acc/p.root");

    // Filenames filtered
    // insert the paths of the filtered files
    vector<TString> filenames_filtered;
    filenames_filtered.push_back("/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST_Acc/e_filt.root");
    filenames_filtered.push_back("/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST_Acc/p_filt.root");


    vector<TFile*> files;
    vector<TTree*> trees;

    vector<TFile*> files_filtered;
    vector<TTree*> trees_filtered;


    // Set the resolution of the detector
    double GAGG_Res = 0.1; // %
    double Plastic_Res = 0.3; // %
    double SD_Res = 6e-3; //MeV


    double Threshold_GAGG = 0.01; // MeV
    double Threshold_Plastic = 0.05; // MeV
    double Threshold_SD = 0.04; // MeV



    for (int i = 0; i < filenames.size(); i++)
    {
        files.push_back(new TFile(filenames[i]));
        trees.push_back((TTree*)files[i]->Get("Edep"));
    }



    for (int i = 0; i < trees.size(); i++)
    {


        trees[i]->SetAlias("gCalo_A1",Form("(Calo_A1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_A1*%g)", GAGG_Res));
        trees[i]->SetAlias("gCalo_A2", Form("(Calo_A2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_A2*%g)", GAGG_Res));
        trees[i]->SetAlias("gCalo_B1", Form("(Calo_B1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_B1*%g)", GAGG_Res));
        trees[i]->SetAlias("gCalo_B2", Form("(Calo_B2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*Calo_B2*%g)", GAGG_Res));

        trees[i]->SetAlias("gSD1", Form("(SD1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g)", SD_Res));
        trees[i]->SetAlias("gSD2", Form("(SD2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g)", SD_Res));

        trees[i]->SetAlias("gVT", Form("(VT + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VT)",Plastic_Res));
        trees[i]->SetAlias("gVB", Form("(VB + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VB)", Plastic_Res));
        trees[i]->SetAlias("gVL0", Form("(VL0 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL0)", Plastic_Res));
        trees[i]->SetAlias("gVL1", Form("(VL1 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL1)", Plastic_Res));
        trees[i]->SetAlias("gVL2", Form("(VL2 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL2)", Plastic_Res));
        trees[i]->SetAlias("gVL3", Form("(VL3 + (TMath::Sin(2 *pi*rndm)*sqrt(-2*log(rndm)))*%g*VL3)", Plastic_Res));
    }


    


    TString GoodEvents_gPID_1 = Form("(gCalo_A1 > %g) && (gSD1 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_PID_1 = Form("(Calo_A1 > %g) && (SD1 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);

    TString GoodEvents_gPID_2 = Form("(gCalo_A2 > %g) && (gSD1 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_PID_2 = Form("(Calo_A2 > %g) && (SD1 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_SD, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);

    TString GoodEvents_gGBM_1 = Form("(gCalo_B1 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_GBM_1 = Form("(Calo_B1 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);

    TString GoodEvents_gGBM_2 = Form("(gCalo_B2 > %g) && (gVT < %g) && (gVB < %g) && (gVL0 < %g) && (gVL1< %g) && (gVL2 < %g) && (gVL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);
    TString GoodEvents_GBM_2 = Form("(Calo_B2 > %g) && (VT < %g) && (VB < %g) && (VL0 < %g) && (VL1< %g) && (VL2 < %g) && (VL3< %g) ", Threshold_GAGG, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic, Threshold_Plastic);





    // 0.1 a 100 MeV

    // ELECTRONS

    double E0 = 0.1;
    double E1 = 200.;
    int N = 50;
    // No filterning
    TH1DLog *h_e = new TH1DLog();
    h_e -> SetXAxis(E0, E1, N);
    h_e -> SetName("h_e_0");
    h_e -> SetTitle("Acceptance");
    h_e -> SetXTitle("Egen");
    h_e -> SetYTitle("Counts");

    h_e -> GenerateHistogram();
    TH1D* h_e_0 = (TH1D*)h_e->GetHistogram();
    h_e_0 -> SetStats(0);

    // Filtering
    h_e = new TH1DLog();
    h_e -> SetXAxis(E0, E1, N);
    h_e -> SetName("h_e_1");
    h_e -> SetTitle("Acceptance");
    h_e -> SetXTitle("Egen");
    h_e -> SetYTitle("Counts");

    h_e -> GenerateHistogram();
    TH1D* h_e_1 = (TH1D*)h_e->GetHistogram();
    h_e_1 -> SetStats(0);

    // PROTONS

    // No filterning
    TH1DLog *h_p = new TH1DLog();
    h_p -> SetXAxis(E0, E1, N);
    h_p -> SetName("h_p_0");
    h_p -> SetTitle("Acceptance");
    h_p -> SetXTitle("Egen");
    h_p -> SetYTitle("Counts");

    h_p -> GenerateHistogram();
    TH1D* h_p_0 = (TH1D*)h_p->GetHistogram();
    h_p_0 -> SetStats(0);

    // Filtering
    h_p = new TH1DLog();
    h_p -> SetXAxis(E0, E1, N);
    h_p -> SetName("h_p_1");
    h_p -> SetTitle("Acceptance");
    h_p -> SetXTitle("Egen");
    h_p -> SetYTitle("Counts");
    h_p -> GenerateHistogram();
    TH1D* h_p_1 = (TH1D*)h_p->GetHistogram();
    h_p_1 -> SetStats(0);


    // Electrons
    trees[0] -> Draw("Egen >> h_e_0", "", "goff"); // No Filtering
    trees[0] -> Draw("Egen >> h_e_1", Form("(%s) || (%s)", GoodEvents_PID_1.Data(), GoodEvents_PID_2.Data()), "goff");

    // Protons
    trees[1] -> Draw("Egen >> h_p_0", "", "goff");
    trees[1] -> Draw("Egen >> h_p_1", Form("(%s) || (%s)", GoodEvents_PID_1.Data(), GoodEvents_PID_2.Data()), "goff");

    // Clone the histograms
    TH1D* h_e_2 = (TH1D*)h_e_1->Clone("h_e_2");
    TH1D* h_p_2 = (TH1D*)h_p_1->Clone("h_p_2");


    double Lx = 10.;
    double Ly = 10.;
    double Area = Lx*Ly;

    double R = 4.;
    //Area = TMath::Pi()*R*R;

    double G_gen = Area*TMath::Pi();


    for(int i = 0; i < h_e_0->GetNbinsX(); ++i)
    {
        double e_gen = h_e_0->GetBinContent(i+1);
        double sigma_e_gen = h_e_0->GetBinError(i+1);

        double e_det = h_e_1->GetBinContent(i+1);
        double sigma_e_det = h_e_1->GetBinError(i+1);

        double p_gen = h_p_0->GetBinContent(i+1);
        double sigma_p_gen = h_p_0->GetBinError(i+1);

        double p_det = h_p_1->GetBinContent(i+1);
        double sigma_p_det = h_p_1->GetBinError(i+1);

        double e_acc = e_det/e_gen*G_gen;
        double sigma_e_acc = sqrt(pow(sigma_e_det/e_det,2.) + pow(sigma_e_gen/e_gen,2.))*e_acc;

        if(e_det == 0. || e_gen == 0.)
        {
            e_acc = 0.;
            sigma_e_acc = 0.;
        }


        h_e_2->SetBinContent(i+1, e_acc);
        h_e_2->SetBinError(i+1, sigma_e_acc);

        double p_acc = p_det/p_gen*G_gen;
        double sigma_p_acc = sqrt(pow(sigma_p_det/p_det, 2.) + pow(sigma_p_gen/p_gen, 2.))*p_acc;

        if(p_det == 0. || p_gen == 0.)
        {
            p_acc = 0.;
            sigma_p_acc = 0.;
        }

        h_p_2->SetBinContent(i+1, p_acc);
        h_p_2->SetBinError(i+1, sigma_p_acc);

    }


    TCanvas* c1 = new TCanvas("c1", "Edep", 800, 600);
    h_e_2 -> Draw();
    h_e_2 -> SetLineColor(kBlue);
    h_e_2 -> SetLineWidth(2);

    h_p_2 -> Draw("same");
    h_p_2 -> SetLineColor(kRed);
    h_p_2 -> SetLineWidth(2);

    TLegend* legend1 = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend1->AddEntry(h_e_2, "Electrons", "l");
    legend1->AddEntry(h_p_2, "Protons", "l");
    legend1->Draw();

    gPad -> SetLogx();
    gPad -> SetGrid();
    gPad -> SetLogy();

    c1 -> SaveAs("Edep.png");
    c1 -> SaveAs("Edep.root");

    h_e_2 -> SaveAs("Edep_Electrons.root");
    h_p_2 -> SaveAs("Edep_Protons.root");


    double Calo_A1, Calo_A2, Calo_B1, Calo_B2;
    double SD1, SD2;
    double VT, VB;
    double VL0, VL1, VL2, VL3;
    double Egen;
    double Xgen, Ygen, Zgen;
    double pDirX, pDirY, pDirZ;

    Calo_A1 = 0.;
    Calo_A2 = 0.;
    Calo_B1 = 0.;
    Calo_B2 = 0.;
    SD1 = 0.;
    SD2 = 0.;
    VT = 0.;
    VB = 0.;
    VL0 = 0.;
    VL1 = 0.;
    VL2 = 0.;
    VL3 = 0.;
    Egen = 0.;
    Xgen = 0.;
    Ygen = 0.;
    Zgen = 0.;
    pDirX = 0.;
    pDirY = 0.;
    pDirZ = 0.;


    for(int i = 0; i < files.size(); ++i)
    {
        files_filtered.push_back(TFile::Open(Form("%s_filtered.root", files[i]->GetName()), "RECREATE"));
        trees_filtered.push_back(new TTree("Edep", "Edep"));

        trees[i]->SetBranchAddress("Calo_A1", &Calo_A1);
        trees[i]->SetBranchAddress("Calo_A2", &Calo_A2);
        trees[i]->SetBranchAddress("Calo_B1", &Calo_B1);
        trees[i]->SetBranchAddress("Calo_B2", &Calo_B2);
        trees[i]->SetBranchAddress("SD1", &SD1);
        trees[i]->SetBranchAddress("SD2", &SD2);
        trees[i]->SetBranchAddress("VT", &VT);
        trees[i]->SetBranchAddress("VB", &VB);
        trees[i]->SetBranchAddress("VL0", &VL0);
        trees[i]->SetBranchAddress("VL1", &VL1);
        trees[i]->SetBranchAddress("VL2", &VL2);
        trees[i]->SetBranchAddress("VL3", &VL3);
        trees[i]->SetBranchAddress("Egen", &Egen);
        trees[i]->SetBranchAddress("Xgen", &Xgen);
        trees[i]->SetBranchAddress("Ygen", &Ygen);
        trees[i]->SetBranchAddress("Zgen", &Zgen);
        trees[i]->SetBranchAddress("pDirX", &pDirX);
        trees[i]->SetBranchAddress("pDirY", &pDirY);
        trees[i]->SetBranchAddress("pDirZ", &pDirZ);

        trees_filtered[i]->Branch("Calo_A1", &Calo_A1, "Calo_A1/D");
        trees_filtered[i]->Branch("Calo_A2", &Calo_A2, "Calo_A2/D");
        trees_filtered[i]->Branch("Calo_B1", &Calo_B1, "Calo_B1/D");
        trees_filtered[i]->Branch("Calo_B2", &Calo_B2, "Calo_B2/D");
        trees_filtered[i]->Branch("SD1", &SD1, "SD1/D");
        trees_filtered[i]->Branch("SD2", &SD2, "SD2/D");
        trees_filtered[i]->Branch("VT", &VT, "VT/D");
        trees_filtered[i]->Branch("VB", &VB, "VB/D");
        trees_filtered[i]->Branch("VL0", &VL0, "VL0/D");
        trees_filtered[i]->Branch("VL1", &VL1, "VL1/D");
        trees_filtered[i]->Branch("VL2", &VL2, "VL2/D");
        trees_filtered[i]->Branch("VL3", &VL3, "VL3/D");
        trees_filtered[i]->Branch("Egen", &Egen, "Egen/D");
        trees_filtered[i]->Branch("Xgen", &Xgen, "Xgen/D");
        trees_filtered[i]->Branch("Ygen", &Ygen, "Ygen/D");
        trees_filtered[i]->Branch("Zgen", &Zgen, "Zgen/D");
        trees_filtered[i]->Branch("pDirX", &pDirX, "pDirX/D");
        trees_filtered[i]->Branch("pDirY", &pDirY, "pDirY/D");
        trees_filtered[i]->Branch("pDirZ", &pDirZ, "pDirZ/D");

        cout << "Processing file: " << files[i]->GetName() << endl;

        for(int j = 0; j < trees[i]->GetEntries(); ++j)
        {
            trees[i]->GetEntry(j);

            bool Good_PID1 = (SD1 > Threshold_SD) && (Calo_A1 > Threshold_GAGG);
            bool Good_PID2 = (SD2 > Threshold_SD) && (Calo_A2 > Threshold_GAGG);
            bool Veto_OFF = (VT < Threshold_Plastic) && (VB < Threshold_Plastic) && (VL0 < Threshold_Plastic) && (VL1 < Threshold_Plastic) && (VL2 < Threshold_Plastic) && (VL3 < Threshold_Plastic);
            bool Good_GBM1 = (Calo_B1 > Threshold_GAGG);
            bool Good_GBM2 = (Calo_B2 > Threshold_GAGG);

            bool Good_Event = (Good_PID1 || Good_PID2) && Veto_OFF && (Good_GBM1 || Good_GBM2);

            if(Good_Event)
            {
                trees_filtered[i]->Fill();
            }

            if(j%10000 == 0)
            {
                cout << "Processing event: " << j << endl;
            }
        }

        trees_filtered[i]->Write();







    }






    return 0;
}