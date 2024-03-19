#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <map>
#include <TH2.h>
#include <TLine.h>
#include <TCanvas.h>
#include <sstream>
#include <TString.h>
#include <TLatex.h>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
using namespace std;

TH2D* GenerateHistogram2(TString Name, TString title, TString XTitle, TString YTitle, 
                        int NbinsX, double xmin, double xmax, int NbinsY, double ymin, double ymax, bool ylog)
{
    std::vector<double> xBins;
    std::vector<double> xBinsLog;

    int NEdgesX = NbinsX + 1;
    double xMaxLog = TMath::Log10(xmax);
    double xMinLog = TMath::Log10(xmin);

    double xBinWidthLog = (xMaxLog - xMinLog) / NbinsX;
    
    for (int i = 0; i < NEdgesX; i++)
    {
        xBinsLog.push_back(xMinLog + i * xBinWidthLog);
        xBins.push_back(TMath::Power(10, xBinsLog.at(i)));
    }

    std::vector<double> yBins;
    std::vector<double> yBinsLog;

    int NEdgesY = NbinsY + 1;
    double yMaxLog = TMath::Log10(ymax);
    double yMinLog = TMath::Log10(ymin);

    double yBinWidthLog = (yMaxLog - yMinLog) / NbinsY;
    double yBinWidth = (ymax - ymin) / NbinsY;

    for (int i = 0; i < NEdgesY; i++)
    {
        yBinsLog.push_back(yMinLog + i * yBinWidthLog);
        if(ylog) yBins.push_back(TMath::Power(10, yBinsLog.at(i)));
        else yBins.push_back(ymin + i * yBinWidth);
    }

    TString Title_tot =  title+";"+XTitle+";"+YTitle;

    TH2D *histogram = new TH2D(Name, Title_tot, NbinsX, &xBins[0], NbinsY, &yBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);

    histogram->Sumw2();

    return histogram;
}



void Emean(){

    gStyle->SetOptStat(0); // no statistics

    //Files to be used
    const int NFiles = 2;
    TString particles[NFiles] = {"p", "e^{-}"};
    TString names[NFiles] = {"../build/GPSproton_run_0_um.root", "../build/GPSelectron_run_0_um.root"};
    TString hNames[NFiles] = {"p", "e^{-}"};


    TH2D *h2_GAGGvsMCEnergy[NFiles];    // Energy deposition in GAGG vs Generation Energy
    TH2D *h2_Si1vsMCEnergy[NFiles];     // Energy deposition in Si1 vs Generation Energy
    TH2D *h2_Si2vsMCEnergy[NFiles];     // Energy deposition in Si2 vs Generation Energy
    TH2D *h2_DrilleVvsMCEnergy[NFiles]; // Energy deposition in drilled veto vs Generation Energy
    TH2D *h2_BottomVvsMCEnergy[NFiles]; // Energy deposition in bottom veto vs Generation Energy

    const int nbins = 400;
    const double Ek_min_p = 1, Ek_max_p = 500; // MeV
    const double Ek_min_e = 0.05, Ek_max_e = 500; // MeV

    // --- initializing histograms --- //
    for(int iFile=0; iFile<NFiles; iFile++){
        if(names[iFile].Contains("proton")) {
            h2_GAGGvsMCEnergy[iFile] = GenerateHistogram2("GAGG_energy_deposition_vs_" + hNames[iFile] + "_energy", "GAGG_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_p, Ek_max_p, nbins, 0.01, Ek_max_p, true);
            h2_Si1vsMCEnergy[iFile] = GenerateHistogram2("Si1_energy_deposition_vs_" + hNames[iFile] + "_energy", "Si1_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_p, Ek_max_p, nbins, 0.01, 40.1, true);
            h2_Si2vsMCEnergy[iFile] = GenerateHistogram2("Si2_energy_deposition_vs_" + hNames[iFile] + "_energy", "Si2_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_p, Ek_max_p, nbins, 0.01, 40.1, true);
            h2_DrilleVvsMCEnergy[iFile] = GenerateHistogram2("DrilledVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "DrilledVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_p, Ek_max_p, nbins, 0.01, Ek_max_p, true);
            h2_BottomVvsMCEnergy[iFile] = GenerateHistogram2("BottomVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "BottomVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_p, Ek_max_p, nbins, 0.01, Ek_max_p, true);
        }
        if(names[iFile].Contains("electron")){
            h2_GAGGvsMCEnergy[iFile] = GenerateHistogram2("GAGG_energy_deposition_vs_" + hNames[iFile] + "_energy", "GAGG_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_e, Ek_max_e, nbins, 0.01, Ek_max_e, true);
            h2_Si1vsMCEnergy[iFile] = GenerateHistogram2("Si1_energy_deposition_vs_" + hNames[iFile] + "_energy", "Si1_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_e, Ek_max_e, nbins, 0.01, 40.1, true);
            h2_Si2vsMCEnergy[iFile] = GenerateHistogram2("Si2_energy_deposition_vs_" + hNames[iFile] + "_energy", "Si2_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_e, Ek_max_e, nbins, 0.01, 40.1, true);
            h2_DrilleVvsMCEnergy[iFile] = GenerateHistogram2("DrilledVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "DrilledVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_e, Ek_max_e, nbins, 0.01, Ek_max_e, true);
            h2_BottomVvsMCEnergy[iFile] = GenerateHistogram2("BottomVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "BottomVeto_energy_deposition_vs_" + hNames[iFile] + "_energy", "E_{k} [MeV]", "E_{dep} [MeV]", nbins, Ek_min_e, Ek_max_e, nbins, 0.01, Ek_max_e, true);
        } 
    }

    //Trigger thresholds and cuts
    const double th_Calo0 = 0.08; // Mev = 80 kev
    const double th_Silicon = 0.;
    const double th_DrilledVeto = 0.;
    const double th_BottomVeto = 0.;


    //Tree variables
    double RandEnergy, Ed_Calo0, Ed_DrilledVeto, Ed_BottomVeto, Thin_x0_y0_ID0, Thin_x0_y1_ID1;

    //Loop on the files
    for(int iFile=0; iFile<NFiles; iFile++){

        TFile *file = new TFile(names[iFile], "OPEN");
        printf("\n Reading file %d/%d: %s\n", iFile+1, NFiles, names[iFile].Data());

        TTree *t = (TTree *)file->Get("Edep");

        //setting the branches
        t->SetBranchAddress("RandEnergy", &RandEnergy);
        t->SetBranchAddress("Ed_Calo0", &Ed_Calo0);
        t->SetBranchAddress("Ed_DrilledVeto", &Ed_DrilledVeto);
        t->SetBranchAddress("Ed_BottomVeto", &Ed_BottomVeto);
        t->SetBranchAddress("Thin_x0_y0_ID0", &Thin_x0_y0_ID0);
        t->SetBranchAddress("Thin_x0_y1_ID1", &Thin_x0_y1_ID1);



        // Loop on the events
        for(int entry=0; entry<t->GetEntries(); entry++){

            t->GetEntry(entry);

            if(Ed_DrilledVeto>0.) h2_DrilleVvsMCEnergy[iFile]->Fill(RandEnergy, Ed_DrilledVeto);
            if(Ed_BottomVeto>0.) h2_BottomVvsMCEnergy[iFile]->Fill(RandEnergy, Ed_BottomVeto);


            // Charged particle selection
            if(Ed_DrilledVeto == th_DrilledVeto &&
               Ed_BottomVeto == th_BottomVeto && 
               Thin_x0_y1_ID1 == th_Silicon &&
               Ed_Calo0 > th_Calo0 &&
               Thin_x0_y0_ID0 > th_Silicon) {

                h2_GAGGvsMCEnergy[iFile]->Fill(RandEnergy, Ed_Calo0);
                h2_Si1vsMCEnergy[iFile]->Fill(RandEnergy, Thin_x0_y0_ID0);
               }
            
            if(Ed_DrilledVeto == th_DrilledVeto &&
               Ed_BottomVeto == th_BottomVeto && 
               Thin_x0_y0_ID0 == th_Silicon &&
               Ed_Calo0 > th_Calo0 &&
               Thin_x0_y1_ID1 > th_Silicon) {

                h2_GAGGvsMCEnergy[iFile]->Fill(RandEnergy, Ed_Calo0);
                h2_Si2vsMCEnergy[iFile]->Fill(RandEnergy, Thin_x0_y1_ID1);
               }
        
        } // end of the event loop

    }

    TFile *f_out = new TFile("E_mean_deposition.root", "RECREATE");
    f_out->cd();
    for(int iFile=0; iFile<NFiles; iFile++){
        h2_GAGGvsMCEnergy[iFile]->Write();
        h2_Si1vsMCEnergy[iFile]->Write();
        h2_Si2vsMCEnergy[iFile]->Write();
        h2_DrilleVvsMCEnergy[iFile]->Write();
        h2_BottomVvsMCEnergy[iFile]->Write();
    }
    f_out->Close();



    return;
}