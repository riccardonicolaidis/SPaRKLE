#include <TFile.h>
#include <TTree.h>
#include <Riostream.h>
#include <TMath.h>
#include <map>
#include <TH1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <sstream>
#include <TString.h>
#include <iostream>
using namespace std;


// TH1D* GenerateHistogram(TString Name, TString title, TString XTitle, TString YTitle, 
//                         int NbinsX, double xmin, double xmax)
// {
//     std::vector<double> xBins;
//     std::vector<double> xBinsLog;

//     int NEdgesX = NbinsX + 1;
//     double xMaxLog = TMath::Log10(xmax);
//     double xMinLog = TMath::Log10(xmin);

//     double xBinWidthLog = (xMaxLog - xMinLog) / NbinsX;
    
//     for (int i = 0; i < NEdgesX; i++)
//     {
//         xBinsLog.push_back(xMinLog + i * xBinWidthLog);
//         xBins.push_back(TMath::Power(10, xBinsLog.at(i)));
//     }

//     TString Title_tot =  title+";"+XTitle+";"+YTitle;

//     TH1D *histogram = new TH1D(Name, Title_tot, NbinsX, &xBins[0]);
//     histogram->GetXaxis()->SetTitle(XTitle);
//     histogram->GetYaxis()->SetTitle(YTitle);
//     histogram->Sumw2();

//     return histogram;
// }

TH2D* GenerateTH2(TString Name, TString title, TString XTitle, TString YTitle, 
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


void PID()
{
    gStyle->SetOptStat(0);
    const int Nfiles = 4;
    TFile *files[Nfiles];
    TString file_names[Nfiles] = {"/home/franz/Desktop/SPaRKLE/build/proton_run_0_um.root", 
                                  "/home/franz/Desktop/SPaRKLE/build/He4_run_0_um.root",
                                  "/home/franz/Desktop/SPaRKLE/build/Li7_run_0_um.root",
                                  "/home/franz/Desktop/SPaRKLE/build/electron_run_0_um.root"};//,
                                //   "/home/franz/Desktop/SPaRKLE/build/He3_run_0_um.root"};

    
    TString histo_names[2*Nfiles] = {"p_thinID0", "p_thinID1",
                                     "He4_thinID0", "He4_thinID1",
                                     "Li7_thinID0", "Li7_thinID1",
                                     "e_thinID0", "e_thinID1"};//,
                                    //  "He3_thinID0", "He3_thinID1"};
    // TString histo_titles[2*Nfiles];
    TH2D *histos[2*Nfiles];
    const int NBinsX = 400, NBinsY = 400;
    const double ETOTmin = 0.05, ETOTmax = 3500; //MeV
    const double E1min = 0.005, E1max = 30; //MeV
    for(int i=0; i<2*Nfiles; i++){
        // histo_titles[i] = histo_names[i] + "; E_{tot} [MeV]; E_{1} [MeV]";
        histos[i] = GenerateTH2(histo_names[i], histo_names[i], "E_{tot} [MeV]", "E_{1} [MeV]", NBinsX, ETOTmin, ETOTmax, NBinsY, E1min, E1max, true);
        //  TH1D(histo_names[i], histo_titles[i], NBins, Emin, Emax);
    }
    TH2D *hComulative = GenerateTH2("hComulative", "Particle Identification", "E_{tot} [MeV]", "E_{1} [MeV]", NBinsX, ETOTmin, ETOTmax, NBinsY, E1min, E1max, true);

    

    TTree *tree;
    double RandEnergy = 0, Ed_DrilledVeto = 0, Ed_BottomVeto = 0, Thin_x0_y0_ID0 = 0, Thin_x0_y1_ID1 = 0, Ed_Calo0 = 0;

    for(int i=0; i<Nfiles; i++){
        files[i] = new TFile(file_names[i]);
        tree = (TTree*)files[i]->Get("Edep");

        //setting branches
        tree->SetBranchAddress("RandEnergy", &RandEnergy);
        tree->SetBranchAddress("Ed_DrilledVeto", &Ed_DrilledVeto);
        tree->SetBranchAddress("Ed_BottomVeto", &Ed_BottomVeto);
        tree->SetBranchAddress("Thin_x0_y0_ID0", &Thin_x0_y0_ID0);
        tree->SetBranchAddress("Thin_x0_y1_ID1", &Thin_x0_y1_ID1);
        tree->SetBranchAddress("Ed_Calo0", &Ed_Calo0);

        // event loop
        for(int iEntry=0; iEntry<tree->GetEntries(); iEntry++){
            tree->GetEntry(iEntry);
            if(Ed_DrilledVeto != 0.) continue;
            if(Ed_BottomVeto != 0.) continue;

            if(file_names[i].Contains("Li7")){

                if(Thin_x0_y1_ID1 == 0. && Ed_Calo0 > 0.08 && !( (Thin_x0_y0_ID0<16) && (Thin_x0_y0_ID0+Ed_Calo0)<32)) histos[i]->Fill(Thin_x0_y0_ID0+Ed_Calo0, Thin_x0_y0_ID0);
                if(Thin_x0_y0_ID0 == 0. && Ed_Calo0 > 0.08 && !( (Thin_x0_y1_ID1<16) && (Thin_x0_y1_ID1+Ed_Calo0)<32)) histos[i+1]->Fill(Thin_x0_y1_ID1+Ed_Calo0, Thin_x0_y1_ID1);

            }
            else{

                if(Thin_x0_y1_ID1 == 0. && Ed_Calo0 > 0.08) histos[i]->Fill(Thin_x0_y0_ID0+Ed_Calo0, Thin_x0_y0_ID0);
                if(Thin_x0_y0_ID0 == 0. && Ed_Calo0 > 0.08) histos[i+1]->Fill(Thin_x0_y1_ID1+Ed_Calo0, Thin_x0_y1_ID1);

            }
        }

        hComulative->Add(histos[i]);
        hComulative->Add(histos[i+1]);
    }

    // style->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    // c1->GetPad(0)->SetTicks(2,2);
    c1->SetLogx();
    c1->SetLogy();
    c1->SetLogz();
    c1->SetGrid();
    hComulative->SetTitleOffset(1.3, "X");
    hComulative->SetTitleOffset(1.2, "Y");
    hComulative->Draw("colz");

    TFile *outFile = new TFile("../SPaRKLE_PID.root", "RECREATE");
    outFile->cd();
    for(int i=0; i<2*Nfiles; i++) histos[i]->Write();
    hComulative->Write();
    c1->Write();
    outFile->Close();


    return 0;

}