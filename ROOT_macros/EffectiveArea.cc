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


TH1D* GenerateHistogram(TString Name, TString title, TString XTitle, TString YTitle, 
                        int NbinsX, double xmin, double xmax)
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

    TString Title_tot =  title+";"+XTitle+";"+YTitle;

    TH1D *histogram = new TH1D(Name, Title_tot, NbinsX, &xBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);
    histogram->Sumw2();

    return histogram;
}

void EffectiveArea()
{
    TFile *file = new TFile("/home/franz/Desktop/SPaRKLE/build/gamma_run_0_um.root");
    TTree *tree = (TTree*)file->Get("Edep");
    //histo boundaries
    int Nbins = 50; 
    double Emin = 0.01, Emax = 100;

    TH1D* h_det = GenerateHistogram("h_det", "Effective area", "E [MeV]", "Area [cm^{2}]", Nbins, Emin, Emax);
    h_det->Print();

    TH1D* h_gen = GenerateHistogram("h_gen", "Events with generated events", "E [MeV]", "Entries", Nbins, Emin, Emax);
    h_gen->Print();

    //Thresholds
    double EminTH = 0., EminCaloTH = 0.02, EminGenTH = 0.03;
    // printf("Entries of the tree: %lld\n", tree->GetEntries());
    long Entries = tree->GetEntries();

    //Generation area in cm^2
    double GenArea = 64.; // cm^2

    //setting branches
    double RandEnergy = 0, Ed_DrilledVeto = 0, Ed_BottomVeto = 0, Thin_x0_y0_ID0 = 0, Thin_x0_y1_ID1 = 0, Ed_Calo0 = 0;
    tree->SetBranchAddress("RandEnergy", &RandEnergy);
    tree->SetBranchAddress("Ed_DrilledVeto", &Ed_DrilledVeto);
    tree->SetBranchAddress("Ed_BottomVeto", &Ed_BottomVeto);
    tree->SetBranchAddress("Thin_x0_y0_ID0", &Thin_x0_y0_ID0);
    tree->SetBranchAddress("Thin_x0_y1_ID1", &Thin_x0_y1_ID1);
    tree->SetBranchAddress("Ed_Calo0", &Ed_Calo0);

    // event loop
    for(int i=0; i<Entries; i++){
        tree->GetEntry(i);
        if(RandEnergy<EminGenTH) continue;
        h_gen->Fill(RandEnergy);

        if( (Ed_DrilledVeto==0.) && (Ed_BottomVeto==0.) && (Thin_x0_y0_ID0==0.) &&  (Thin_x0_y1_ID1==0.) && (Ed_Calo0>EminCaloTH) ) h_det->Fill(RandEnergy);
    }

    printf("Entries of the histos: %f, %f\n", h_det->GetEntries(), h_gen->GetEntries());

    h_det->Divide(h_gen);
    h_det->Scale(GenArea);
    h_det->SetLineColor(kRed);
    h_det->SetLineWidth(3);
    h_det->SetTitleOffset(1.2, "X");
    h_det->SetTitleOffset(0.9, "Y");
    gStyle->SetOptStat(0);
    
    // style->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
    // c1->GetPad(0)->SetTicks(2,2);
    c1->SetLogx();
    c1->SetLogy();
    c1->SetGrid();
    h_det->Draw();

    TFile *outFile = new TFile("../SPaRKLE_EffectiveArea.root", "RECREATE");
    outFile->cd();
    h_det->Write();
    h_gen->Write();
    c1->Write();
    outFile->Close();

    return 0;
}
