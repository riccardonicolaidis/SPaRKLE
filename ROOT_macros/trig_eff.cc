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








TH1D* GenerateHistogram(TString Name, TString title, TString XTitle, TString YTitle, 
                        int NbinsX, double xmin, double xmax)
{
    std::vector<double> xBins;
    std::vector<double> xBinsLog;

    int NEdgesX = NbinsX + 1;
    double xMaxLog = TMath::Log10(xmax);
    double xMinLog = TMath::Log10(xmin);

    double xBinWidthLog_p = (xMaxLog - xMinLog) / NbinsX;
    
    for (int i = 0; i < NEdgesX; i++)
    {
        xBinsLog.push_back(xMinLog + i * xBinWidthLog_p);
        xBins.push_back(TMath::Power(10, xBinsLog.at(i)));
    }

    TString Title_tot =  title+";"+XTitle+";"+YTitle;

    TH1D *histogram = new TH1D(Name, Title_tot, NbinsX, &xBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);
    histogram->Sumw2();

    return histogram;
}






void trig_eff(){
    gStyle->SetOptStat(0); // no statistics

    //Files to be used
    const int NFiles = 2;
    TString particles[NFiles] = {"p", "e^{-}"};
    TString names[NFiles] = {"../build/GPSproton_run_0_um.root", "../build/GPSelectron_run_0_um.root"};
    TString h1_names[NFiles] = {"Gen_spectra_P", "Gen_spectra_e"};


    //Histos for efficiencies
    TH1D *h1_eff[NFiles];
    const int LineWidth = 3;
    Color_t colors[NFiles] = {kBlue, kRed};
    TPaveText* title;
    const int nbins = 100;
    const double Ek_min_p = 1, Ek_max_p = 500; // MeV
    const double Ek_min_e = 0.05, Ek_max_e = 500; // MeV
    const double BinWidthLog_p = (TMath::Log10(Ek_max_p) - TMath::Log10(Ek_min_p)) / nbins; // bin width in log scale
    const double BinWidthLog_e = (TMath::Log10(Ek_max_e) - TMath::Log10(Ek_min_e)) / nbins; // bin width in log scale

    //Generation events info and histos
    const double GenArea = 6400; // m^2 = 64 cm^2 probabilmente troppo poco per dire che il flusso e' davvero isotropo ma per il momento AMEN
    const double PlaneAcceptance = TMath::Pi() * GenArea; // m^2 sr
    TH1D *h1_genEvents[NFiles];

    //Intializing the histograms

    for(int i=0; i<NFiles; i++) {
        if(names[i].Contains("proton")) {
            printf("Creating histogram for protons\n");
            h1_eff[i] = GenerateHistogram("h1_eff_"+particles[i], "Trigger_efficiency_"+particles[i], "E_{k} [MeV/c]", "#epsilon", nbins, Ek_min_p, Ek_max_p);
        }
        if(names[i].Contains("electron")) {
            printf("Creating histogram for electrons\n");
            h1_eff[i] = GenerateHistogram("h1_eff_"+particles[i], "Trigger_efficiency_"+particles[i], "E_{k} [MeV/c]", "#epsilon", nbins, Ek_min_e, Ek_max_e);
        }
        
        h1_eff[i]->SetLineWidth(LineWidth);
        h1_eff[i]->SetLineColor(colors[i]);
        h1_eff[i]->GetYaxis()->SetTitleOffset(1.0);
        h1_eff[i]->GetXaxis()->SetTitleOffset(1.2);

    }


    //Trigger thresholds and cuts
    const double th_Calo0 = 0.08; // Mev = 80 kev
    const double th_Silicon = 0.;
    const double th_DrilledVeto = 0.;
    const double th_BottomVeto = 0.;

    const int NCut = 5;
    TString charged_Sel[NCut] = {"Ed_DrilledVeto", "Ed_BottomVeto==0", "Thin_x0_y1_ID1", "Ed_Calo0", "Thin_x0_y0_ID0"};

     // printf("\n ------- CURRENT CHARGED PARTICLE SELECTION -------\n");
    // for(int j=0;j<NCut;j++){
    //     switch(j){
    //         case 0: case 1:
    //             printf(" %s \t > %.2f \n", charged_Sel[j].Data(), th_Calo0);
    //             break;
    //         case 2:
    //             printf(" %s \t > %.2f \n", charged_Sel[j].Data(), th_DrilledVeto);
    //             break;
    //         case 3: case 4:
    //             printf(" %s \t < %.2f \n", charged_Sel[j].Data(), th_BottomVeto);
    //             printf(" %s \t > %.2f \n", charged_Sel[j].Data(), 0.0);
    //             if(j==4) {
    //             }
    //             break;
    //     }        
    // }



    //Tree variables
    double RandEnergy, Ed_Calo0, Ed_DrilledVeto, Ed_BottomVeto, Thin_x0_y0_ID0, Thin_x0_y1_ID1;


    //Loop on the files
    for(int iFile=0; iFile<NFiles; iFile++){

        TFile *file = new TFile(names[iFile], "OPEN");
        printf("\n Reading file %d/%d: %s\n", iFile+1, NFiles, names[iFile].Data());

        TTree *t = (TTree *)file->Get("Edep");
        h1_genEvents[iFile] = (TH1D *)file->Get(h1_names[iFile]);
        h1_genEvents[iFile]->SetLineColor(colors[iFile]);


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


            // Charged particle selection
            if(Ed_DrilledVeto == th_DrilledVeto &&
               Ed_BottomVeto == th_BottomVeto && 
               Thin_x0_y1_ID1 == th_Silicon &&
               Ed_Calo0 > th_Calo0 &&
               Thin_x0_y0_ID0 > th_Silicon) h1_eff[iFile]->Fill(RandEnergy);
            
            if(Ed_DrilledVeto == th_DrilledVeto &&
               Ed_BottomVeto == th_BottomVeto && 
               Thin_x0_y0_ID0 == th_Silicon &&
               Ed_Calo0 > th_Calo0 &&
               Thin_x0_y1_ID1 > th_Silicon) h1_eff[iFile]->Fill(RandEnergy);
        
        } // end of the event loop

        //Getting the efficiencies
        h1_eff[iFile]->Divide(h1_genEvents[iFile]);

    } // end of the file loop



    // Test to visualize the generated events histo
    // -----------------------------------------------------------------
    TH1D *h1_dummy = new TH1D("h1_dummy", "Generated events in each bin for protons and electrons; E_{k} [MeV/c]; Entries", 100, Ek_min_e, Ek_max_p); //dummy histo to set the axis
    TCanvas *c_test = new TCanvas("c_test", "c_test", 800, 600);
    c_test->SetLogx();
    double max = 0;
    if(h1_genEvents[0]->GetMaximum() > h1_genEvents[1]->GetMaximum()) max = h1_genEvents[0]->GetMaximum()*1.1;
    else max = h1_genEvents[1]->GetMaximum() * 1.1;
    h1_dummy->GetYaxis()->SetRangeUser(1,max);
    h1_dummy->DrawCopy("hist");
    h1_genEvents[0]->Draw("e same");
    h1_genEvents[1]->Draw("e same");
    // -----------------------------------------------------------------



   

    TCanvas *c = new TCanvas("c", "Efficiencies", 800, 800);
    c->SetLogx();
    c->SetLogy();
    c->SetGrid();
    h1_dummy->SetTitle("Efficiencies");
    h1_dummy->GetYaxis()->SetTitle("#epsilon");
    max = 0.001;
    h1_dummy->GetYaxis()->SetRangeUser(1e-8, max);
    h1_dummy->GetXaxis()->SetTitle("E_{k} [MeV/c]");
    h1_dummy->DrawCopy();
    for(int i=0; i<NFiles; i++) {
        h1_eff[i]->SetTitle(particles[i].Data());
        h1_eff[i]->Draw("e same");

        if(i!=NFiles-1) continue;
        c->BuildLegend();
    }

    



    // ----- RATE ESTIMATION -----

    //histo prompt
    TH1D *h_pRate = GenerateHistogram("h_pRate", "p", "E_{k} [MeV/c]", "Rate [s^{-1}]", nbins, Ek_min_p, Ek_max_p);
    h_pRate->SetLineColor(colors[0]);
    h_pRate->SetLineWidth(LineWidth);

    TH1D *h_eRate = GenerateHistogram("h_eRate", "e^{-}", "E_{k} [MeV/c]", "Rate [s^{-1}]", nbins, Ek_min_e, Ek_max_e);
    h_eRate->SetLineColor(colors[1]);
    h_eRate->SetLineWidth(LineWidth);



    //getting the flux TGraph
    // TGraphErrors *g_Pflux = new TGraphErrors("protoni_conErrori.txt", "%lg %lg %lg %lg");
    TGraph *g_Pflux = new TGraph("txt_file/AP9_fluxIntvAVG_95.txt", "%lg %lg");
    TGraph *g_eflux = new TGraph("txt_file/AE9_fluxIntvAVG_95.txt", "%lg %lg");
    g_Pflux->Sort();
    g_eflux->Sort();

    g_Pflux->SetLineWidth(LineWidth);
    g_eflux->SetLineWidth(LineWidth);

    g_Pflux->SetLineColor(colors[0]);
    g_eflux->SetLineColor(colors[1]);

    g_Pflux->SetLineStyle(2);
    g_eflux->SetLineStyle(2);

    g_Pflux->SetMarkerStyle(20);
    g_eflux->SetMarkerStyle(20);

    g_Pflux->SetMarkerSize(1.2);
    g_eflux->SetMarkerSize(1.2);

    g_Pflux->SetMarkerColor(kMagenta);
    g_eflux->SetMarkerColor(kOrange);

    TCanvas *c2 = new TCanvas("c2", "Fluxes", 800, 600);
    c2->SetLogx();
    c2->SetLogy();
    g_Pflux->SetTitle("p");
    g_eflux->SetTitle("e^{-}");
    // g_eflux->GetYaxis()->SetTitle("Flux [s^{-1} cm^{-2} sr^{-1} (MeV)^{-1}]");
    h1_dummy->GetXaxis()->SetTitle("E_{k} [MeV/c]");
    h1_dummy->GetYaxis()->SetRangeUser(1e-6, 3e4);
    h1_dummy->GetYaxis()->SetTitle("Flux [s^{-1} cm^{-2} sr^{-1} (MeV)^{-1}]");
    h1_dummy->SetTitle("Fluxes");
    h1_dummy->DrawCopy();
    g_eflux->Draw("pl same");
    g_Pflux->Draw("pl same");
    c2->BuildLegend();


    TCanvas *c3 = new TCanvas("c3", "Rate", 800, 600);
    c3->cd();
    c3->SetLogx();
    c3->SetLogy();
    c3->SetGridx(); 
    c3->SetGridy();
    for(int i=1; i<=h_pRate->GetNbinsX(); i++){
        // if(i<=h1_eff[0]->GetNbinsX()){

            //prompt proton
            double binCenter = h1_eff[0]->GetBinCenter(i);        //MeV
            double binWidth = h1_eff[0]->GetBinWidth(i);  
            double fluxValue = g_Pflux->Eval(binCenter);
            double effValue = h1_eff[0]->GetBinContent(i);
            double rate = (binWidth * PlaneAcceptance * effValue) * fluxValue;
            double rate_err = rate * sqrt( pow(h1_eff[0]->GetBinError(i)/h1_eff[0]->GetBinContent(i), 2) + pow(0.1, 2)); // 10% error on flux
            // printf("Rate: %.9e ; binCenter: %.2e ; binWidth: %.2e ; efficiency: %.9e ; AcceptancePlane: %.2e ; flux: %.2e  \n", rate, binCenter, binWidth, effValue, PlaneAcceptance, fluxValue);
            h_pRate->SetBinContent(i, rate);
            h_pRate->SetBinError(i, rate_err);

            //prompt electron
            binCenter = h1_eff[1]->GetBinCenter(i);        //MeV
            binWidth = h1_eff[1]->GetBinWidth(i);    //GeV
            fluxValue = g_eflux->Eval(binCenter);
            effValue = h1_eff[1]->GetBinContent(i);
            rate = (binWidth * PlaneAcceptance * effValue) * fluxValue;
            rate_err = rate * sqrt( pow(h1_eff[1]->GetBinError(i)/h1_eff[1]->GetBinContent(i), 2) + pow(0.1, 2)); // 10% error on flux
            h_eRate->SetBinContent(i, rate);
            h_eRate->SetBinError(i, rate_err);
        // } 
    }
    h_eRate->Draw("e");
    h_pRate->Draw("e same");
    c3->BuildLegend();
    h_eRate->SetTitle("Expected rates");
    printf("Electron prompt  ----- Total rate: %.5e ; \n", h_eRate->Integral());
    printf("Proton prompt     ----- Total rate: %.5e ; \n", h_pRate->Integral());

    TFile *f_out = new TFile("rate_proton_electron.root", "RECREATE");
    h_pRate->Write();
    h_eRate->Write();
    f_out->Close();

    return;

}







