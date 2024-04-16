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

enum detectors_order{GAGG_i, Si1_i, Si2_i};
enum particles_order{proton, electron};




void GeometricAcceptance(){
    // gStyle->SetOptStat(0); // no statistics

    //Files to be used
    const int NFiles = 2;
    const int NDetectors = 3;
    TString particles[NFiles] = {"p", "e^{-}"};
    TString names[NFiles] = {"../build/GPSproton_8cm_run_0_um.root", "../build/GPSelectron_8cm_run_0_um.root"};
    TString h1_names[NFiles] = {"Gen_spectra_P", "Gen_spectra_e"};
    TString h1_detNames[NDetectors] = {"GAGG", "Si1", "Si2"};


    //Histos for Geometric acceptances
    TH1D *h1_acc[NFiles];
    TH1D *h1_accDetectors[NFiles*NDetectors];

    const int LineWidth = 3;
    Color_t part_colors[NFiles] = {kBlue, kRed};
    Color_t det_colors[NDetectors] = {kGreen+1, kMagenta+1, kOrange+1};
    TPaveText* title;
    const int nbins = 100;
    const double Ek_min_p = 1, Ek_max_p = 500; // MeV
    const double Ek_min_e = 0.05, Ek_max_e = 500; // MeV
    const double BinWidthLog_p = (TMath::Log10(Ek_max_p) - TMath::Log10(Ek_min_p)) / nbins; // bin width in log scale
    const double BinWidthLog_e = (TMath::Log10(Ek_max_e) - TMath::Log10(Ek_min_e)) / nbins; // bin width in log scale

    //Generation events info and histos
    const double GenArea = 64; // cm^2 
    const double PlaneAcceptance = TMath::Pi() * GenArea; // cm^2 sr
    TH1D *h1_genEvents[NFiles];

    //Intializing the histograms

    for(int i=0; i<NFiles; i++) {
        if(names[i].Contains("proton")) {
            printf("Creating histogram for protons\n");
            h1_acc[i] = GenerateHistogram("h1_acc_"+particles[i], "Geometric_acceptance_"+particles[i], "E_{k} [MeV/c]", "Acceptance [cm^{2} sr]", nbins, Ek_min_p, Ek_max_p);
            for(int index = 0; index<NDetectors; index++){
                h1_accDetectors[(NDetectors*i) +index] = GenerateHistogram("h1_"+h1_detNames[index]+"_"+particles[i], "Geometric_acceptance_"+h1_detNames[index]+"_"+particles[i], "E_{k} [MeV/c]", "Acceptance [cm^{2} sr]", nbins, Ek_min_p, Ek_max_p);
                h1_accDetectors[(NDetectors*i) +index]->SetLineWidth(LineWidth);
                h1_accDetectors[(NDetectors*i) +index]->SetLineColor(det_colors[index]);
                h1_accDetectors[(NDetectors*i) +index]->GetYaxis()->SetTitleOffset(1.0);
                h1_accDetectors[(NDetectors*i) +index]->GetXaxis()->SetTitleOffset(1.2);
            }

        }
        if(names[i].Contains("electron")) {
            printf("Creating histogram for electrons\n");
            h1_acc[i] = GenerateHistogram("h1_acc_"+particles[i], "Geometric_acceptance_"+particles[i], "E_{k} [MeV/c]", "Acceptance [cm^{2} sr]", nbins, Ek_min_e, Ek_max_e);
            for(int index = 0; index<NDetectors; index++){
                h1_accDetectors[(NDetectors*i) + index] = GenerateHistogram("h1_"+h1_detNames[index]+"_"+particles[i], "Geometric_acceptance_"+h1_detNames[index]+"_"+particles[i], "E_{k} [MeV/c]", "Acceptance [cm^{2} sr]", nbins, Ek_min_e, Ek_max_e);
                h1_accDetectors[(NDetectors*i) +index]->SetLineWidth(LineWidth);
                h1_accDetectors[(NDetectors*i) +index]->SetLineColor(det_colors[index]);
                h1_accDetectors[(NDetectors*i) +index]->GetYaxis()->SetTitleOffset(1.0);
                h1_accDetectors[(NDetectors*i) +index]->GetXaxis()->SetTitleOffset(1.2);
            }
        }
        
        h1_acc[i]->SetLineWidth(LineWidth);
        h1_acc[i]->SetLineColor(part_colors[i]);
        h1_acc[i]->GetYaxis()->SetTitleOffset(1.0);
        h1_acc[i]->GetXaxis()->SetTitleOffset(1.2);
        }


    //Trigger thresholds and cuts
    const double th_Calo0 = 0.08; // Mev = 80 kev
    const double th_Silicon = 0.;
    const double th_DrilledVeto = 0.;
    const double th_BottomVeto = 0.;

    const int NCut = 5;
    TString charged_Sel[NCut] = {"Ed_DrilledVeto", "Ed_BottomVeto==0", "Thin_x0_y1_ID1", "Ed_Calo0", "Thin_x0_y0_ID0"};


    //Tree variables
    double RandEnergy, Ed_Calo0, Ed_DrilledVeto, Ed_BottomVeto, Thin_x0_y0_ID0, Thin_x0_y1_ID1;


    //Loop on the files
    for(int iFile=0; iFile<NFiles; iFile++){

        TFile *file = new TFile(names[iFile], "OPEN");
        printf("\n Reading file %d/%d: %s\n", iFile+1, NFiles, names[iFile].Data());

        TTree *t = (TTree *)file->Get("Edep");
        h1_genEvents[iFile] = (TH1D *)file->Get(h1_names[iFile]);
        h1_genEvents[iFile]->SetLineColor(part_colors[iFile]);


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
               Thin_x0_y0_ID0 > th_Silicon) {

                h1_acc[iFile]->Fill(RandEnergy);
                for(int index = 0; index<NDetectors; index++){
                    if(index == Si2_i) continue;
                    h1_accDetectors[(NDetectors*iFile) + index]->Fill(RandEnergy);
                }

               }
            
            if(Ed_DrilledVeto == th_DrilledVeto &&
               Ed_BottomVeto == th_BottomVeto && 
               Thin_x0_y0_ID0 == th_Silicon &&
               Ed_Calo0 > th_Calo0 &&
               Thin_x0_y1_ID1 > th_Silicon) {

                h1_acc[iFile]->Fill(RandEnergy);
                for(int index = 0; index<NDetectors; index++){
                    if(index == Si1_i ) continue;
                    h1_accDetectors[(NDetectors*iFile) + index]->Fill(RandEnergy);
                }
            }
        
        } // end of the event loop

        //Getting the Geometric acceptances
        h1_acc[iFile]->Divide(h1_genEvents[iFile]);
        h1_acc[iFile]->Scale(PlaneAcceptance);

        for(int index = 0; index<NDetectors; index++){
            h1_accDetectors[(NDetectors*iFile) + index]->Divide(h1_genEvents[iFile]);
            h1_accDetectors[(NDetectors*iFile) + index]->Scale(PlaneAcceptance);
        }

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

   

    TCanvas *c = new TCanvas("c", "Geometric acceptances", 800, 800);
    c->SetLogx();
    c->SetLogy();
    c->SetGrid();
    h1_dummy->SetTitle("Geometric acceptances");
    h1_dummy->GetYaxis()->SetTitle("Acceptance [cm^{2} sr]");
    max = 1;
    h1_dummy->GetYaxis()->SetRangeUser(1e-3, max);
    h1_dummy->GetXaxis()->SetTitle("E_{k} [MeV/c]");
    h1_dummy->DrawCopy();
    for(int i=0; i<NFiles; i++) {
        h1_acc[i]->SetTitle(particles[i].Data());
        h1_acc[i]->Draw("e same");

        // for(int index = 0; index<NDetectors; index++){
        //     h1_accDetectors[(NDetectors*i) + index]->Draw("e same");
        // }

        if(i!=NFiles-1) continue;
        c->BuildLegend();
    }

    
    // ----- RATE ESTIMATION -----

    //histo particles
    TH1D *h1_pRate[NFiles];
    TString rate_names[NFiles] = {"h1_pRate", "h1_eRate"};

    //Detectors histos
    TH1D* h1_detRate[NDetectors*NFiles];


    //Declaration of the histos
    for(int i_hist=0;i_hist<NFiles;i_hist++){

        if(i_hist==0){
            h1_pRate[i_hist] = GenerateHistogram(rate_names[i_hist]+"_"+particles[i_hist], particles[i_hist], "E_{k} [MeV/c]", "Rate [s^{-1}]", nbins, Ek_min_p, Ek_max_p);

            for(int i_det=0; i_det<NDetectors; i_det++){
                h1_detRate[(NDetectors*i_hist) + i_det] = GenerateHistogram("h1_Rate_"+h1_detNames[i_det]+"_"+particles[i_hist], h1_detNames[i_det]+"_"+particles[i_hist], "E_{k} [MeV/c]", "Rate [s^{-1}]", nbins, Ek_min_p, Ek_max_p);
                h1_detRate[(NDetectors*i_hist) + i_det]->SetLineColor(det_colors[i_det]);
                h1_detRate[(NDetectors*i_hist) + i_det]->SetLineWidth(LineWidth);
            }
        }
        else{
            h1_pRate[i_hist] = GenerateHistogram(rate_names[i_hist]+"_"+particles[i_hist], particles[i_hist], "E_{k} [MeV/c]", "Rate [s^{-1}]", nbins, Ek_min_e, Ek_max_e);

            for(int i_det=0; i_det<NDetectors; i_det++){
                h1_detRate[(NDetectors*i_hist) + i_det] = GenerateHistogram("h1_Rate_"+h1_detNames[i_det]+"_"+particles[i_hist], h1_detNames[i_det]+"_"+particles[i_hist], "E_{k} [MeV/c]", "Rate [s^{-1}]", nbins, Ek_min_e, Ek_max_e);
                h1_detRate[(NDetectors*i_hist) + i_det]->SetLineColor(det_colors[i_det]);
                h1_detRate[(NDetectors*i_hist) + i_det]->SetLineWidth(LineWidth);
            }
        } 

        h1_pRate[i_hist]->SetLineColor(part_colors[i_hist]);
        h1_pRate[i_hist]->SetLineWidth(LineWidth);       
    }

    // getting the flux TGraph
    TGraph *g_fluxes[NFiles];
    TString txt_files[NFiles] = {"txt_file/AP9_fluxIntvAVG_95.txt", "txt_file/AE9_fluxIntvAVG_95.txt"};
    for(int i_flux=0; i_flux<NFiles; i_flux++){
        g_fluxes[i_flux] = new TGraph(txt_files[i_flux], "%lg %lg");
        g_fluxes[i_flux]->Sort();
        g_fluxes[i_flux]->SetLineWidth(LineWidth);
        g_fluxes[i_flux]->SetLineColor(part_colors[i_flux]);
        g_fluxes[i_flux]->SetMarkerColor(part_colors[i_flux]);
        g_fluxes[i_flux]->SetLineStyle(2);
        g_fluxes[i_flux]->SetMarkerStyle(20);
        g_fluxes[i_flux]->SetMarkerSize(1.2);
        g_fluxes[i_flux]->SetTitle(particles[i_flux]);

    }

    TCanvas *c2 = new TCanvas("c2", "Fluxes", 800, 600);
    c2->SetLogx();
    c2->SetLogy();
    h1_dummy->GetXaxis()->SetTitle("E_{k} [MeV/c]");
    h1_dummy->GetYaxis()->SetRangeUser(1e-6, 3e4);
    h1_dummy->GetYaxis()->SetTitle("Flux [s^{-1} cm^{-2} sr^{-1} (MeV)^{-1}]");
    h1_dummy->SetTitle("Fluxes");
    h1_dummy->DrawCopy();
    for(int i=0; i<NFiles; i++) g_fluxes[i]->Draw("pl same");
    c2->BuildLegend();


    
    for(int i=1; i<=h1_pRate[proton]->GetNbinsX(); i++){

        // proton
        double binCenter = h1_acc[proton]  ->GetBinCenter(i); //MeV
        double binWidth = h1_acc[proton]   ->GetBinWidth(i);  
        double fluxValue = g_fluxes[proton]->Eval(binCenter);
        double accValue = h1_acc[proton]   ->GetBinContent(i);
        double rate = (binWidth * accValue) * fluxValue;
        double rate_err = rate * sqrt( pow(h1_acc[proton]->GetBinError(i)/h1_acc[proton]->GetBinContent(i), 2) + pow(0.1, 2)); // 10% error on flux
        if(rate == 0){
            h1_pRate[proton]->SetBinContent(i, 0);
            h1_pRate[proton]->SetBinError(i, 0);
        }
        else{
            h1_pRate[proton]->SetBinContent(i, rate);
            h1_pRate[proton]->SetBinError(i, rate_err);
        }
        

        for(int i_det=0; i_det<NDetectors; i_det++){
            binCenter = h1_accDetectors[(proton*NDetectors) + i_det] ->GetBinCenter(i);  //MeV
            binWidth = h1_accDetectors[(proton*NDetectors) + i_det]  ->GetBinWidth(i);   //MeV
            fluxValue = g_fluxes[proton]                             ->Eval(binCenter);
            accValue = h1_accDetectors[(proton*NDetectors) + i_det]  ->GetBinContent(i);
            rate = (binWidth * accValue) * fluxValue;
            if(rate==0){
                h1_detRate[(proton*NDetectors) + i_det]->SetBinContent(i, 0);
                h1_detRate[(proton*NDetectors) + i_det]->SetBinError(i, 0);
            }
            else{
                rate_err = rate * sqrt( pow(h1_accDetectors[(proton*NDetectors) + i_det]->GetBinError(i)/h1_accDetectors[(proton*NDetectors) + i_det]->GetBinContent(i), 2) + pow(0.1, 2)); // 10% error on flux
                h1_detRate[(proton*NDetectors) + i_det]->SetBinContent(i, rate);
                h1_detRate[(proton*NDetectors) + i_det]->SetBinError(i, rate_err);
            }
        } // detectors rate for protons


        // electron
        binCenter = h1_acc[electron]  ->GetBinCenter(i);  //MeV
        binWidth = h1_acc[electron]   ->GetBinWidth(i);   //MeV
        fluxValue = g_fluxes[electron]->Eval(binCenter);
        accValue = h1_acc[electron]   ->GetBinContent(i);
        rate = (binWidth * accValue) * fluxValue;
        rate_err = rate * sqrt( pow(h1_acc[electron]->GetBinError(i)/h1_acc[electron]->GetBinContent(i), 2) + pow(0.1, 2)); // 10% error on flux
        if(rate==0){
            h1_pRate[electron]->SetBinContent(i, 0.);
            h1_pRate[electron]->SetBinError(i, 0.);
        }
        else{
            h1_pRate[electron]->SetBinContent(i, rate);
            h1_pRate[electron]->SetBinError(i, rate_err);
        }
        

        for(int i_det=0; i_det<NDetectors; i_det++){
            binCenter = h1_accDetectors[(electron*NDetectors) + i_det] ->GetBinCenter(i);  //MeV
            binWidth = h1_accDetectors[(electron*NDetectors) + i_det]  ->GetBinWidth(i);   //MeV
            fluxValue = g_fluxes[electron]                             ->Eval(binCenter);
            accValue = h1_accDetectors[(electron*NDetectors) + i_det]  ->GetBinContent(i);
            rate = (binWidth * accValue) * fluxValue;
            if(rate==0){
                h1_detRate[(electron*NDetectors) + i_det]->SetBinContent(i, 0);
                h1_detRate[(electron*NDetectors) + i_det]->SetBinError(i, 0);
            }
            else{
                rate_err = rate * sqrt( pow(h1_accDetectors[(electron*NDetectors) + i_det]->GetBinError(i)/h1_accDetectors[(electron*NDetectors) + i_det]->GetBinContent(i), 2) + pow(0.1, 2)); // 10% error on flux
                h1_detRate[(electron*NDetectors) + i_det]->SetBinContent(i, rate);
                h1_detRate[(electron*NDetectors) + i_det]->SetBinError(i, rate_err);
            }
        } //detectors rate for electrons

    } //end of the bin loop



    TCanvas *c3 = new TCanvas("c3", "Rate", 800, 600);
    c3->cd();
    c3->SetLogx();
    c3->SetLogy();
    c3->SetGridx(); 
    c3->SetGridy();
    h1_pRate[electron]->GetYaxis()->SetRangeUser(1e-4, 1e3);

    h1_pRate[electron]->Draw("e");
    h1_pRate[proton]->Draw("e same");

    c3->BuildLegend();
    h1_pRate[electron]->SetTitle("Expected rates");

    double ElectronCheck = 0, ProtonCheck = 0;
    ProtonCheck = (h1_detRate[GAGG_i]->Integral() +  h1_detRate[Si1_i]->Integral() + h1_detRate[Si2_i]->Integral()) - h1_pRate[proton]->Integral();

    ElectronCheck = (h1_detRate[NDetectors+GAGG_i]->Integral() + h1_detRate[NDetectors+Si1_i]->Integral() + h1_detRate[NDetectors+Si2_i]->Integral()) - h1_pRate[electron]->Integral();

    // La somma degli integrali dei rate attesi sui detectors deve essere il doppio del rate atteso per una data particella. 
    // Di conseguenza le variabili ElectronCheck e ProtonCheck devono essere uguali agli integrali di h1_pRate[electron] e h1_pRate[proton] (rispettivamente)
    printf("Electron prompt  ----- Total rate: %.5e  ElectronCheck: %.5e; \n", h1_pRate[electron]->Integral(), ElectronCheck);
    printf("Proton prompt     ----- Total rate: %.5e ProtonCheck: %.5e; \n", h1_pRate[proton]->Integral(), ProtonCheck);

    TFile *f_out = new TFile("SPaRKLE_AcquisitionRate.root", "RECREATE");
    c->Write();
    for(int i_hist=0; i_hist<NFiles; i_hist++){
        h1_acc[i_hist]->Write();
        h1_genEvents[i_hist]->Write();
        for(int i_det=0; i_det<NDetectors; i_det++){
            h1_accDetectors[(NDetectors*i_hist) + i_det]->Write();
        }
        for(int i_det=0; i_det<NDetectors; i_det++){
            h1_detRate[(NDetectors*i_hist) + i_det]->Write();
        }
        h1_pRate[i_hist]->Write();
    }

    f_out->Close();

    return;

}







