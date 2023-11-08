#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "TCanvas.h"
#include "TMath.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TBranch.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TColor.h"
#include "THStack.h"
#include "TColor.h"

#include "./Custom/TH1DLog.h"

using namespace std;



void GeometricalFactor_CZT()
{
    int Nfiles = 3;
    int Nx = 4;
    int Ny = 4;
    int Ntot = Nx * Ny; 
    int j;
    int CopyNumber;


    double ResSilicon;
    double ResPlastic;
    double E_min_thin;   // Thick layer
    double E_min_thick;  // Thin layer
    double E_th_Vetoed;  // Energy dispersion (Veto threshold)
    double E_th_plastic;
    
    double ERange[3][2];
    double LActive;
    
    // ********************************************
    // ********************************************
    //                 SETTINGS
    // ********************************************
    // ********************************************
    Nfiles       = 3;
    // Geometric parameter
    Nx           = 4;
    Ny           = 4;
    Ntot         = Nx * Ny;
    // Smearing parameters
    ResSilicon   = 0.01;
    ResPlastic   = 0.3;
    // Trigger and veto thresholds
    E_min_thin   = 0.02;
    E_min_thick  = 0.04;
    E_th_Vetoed  = 0.1;
    E_th_plastic = 0.2;
    
    // ********************************************
    //           GEOMETRICAL FACTOR SETTINGS
    // ********************************************
    
    ERange[0][0] = 0.08;
    ERange[0][1] = 14.;
    ERange[1][0] = 0.01;
    ERange[1][1] = 100.;
    ERange[2][0] = 2.;
    ERange[2][1] = 400.;
    // Number of bins
    int Nbins[3] = {100, 100, 100};
    LActive = 2 * 6.0;
    
    



    TString FileName[Nfiles];

    // Read the FileNames.txt and for each line of the file fill the TString FileName[]
    ifstream FileNames("FileNames.txt");
    if(!FileNames.is_open())
    {
        cout << "FileNames.txt not found" << endl;
        return;
    }
    else
    {
        for(int i = 0; i < Nfiles; i++)
        {
            FileNames >> FileName[i];
        }
    }

    // Definitions of the Ttrees and Tfiles
    TFile *file[Nfiles];
    TTree *Edep[Nfiles];

    for(int i = 0; i < Nfiles ; ++i)
    {
        file[i] = TFile::Open(FileName[i]);
        Edep[i] = (TTree*) file[i] -> Get("Edep");
        Edep[i] -> Print();
    }

    // Definitions of the initial quantities from the data of the simulation
    TString BranchName[47];
    TString TotEnergyName[16];
    TString TotEnergyCondition[16];
    TString PIDName[16];
    TString DirName[3];
    TString PolarAngle[2];
    TString NewPolarAngle[2];
    TString AliasETot;
    TString ThinTot;
    TString ThickTot;




    AliasETot        = "ETot";
    DirName[0]       = "DirX";
    DirName[1]       = "DirY";
    DirName[2]       = "DirZ";
    PolarAngle[0]    = "dirTheta";
    PolarAngle[1]    = "dirPhi";
    NewPolarAngle[0] = "dirThetaNew";
    NewPolarAngle[1] = "dirPhiNew";



    // Building the branch names
    BranchName[0]  = "NumberID";
    BranchName[1]  = "RandEnergy";
    BranchName[2]  = "RandNumber";
    BranchName[3]  = "Xgen";
    BranchName[4]  = "Ygen";
    BranchName[5]  = "Zgen";
    BranchName[6]  = "pDirX";
    BranchName[7]  = "pDirY";
    BranchName[8]  = "pDirZ";
    BranchName[9]  = "Ed_Veto0";
    BranchName[10] = "Ed_DrilledVeto";
    BranchName[43] = "Ed_Lateral0";
    BranchName[44] = "Ed_Lateral1";
    BranchName[45] = "Ed_Lateral2";
    BranchName[46] = "Ed_Lateral3";


    j = 11;
    CopyNumber = 0;
    
    for(int ix = 0; ix < Nx; ++ix)
    {
        for(int iy = 0; iy < Ny; ++iy)
        {   
            BranchName[j]        = Form("Thin_x%d_y%d_ID%d",ix,iy,CopyNumber);
            BranchName[j + Ntot] = Form("Thick_x%d_y%d_ID%d",ix,iy,CopyNumber);
            if (CopyNumber == 0)
            {
                ThinTot  = Form("(%s)", BranchName[j].Data());
                ThickTot = Form("(%s)", BranchName[j + Ntot].Data());
            }
            else
            {
                ThinTot  = Form("(%s + %s)",ThinTot.Data(), BranchName[j].Data());
                ThickTot = Form("(%s + %s)",ThickTot.Data(), BranchName[j + Ntot].Data());
            }

            ++j;
            ++CopyNumber;
        }
    }

    cout << "Branches Names : \n############################" << endl;
    for(int i = 0; i < 47; ++i)
    {
        cout << i << " " << BranchName[i].Data() << endl;
    }

    /* -------------------------------------------------------------------------- */
    /*                              Alias definitions                             */
    /* -------------------------------------------------------------------------- */

    /* -------------------------------------------------------------------------- */
    /*                              ALIAS DEFINITIONS                             */
    /* -------------------------------------------------------------------------- */

    for(int i = 0; i < Nfiles ; ++i)
    {
        cout << "Setting alias in File " << i << " : " << FileName[i].Data() << endl;
        /* -------------------------------- SMEARING -------------------------------- */
        /* ------------------------ YOU ONLY NEED TO ADD A g ------------------------ */
        for(int k = 9; k < 47; ++k)
        {
            Edep[i] -> SetAlias(Form("wnorm%d",k),"(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))");
            if(k == 9 || k == 10 || (k >= 43 && k <= 46))
            {
                Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResPlastic));
            }
            else 
            {
                Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
            }
        }

        /* ------------------- Incident direction of the particle ------------------- */

        Edep[i] -> SetAlias("NormP", Form("(TMath::Sqrt(TMath::Power(%s,2) + TMath::Power(%s,2) + TMath::Power(%s,2)))", BranchName[6].Data(), BranchName[7].Data(), BranchName[8].Data()));
        for(int k = 6; k <= 8; ++k)
        {
            Edep[i] -> SetAlias(DirName[k-6].Data(), Form("((%s)/(NormP))", BranchName[k].Data()));
        }

        // Legend
        // [0] : Theta
        // [1] : Phi
        Edep[i] -> SetAlias(PolarAngle[0].Data()   , Form("((TMath::ACos(%s)))", DirName[2].Data())); // *(180/3.415927))
        Edep[i] -> SetAlias(PolarAngle[1].Data()   , Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data())); // *(180/3.415927))
        Edep[i] -> SetAlias(NewPolarAngle[1].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(TMath::Cos(%s))))",              PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() )); // *(180/3.415927))
        Edep[i] -> SetAlias(NewPolarAngle[0].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)),TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data())); // *(180/3.415927))
        
        /* -------------------- Particle identification parameter ------------------- */

        j = 11;
        CopyNumber = 0;
        cout << "Defining total energy in File " << i << " : " << FileName[i].Data() << endl;
        for(int ix = 0; ix < Nx; ++ix)
        {
            for(int iy = 0; iy < Ny; ++iy)
            {   
                TotEnergyName[CopyNumber] = Form("Tot_x%d_y%d_ID%d",ix,iy,CopyNumber);
                PIDName[CopyNumber] = Form("PID_x%d_y%d_ID%d",ix,iy,CopyNumber);
                TotEnergyCondition[CopyNumber]  = Form("(%s + %s",BranchName[j].Data(), BranchName[j + Ntot].Data());
                //TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[9].Data(), BranchName[9].Data(),E_th_plastic);
                //TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[43].Data(), BranchName[43].Data(),E_th_plastic);
                //TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[44].Data(), BranchName[44].Data(),E_th_plastic);
                //TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[45].Data(), BranchName[45].Data(),E_th_plastic);
                TotEnergyCondition[CopyNumber] += ")";
                Edep[i] -> SetAlias(TotEnergyName[CopyNumber], TotEnergyCondition[CopyNumber]);
                Edep[i] -> SetAlias(PIDName[CopyNumber], Form("(TMath::Log10(g%s * %s))",BranchName[j].Data(), TotEnergyName[CopyNumber].Data()));
                ++j;
                ++CopyNumber;
                cout << Form("ix= %d iy= %d CopyNumber= %d",ix, iy, CopyNumber) << endl;
            }
        }
        /* ------------------------------- Total Thin ------------------------------- */
        /* ------------------------------- Total Thick ------------------------------ */
    }


    /* -------------------------------------------------------------------------- */
    /*                   Definition of the SELECTION conditions                   */
    /* -------------------------------------------------------------------------- */

    TString ConditionPairSilicon[Ntot];
    TString ConditionPairSiliconAll;
    TString ConditionEnergyDispersion;
    TString ConditionDrilledVeto;
    TString ConditionGoodEvents;
    TString ConditionGoodEventsSinglePair[Ntot];


    for(int i = 0; i< Ntot; ++i)
    {
        cout << "Defining good events for pair" << BranchName[11+i].Data() << " & " <<BranchName[11+ Ntot+i].Data() << endl;
        ConditionPairSilicon[i] = Form("((%s > %g) && (%s > %g))", BranchName[11+i].Data(), E_min_thin, BranchName[11 + Ntot + i].Data(), E_min_thick);
        cout << "ConditionPairSilicon[" << i << "] = " << ConditionPairSilicon[i].Data() << endl;
        if(i == 0)
        {
            ConditionPairSiliconAll = Form("((%s > %g) && (%s > %g))", BranchName[11+i].Data(), E_min_thin, BranchName[11 + Ntot + i].Data(), E_min_thick);
        }
        else 
        {
            ConditionPairSiliconAll += Form("|| ((%s > %g) && (%s > %g))", BranchName[11+i].Data(), E_min_thin, BranchName[11 + Ntot + i].Data(), E_min_thick);
        }
    }

    ConditionDrilledVeto       = Form("(%s < %g)", BranchName[10].Data(), E_th_Vetoed);
    //ConditionEnergyDispersion  = Form("((%s) - (%s) - (%s) - (%s) - (%s) - (%s) - (%s)) < %g", BranchName[1].Data(), ThinTot.Data(), ThickTot.Data(), BranchName[9].Data(),BranchName[43].Data(), BranchName[44].Data(), BranchName[45].Data(), E_th_Vetoed);
    ConditionEnergyDispersion  = Form("((%s)" , BranchName[1].Data());
    ConditionEnergyDispersion += Form("- (%s)", ThinTot.Data());
    ConditionEnergyDispersion += Form("- (%s)", ThickTot.Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[9].Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[43].Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[44].Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[45].Data());
    ConditionEnergyDispersion += Form(") < %g", E_th_Vetoed);
    
    ConditionGoodEvents        = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
    for(int i = 0; i < Ntot; ++i)
    {
        cout << "Defining good events for pair" << BranchName[11+i].Data() << " & " <<BranchName[11+ Ntot+i].Data() << endl;
        ConditionGoodEventsSinglePair[i] = Form("(%s) && (%s) && (%s)", ConditionPairSilicon[i].Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
    }


    /* -------------------------------------------------------------------------- */
    /*                    Specific part for GEOMETRICAL FACTOR                    */
    /* -------------------------------------------------------------------------- */

    
    TH1DLog *hGenerated[3];
    TH1DLog *hAccepted[3];
    TH1D *hGen[3];
    TH1D *hAcc[3];
    
    int Colors [3] = {kOrange+6, kMagenta+1, kAzure+10};
    
    TCanvas *cAcceptance[3];
    TString ParticleName[3] = {"electrons", "protons", "alpha"};

    int NBins = 70;

/*
    for(int i = 0; i < 3; ++i)
    {
        cout << "Entering loop for acceptance plot histogram 1 : i = " <<   i << endl;

        hGenerated[i] = new TH1DLog();
        hGenerated[i] -> SetName(Form("h_Gen%d", i));
        hGenerated[i] -> SetTitle(Form("Generated %s", ParticleName[i].Data()));
        hGenerated[i] -> SetXTitle("Energy [MeV]");
        hGenerated[i] -> SetYTitle("Counts");
        hGenerated[i] -> SetXAxis(ERange[i][0], ERange[i][1], NBins);
        hGenerated[i] -> GenerateHistogram();

        hGen[i] = hGenerated[i] -> GetHistogram();

        hAccepted[i] = new TH1DLog();
        hAccepted[i] -> SetName(Form("h_Acc%d", i));
        hAccepted[i] -> SetTitle(Form("Acceptance %s", ParticleName[i].Data()));
        hAccepted[i] -> SetXTitle("Energy [MeV]");
        hAccepted[i] -> SetYTitle("Acceptance [cm^{2} sr]");
        hAccepted[i] -> SetXAxis(ERange[i][0], ERange[i][1], NBins);
        hAccepted[i] -> GenerateHistogram();

        hAcc[i] = hAccepted[i] -> GetHistogram();


        Edep[i] -> Draw(Form("%s >> h_Gen%d", BranchName[1].Data(), i), "", "goff");

        Edep[i] -> Draw(Form("%s >> h_Acc%d", BranchName[1].Data(), i), ConditionGoodEvents.Data(), "goff");

        hAcc[i] -> Sumw2(1);
        hGen[i] -> Sumw2(1);

        hAcc[i] -> Divide(hGen[i]);
        
        hAcc[i] -> Scale((LActive*LActive * TMath::Pi()));

        cAcceptance[i] = new TCanvas(Form("cAcceptance%d", i), Form("cAcceptance%d", i), 800, 600);

        hAcc[i] -> Draw("lp");

        hAcc[i] -> SetLineWidth(2);
        hAcc[i] -> SetLineColor(kBlue);
        hAcc[i] -> SetMarkerStyle(20);
        hAcc[i] -> SetMarkerSize(0.5);
        hAcc[i] -> SetMarkerColor(kBlue);

        gPad -> SetGrid();
        gPad -> SetLogx();

        cAcceptance[i] -> SaveAs(Form("Acceptance_%s.pdf", ParticleName[i].Data()));
        cAcceptance[i] -> SaveAs(Form("Acceptance_%s.png", ParticleName[i].Data()));
    }

    THStack *hsAcceptance = new THStack("hsAcceptance", "hsAcceptance");


    

    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);

    for(int i = 0; i < 3; ++i)
    {
        hAcc[i] -> SetLineColor(Colors[i]);
        hAcc[i] -> SetMarkerColor(Colors[i]);
        leg -> AddEntry(hAcc[i], ParticleName[i].Data(), "lp");
        hsAcceptance -> Add(hAcc[i]);
    }

    */
    

    TH1DLog *hGenerated_2[3];
    TH1DLog *hAccepted_2[3];
    TH1D *hGen_2[3];
    TH1D *hAcc_2[3];
    THStack *hsAcceptance_2 = new THStack("hsAcceptance_2", "hsAcceptance_2");
    TLegend *leg_2 = new TLegend(0.1, 0.1, 0.4, 0.3);

    int NBins_2 = 150;
    double ERange_2[2] = {0.08,100.};

    for(int i = 0; i < Nfiles; ++i)
    {
        cout << "Entering loop for acceptance plot histogram 2 : i = " <<   i << endl;
        hGenerated_2[i] = new TH1DLog();
        hGenerated_2[i] -> SetName(Form("h_Gen_2_%d", i));
        hGenerated_2[i] -> SetTitle(Form("Generated %s", ParticleName[i].Data()));
        hGenerated_2[i] -> SetXTitle("Energy [MeV]");
        hGenerated_2[i] -> SetYTitle("Counts");
        hGenerated_2[i] -> SetXAxis(ERange_2[0], ERange_2[1], NBins_2);
        hGenerated_2[i] -> GenerateHistogram();

        hGen_2[i] = hGenerated_2[i] -> GetHistogram();

        hAccepted_2[i] = new TH1DLog();
        hAccepted_2[i] -> SetName(Form("h_Acc_2_%d", i));
        hAccepted_2[i] -> SetTitle(Form("Acceptance %s", ParticleName[i].Data()));
        hAccepted_2[i] -> SetXTitle("Energy [MeV]");
        hAccepted_2[i] -> SetYTitle("Acceptance [cm^{2} sr]");
        hAccepted_2[i] -> SetXAxis(ERange_2[0], ERange_2[1], NBins_2);
        hAccepted_2[i] -> GenerateHistogram();

        hAcc_2[i] = hAccepted_2[i] -> GetHistogram();

        Edep[i] -> Draw(Form("%s >> h_Gen_2_%d", BranchName[1].Data(), i), "", "goff");
        Edep[i] -> Draw(Form("%s >> h_Acc_2_%d", BranchName[1].Data(), i), ConditionGoodEvents.Data(), "goff");
    
        hAcc_2[i] -> Sumw2(1);
        hGen_2[i] -> Sumw2(1);

        hAcc_2[i] -> Divide(hGen_2[i]);

        hAcc_2[i] -> Scale((LActive*LActive * TMath::Pi()));

        hAcc_2[i] -> SetLineColor(Colors[i]);
        hAcc_2[i] -> SetMarkerColor(Colors[i]);
        hAcc_2[i] -> SetLineWidth(3);
        hAcc_2[i] -> SetMarkerStyle(20);
        hAcc_2[i] -> SetMarkerSize(0.8);

        leg_2 -> AddEntry(hAcc_2[i], ParticleName[i].Data(), "lp");
        hsAcceptance_2 -> Add(hAcc_2[i]);

    }


    TCanvas *c_Acc_all = new TCanvas("c_Acc_all", "c_Acc_all", 1000, 600);; 
    hsAcceptance_2 -> Draw("nostack");
    hsAcceptance_2 -> GetXaxis() -> SetTitle("Energy [MeV]");
    hsAcceptance_2 -> GetYaxis() -> SetTitle("Acceptance [cm^{2} sr]");
    hsAcceptance_2 -> SetTitle("");
    leg_2 -> Draw("same");
    gPad -> SetGrid();
    gPad -> SetLogx();
    gPad -> SetLogy();
    c_Acc_all -> SaveAs("Acceptance_all.pdf");
    c_Acc_all -> SaveAs("Acceptance_all.root");




    return;
}