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

using namespace std;


void PID()
{

    // Some initial definitions
    bool DebugMuons = false;
    bool Plastic = false;
    int CopyNoToPlot;

    int Nfiles = 3;
    int Nx = 4;
    int Ny = 4;
    int Ntot = Nx * Ny; 
    int j;
    int CopyNumber;
    double LThin;
    double xDelta;
    double yDelta;
    double thetaMax;
    double RMax;
    double distanceR;
    double ResSilicon;
    double ResPlastic;
    double E_min_thin;   // Thick layer
    double E_min_thick;  // Thin layer
    double E_th_Vetoed;  // Energy dispersion (Veto threshold)
    double E_th_plastic;
    double Emaxx;
    double Emaxy;
    double Eminx;
    double Eminy;
    double Nbinsx;
    double Nbinsy;


    // ********************************************
    // ********************************************
    //                 SETTINGS
    // ********************************************
    // ********************************************
    DebugMuons = false;
    Nfiles       = 3;
    // Geometric parameter
    Nx           = 4;
    Ny           = 4;
    Ntot         = Nx * Ny;
    LThin        = 10.;
    xDelta       = LThin / (Nx);
    yDelta       = LThin / (Ny);
    thetaMax     = 35. * TMath::Pi() / 180.;
    RMax         = std::sqrt(2) * 1.5;
    // Smearing parameters
    ResSilicon   = 0.01;
    ResPlastic   = 0.3;
    // Trigger and veto thresholds
    E_min_thin   = 0.02;
    E_min_thick  = 0.04;
    E_th_Vetoed  = 0.1;
    E_th_plastic = 0.2;
    // Plotting parameters (canvas)
    // Energy is in log10 scale
    Emaxx       = 3;
    Emaxy       = 3.4;
    Eminx       = -1.4;
    Eminy       = -3.2;
    Nbinsx      = 300;
    Nbinsy      = 250;

    CopyNoToPlot = 6;
    // ********************************************
    // ********************************************


    if(DebugMuons)
    Nfiles = 4;
    cout << "Nfiles = " << Nfiles << endl;


    TString FileName[Nfiles];
    //FileName[0] = "Edep_05_08_2022_ELECTRON.root";
    //FileName[1] = "Edep_05_08_2022_PROTON.root";
    //FileName[2] = "Edep_05_08_2022_ALPHA.root";
    FileName[0] = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST/ELECTRON_2cm_0_um.root";
    FileName[1] = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST/PROTON_2cm_0_um.root";
    FileName[2] = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/DST/ALPHA_2cm_0_um.root";

    if(DebugMuons)
    FileName[3] = "01_08_2022_MUON.root";

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

    double ThetaDir[Ntot];
    double PhiDir[Ntot];
    double xHole[Ntot];
    double yHole[Ntot];

    TMarker *markerPositionHole[Ntot];
    TMarker *markerProjDir[Ntot];


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
            xHole[CopyNumber]              = -(LThin/2.) + (xDelta/2.) + ix * xDelta;
            yHole[CopyNumber]              = -(LThin/2.) + (yDelta/2.) + iy * yDelta;
            PhiDir[CopyNumber]             = std::atan2((iy-1.5),(ix-1.5));
            distanceR                      = std::sqrt((iy-1.5)*(iy-1.5) + (ix-1.5)*(ix-1.5));
            ThetaDir[CopyNumber]           = distanceR * thetaMax / RMax;
            markerPositionHole[CopyNumber] = new TMarker(xHole[CopyNumber],yHole[CopyNumber],20);
            markerPositionHole[CopyNumber] -> SetMarkerColor(kRed);
            markerPositionHole[CopyNumber] -> SetMarkerSize(1.2);
            markerProjDir[CopyNumber]      = new TMarker(TMath::Cos(PhiDir[CopyNumber])*ThetaDir[CopyNumber],TMath::Sin(PhiDir[CopyNumber])*ThetaDir[CopyNumber],20);
            markerProjDir[CopyNumber]      -> SetMarkerColor(kRed);
            markerProjDir[CopyNumber]      -> SetMarkerSize(1.2);

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
                TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[9].Data(), BranchName[9].Data(),E_th_plastic);
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
    TString ConditionEnergyDispersionNoPS;
    TString ConditionDrilledVeto;
    TString ConditionGoodEvents;
    TString ConditionGoodEventsNoPS;
    TString ConditionGoodEventsSinglePair[Ntot];
    TString ConditionGoodEventsSinglePairNoPS[Ntot];


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
    ConditionEnergyDispersion += Form("- (%s)", BranchName[9].Data());
    ConditionEnergyDispersion += Form(") < %g", E_th_Vetoed);
    

    ConditionEnergyDispersionNoPS  = Form("((%s)" , BranchName[1].Data());
    ConditionEnergyDispersionNoPS += Form("- (%s)", ThinTot.Data());
    ConditionEnergyDispersionNoPS += Form("- (%s)", ThickTot.Data());
    ConditionEnergyDispersionNoPS += Form(") < %g", E_th_Vetoed);
    


    ConditionGoodEvents        = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionDrilledVeto.Data(),     "(Ed_Lateral0 < 0.1) && (Ed_Lateral1 < 0.1) && (Ed_Lateral2 < 0.1) && (Ed_Lateral3 < 0.1) && (Ed_DrilledVeto < 0.1)");
    ConditionGoodEventsNoPS    = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionDrilledVeto.Data(),     "(Ed_Lateral0 < 0.1) && (Ed_Lateral1 < 0.1) && (Ed_Lateral2 < 0.1) && (Ed_Lateral3 < 0.1) && (Ed_DrilledVeto < 0.1)");
    
    for(int i = 0; i < Ntot; ++i)
    {
        cout << "Defining good events for pair" << BranchName[11+i].Data() << " & " <<BranchName[11+ Ntot+i].Data() << endl;
        ConditionGoodEventsSinglePair[i] = Form("(%s) && (%s) && (%s)", ConditionPairSilicon[i].Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
        ConditionGoodEventsSinglePairNoPS[i] = Form("(%s) && (%s) && (%s)", ConditionPairSilicon[i].Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersionNoPS.Data());
    }


    TCanvas *c1[2]; 
    TString HistName;
    TString HistTitle;
    c1[0] = new TCanvas("","PID Vs E_{tot} with Plastic Scintillator", 2480, 1748);
    cout << "Making Canvas c1" << endl;
    
    Emaxx       = 3;
    Emaxy       = 3.4;
    Eminx       = -1.4;
    Eminy       = -3.2;
    Nbinsx      = 400;
    Nbinsy      = 400;

    TH2D *hPID = new TH2D("hPID", "PID Vs E_{tot}", Nbinsx, Eminx, Emaxx, Nbinsy, Eminy, Emaxy);
    hPID -> SetTitle(" ; log_{10} (E_{tot} / 1 MeV)  [];PID_{proxy} []");
    hPID -> SetStats(0);

    for(int k = 0; k < Nfiles; ++k)
    {
        if(k == 0)
        {
            Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>%s", PIDName[CopyNoToPlot].Data(), TotEnergyName[CopyNoToPlot].Data(), hPID->GetName()), ConditionGoodEventsSinglePair[CopyNoToPlot].Data(),"colz");
            gPad -> SetGrid();
            gPad -> SetLogz();
        }
        else 
        {
            Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>+%s", PIDName[CopyNoToPlot].Data(), TotEnergyName[CopyNoToPlot].Data(), hPID->GetName()), ConditionGoodEventsSinglePair[CopyNoToPlot].Data(),"colz");
        }
    }
    c1[0] -> SaveAs("../docs/assets/PID/EtotPID_PS.png");
    c1[0] -> SaveAs("../docs/assets/PID/EtotPID_PS.pdf");
    c1[0] -> SaveAs("../docs/assets/PID/EtotPID_PS.root");

    c1[1] = new TCanvas("","PID Vs E_{tot} no Plastic Scintillator", 2480, 1748);
    cout << "Making Canvas c1" << endl;

    Emaxx       = 1.8;
    Emaxy       = 3.4;
    Eminx       = -1.4;
    Eminy       = -3.2;
    Nbinsx      = 400;
    Nbinsy      = 400;

    hPID = new TH2D("hPID", "PID Vs E_{tot}", Nbinsx, Eminx, Emaxx, Nbinsy, Eminy, Emaxy);
    hPID -> SetTitle(" ; log_{10} (E_{tot} / 1 MeV)  [];PID_{proxy} []");
    hPID -> SetStats(0);

    for(int k = 0; k < Nfiles; ++k)
    {
        if(k == 0)
        {
            Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>%s", PIDName[CopyNoToPlot].Data(), TotEnergyName[CopyNoToPlot].Data(), hPID->GetName()), ConditionGoodEventsSinglePairNoPS[CopyNoToPlot].Data(),"colz");
            gPad -> SetGrid();
            gPad -> SetLogz();
        }
        else 
        {
            Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>+%s", PIDName[CopyNoToPlot].Data(), TotEnergyName[CopyNoToPlot].Data(), hPID->GetName()), ConditionGoodEventsSinglePairNoPS[CopyNoToPlot].Data(),"colz");
        }
    }

    c1[1] -> SaveAs("../docs/assets/PID/EtotPID_noPS.png");
    c1[1] -> SaveAs("../docs/assets/PID/EtotPID_noPS.pdf");
    c1[1] -> SaveAs("../docs/assets/PID/EtotPID_noPS.root");


    return;
}