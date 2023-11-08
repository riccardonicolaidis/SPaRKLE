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
#include "TGraphErrors.h"


#define NSIZEMAX 100
//#define FINALLINE 241908

using namespace std;


void AnalysisSimulation3()
{
    /* -------------------------------------------------------------------------- */
    /*                                   OPTIONS                                  */
    /* -------------------------------------------------------------------------- */

    bool NoEnergy = true;
    int Nfiles = 3;


    cout << "Nfiles = " << Nfiles << endl;


    /* -------------------------------------------------------------------------- */
    /*                                 /FILE NAMES                                */
    /* -------------------------------------------------------------------------- */

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


    for(int i = 0; i < Nfiles; i++)
    {
        cout << "FileName[" << i << "] = " << FileName[i] << endl;
    }

    TFile *file[Nfiles];
    TTree *Edep[Nfiles];


    for(int i = 0; i < Nfiles ; ++i)
    {
        file[i] = TFile::Open(FileName[i]);
        Edep[i] = (TTree*) file[i] -> Get("Edep");
        Edep[i] -> Print();
    }



    /* -------------------------------------------------------------------------- */
    /*                                BRANCH NAMES                                */
    /* -------------------------------------------------------------------------- */

    TString BranchName[47];
    TString TotEnergyName[16];
    TString TotEnergyCondition[16];
    TString PIDName[16];
    TString DirName[3];
    TString PolarAngle[2];
    TString NewPolarAngle[2];
    TString AliasETot = "ETot";

    TString ThinTot;
    TString ThickTot;

    DirName[0] = "DirX";
    DirName[1] = "DirY";
    DirName[2] = "DirZ";

    PolarAngle[0] = "dirTheta";
    PolarAngle[1] = "dirPhi";

    NewPolarAngle[0] = "dirThetaNew";
    NewPolarAngle[1] = "dirPhiNew";

    int Nx;
    int Ny;
    int Ntot; 
    int j;
    int CopyNumber;
    int VERBOSE = 7;


    /* -------------------------------------------------------------------------- */
    /*                           GEOMETRICAL DEFINITIONS                          */
    /* -------------------------------------------------------------------------- */
    Nx = 4;
    Ny = 4;
    Ntot = Nx*Ny; 

    double ThetaDir[Ntot];
    double PhiDir[Ntot];
    double xHole[Ntot];
    double yHole[Ntot];
    TMarker *markerPositionHole[Ntot];
    TMarker *markerProjDir[Ntot];

    double LThin = 10.;
    double xDelta = LThin / (Nx);
    double yDelta = LThin / (Ny);
    double thetaMax = 35. * TMath::Pi() / 180.;
    double RMax = std::sqrt(2) * 1.5;
   
    double distanceR;

    /* -------------------------------------------------------------------------- */
    /*                             PLOTTING PARAMETERS                            */
    /* -------------------------------------------------------------------------- */

    double ResSilicon  = 0.01;
    double ResCZT      = 0.05;
    double ResPlastic  = 0.3;
    double E_min_thin  = 0.01;  // Thick layer
    double E_min_thick = 0.03; // Thin layer
    double E_th_Vetoed = 0.1;  // Energy dispersion (Veto threshold)
    double E_th_plastic = 0.2;


    double Emaxx = 3;
    double Emaxy = 3.4;
    double Eminx = -1.4;
    double Eminy = -3.2;
    double Nbinsx = 600;
    double Nbinsy = 600;


/* -------------------------------------------------------------------------- */
/*               Constructiong the Branch names to read the tree              */
/* -------------------------------------------------------------------------- */


    BranchName[0] = "NumberID";
    BranchName[1] = "RandEnergy";
    BranchName[2] = "RandNumber";
    BranchName[3] = "Xgen";
    BranchName[4] = "Ygen";
    BranchName[5] = "Zgen";
    BranchName[6] = "pDirX";
    BranchName[7] = "pDirY";
    BranchName[8] = "pDirZ";

    BranchName[9] = "Ed_Veto0";
    
    BranchName[10] = "Ed_DrilledVeto";
    BranchName[43] = "Ed_Lateral0";
    BranchName[44] = "Ed_lateral1";
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
        cout << i << " " << BranchName[i] << endl;
    }


    /* -------------------------------------------------------------------------- */
    /*                              ALIAS DEFINITIONS                             */
    /* -------------------------------------------------------------------------- */

    for(int i = 0; i < Nfiles ; ++i)
    {
        cout << "Setting alias in File " << i << " : " << FileName[i] << endl;
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
                Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResCZT));
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
        cout << "Defining total energy in File " << i << " : " << FileName[i] << endl;
        for(int ix = 0; ix < Nx; ++ix)
        {
            for(int iy = 0; iy < Ny; ++iy)
            {   
                TotEnergyName[CopyNumber] = Form("Tot_x%d_y%d_ID%d",ix,iy,CopyNumber);
                PIDName[CopyNumber] = Form("PID_x%d_y%d_ID%d",ix,iy,CopyNumber);
                TotEnergyCondition[CopyNumber]  = Form("(%s + %s",BranchName[j].Data(), BranchName[j + Ntot].Data());
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
/*                         Figure definition and plots                        */
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */
/*                   Generation of the selection mask string                  */
/* -------------------------------------------------------------------------- */

    cout << "Image generation" << endl;    
    if(NoEnergy)
    {
        cout << "Skipping energy plots" << endl;
    }

    TCanvas *c1;
    if(!NoEnergy)
    {
        c1 = new TCanvas("","Si 100 um Si 500 um", 1600, 900); 
    }
    TCanvas *call[Ntot];
    TCanvas *call2[Ntot];
    TH2D *h2[Ntot];
    TH2D *h2lin[Ntot];
    TString HistName;
    TString HistTitle;
    TString ConditionPairSilicon[Ntot];
    TString ConditionPairSiliconAll;
    TString ConditionEnergyDispersion;
    TString ConditionDrilledVeto;
    TString ConditionGoodEvents;
    TString ConditionGoodEventsSinglePair[Ntot];




    if(!NoEnergy)
    {
        c1 -> Divide(Ny, Nx);
    }


    /* -------------------------------------------------------------------------- */
    /*                          DEFINITION OF GOOD EVENTS                         */
    /* -------------------------------------------------------------------------- */
    
    cout << "Defining good events" << endl;

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
    ConditionEnergyDispersion  = Form("((%s)" , BranchName[1].Data());
    ConditionEnergyDispersion += Form("- (%s)", ThinTot.Data());
    ConditionEnergyDispersion += Form("- (%s)", ThickTot.Data());
    ConditionEnergyDispersion += Form(") < %g", E_th_Vetoed);
    
    ConditionGoodEvents        = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
    for(int i = 0; i < Ntot; ++i)
    {
        cout << "Defining good events for pair" << BranchName[11+i].Data() << " & " <<BranchName[11+ Ntot+i].Data() << endl;
        ConditionGoodEventsSinglePair[i] = Form("(%s) && (%s) && (%s)", ConditionPairSilicon[i].Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
    }



    TCanvas *cPIDAll = new TCanvas("cPIDAll", "PID All", 1600, 900);
    TH2D *hPIDAll = new TH2D("hPIDAll", "PID All; Log 10 Energy [Log10(MeV)]; PID [] ",Nbinsx, Eminx ,Emaxx ,Nbinsy,Eminy, Emaxy);
    for(int i = 0; i<3; ++i)
    {
        for(int j = 0; j < Ntot; ++j)
        {
            if(i == 0 && j == 0)
            {
                Edep[i] -> Draw(Form("%s:TMath::Log10(%s)>>hPIDAll", PIDName[j].Data() , TotEnergyName[j].Data()), ConditionGoodEventsSinglePair[j].Data(), "colz");
            }
            else
            {
                Edep[i] -> Draw(Form("%s:TMath::Log10(%s)>>+hPIDAll", PIDName[j].Data() , TotEnergyName[j].Data()), ConditionGoodEventsSinglePair[j].Data(), "colz");
            }

            cout << "Processing file " << i << "Detector " << j << endl;
        }
        
    }

    gPad -> SetGrid();
    gPad -> SetLogz();



    if(!NoEnergy)
    {
        for(int i = 0; i < 16; ++i)
        {
            if(!NoEnergy)
            {
                c1 -> cd(i+1);
            }
            
            HistName = Form("h%d",i);
            HistTitle = Form("Copy number %d; Log 10 Energy [Log10(MeV)]; PID [] ",i);
            h2[i] = new TH2D(HistName, HistTitle ,Nbinsx, Eminx ,Emaxx ,Nbinsy,Eminy, Emaxy);
            for(int k = 0; k < Nfiles; ++k)
            {
                if(k == 0)
                {
                    Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>%s", PIDName[i].Data() , TotEnergyName[i].Data(), HistName.Data()), ConditionGoodEventsSinglePair[i].Data(), "colz");
                    gPad -> SetGrid();
                    //gPad -> SetLogz();
                    h2[i] -> SetStats(0);
                }
                else
                {
                    Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>+%s", PIDName[i].Data() , TotEnergyName[i].Data(), HistName.Data()), ConditionGoodEventsSinglePair[i].Data(), "colz");
                
                }
            }

            if(!NoEnergy)
            {
                cout << "Drawing  " << i << endl;
                call[i] = new TCanvas("","Si 100 um Si 500 um", 1000, 900);
                call[i] -> cd();
                h2[i] -> Draw("colz");        
                gPad -> SetGrid();
                gPad -> SetLogz();

                call2[i] = new TCanvas("","Si 100 um Si 500 um", 1000, 900);
                call2[i] -> cd();
                HistName = Form("hlin%d",i);
                HistTitle = Form("Copy number %d; Energy [MeV]; Delta E ",i);
                h2lin[i] = new TH2D(HistName, HistTitle ,200, -2.5 ,3.5 ,200,-2.,3.);
                for(int k = 0; k < Nfiles; ++k)
                {
                    if(k == 0)
                    {
                        Edep[k] -> Draw(Form("TMath::Log10(%s):TMath::Log10(%s)>>%s", ThinTot.Data() , TotEnergyName[i].Data(), HistName.Data()), ConditionGoodEventsSinglePair[i].Data(), "colz");
                        gPad -> SetGrid();
                        gPad -> SetLogz();
                        h2lin[i] -> SetStats(0);
                    }
                    else
                    {
                        
                        Edep[k] -> Draw(Form("TMath::Log10(%s):TMath::Log10(%s)>>+%s", ThinTot.Data() , TotEnergyName[i].Data(), HistName.Data()), ConditionGoodEventsSinglePair[i].Data(), "colz");
                    
                    }
                }
                
                call[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PID_%d.pdf", i));
                call[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PID_%d.png", i));
                call[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PID_%d.C", i));
                call[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PID_%d.root", i));
                call2[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PIDlin_%d.pdf", i));
                call2[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PIDlin_%d.png", i));
                call2[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PIDlin_%d.C", i));
                call2[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/PIDlin_%d.root", i));

            }
        }

        if((!NoEnergy))
        {
            c1 -> SaveAs("../docs/assets/AnalysisSimulation2/PID.pdf");
            c1 -> SaveAs("../docs/assets/AnalysisSimulation2/PID.png");
            c1 -> SaveAs("../docs/assets/AnalysisSimulation2/PID.C");
            c1 -> SaveAs("../docs/assets/AnalysisSimulation2/PID.root");
        }
    }

    

    return;
}