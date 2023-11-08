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

#define NSIZEMAX 100
//#define FINALLINE 241908

using namespace std;


void AnalysisSimulation2()
{
    /* -------------------------------------------------------------------------- */
    /*                                   OPTIONS                                  */
    /* -------------------------------------------------------------------------- */

    bool OnlyEnergy = true;
    bool NoEnergy = false;
    bool AnglePicture = false;
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
    double ResPlastic  = 0.3;
    double E_min_thin  = 0.01;  // Thick layer
    double E_min_thick = 0.03; // Thin layer
    double E_th_Vetoed = 0.1;  // Energy dispersion (Veto threshold)
    double E_th_plastic = 0.2;


    double Emaxx = 3;
    double Emaxy = 3.4;
    double Eminx = -1.4;
    double Eminy = -3.2;
    double Nbinsx = 400;
    double Nbinsy = 400;


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
    BranchName[43] = "Ed_Veto1";
    BranchName[44] = "Ed_Veto2";
    BranchName[45] = "Ed_Veto3";
    BranchName[46] = "Ed_Veto4";

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
        cout << "Defining total energy in File " << i << " : " << FileName[i] << endl;
        for(int ix = 0; ix < Nx; ++ix)
        {
            for(int iy = 0; iy < Ny; ++iy)
            {   
                TotEnergyName[CopyNumber] = Form("Tot_x%d_y%d_ID%d",ix,iy,CopyNumber);
                PIDName[CopyNumber] = Form("PID_x%d_y%d_ID%d",ix,iy,CopyNumber);
                TotEnergyCondition[CopyNumber]  = Form("(%s + %s",BranchName[j].Data(), BranchName[j + Ntot].Data());
                TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[9].Data(), BranchName[9].Data(),E_th_plastic);
                TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[43].Data(), BranchName[43].Data(),E_th_plastic);
                TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[44].Data(), BranchName[44].Data(),E_th_plastic);
                TotEnergyCondition[CopyNumber] += Form("+ (g%s)*(g%s > %g)", BranchName[45].Data(), BranchName[45].Data(),E_th_plastic);
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

    TString SumE;

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
    //ConditionEnergyDispersion  = Form("((%s) - (%s) - (%s) - (%s) - (%s) - (%s) - (%s)) < %g", BranchName[1].Data(), ThinTot.Data(), ThickTot.Data(), BranchName[9].Data(),BranchName[43].Data(), BranchName[44].Data(), BranchName[45].Data(), E_th_Vetoed);
    ConditionEnergyDispersion  = Form("((%s)" , BranchName[1].Data());
    ConditionEnergyDispersion += Form("- (%s)", ThinTot.Data());
    ConditionEnergyDispersion += Form("- (%s)", ThickTot.Data());
    ConditionEnergyDispersion += Form("- (%s)", BranchName[9].Data());
    ConditionEnergyDispersion += Form("- (%s)", BranchName[43].Data());
    ConditionEnergyDispersion += Form("- (%s)", BranchName[44].Data());
    ConditionEnergyDispersion += Form("- (%s)", BranchName[45].Data());
    ConditionEnergyDispersion += Form(") < %g", E_th_Vetoed);
    
    ConditionGoodEvents        = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
    for(int i = 0; i < Ntot; ++i)
    {
        cout << "Defining good events for pair" << BranchName[11+i].Data() << " & " <<BranchName[11+ Ntot+i].Data() << endl;
        ConditionGoodEventsSinglePair[i] = Form("(%s) && (%s) && (%s)", ConditionPairSilicon[i].Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
    }




    for(int i = 0; i < Ntot; ++i)
    {
        
        if(i == 0)
        {
            SumE = Form("%s", TotEnergyName[i].Data());
        }
        else
        {
            SumE = Form("%s + %s", SumE.Data() ,TotEnergyName[i].Data());
        }
    }

    SumE = Form("(%s + %s + %s + %s + %s)", SumE.Data(), BranchName[9].Data(),BranchName[43].Data(),BranchName[44].Data(), BranchName[45].Data() );
    


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

        if(OnlyEnergy)
        {
            return;
        }
    }

    
    /* -------------------------------------------------------------------------- */
    /*                    Acceptance histogram of the detector                    */
    /* -------------------------------------------------------------------------- */

    TH1D    *hAcceptance[3];
    TH1D    *hFinal[3];
    TH1D    *hIncident[3];

    TCanvas *cAcceptance[3];
    TCanvas *cIncident[3];
    TCanvas *cFinal[3];
    double ERange[3][2];
    double ERangeGen[3][2];
    double dNdEGen;
    double NgenCalibrated;


    ERangeGen[0][0] = 0.001;
    ERangeGen[0][1] = 1.;
    ERangeGen[1][0] = 0.01;
    ERangeGen[1][1] = 10.;
    ERangeGen[2][0] = 0.01;
    ERangeGen[2][1] = 50.;


    int Nbins[3] = {100, 100, 100};

    ERange[0][0] = 0.001;
    ERange[0][1] = 1.;
    ERange[1][0] = 0.01;
    ERange[1][1] = 10.;
    ERange[2][0] = 0.01;
    ERange[2][1] = 50.;

    double LActive = 2 * 6.1;
    double Radius = 0.005;
    int NGen = 10000000;
    //         100000000


    /*
    for(int i = 0 ; i < 3 ; ++i)
    {
        hIncident[i]   = new TH1D(Form("hIncident_%d", i), Form("Incident energy %d", i), Nbins[i], ERange[i][0], ERange[i][1]);
        hAcceptance[i] = new TH1D(Form("hAcceptance_%d", i), Form("Acceptance %d", i), Nbins[i], ERange[i][0], ERange[i][1]);
        
        cIncident[i] = new TCanvas(Form("cIncident_%d", i), Form("Incident energy %d", i), 1000, 900);
        gPad -> SetGrid();
        gPad -> SetLogy();
        hIncident[i] -> SetTitle(Form("Incident energy %d: Log10(E[MeV]): Counts", i));
        Edep[i] -> Draw(Form("%s>>hIncident_%d", BranchName[1].Data(), i), "", "");
        

        cAcceptance[i] = new TCanvas(Form("Acceptance %d", i), Form("Accepted events %d", i), 900, 700);
        gPad -> SetGrid();
        gPad -> SetLogy();
        hAcceptance[i] -> SetTitle(Form("Accepted events %d: Log10(E[MeV]): Counts", i));
        Edep[i] -> Draw(Form("%s>>hAcceptance_%d", BranchName[1].Data(), i), ConditionGoodEvents.Data(), "");
        //hAcceptance[i] -> Sumw2();

        NgenCalibrated = NGen * (ERange[i][1] - ERange[i][0])/ (ERangeGen[i][1] - ERangeGen[i][0]);;
        dNdEGen = NGen * (hAcceptance[i] -> GetBinWidth(1))/ (ERangeGen[i][1] - ERangeGen[i][0]);
        
        cFinal[i] = new TCanvas(Form("cFinal_%d", i), Form("Final energy %d", i), 900, 700);
        hFinal[i] = (TH1D*) hAcceptance[i] -> Clone(Form("hFinal_%d", i));
        hFinal[i] -> Scale(TMath::Pi() * LActive * LActive / dNdEGen);
        hFinal[i] -> SetTitle(Form("Geometrical factor %d; Energy [MeV]; Geometrical Factor [cm^2 rad]",i));
        hFinal[i] -> SetStats(0);
        hFinal[i] -> SetMarkerStyle(8);
        hFinal[i] -> SetMarkerSize(0.8);
        hFinal[i] -> SetMarkerColor(kBlue);
        hFinal[i] -> SetLineWidth(2);
        hFinal[i] -> Draw();
        gPad -> SetGrid();
        gPad -> SetLogy();
        cout << "#############################" << endl;
        cout << "Geometrical factor " << i << endl;
        cFinal[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Acceptance%d.png",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Acceptance%d.pdf",i));
    }
    */


    TCanvas *cAngles[Nfiles];
    TH2D *hAngles[Nfiles];

/* -------------------------------------------------------------------------- */
/*     Definisco i TGraph per gli angoli in modo da avere diversi grafici     */
/* -------------------------------------------------------------------------- */

    TGraph *gAngles[Nfiles][16];
    TMultiGraph *mgAngles[Nfiles];
    TMultiGraph *mgAnglesProtonsAlpha;
    TCanvas *CAnglesProtonsAlpha;
    TFile *tProjFile[Nfiles];
    TTree *tProj[Nfiles];
    TCanvas *cAngleGraphs[3];
    int ColorsPlot[16] = {1, 2, 4, 5, 4, 5, 1, 2, 1, 2, 4, 5, 4, 5, 1, 2 };
    TLine *lineAngles[10];

    for(int i = 0; i < 5; ++i)
    {
        lineAngles[i] = new TLine(-30 + i*15, -40 , -30 + i*15, 40);
        lineAngles[i] -> SetLineColor(1);
        lineAngles[i] -> SetLineWidth(3);
        lineAngles[i] -> SetLineStyle(4);
        lineAngles[i] -> SetLineColor(12);
    }

    for(int i = 5; i < 10; ++i)
    {
        lineAngles[i] = new TLine(-40, -30 + (i-5)*15 ,40, -30 + (i-5)*15);
        lineAngles[i] -> SetLineColor(1);
        lineAngles[i] -> SetLineWidth(3);
        lineAngles[i] -> SetLineStyle(4);
        lineAngles[i] -> SetLineColor(12);
    }

    double projXAngle, projYAngle;
    double pDirX, pDirY, pDirZ;
    double Norm;
    double nDirX, nDirY, nDirZ;
    double Theta, Phi;

    mgAnglesProtonsAlpha = new TMultiGraph();
  

    if(AnglePicture)
    {
        for(int i = 0; i < Nfiles; ++i)
        {
            cAngleGraphs[i] = new TCanvas(Form("Directions: %d", i),Form("Directions: %d", i),900,700);
            cAngleGraphs[i] -> cd();

            tProjFile[i] = new TFile(Form("../ROOT_macros/Tree%d.root",i),"RECREATE");
            tProj[i] = new TTree(Form("Proj_%d", i),Form("Proj_%d", i));

            tProj[i] -> Branch("ProjX", &projXAngle, "ProjX/D");
            tProj[i] -> Branch("ProjY", &projYAngle, "ProjY/D");

            Edep[i] -> SetBranchAddress(BranchName[6].Data(), &pDirX);
            Edep[i] -> SetBranchAddress(BranchName[7].Data(), &pDirY);
            Edep[i] -> SetBranchAddress(BranchName[8].Data(), &pDirZ);

            mgAngles[i] = new TMultiGraph();

            for(int k = 0; k < (Edep[i] -> GetEntries()); ++k)
            {
                Edep[i] -> GetEntry(k);
                Norm= TMath::Sqrt(pDirX*pDirX + pDirY*pDirY + pDirZ*pDirZ);
                nDirX = pDirX/Norm;
                nDirY = pDirY/Norm;
                nDirZ = pDirZ/Norm;

                Theta = TMath::ACos(nDirZ);
                Phi   = TMath::ATan2(nDirY, nDirX);

                projXAngle = Theta * TMath::Cos(Phi) * (180/TMath::Pi());
                projYAngle = Theta * TMath::Sin(Phi) * (180/TMath::Pi());

                tProj[i] -> Fill();
            }

            Edep[i] -> AddFriend(tProj[i]);
            tProjFile[i] -> Write();
            

            for(CopyNumber = 0 ; CopyNumber < Ntot; ++CopyNumber)
            {
                
                Edep[i] -> Draw("ProjY:ProjX", ConditionGoodEventsSinglePair[CopyNumber].Data(), "");
                gAngles[i][CopyNumber] = new TGraph(Edep[i]->GetSelectedRows(), Edep[i]->GetV2(), Edep[i]->GetV1());
                gAngles[i][CopyNumber] -> SetMarkerColor(ColorsPlot[CopyNumber]);
                gAngles[i][CopyNumber] -> SetMarkerStyle(8);
                gAngles[i][CopyNumber] -> SetMarkerSize(0.5);
                mgAngles[i] -> Add(gAngles[i][CopyNumber]);
                if(i > 0)
                {
                    mgAnglesProtonsAlpha -> Add(gAngles[i][CopyNumber]);
                }
            }

            //gPad -> SetGrid();
            mgAngles[i] -> GetXaxis() -> SetTitle("Angle projection X [deg]");
            mgAngles[i] -> GetYaxis() -> SetTitle("Angle projection Y [deg]");

            mgAngles[i] -> GetXaxis() -> SetRangeUser(-60, 60);
            mgAngles[i] -> GetYaxis() -> SetRangeUser(-60, 60);

            mgAngles[i] -> Draw("AP");

            for(int q = 0; q < 10; ++q)
            {
                lineAngles[q] -> Draw("same");
            }

            cAngleGraphs[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Angles_%d.pdf", i));
            cAngleGraphs[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Angles_%d.svg", i));
            cAngleGraphs[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Angles_%d.png", i));
        }
        

        CAnglesProtonsAlpha = new TCanvas("AnglesProtonsAlpha", "AnglesProtonsAlpha", 900, 700);
        CAnglesProtonsAlpha -> cd();
        mgAnglesProtonsAlpha -> GetXaxis() -> SetTitle("Angle projection X [deg]");
        mgAnglesProtonsAlpha -> GetYaxis() -> SetTitle("Angle projection Y [deg]");

        mgAnglesProtonsAlpha -> GetXaxis() -> SetRangeUser(-60, 60);
        mgAnglesProtonsAlpha -> GetYaxis() -> SetRangeUser(-60, 60);

        mgAnglesProtonsAlpha -> Draw("AP");
        //gPad -> SetGrid();

        for(int q = 0; q < 10; ++q)
        {
            lineAngles[q] -> Draw("same");
        }

        CAnglesProtonsAlpha -> SaveAs("../docs/assets/AnalysisSimulation2/AnglesProtonsAlpha.pdf");
        CAnglesProtonsAlpha -> SaveAs("../docs/assets/AnalysisSimulation2/AnglesProtonsAlpha.svg");
        CAnglesProtonsAlpha -> SaveAs("../docs/assets/AnalysisSimulation2/AnglesProtonsAlpha.png");
    }
    


    /*
    for(int i = 0 ; i < Nfiles; ++i)
    {
        cAngles[i] = new TCanvas(Form("Angles %d", i), Form("Angles %d", i), 900, 700);
        hAngles[i] = new TH2D(Form("hAngles_%d", i), Form("Angles %d", i), 60, -80, 80, 60, -80, 80);
        Edep[i] -> Draw(Form("((%s)*TMath::Sin(%s)):((%s)*TMath::Cos(%s))>>%s", PolarAngle[0].Data(),PolarAngle[1].Data(),PolarAngle[0].Data(),PolarAngle[1].Data(), Form("hAngles_%d", i)), MaskTotalChip.Data(), "colz");
        gPad -> SetGrid();
        // Draw all the markers Dir
        for(CopyNumber = 0; CopyNumber < Ntot ; ++CopyNumber)
        {
            markerProjDir[CopyNumber] -> Draw("same");
        }
        //gPad -> SetLogz();
    // && ((%s ==%d) || (%s ==%d)) BranchName[0].Data(),0, BranchName[0].Data(),15
    }



    TH2D *hmap[Nfiles];
    TCanvas *c2;
    c2 = new TCanvas("c2", "Acceptance direction map No filtering", 600, 800);
    c2 -> Divide(1,Nfiles);

    for(int i = 0; i < Nfiles ; ++i)
    {
        c2 -> cd(i+1);
        //hmap[i] = new TH2D(Form("hmap%d",i), Form("Acceptance direction %d; Phi; Theta",i), 30, -200, 200, 60 ,-10.0, 60.0 );   
        //Edep[i] -> Draw(Form("%s:%s>>hmap%d",PolarAngle[0].Data(), PolarAngle[1].Data(), i), MaskTotal ,"colz aitoff");
        hmap[i] = new TH2D(Form("hmap%d",i), Form("Acceptance direction X Y No filtering%d; X generation; Y generation",i), 200, -100, 100, 200 ,-100, 100);   
        Edep[i] -> Draw(Form("%s:%s>>hmap%d",BranchName[4].Data(), BranchName[3].Data(), i), "" ,"colz");
        gPad -> SetGrid();
    }

    
    
    c2 -> SaveAs("../docs/assets/AnalysisSimulation2/SkyMapPositionNoMask.png");
    c2 -> SaveAs("../docs/assets/AnalysisSimulation2/SkyMapPositionNoMask.pdf");
    
    


    TH2D *hmapM[Nfiles];
    TCanvas *c3;
    c3 = new TCanvas("c3", "Acceptance direction map with filtering", 600, 800);
    c3 -> Divide(1,Nfiles);



    for(int i = 0; i < Nfiles ; ++i)
    {
        c3 -> cd(i+1);
        //hmap[i] = new TH2D(Form("hmap%d",i), Form("Acceptance direction %d; Phi; Theta",i), 30, -200, 200, 60 ,-10.0, 60.0 );   
        //Edep[i] -> Draw(Form("%s:%s>>hmap%d",PolarAngle[0].Data(), PolarAngle[1].Data(), i), MaskTotal ,"colz aitoff");
        hmapM[i] = new TH2D(Form("hmapM%d",i), Form("Acceptance direction X Y With filtering %d; X generation; Y generation",i), 200, -100, 100, 200 ,-100, 100);   
        Edep[i] -> Draw(Form("%s:%s>>hmapM%d",BranchName[4].Data(), BranchName[3].Data(), i), MaskTotal ,"colz");
        gPad -> SetGrid();

    }

    if(DebugMuons)
    {
        c3 -> SaveAs("../docs/assets/AnalysisSimulation2/SkyMapPosition_muons.png");
        c3 -> SaveAs("../docs/assets/AnalysisSimulation2/SkyMapPosition_muons.pdf");
    } 
    else 
    {
        c3 -> SaveAs("../docs/assets/AnalysisSimulation2/SkyMapPosition.png");
        c3 -> SaveAs("../docs/assets/AnalysisSimulation2/SkyMapPosition.pdf");
    }
    */


    TH1D *htheta [Nfiles][Ntot];
    TCanvas *ctheta[Nfiles];

    TH1D *hphi [Nfiles][Ntot];
    TCanvas *cphi[Nfiles];

    

    for(int i = 0; i < Nfiles ; ++i)
    {
        cphi[i] = new TCanvas(Form("cphi%d",i), Form("Phi accepted %d", i), 3508 ,2480);
        cphi[i] -> Divide(Ny, Nx);
        ctheta[i] = new TCanvas(Form("ctheta%d",i), Form("Theta accepted %d", i), 3508 ,2480);
        ctheta[i] -> Divide(Ny, Nx);
        for(int k = 0 ; k < Ntot ; ++k)
        {
            cphi[i] -> cd(k+1);
            HistName = Form("hphi%d%d",i,k);
            HistTitle = Form("#phi Accepted Number: %d; #phi [deg]; counts []", k);
            hphi[i][k] = new TH1D(HistName, HistTitle, 160, -180, 180);
            hphi[i][k] -> SetStats(0);
            hphi[i][k] -> SetLineColor(kBlue);
            hphi[i][k] -> SetLineWidth(1.8);
            // Enlarge font size
            hphi[i][k] -> SetTitleSize(0.08, "XYZ");
            hphi[i][k] -> SetLabelSize(0.08, "XYZ");
            hphi[i][k] -> SetTitleOffset(1.2, "Y");
            hphi[i][k] -> SetTitleOffset(1.2, "X");

            
            if(i == 3)
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[1].Data(), HistName.Data()), ConditionGoodEventsSinglePair[k].Data(), "");
            }
            else
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[1].Data(), HistName.Data()), ConditionGoodEventsSinglePair[k].Data(), "");
            }

            ctheta[i] -> cd(k+1);
            HistName = Form("htheta%d%d",i,k);
            HistTitle = Form("#theta Accepted Number: %d; #theta [deg]; counts []", k);
            htheta[i][k] = new TH1D(HistName, HistTitle, 70, 0, 70);
            htheta[i][k] -> SetStats(0);
            htheta[i][k] -> SetLineColor(kRed);
            htheta[i][k] -> SetLineWidth(1.8);
            // Enlarge font size
            htheta[i][k] -> SetTitleSize(0.08, "XYZ");
            htheta[i][k] -> SetLabelSize(0.08, "XYZ");
            htheta[i][k] -> SetTitleOffset(1.2, "Y");
            htheta[i][k] -> SetTitleOffset(1.2, "X");
            if(i == 3)
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[0].Data(), HistName.Data()), ConditionGoodEventsSinglePair[k].Data(), "");
            }
            else
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[0].Data(), HistName.Data()), ConditionGoodEventsSinglePair[k].Data(), "");
            }
            
        }
        cphi[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Phi%d.png", i));
        cphi[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Phi%d.pdf", i));

        ctheta[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Theta%d.png", i));
        ctheta[i] -> SaveAs(Form("../docs/assets/AnalysisSimulation2/Theta%d.pdf", i));
    }


    return;
}