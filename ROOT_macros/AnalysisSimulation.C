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


void AnalysisSimulation()
{
    bool OnlyEnergy = false;
    bool DebugMuons = false;
    bool NoEnergy = true;

    int Nfiles = 3;

    if(DebugMuons)
    {
        Nfiles = 4;
    }
    /* -------------------------------------------------------------------------- */
    /*                               Style settings                               */
    /* -------------------------------------------------------------------------- */

    /*
    TStyle *myStyle = new TStyle("myStyle", "Personal style");
    myStyle -> SetPalette(1,0);
    myStyle -> SetOptStat(0);
    myStyle -> SetOptDate(0);
    myStyle -> SetLabelSize(0.03, "xyz");
    myStyle -> SetTitleFont(22,"xyz");
    myStyle -> SetLabelFont(22,"xyz");
    myStyle -> SetTitleOffset(1.2,"y");
    myStyle -> SetCanvasDefW(600);
    myStyle -> SetCanvasDefH(600);
    myStyle->SetCanvasColor(0); // canvas...
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetCanvasBorderSize(0);    
    myStyle->SetPadBottomMargin(0.1); //margins...
    myStyle->SetPadTopMargin(0.1);
    myStyle->SetPadLeftMargin(0.1);
    myStyle->SetPadRightMargin(0.1);
    myStyle->SetPadGridX(0); // grids, tickmarks
    myStyle->SetPadGridY(0);
    myStyle->SetPadTickX(1);
    myStyle->SetPadTickY(1);
    myStyle->SetFrameBorderMode(0);
    myStyle->SetPaperSize(20,24); // US letter size
    */
    
    //myStyle -> cd();


    //defaultS -> cd();


    /* ------------------------------- File names ------------------------------- */
    TString FileName[Nfiles];
    FileName[0] = "100um_500um_e.root";
    FileName[1] = "100um_500um_p.root";
    FileName[2] = "100um_500um_He.root";
    if(DebugMuons)
    {
        FileName[3] = "100um_500um_mu.root";
    }

    TFile *file[Nfiles];
    TTree *Edep[Nfiles];

    /* ------------------ Tree Branches names and labels names ------------------ */

    TString BranchName[43];
    TString TotEnergyName[16];
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

    /* ------------- Physical parameters and parameters for plotting ------------ */

    double ResSilicon  = 0.01;
    double ResPlastic  = 0.3;
    double E_min_thin  = 0.03;  // Thick layer
    double E_min_thick = 0.03; // Thin layer
    double E_th_Vetoed = 0.1;  // Energy dispersion (Veto threshold)


    double Emaxx = 2;
    double Emaxy = 3.4;
    double Eminx = -1.4;
    double Eminy = -3.2;
    double Nbinsx = 600;
    double Nbinsy = 400;


    for(int i = 0; i < Nfiles ; ++i)
    {
        file[i] = TFile::Open(FileName[i]);
        Edep[i] = (TTree*) file[i] -> Get("Edep");
        Edep[i] -> Print();
    }

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
    BranchName[9] = "Ed_Veto";
    BranchName[10] = "Ed_DrilledVeto";
    
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
    
    j = 11;
    CopyNumber = 0;
    
    for(int ix = 0; ix < Nx; ++ix)
    {
        for(int iy = 0; iy < Ny; ++iy)
        {   
            BranchName[j] = Form("Thin_x%d_y%d_ID%d",ix,iy,CopyNumber);
            BranchName[j + Ntot] = Form("Thick_x%d_y%d_ID%d",ix,iy,CopyNumber);

            if (CopyNumber == 0)
            {
                ThinTot = Form("(%s)", BranchName[j].Data());
                ThickTot = Form("(%s)", BranchName[j + Ntot].Data());
            }
            else
            {
                ThinTot = Form("(%s + %s)",ThinTot.Data(), BranchName[j].Data());
                ThickTot = Form("(%s + %s)",ThickTot.Data(), BranchName[j + Ntot].Data());
            }
            
            xHole[CopyNumber] = -(LThin/2.) + (xDelta/2.) + ix * xDelta;
            yHole[CopyNumber] = -(LThin/2.) + (yDelta/2.) + iy * yDelta;
            PhiDir[CopyNumber] = std::atan2((iy-1.5),(ix-1.5));
            distanceR = std::sqrt((iy-1.5)*(iy-1.5) + (ix-1.5)*(ix-1.5));
            ThetaDir[CopyNumber] = distanceR * thetaMax / RMax;

            markerPositionHole[CopyNumber] = new TMarker(xHole[CopyNumber],yHole[CopyNumber],20);
            markerPositionHole[CopyNumber] -> SetMarkerColor(kRed);
            markerPositionHole[CopyNumber] -> SetMarkerSize(1.2);

            markerProjDir[CopyNumber] = new TMarker(TMath::Cos(PhiDir[CopyNumber])*ThetaDir[CopyNumber]*(180/TMath::Pi()),TMath::Sin(PhiDir[CopyNumber])*ThetaDir[CopyNumber]*(180/TMath::Pi()),20);
            markerProjDir[CopyNumber] -> SetMarkerColor(kRed);
            markerProjDir[CopyNumber] -> SetMarkerSize(1.2);

            ++j;
            ++CopyNumber;
        }
    }
    

/* -------------------------------------------------------------------------- */
/*                      Verbosity levelo of the software                      */
/* -------------------------------------------------------------------------- */
    
    if(VERBOSE >= 6)
    {
        for(int i = 0; i <= 42; ++i)
        {
            cout << i << " " << BranchName[i] << endl;
        }
    }
    





    /* -------------------------------------------------------------------------- */
    /*                         Alias declaration and plots                        */
    /* -------------------------------------------------------------------------- */

    for(int i = 0; i < Nfiles ; ++i)
    {
        /* -------------------------------- Smearing -------------------------------- */
        
        for(int k = 9; k < 43; ++k)
        {
            Edep[i] -> SetAlias(Form("wnorm%d",k),"(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))");
            if(k == 9 || k == 10)
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
        Edep[i] -> SetAlias(PolarAngle[0].Data(), Form("((TMath::ACos(%s)))", DirName[2].Data())); // *(180/3.415927))
        Edep[i] -> SetAlias(PolarAngle[1].Data(), Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data())); // *(180/3.415927))
        Edep[i] -> SetAlias(NewPolarAngle[1].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(- TMath::Cos(%s))))",PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() )); // *(180/3.415927))
        Edep[i] -> SetAlias(NewPolarAngle[0].Data(),Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)), -TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data())); // *(180/3.415927))

        
        /* -------------------- Particle identification parameter ------------------- */

        j = 11;
        CopyNumber = 0;
        for(int ix = 0; ix < Nx; ++ix)
        {
            for(int iy = 0; iy < Ny; ++iy)
            {   
                TotEnergyName[CopyNumber] = Form("Tot_x%d_y%d_ID%d",ix,iy,CopyNumber);
                PIDName[CopyNumber] = Form("PID_x%d_y%d_ID%d",ix,iy,CopyNumber);
                Edep[i] -> SetAlias(TotEnergyName[CopyNumber], Form("(%s + %s)",Form("g%s",BranchName[j].Data()), Form("g%s",BranchName[j + Ntot].Data())));
                Edep[i] -> SetAlias(PIDName[CopyNumber], Form("(TMath::Log10(%s * %s))",Form("g%s",BranchName[j].Data()), TotEnergyName[CopyNumber].Data()));
                ++j;
                ++CopyNumber;
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

    
    TCanvas *c1;
    if(!NoEnergy)
    {
        c1 = new TCanvas("","Si 100 um Si 500 um", 1600, 900); 
    }
    TCanvas *call[Ntot];
    TH2D *h2[Ntot];    
    TString HistName;
    TString HistTitle;
    TString MaskVeto[Ntot];
    TString MaskTotal;
    TString MaskTotalAndThreshold;

    if(!NoEnergy)
    {
        c1 -> Divide(Ny, Nx);
    }
    

    for(int i = 0; i < Ntot; ++i)
    {
        if(!NoEnergy)
        {
            c1 -> cd(i+1);
        }
        
        HistName = Form("h%d",i);
        HistTitle = Form("Copy number %d; Log 10 Energy [Log10(MeV)]; PID []",i);
        h2[i] = new TH2D(HistName, HistTitle ,Nbinsx, Eminx ,Emaxx ,Nbinsy,Eminy, Emaxy);

        MaskVeto[i] = Form("((%s - %s - %s) <= %f) && (%s >= %f) && (%s >= %f)", BranchName[1].Data(), BranchName[11+i].Data(), BranchName[11+Ntot+i].Data(), E_th_Vetoed, BranchName[11+i].Data(),E_min_thin, BranchName[11+Ntot+i].Data(), E_min_thick);

        for(int k = 0; k < Nfiles; ++k)
        {
            if(k == 0)
            {
                Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>%s", PIDName[i].Data() , TotEnergyName[i].Data(), HistName.Data()), MaskVeto[i], "colz");
                //gPad -> SetLogx();
                gPad -> SetGrid();
                gPad -> SetLogz();
                h2[i] -> SetStats(0);
            }
            else
            {
                if((k == 3) & DebugMuons)
                {
                    Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>+%s", PIDName[i].Data() , TotEnergyName[i].Data(), HistName.Data()), "" , "colz");
                }
                else
                {
                    Edep[k] -> Draw(Form("%s:TMath::Log10(%s)>>+%s", PIDName[i].Data() , TotEnergyName[i].Data(), HistName.Data()), MaskVeto[i], "colz");
                }
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
            if(DebugMuons)
            {
                call[i] -> SaveAs(Form("../docs/assets/images/PID_%d_muons.pdf", i));
                call[i] -> SaveAs(Form("../docs/assets/images/PID_%d_muons.svg", i));
                call[i] -> SaveAs(Form("../docs/assets/images/PID_%d_muons.png", i));
            }
            else
            {
                call[i] -> SaveAs(Form("../docs/assets/images/PID_%d.pdf", i));
                call[i] -> SaveAs(Form("../docs/assets/images/PID_%d.svg", i));
                call[i] -> SaveAs(Form("../docs/assets/images/PID_%d.png", i));
            }
        }
        
        

    }

    if(DebugMuons && (!NoEnergy))
    {
        c1 -> SaveAs("../docs/assets/images/PID_total_muons.pdf");
        c1 -> SaveAs("../docs/assets/images/PID_total_muons.svg");
    }
    else if(!NoEnergy)
    {
        c1 -> SaveAs("../docs/assets/images/PID_total.pdf");
        c1 -> SaveAs("../docs/assets/images/PID_total.svg");
    }

    if(OnlyEnergy)
    {
        return;
    }

    TString SumE;
    for(int i = 0; i < Ntot ; ++i)
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
    SumE = Form("(%s)", SumE.Data());
    cout << SumE << endl;

    MaskTotal = Form("(%s - %s) <= %f", BranchName[1].Data(), SumE.Data(), E_th_Vetoed);
    MaskTotalAndThreshold = Form("(%s) && (%s >= %f) && (%s >= %f)", MaskTotal.Data(), ThinTot.Data(), E_min_thin, ThickTot.Data(), E_min_thick);

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

    ERangeGen[0][0] = 0.0001;
    ERangeGen[0][1] = 10.;
    ERangeGen[1][0] = 0.0001;
    ERangeGen[1][1] = 40.;
    ERangeGen[2][0] = 0.0001;
    ERangeGen[2][1] = 120.;


    int Nbins[3] = {40, 40, 40};

    ERange[0][0] = 0.05;
    ERange[0][1] = 2.;
    ERange[1][0] = 1.;
    ERange[1][1] = 20.;
    ERange[2][0] = 9.;
    ERange[2][1] = 50.;

    double LActive = 15.;
    double Radius = 0.005;
    int NGen = 30000000;


    for(int i = 0 ; i < 3 ; ++i)
    {
        
        hIncident[i] = new TH1D(Form("hIncident_%d", i), Form("Incident energy %d", i), Nbins[i], ERange[i][0], ERange[i][1]);
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
        Edep[i] -> Draw(Form("%s>>hAcceptance_%d", BranchName[1].Data(), i), MaskTotalAndThreshold.Data(), "");
        //hAcceptance[i] -> Sumw2();

        NgenCalibrated = NGen * (ERange[i][1] - ERange[i][0])/ (ERangeGen[i][1] - ERangeGen[i][0]);;
        dNdEGen = NGen * (hAcceptance[i] -> GetBinWidth(1))/ (ERangeGen[i][1] - ERangeGen[i][0]);
        
        cFinal[i] = new TCanvas(Form("cFinal_%d", i), Form("Final energy %d", i), 900, 700);
        hFinal[i] = (TH1D*) hAcceptance[i] -> Clone(Form("hFinal_%d", i));
        hFinal[i] -> Scale(TMath::Pi() * LActive * LActive / dNdEGen);
        hFinal[i] -> SetTitle(Form("Geometrical factor %d; Energy [MeV]; Geometrical Factor [cm^2 rad]",i));
        hFinal[i] -> Draw();
        gPad -> SetGrid();
        gPad -> SetLogy();
        cout << "#############################" << endl;
        cout << "Geometrical factor " << i << ": " << endl;
        cout << TMath::Pi() * LActive * LActive * (hAcceptance[i] -> GetEntries()) / NgenCalibrated << endl;


    }

    return;



    TCanvas *cAngles[Nfiles];
    TH2D *hAngles[Nfiles];

    TString ChipSelectionMask;
    TString MaskTotalChip;

    
    int ChipList[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

    for(int i = 0; i < 16; ++i)
    {
        if(i == 0)
        {
            ChipSelectionMask = Form("(%s == %i)", BranchName[0].Data(), ChipList[i]);
            cout << ChipSelectionMask.Data() << endl;
        } 
        else
        {
            ChipSelectionMask = Form("%s || (%s == %i)", ChipSelectionMask.Data(), BranchName[0].Data(), ChipList[i]);
            cout << ChipSelectionMask.Data() << endl;
        }
    }
    MaskTotalChip = Form("(%s) && (%s)", ChipSelectionMask.Data(), MaskTotal.Data());

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
        lineAngles[i] = new TLine(-30 + i*15, -50 , -30 + i*15, 50);
        lineAngles[i] -> SetLineColor(1);
        lineAngles[i] -> SetLineWidth(3);
        lineAngles[i] -> SetLineStyle(4);
        lineAngles[i] -> SetLineColor(12);
    }

    for(int i = 5; i < 10; ++i)
    {
        lineAngles[i] = new TLine(-50, -30 + (i-5)*15 ,50, -30 + (i-5)*15);
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
            Norm = TMath::Sqrt(pDirX*pDirX + pDirY*pDirY + pDirZ*pDirZ);
            nDirX = pDirX/Norm;
            nDirY = pDirY/Norm;
            nDirZ = pDirZ/Norm;

            Theta = TMath::ACos(nDirZ);
            Phi = TMath::ATan2(nDirY, nDirX);

            projXAngle = Theta * TMath::Cos(Phi) * (180/TMath::Pi());
            projYAngle = Theta * TMath::Sin(Phi) * (180/TMath::Pi());

            tProj[i] -> Fill();
        }

        Edep[i] -> AddFriend(tProj[i]);
        tProjFile[i] -> Write();
        

        for(CopyNumber = 0 ; CopyNumber < Ntot; ++CopyNumber)
        {
            
            Edep[i] -> Draw("ProjY:ProjX", Form("(%s == %d) && (%s)", BranchName[0].Data(), CopyNumber, MaskTotal.Data()), "");
            gAngles[i][CopyNumber] = new TGraph(Edep[i]->GetSelectedRows(), Edep[i]->GetV2(), Edep[i]->GetV1());
            gAngles[i][CopyNumber] -> SetMarkerColor(ColorsPlot[CopyNumber]);
            gAngles[i][CopyNumber] -> SetMarkerStyle(7);
            mgAngles[i] -> Add(gAngles[i][CopyNumber]);
            if(i > 0)
            {
                mgAnglesProtonsAlpha -> Add(gAngles[i][CopyNumber]);
            }
        }

        //gPad -> SetGrid();
        mgAngles[i] -> GetXaxis() -> SetTitle("Angle projection X [deg]");
        mgAngles[i] -> GetYaxis() -> SetTitle("Angle projection Y [deg]");

        mgAngles[i] -> GetXaxis() -> SetRangeUser(-50, 50);
        mgAngles[i] -> GetYaxis() -> SetRangeUser(-50, 50);

        mgAngles[i] -> Draw("AP");

        for(int q = 0; q < 10; ++q)
        {
            lineAngles[q] -> Draw("same");
        }

        cAngleGraphs[i] -> SaveAs(Form("../docs/assets/images/Angles_%d.pdf", i));
        cAngleGraphs[i] -> SaveAs(Form("../docs/assets/images/Angles_%d.svg", i));
        cAngleGraphs[i] -> SaveAs(Form("../docs/assets/images/Angles_%d.png", i));


    }
    

    CAnglesProtonsAlpha = new TCanvas("AnglesProtonsAlpha", "AnglesProtonsAlpha", 900, 700);
    CAnglesProtonsAlpha -> cd();
    mgAnglesProtonsAlpha -> GetXaxis() -> SetTitle("Angle projection X [deg]");
    mgAnglesProtonsAlpha -> GetYaxis() -> SetTitle("Angle projection Y [deg]");

    mgAnglesProtonsAlpha -> GetXaxis() -> SetRangeUser(-50, 50);
    mgAnglesProtonsAlpha -> GetYaxis() -> SetRangeUser(-50, 50);

    mgAnglesProtonsAlpha -> Draw("AP");
    //gPad -> SetGrid();

    for(int q = 0; q < 10; ++q)
    {
        lineAngles[q] -> Draw("same");
    }

    

    CAnglesProtonsAlpha -> SaveAs("../docs/assets/images/AnglesProtonsAlpha.pdf");
    CAnglesProtonsAlpha -> SaveAs("../docs/assets/images/AnglesProtonsAlpha.svg");
    CAnglesProtonsAlpha -> SaveAs("../docs/assets/images/AnglesProtonsAlpha.png");


    for(int i = 0 ; i < Nfiles; ++i)
    {
        cAngles[i] = new TCanvas(Form("Angles %d", i), Form("Angles %d", i), 900, 700);
        hAngles[i] = new TH2D(Form("hAngles_%d", i), Form("Angles %d", i), 60, -55, 55, 60, -55, 55);
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

    if(DebugMuons)
    {
        c2 -> SaveAs("../docs/assets/images/SkyMapPositionNoMask_muons.png");
        c2 -> SaveAs("../docs/assets/images/SkyMapPositionNoMask_muons.pdf");
    }
    else
    {
        c2 -> SaveAs("../docs/assets/images/SkyMapPositionNoMask.png");
        c2 -> SaveAs("../docs/assets/images/SkyMapPositionNoMask.pdf");
    }
    


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
        c3 -> SaveAs("../docs/assets/images/SkyMapPosition_muons.png");
        c3 -> SaveAs("../docs/assets/images/SkyMapPosition_muons.pdf");
    } 
    else 
    {
        c3 -> SaveAs("../docs/assets/images/SkyMapPosition.png");
        c3 -> SaveAs("../docs/assets/images/SkyMapPosition.pdf");
    }
    

    



    TH1D *htheta [Nfiles][Ntot];
    TCanvas *ctheta[Nfiles];
    TString MaskID[Ntot];

    TH1D *hphi [Nfiles][Ntot];
    TCanvas *cphi[Nfiles];

    for(int i = 0; i < Nfiles ; ++i)
    {
        cphi[i] = new TCanvas(Form("cphi%d",i), Form("Phi accepted %d", i), 780, 576);
        cphi[i] -> Divide(Ny, Nx);
        ctheta[i] = new TCanvas(Form("ctheta%d",i), Form("Theta accepted %d", i), 780, 576);
        ctheta[i] -> Divide(Ny, Nx);
        for(int k = 0 ; k < Ntot ; ++k)
        {
            cphi[i] -> cd(k+1);
            HistName = Form("hphi%d%d",i,k);
            HistTitle = Form("Phi accepted Copy No. %d; phi [deg]; counts []", k);
            hphi[i][k] = new TH1D(HistName, HistTitle, 50, -180, 180);
            MaskID[k] = Form("(%s == %d) && %s", BranchName[0].Data(), k, MaskTotal.Data());
            
            if(i == 3)
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[1].Data(), HistName.Data()), Form("(%s == %d)", BranchName[0].Data(), k), "");
            }
            else
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[1].Data(), HistName.Data()), MaskID[k], "");
            }

            ctheta[i] -> cd(k+1);
            HistName = Form("htheta%d%d",i,k);
            HistTitle = Form("Theta accepted Copy No. %d; phi [deg]; counts []", k);
            hphi[i][k] = new TH1D(HistName, HistTitle, 20, 0, 70);
            //MaskID[k] = Form("(%s == %d) && %s", BranchName[0].Data(), k, MaskTotal.Data());

            if(i == 3)
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[0].Data(), HistName.Data()), Form("(%s == %d)", BranchName[0].Data(), k), "");
            }
            else
            {
                Edep[i] -> Draw(Form("(%s*(180/3.415927))>>%s",PolarAngle[0].Data(), HistName.Data()), MaskID[k], "");
            }
            
        }
        cphi[i] -> SaveAs(Form("../docs/assets/images/Phi%d.png", i));
        cphi[i] -> SaveAs(Form("../docs/assets/images/Phi%d.pdf", i));

        ctheta[i] -> SaveAs(Form("../docs/assets/images/Theta%d.png", i));
        ctheta[i] -> SaveAs(Form("../docs/assets/images/Theta%d.pdf", i));
    }


    return;
}