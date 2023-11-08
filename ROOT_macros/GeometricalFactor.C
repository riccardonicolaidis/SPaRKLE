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

using namespace std;



void GeometricalFactor()
{
    // Some initial definitions
    bool DebugMuons = false;
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
    
    // Definitions for the Geometrical Factor part
    TH1D    *hAcceptance[3];
    TH1D    *hFinal[3];
    TH1D    *hIncident[3];
    TCanvas *cFinal[3];
    double ERange[3][2];
    double ERangeGen[3][2];
    double dNdEGen;
    double NgenCalibrated;
    double LActive;
    double Radius;
    int NGen;
    
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
    
    // ********************************************
    //           GEOMETRICAL FACTOR SETTINGS
    // ********************************************
    // Range of energy generated (see in run_..._.mac)
    ERangeGen[0][0] = 0.001;
    ERangeGen[0][1] = 12.;
    ERangeGen[1][0] = 0.01;
    ERangeGen[1][1] = 100.;
    ERangeGen[2][0] = 0.01;
    ERangeGen[2][1] = 400.;
    // Range of energy of the plots
    ERange[0][0] = 0.001;
    ERange[0][1] = 12.;
    ERange[1][0] = 0.01;
    ERange[1][1] = 100.;
    ERange[2][0] = 0.01;
    ERange[2][1] = 400.;
    // Number of bins
    int Nbins[3] = {90, 100, 100};
    LActive = 2 * 6.1;
    Radius = 0.005;
    NGen = 100000000;
    // ********************************************
    // ********************************************

    


    if(DebugMuons)
    Nfiles = 4;
    cout << "Nfiles = " << Nfiles << endl;

    // Definitions of the files 
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


    /* -------------------------------------------------------------------------- */
    /*                    Specific part for GEOMETRICAL FACTOR                    */
    /* -------------------------------------------------------------------------- */

    

    for(int i = 0 ; i < 3 ; ++i)
    {
        cFinal[i] = new TCanvas(Form("cFinal_%d", i), Form("Final energy %d", i), 900, 700);
        
        hIncident[i]   = new TH1D(Form("hIncident_%d", i), Form("Incident energy %d", i), Nbins[i], ERange[i][0], ERange[i][1]);
        hIncident[i] -> SetTitle(Form("Incident energy %d: Log10(E[MeV]): Counts", i));
        hAcceptance[i] = new TH1D(Form("hAcceptance_%d", i), Form("Acceptance %d", i), Nbins[i], ERange[i][0], ERange[i][1]);        
        hAcceptance[i] -> SetTitle(Form("Accepted events %d: Log10(E[MeV]): Counts", i));
        Edep[i] -> Draw(Form("%s>>hIncident_%d", BranchName[1].Data(), i), "", "");
        Edep[i] -> Draw(Form("%s>>hAcceptance_%d", BranchName[1].Data(), i), ConditionGoodEvents.Data(), "");
        
        NgenCalibrated = NGen * (ERange[i][1] - ERange[i][0])/ (ERangeGen[i][1] - ERangeGen[i][0]);;
        dNdEGen = NGen * (hAcceptance[i] -> GetBinWidth(1))/ (ERangeGen[i][1] - ERangeGen[i][0]);
        
        hFinal[i] = (TH1D*) hAcceptance[i] -> Clone(Form("hFinal_%d", i));
        hFinal[i] -> Scale(TMath::Pi() * LActive * LActive / dNdEGen);
        hFinal[i] -> SetTitle(Form(" ; Energy [MeV]; Geometrical Factor [cm^2 rad]"));
        hFinal[i] -> SetStats(0);
        hFinal[i] -> SetMarkerStyle(8);
        hFinal[i] -> SetMarkerSize(0.8);
        hFinal[i] -> SetMarkerColor(kBlue);
        hFinal[i] -> SetLineWidth(2);
        hFinal[i] -> Draw();
        gPad -> SetGrid();
        gPad -> SetLogy();
        cFinal[i] -> Update();
        cFinal[i] -> Draw();
        cout << "#############################" << endl;
        cout << "Geometrical factor " << i << ": " << endl;
        cout << TMath::Pi() * LActive * LActive * (hAcceptance[i] -> GetEntries()) / NgenCalibrated << endl;

        // Save image in all possible formats
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.png",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.pdf",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.root",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.C",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.eps",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.ps",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.jpg",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.gif",i));
        cFinal[i] -> SaveAs(Form("../docs/assets/GeometricalFactor/Acceptance%d.jpeg",i));
    }


    /* -------------------------------------------------------------------------- */
    /*                                New Settings                                */
    /* -------------------------------------------------------------------------- */

    TH1D    *hAcceptance2[3];
    TH1D    *hFinal2[3];
    TH1D    *hIncident2[3];
    double Xmin = 0.09;
    double Xmax = 400.0;
    // Number of bins
    int Nbins2[3] = {3000, 800, 600};

    
    TCanvas *c1 = new TCanvas("c1", "Geometrical factor", 1000, 600);

    


    for(int i = 0 ; i < 3 ; ++i)
    {
        c1 -> cd();
        hIncident2[i]   = new TH1D(Form("hIncident2_%d", i), Form("Incident energy %d", i), Nbins2[i], Xmin, Xmax);
        hIncident2[i] -> SetTitle(Form("Incident energy %d: Log10(E[MeV]): Counts", i));
        hAcceptance2[i] = new TH1D(Form("hAcceptance2_%d", i), Form("Acceptance %d", i), Nbins2[i], Xmin, Xmax);        
        hAcceptance2[i] -> SetTitle(Form("Accepted events %d: Log10(E[MeV]): Counts", i));
        Edep[i] -> Draw(Form("%s>>hIncident2_%d", BranchName[1].Data(), i), "", "");
        Edep[i] -> Draw(Form("%s>>hAcceptance2_%d", BranchName[1].Data(), i), ConditionGoodEvents.Data(), "");
        
        dNdEGen = NGen * (hAcceptance2[i] -> GetBinWidth(1))/ (ERangeGen[i][1] - ERangeGen[i][0]);
        
        hFinal2[i] = (TH1D*) hAcceptance2[i] -> Clone(Form("hFinal2_%d", i));
        hFinal2[i] -> Scale(TMath::Pi() * LActive * LActive / dNdEGen);
        hFinal2[i] -> SetTitle(Form(" ; Energy [MeV]; Geometrical Factor [cm^2 sr]"));
        hFinal2[i] -> SetStats(0);
        hFinal2[i] -> SetMarkerStyle(8);
        hFinal2[i] -> SetMarkerSize(0.8);
        hFinal2[i] -> SetLineWidth(2);
    }
    delete c1;

    TCanvas *cGeomAll = new TCanvas("cGeomAll", "Geometrical factor", 1200, 600);



    TColor *color[3];
    color[0] = new TColor(1756, 242./255., 53./255, 141./255);
    color[1] = new TColor(1757, 4./255., 104./255, 191./255);
    color[2] = new TColor(1758, 242./255., 152./255, 73./255);



    // Create TGraph from file ../ROOT_macros/DEMETER.csv
    
    //mg -> Add(gDEMETER, "PLsame");
    /*
    Double_t EnHEPPL[2] = {0.1, 3.0};
    Double_t GeomHEPPL[2] = {(0.54+0.645)/2., (0.54+0.645)/2.};
    Double_t GeomHEPPLErr[2] = {(0.645-0.54)/2., (0.645-0.54)/2.};
    Double_t GeomHEPPL2[2] = {(3.05+2.03)/2., (3.05+2.03)/2.};
    Double_t GeomHEPPLErr2[2] = {(3.05-2.03)/2., (3.05-2.03)/2.}; 
    TGraphErrors *gHEPPL = new TGraphErrors(2, EnHEPPL, GeomHEPPL, 0, GeomHEPPLErr);
    gHEPPL -> SetFillColor(kGray);
    gHEPPL -> SetLineColor(kGray);
    gHEPPL -> SetFillStyle(3010);
    //mg -> Add(gHEPPL, "a3same");

    TGraphErrors *gHEPPL2 = new TGraphErrors(2, EnHEPPL, GeomHEPPL2, 0, GeomHEPPLErr2);
    gHEPPL2 -> SetFillColor(4);
    gHEPPL2 -> SetLineColor(4);
    gHEPPL2 -> SetFillStyle(3010);
    //mg -> Add(gHEPPL2, "a3same");
    //mg -> Draw();

    //mg -> GetXaxis() -> SetLimits(1e-2,300.);
    //mg -> GetYaxis() -> SetRangeUser(1e-4,3.2);
    */



    for(int i = 0 ; i < 3 ; ++i)
    {
        hFinal2[i] -> GetYaxis() -> SetRangeUser(1e-4,1.3);
        hFinal2[i] -> GetXaxis() -> SetRangeUser(1e-2,300.);
        hFinal2[i] -> SetLineColor(color[i] -> GetNumber());
        hFinal2[i] -> SetMarkerColor(color[i] -> GetNumber());

        hFinal2[i] -> Draw("Esame");        
    }

    TGraph *gDEMETER = new TGraph("../ROOT_macros/DEMETER.csv", "%lg %lg");
    gDEMETER -> SetLineColor(kBlack);
    gDEMETER -> SetLineWidth(2);
    gDEMETER -> SetMarkerColor(kBlack);
    gDEMETER -> SetMarkerStyle(8);
    gDEMETER -> SetMarkerSize(0.8);
    //gDEMETER -> Draw("PLsame");

    TLegend  *legend = new TLegend(0.1, 0.1, 0.3, 0.3);
    legend -> AddEntry(hFinal2[0], "electrons", "lep");
    legend -> AddEntry(hFinal2[1], "protons", "lep");
    legend -> AddEntry(hFinal2[2], "alpha", "lep");
    legend -> AddEntry(gDEMETER, "DEMETER", "lp");
    //legend -> AddEntry(gHEPPL, "HEPP-L : narrow", "f");
    //legend -> AddEntry(gHEPPL2, "HEPP-L : wide", "f");
    legend -> Draw();

    gPad -> SetGrid();
    gPad -> SetLogy();
    gPad -> SetLogx();

    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.png");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.pdf");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.root");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.C");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.eps");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.ps");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.jpg");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.gif");
    cGeomAll -> SaveAs("../docs/assets/GeometricalFactor/GeometricalFactorAll.jpeg");



    return;
}