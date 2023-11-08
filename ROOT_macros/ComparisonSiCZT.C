#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include "TMarker.h"
#include "TH2D.h"
#include "TLegend.h"



using namespace std;

void ComparisonSiCZT()
{
    int Nfiles = 6;
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
    double ResCZT;
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
    Nfiles       = 8;
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
    ResCZT       = 0.02;
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
    Nbinsx      = 400;
    Nbinsy      = 400;
    // ********************************************
    // ********************************************


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
            cout << "FileName[" << i << "] = " << FileName[i].Data() << endl;
        }
    }



    // Definitions of the Ttrees and Tfiles
    TFile *file[Nfiles];
    TTree *Edep[Nfiles];

    for(int i = 0; i < Nfiles ; ++i)
    {
        file[i] = TFile::Open(FileName[i]);
        Edep[i] = (TTree*) file[i] -> Get("Edep");
        //Edep[i] -> Print();
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
    TString gThinTot;
    TString gThickTot;


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
                ThinTot   = Form("(%s", BranchName[j].Data());
                ThickTot  = Form("(%s", BranchName[j + Ntot].Data());
                gThinTot  = Form("(g%s", BranchName[j].Data());
                gThickTot = Form("(g%s", BranchName[j + Ntot].Data());
            }
            else
            {
                ThinTot   += Form("+%s", BranchName[j].Data());
                ThickTot  += Form("+%s", BranchName[j + Ntot].Data());
                gThinTot  += Form("+g%s", BranchName[j].Data());
                gThickTot += Form("+g%s", BranchName[j + Ntot].Data());
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

    ThinTot  += ")";
    ThickTot += ")";
    gThinTot += ")";
    gThickTot += ")";

    for(int i = 0; i < Nfiles ; ++i)
    {
        Edep[i] -> SetAlias("ThinTot", ThinTot.Data());
        Edep[i] -> SetAlias("ThickTot", ThickTot.Data());
        Edep[i] -> SetAlias("gThinTot", gThinTot.Data());
        Edep[i] -> SetAlias("gThickTot", gThickTot.Data());
        Edep[i] -> SetAlias("PID", "TMath::Log10(gThinTot * (gThinTot + gThickTot))");
        Edep[i] -> SetAlias("gETOT", "(gThinTot + gThickTot)");
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
                if(k >= 11 && k <= 42)
                {
                    Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
                }
                else
                {
                    if(i < (Nfiles/2))
                    {
                        Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
                    }
                    else
                    {
                        Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResCZT));
                    }
                    
                }
                // Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
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
    /*                         COMPARISON BETWEEN THE DATA                        */
    /* -------------------------------------------------------------------------- */

    TString Condition = Form("(%s) && ( (%s - ThinTot - ThickTot) < %g)", ConditionPairSiliconAll.Data(), BranchName[1].Data(), E_th_Vetoed);

    

    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    TH2D *h2Si  = new TH2D("h2Si","h2Si; Energy",Nbinsx,Eminx,Emaxx,Nbinsy,Eminy,Emaxy);
    TH2D *h2CZT = new TH2D("h2CZT","h2CZT",Nbinsx,Eminx,Emaxx,Nbinsy,Eminy,Emaxy);
    c1 -> Divide(2,1);
    c1 -> cd(1);
    for(int i = 0; i <=2 ; ++i)
    {
        cout << "i = " << i << endl;
        if(i == 0)
        {
            Edep[i] -> Draw("PID:TMath::Log10(gETOT)>>h2Si", Condition.Data(), "colz");
        }
        else
        {
            Edep[i] -> Draw("PID:TMath::Log10(gETOT)>>+h2Si", Condition.Data(), "colz");
        }
    }
    h2Si -> SetStats(0);
    h2Si -> SetTitle("; log_{10}(E/ 1 MeV); PID []");
    gPad -> SetGrid();
    gPad -> SetLogz();

    c1 -> cd(2);

    for(int i = 3; i <= 5; ++i)
    {
        cout << "i = " << i << endl;
        if(i == 3)
        {
            Edep[i] -> Draw("PID:TMath::Log10(gETOT)>>h2CZT", Condition.Data(), "colz");
        }
        else
        {
            Edep[i] -> Draw("PID:TMath::Log10(gETOT)>>+h2CZT", Condition.Data(), "colz");
        }
    }
    h2CZT -> SetStats(0);
    h2CZT -> SetTitle("; log_{10}(E/ 1 MeV); PID");
    gPad -> SetGrid();
    gPad -> SetLogz();

    c1 -> SaveAs("../docs/assets/CZT/SiCZTComparison.pdf");
    c1 -> SaveAs("../docs/assets/CZT/SiCZTComparison.png");
    c1 -> SaveAs("../docs/assets/CZT/SiCZTComparison.C");
    c1 -> SaveAs("../docs/assets/CZT/SiCZTComparison.root");

    

    /* -------------------------------------------------------------------------- */
    /*                           Section for Gamma yield                          */
    /* -------------------------------------------------------------------------- */

    int Nbins = 600;
    double Emingen = 0.05;
    double Emaxgen = 2.;
    double Emin = 0.05;
    double Emax = 2.;
    int Ngen = 5000000;
    E_th_Vetoed = 0.1;

    double CorrectionFactor;

    TH1D *hSi = new TH1D("hSi","hSi",Nbins,log10(Emin),log10(Emax));
    TH1D *hCZT = new TH1D("hCZT","hCZT",Nbins,log10(Emin),log10(Emax));

    TCanvas *c2 = new TCanvas("c2","c2",1200,600);

    Condition = Form(" ( (%s - ThickTot) < %g)", BranchName[1].Data(), E_th_Vetoed);
    cout << "Drawing 1" << endl;   
    Edep[6] -> Draw("log10(ThickTot)>>hSi", Condition.Data(), "");
    cout << "Drawing 2" << endl;
    Edep[7] -> Draw("log10(ThickTot)>>hCZT", Condition.Data(), "same");

    CorrectionFactor = (log10(Emaxgen) - log10(Emingen))/(Ngen*(hSi -> GetBinWidth(1)));
    hSi -> Scale(CorrectionFactor);
    CorrectionFactor = (log10(Emaxgen) - log10(Emingen))/(Ngen*(hCZT -> GetBinWidth(1)));
    hCZT -> Scale(CorrectionFactor);
    hSi -> SetStats(0);
    hSi -> SetTitle("; Log_{10}( E_{dep}[MeV] ); Counts");
    hSi -> SetLineColor(kRed);
    hSi -> SetMarkerColor(kRed);
    hSi -> SetMarkerStyle(8);
    hSi -> SetMarkerSize(0.5);
    hSi -> SetLineWidth(2);
    hCZT -> SetStats(0);
    hCZT -> SetTitle("; Log_{10}( E_{dep}[MeV] ); Counts");
    hCZT -> SetLineColor(kBlue);
    hCZT -> SetMarkerColor(kBlue);
    hCZT -> SetMarkerStyle(8);
    hCZT -> SetMarkerSize(0.4);
    hCZT -> SetLineWidth(2);

    hSi  -> GetYaxis() -> SetRangeUser(0.0001, 1.1);
    //hCZT -> GetYaxis() -> SetRangeUser(0.0001, 1.1);


    gPad -> SetGrid();
    gPad -> SetLogy();
    //gPad -> SetLogx();

    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg -> AddEntry(hSi, "500 um Si", "l");
    leg -> AddEntry(hCZT, "1 mm CZT", "l");
    leg -> Draw();

    c2 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGamma.pdf");
    c2 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGamma.png");
    c2 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGamma.C");
    c2 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGamma.root");





    /* -------------------------------------------------------------------------- */
    double AreaSilicon = 50. * 16 * 1e-2; // cm^2
    double AreaCZT = 10. * 10. * 16 * 1e-2; // cm^2


    TCanvas *c3 = new TCanvas("c3","c3",1200,600);

    TH1D *hSi2 = new TH1D("hSi2","hSi2",Nbins,log10(Emin),log10(Emax));
    TH1D *hCZT2 = new TH1D("hCZT2","hCZT2",Nbins,log10(Emin),log10(Emax));


    Edep[6] -> Draw("TMath::Log10(ThickTot)>>hSi2", Condition.Data(), "");
    Edep[7] -> Draw("TMath::Log10(ThickTot)>>hCZT2", Condition.Data(), "same");

    CorrectionFactor = (log10(Emaxgen) - log10(Emingen))/(Ngen*(hSi2 -> GetBinWidth(1))) * AreaSilicon;
    hSi2 -> Scale(CorrectionFactor);
    CorrectionFactor = (log10(Emaxgen) - log10(Emingen))/(Ngen*(hCZT2 -> GetBinWidth(1))) * AreaCZT;
    hCZT2 -> Scale(CorrectionFactor);
    hSi2 -> SetStats(0);
    hSi2 -> SetTitle("; Log_{10}( E_{dep}[MeV] );  Effective Area [cm^{2}]");
    hSi2 -> SetLineColor(kRed);
    hSi2 -> SetMarkerColor(kRed);
    hSi2 -> SetMarkerStyle(8);
    hSi2 -> SetMarkerSize(0.5);
    hSi2 -> SetLineWidth(2);
    hCZT2 -> SetStats(0);
    hCZT2 -> SetTitle("; Log_{10}( E_{dep}[MeV] ); Effective Area [cm^{2}]");
    hCZT2 -> SetLineColor(kBlue);
    hCZT2 -> SetMarkerColor(kBlue);
    hCZT2 -> SetMarkerStyle(8);
    hCZT2 -> SetMarkerSize(0.4);
    hCZT2 -> SetLineWidth(2);

    hSi2  -> GetYaxis() -> SetRangeUser(0.001, 30.0);


    gPad -> SetGrid();
    gPad -> SetLogy();
    //gPad -> SetLogx();

    TLegend *leg2 = new TLegend(0.7,0.7,0.9,0.9);
    leg2 -> AddEntry(hSi2, "500 um Si", "l");
    leg2 -> AddEntry(hCZT2, "1 mm CZT", "l");
    leg2 -> Draw();

    c3 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGammaArea.pdf");
    c3 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGammaArea.png");
    c3 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGammaArea.C");
    c3 -> SaveAs("../docs/assets/CZT/SiCZTComparisonGammaArea.root");




    return;
}
