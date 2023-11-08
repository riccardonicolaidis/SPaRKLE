#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <regex>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"

using namespace std;


int AnalysisCZTThickness() 
{
    // Open FileNames.txt and store the file names in a vector

    ifstream file("FileNames.txt");
    string str;
    vector<string> fileNames;
    vector<string> fileNameElectron;
    vector<string> fileNameProton;
    vector<string> fileNameAlpha;
    vector<string> fileNameMuon;



    for(int i = 0 ; i < 76; ++i)
    {
        getline(file, str);
        fileNames.push_back(str);

        if(fileNames[i].find("LECTRON") != string::npos)
        {
            fileNameElectron.push_back(fileNames[i]);
            cout << "LECTRON: " << fileNames[i] << endl;
        }
        else if(fileNames[i].find("PROTON") != string::npos)
        {
            fileNameProton.push_back(fileNames[i]);
            cout << "PROTON: " << fileNames[i] << endl;
        }
        else if(fileNames[i].find("ALPHA") != string::npos)
        {
            fileNameAlpha.push_back(fileNames[i]);
            cout << "ALPHA: " << fileNames[i] << endl;
        }
        else if(fileNames[i].find("MUON") != string::npos)
        {
            fileNameMuon.push_back(fileNames[i]);
            cout << "MUON: " << fileNames[i] << endl;
        }
    }

    vector<TFile*> TFileElectron;
    vector<TFile*> TFileProton;
    vector<TFile*> TFileAlpha;
    vector<TFile*> TFileMuon;

    vector<TTree*> TTreeElectron;
    vector<TTree*> TTreeProton;
    vector<TTree*> TTreeAlpha;
    vector<TTree*> TTreeMuon;

    for(int i = 0; i < fileNameElectron.size(); ++i)
    {
        TFileElectron.push_back(new TFile(fileNameElectron[i].c_str()));
        TTreeElectron.push_back((TTree*)TFileElectron[i]->Get("Edep"));
    }

    for(int i = 0; i < fileNameProton.size(); ++i)
    {
        TFileProton.push_back(new TFile(fileNameProton[i].c_str()));
        TTreeProton.push_back((TTree*)TFileProton[i]->Get("Edep"));
    }

    for(int i = 0; i < fileNameAlpha.size(); ++i)
    {
        TFileAlpha.push_back(new TFile(fileNameAlpha[i].c_str()));
        TTreeAlpha.push_back((TTree*)TFileAlpha[i]->Get("Edep"));
    }

    for(int i = 0; i < fileNameMuon.size(); ++i)
    {
        TFileMuon.push_back(new TFile(fileNameMuon[i].c_str()));
        TTreeMuon.push_back((TTree*)TFileMuon[i]->Get("Edep"));
    }


    /* -------------------------------------------------------------------------- */
    /*                       DEFINISCO TUTTI I NOMI DEI TREE                      */
    /* -------------------------------------------------------------------------- */

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
    ResCZT       = 0.05;
    // Trigger and veto thresholds
    E_min_thin   = 0.01;
    E_min_thick  = 0.04;
    E_th_Vetoed  = 0.1;
    E_th_plastic = 0.2;

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
                gThinTot  = Form("(g%s)", BranchName[j].Data());
                gThickTot = Form("(g%s)", BranchName[j + Ntot].Data());
            }
            else
            {
                ThinTot  = Form("(%s + %s)",ThinTot.Data(), BranchName[j].Data());
                ThickTot = Form("(%s + %s)",ThickTot.Data(), BranchName[j + Ntot].Data());
                gThinTot  = Form("(%s + g%s)",gThinTot.Data(), BranchName[j].Data());
                gThickTot = Form("(%s + g%s)",gThickTot.Data(), BranchName[j + Ntot].Data());
            }
            xHole[CopyNumber]              = -(LThin/2.) + (xDelta/2.) + ix * xDelta;
            yHole[CopyNumber]              = -(LThin/2.) + (yDelta/2.) + iy * yDelta;
            PhiDir[CopyNumber]             = std::atan2((iy-1.5),(ix-1.5));
            distanceR                      = std::sqrt((iy-1.5)*(iy-1.5) + (ix-1.5)*(ix-1.5));
            ThetaDir[CopyNumber]           = distanceR * thetaMax / RMax;

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

    for(int i = 0; i < fileNameElectron.size() ; ++i)
    {
        cout << "Setting alias in File " << i << " : " << fileNameElectron[i] << endl;
        cout << "Setting alias in File " << i << " : " << fileNameProton[i] << endl;
        cout << "Setting alias in File " << i << " : " << fileNameMuon[i] << endl;
        cout << "Setting alias in File " << i << " : " << fileNameAlpha[i] << endl;

        TTreeElectron[i] -> SetAlias("ThinTot", ThinTot.Data());
        TTreeElectron[i] -> SetAlias("ThickTot", ThickTot.Data());
        TTreeElectron[i] -> SetAlias("gThinTot", gThinTot.Data());
        TTreeElectron[i] -> SetAlias("gThickTot", gThickTot.Data());
        TTreeElectron[i] -> SetAlias("PID", "TMath::Log10(ThinTot * (ThinTot + ThickTot))");
        TTreeElectron[i] -> SetAlias("gPID", "TMath::Log10(gThinTot * (gThinTot + gThickTot))");
        TTreeElectron[i] -> SetAlias("ETOT", "(ThinTot + ThickTot)");
        TTreeElectron[i] -> SetAlias("gETOT", "(gThinTot + gThickTot)");

        TTreeProton[i]   -> SetAlias("ThinTot", ThinTot.Data());
        TTreeProton[i]   -> SetAlias("ThickTot", ThickTot.Data());
        TTreeProton[i]   -> SetAlias("gThinTot", gThinTot.Data());
        TTreeProton[i]   -> SetAlias("gThickTot", gThickTot.Data());
        TTreeProton[i]   -> SetAlias("PID", "TMath::Log10(ThinTot * (ThinTot + ThickTot))");
        TTreeProton[i]   -> SetAlias("gPID", "TMath::Log10(gThinTot * (gThinTot + gThickTot))");
        TTreeProton[i]   -> SetAlias("ETOT", "(ThinTot + ThickTot)");
        TTreeProton[i]   -> SetAlias("gETOT", "(gThinTot + gThickTot)");

        TTreeMuon[i]     -> SetAlias("ThinTot", ThinTot.Data());
        TTreeMuon[i]     -> SetAlias("ThickTot", ThickTot.Data());
        TTreeMuon[i]     -> SetAlias("gThinTot", gThinTot.Data());
        TTreeMuon[i]     -> SetAlias("gThickTot", gThickTot.Data());
        TTreeMuon[i]     -> SetAlias("PID", "TMath::Log10(ThinTot * (ThinTot + ThickTot))");
        TTreeMuon[i]     -> SetAlias("gPID", "TMath::Log10(gThinTot * (gThinTot + gThickTot))");
        TTreeMuon[i]     -> SetAlias("ETOT", "(ThinTot + ThickTot)");
        TTreeMuon[i]     -> SetAlias("gETOT", "(gThinTot + gThickTot)");

        TTreeAlpha[i]    -> SetAlias("ThinTot", ThinTot.Data());
        TTreeAlpha[i]    -> SetAlias("ThickTot", ThickTot.Data());
        TTreeAlpha[i]    -> SetAlias("gThinTot", gThinTot.Data());
        TTreeAlpha[i]    -> SetAlias("gThickTot", gThickTot.Data());
        TTreeAlpha[i]    -> SetAlias("PID", "TMath::Log10(ThinTot * (ThinTot + ThickTot))");
        TTreeAlpha[i]    -> SetAlias("gPID", "TMath::Log10(gThinTot * (gThinTot + gThickTot))");
        TTreeAlpha[i]    -> SetAlias("ETOT", "(ThinTot + ThickTot)");
        TTreeAlpha[i]    -> SetAlias("gETOT", "(gThinTot + gThickTot)");

        /* -------------------------------- SMEARING -------------------------------- */
        /* ------------------------ YOU ONLY NEED TO ADD A g ------------------------ */
        for(int k = 9; k < 47; ++k)
        {
            TTreeElectron[i] -> SetAlias(Form("wnorm%d",k),"(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))");
            TTreeProton[i]   -> SetAlias(Form("wnorm%d",k),"(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))");
            TTreeMuon[i]     -> SetAlias(Form("wnorm%d",k),"(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))");
            TTreeAlpha[i]    -> SetAlias(Form("wnorm%d",k),"(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))");
            if(k == 9 || k == 10 || (k >= 43 && k <= 46))
            {
                TTreeElectron[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResPlastic));
                TTreeProton[i]   -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResPlastic));
                TTreeMuon[i]     -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResPlastic));
                TTreeAlpha[i]    -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResPlastic));
            }
            else 
            {
                TTreeElectron[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
                TTreeProton[i]   -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
                TTreeMuon[i]     -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
                TTreeAlpha[i]    -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
            }
        }

        /* ------------------- Incident direction of the particle ------------------- */

        TTreeElectron[i] -> SetAlias("NormP", Form("(TMath::Sqrt(TMath::Power(%s,2) + TMath::Power(%s,2) + TMath::Power(%s,2)))", BranchName[6].Data(), BranchName[7].Data(), BranchName[8].Data()));
        TTreeProton[i]   -> SetAlias("NormP", Form("(TMath::Sqrt(TMath::Power(%s,2) + TMath::Power(%s,2) + TMath::Power(%s,2)))", BranchName[6].Data(), BranchName[7].Data(), BranchName[8].Data()));
        TTreeMuon[i]     -> SetAlias("NormP", Form("(TMath::Sqrt(TMath::Power(%s,2) + TMath::Power(%s,2) + TMath::Power(%s,2)))", BranchName[6].Data(), BranchName[7].Data(), BranchName[8].Data()));
        TTreeAlpha[i]    -> SetAlias("NormP", Form("(TMath::Sqrt(TMath::Power(%s,2) + TMath::Power(%s,2) + TMath::Power(%s,2)))", BranchName[6].Data(), BranchName[7].Data(), BranchName[8].Data()));
        for(int k = 6; k <= 8; ++k)
        {
            TTreeElectron[i] -> SetAlias(DirName[k-6].Data(), Form("((%s)/(NormP))", BranchName[k].Data()));
            TTreeProton[i]   -> SetAlias(DirName[k-6].Data(), Form("((%s)/(NormP))", BranchName[k].Data()));
            TTreeMuon[i]     -> SetAlias(DirName[k-6].Data(), Form("((%s)/(NormP))", BranchName[k].Data()));
            TTreeAlpha[i]    -> SetAlias(DirName[k-6].Data(), Form("((%s)/(NormP))", BranchName[k].Data()));
        }

        // Legend
        // [0] : Theta
        // [1] : Phi
        TTreeElectron[i] -> SetAlias(PolarAngle[0].Data()   , Form("((TMath::ACos(%s)))", DirName[2].Data())); // *(180/3.415927))
        TTreeElectron[i] -> SetAlias(PolarAngle[1].Data()   , Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data())); // *(180/3.415927))
        TTreeElectron[i] -> SetAlias(NewPolarAngle[1].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(TMath::Cos(%s))))",              PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() )); // *(180/3.415927))
        TTreeElectron[i] -> SetAlias(NewPolarAngle[0].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)),TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data())); // *(180/3.415927))

        TTreeProton[i]   -> SetAlias(PolarAngle[0].Data()   , Form("((TMath::ACos(%s)))", DirName[2].Data())); // *(180/3.415927))
        TTreeProton[i]   -> SetAlias(PolarAngle[1].Data()   , Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data())); // *(180/3.415927))
        TTreeProton[i]   -> SetAlias(NewPolarAngle[1].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(TMath::Cos(%s))))",              PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() )); // *(180/3.415927))
        TTreeProton[i]   -> SetAlias(NewPolarAngle[0].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)),TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data())); // *(180/3.415927))

        TTreeMuon[i]     -> SetAlias(PolarAngle[0].Data()   , Form("((TMath::ACos(%s)))", DirName[2].Data())); // *(180/3.415927))
        TTreeMuon[i]     -> SetAlias(PolarAngle[1].Data()   , Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data())); // *(180/3.415927))
        TTreeMuon[i]     -> SetAlias(NewPolarAngle[1].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(TMath::Cos(%s))))",              PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() )); // *(180/3.415927))
        TTreeMuon[i]     -> SetAlias(NewPolarAngle[0].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)),TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data())); // *(180/3.415927))

        TTreeAlpha[i]    -> SetAlias(PolarAngle[0].Data()   , Form("((TMath::ACos(%s)))", DirName[2].Data())); // *(180/3.415927))
        TTreeAlpha[i]    -> SetAlias(PolarAngle[1].Data()   , Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data())); // *(180/3.415927))
        TTreeAlpha[i]    -> SetAlias(NewPolarAngle[1].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(TMath::Cos(%s))))",              PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() )); // *(180/3.415927))
        TTreeAlpha[i]    -> SetAlias(NewPolarAngle[0].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)),TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data())); // *(180/3.415927))


        /* -------------------- Particle identification parameter ------------------- */

        j = 11;
        CopyNumber = 0;
        cout << "Defining total energy in File " << i << " : " << fileNameElectron[i] << endl;
        cout << "Defining total energy in File " << i << " : " << fileNameProton[i] << endl;
        cout << "Defining total energy in File " << i << " : " << fileNameMuon[i] << endl;
        cout << "Defining total energy in File " << i << " : " << fileNameAlpha[i] << endl;


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

                TTreeElectron[i] -> SetAlias(PIDName[CopyNumber], Form("(TMath::Log10(g%s * %s))",BranchName[j].Data(), TotEnergyName[CopyNumber].Data()));

                TTreeProton[i] -> SetAlias(TotEnergyName[CopyNumber], TotEnergyCondition[CopyNumber]);
                TTreeProton[i] -> SetAlias(PIDName[CopyNumber], Form("(TMath::Log10(g%s * %s))",BranchName[j].Data(), TotEnergyName[CopyNumber].Data()));

                TTreeMuon[i] -> SetAlias(TotEnergyName[CopyNumber], TotEnergyCondition[CopyNumber]);
                TTreeMuon[i] -> SetAlias(PIDName[CopyNumber], Form("(TMath::Log10(g%s * %s))",BranchName[j].Data(), TotEnergyName[CopyNumber].Data()));

                TTreeAlpha[i] -> SetAlias(TotEnergyName[CopyNumber], TotEnergyCondition[CopyNumber]);
                TTreeAlpha[i] -> SetAlias(PIDName[CopyNumber], Form("(TMath::Log10(g%s * %s))",BranchName[j].Data(), TotEnergyName[CopyNumber].Data()));


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

    TString Condition = Form("(%s) && ( (%s - ThinTot - ThickTot) < %g)", ConditionPairSiliconAll.Data(), BranchName[1].Data(), E_th_Vetoed);

    /* -------------------------------------------------------------------------- */
    /*                      VETTORE CON GLI SPESSORI DEL CZT                      */
    /* -------------------------------------------------------------------------- */

    vector<double> Thickness;
    // Information is in the filenames in the form "_1000_um"
    regex pattern = regex("\\d{4}_um");
    regex pattern2 = regex("\\d{5}_um");


    for(int i = 0; i < fileNameElectron.size(); ++i)
    {
        smatch match;
        smatch match2;
        if(regex_search(fileNameElectron[i], match, pattern))
        {
            if(regex_search(fileNameElectron[i], match2, pattern2))
            {
                string match_str = match2.str();
                Thickness.push_back(stod(match_str));
                cout << "Thickness[" << i << "] = " << Thickness[i] << endl;
            }
            else
            {
                string match_str = match.str();
                Thickness.push_back(stod(match_str));
                cout << "Thickness[" << i << "] = " << Thickness[i] << endl;
            }
        }
    }


    vector<TCanvas*> CanvasPID;
    vector<TH2D*> Histogram_PID;

    for(int i = 0; i < TTreeMuon.size(); ++i)
    {
        CanvasPID.push_back(new TCanvas(Form("CanvasPID_%d", i), Form("CanvasPID_%d", i), 800, 400));
        Histogram_PID.push_back(new TH2D(Form("Histogram_PID_%d", i), Form("Histogram_PID_%d", i), 200, -1.5, 3.4, 200, -3.4, 3.4));
        TTreeElectron[i] -> Draw(Form("gPID:TMath::Log10(gETOT)>>%s", Histogram_PID[i]->GetName()), Condition.Data(), "colz");
        TTreeProton[i] -> Draw(Form("gPID:TMath::Log10(gETOT)>>+%s", Histogram_PID[i]->GetName()), Condition.Data(), "colz");
        TTreeAlpha[i] -> Draw(Form("gPID:TMath::Log10(gETOT)>>+%s", Histogram_PID[i]->GetName()), Condition.Data(), "colz");
        gPad -> SetGrid();
        gPad -> SetLogz();

        Histogram_PID[i] -> GetXaxis() -> SetTitle("Log10(Energy / 1 MeV)");
        Histogram_PID[i] -> GetYaxis() -> SetTitle("Log10(PID)");
        Histogram_PID[i] -> SetTitle(Form("CZT Thickness : %g mm", Thickness[i]*1e-3));
        Histogram_PID[i] -> SetStats(0);

        CanvasPID[i] -> SaveAs(Form("../docs/PIDCZT/%g.png", Thickness[i]));
        CanvasPID[i] -> SaveAs(Form("../docs/PIDCZT/%g.pdf", Thickness[i]));

    }








    vector<TCanvas*> Canvas;
    vector<TH1D*> Histogram_CZT;
    vector<TF1*> Fit_Landau;
    vector<double> MPV;
    vector<double> MPV_error;

    


    for(int i = 0; i < TTreeMuon.size(); ++i)
    {
        Canvas.push_back(new TCanvas(Form("Canvas_%d", i), Form("Canvas_%d", i), 800, 400));
        Histogram_CZT.push_back(new TH1D(Form("Histogram_CZT_%d", i), Form("Histogram_CZT_%d", i), 400, 0.001, 10));
        TTreeMuon[i] -> Draw(Form("(%s +0*TMath::Cos(%s)) >> %s", BranchName[27].Data(), PolarAngle[0].Data(), Histogram_CZT[i] -> GetName()), "", "");
    
        Fit_Landau.push_back(new TF1(Form("Fit_Landau_%d", i), "landau(0) + [3] + [4]*x", 0.001, 10));
        Fit_Landau[i] -> SetParameter(1, Histogram_CZT[i] -> GetMean());
        Fit_Landau[i] -> SetParameter(2, Histogram_CZT[i] -> GetRMS());
        Fit_Landau[i] -> FixParameter(3, 0);
        Fit_Landau[i] -> FixParameter(4, 0);

        double mean = Histogram_CZT[i] -> GetMean();
        double rms = Histogram_CZT[i] -> GetRMS();

        Histogram_CZT[i] -> Fit(Fit_Landau[i], "", "", mean-rms, mean+3*rms);

        Fit_Landau[i] -> ReleaseParameter(3);
        Fit_Landau[i] -> ReleaseParameter(4);

        Histogram_CZT[i] -> Fit(Fit_Landau[i], "", "", mean-rms, mean+3*rms);
        
        MPV.push_back(Fit_Landau[i] -> GetParameter(1));
        MPV_error.push_back(Fit_Landau[i] -> GetParError(1));

        //Canvas[i] -> Close();
    }

    TGraphErrors *Graph = new TGraphErrors(Thickness.size(), &Thickness[0], &MPV[0], 0, &MPV_error[0]);

    TCanvas *CanvasGraph = new TCanvas("Canvas", "Canvas", 100, 50);
    CanvasGraph -> Divide(1, 2);
    // Second subfigure with smaller height

    CanvasGraph -> cd(1);


    Graph -> Draw("AP");
    //Graph -> GetXaxis() -> SetTitle("Thickness [#mum]");
    Graph -> GetYaxis() -> SetTitle("MPV [MeV]");

    Graph -> SetMarkerStyle(20);
    Graph -> SetMarkerSize(1.5);
    Graph -> SetMarkerColor(kBlack);


    TF1 *LinFit = new TF1("Fit", "[0] + [1]*x", 0, 10000);

    Graph -> Fit(LinFit, "", "", 0, 10000);
    LinFit -> Draw("same");
    
    CanvasGraph -> cd(2);
    // Residual plot
    vector<double> Residual;
    for(int i = 0; i < MPV.size(); ++i)
    {
        Residual.push_back(MPV[i] - LinFit -> Eval(Thickness[i]));
    }

    TGraphErrors *GraphResidual = new TGraphErrors(Thickness.size(), &Thickness[0], &Residual[0], 0, &MPV_error[0]);
    
    GraphResidual -> Draw("AP");
    GraphResidual -> GetXaxis() -> SetTitle("Thickness [#mum]");
    GraphResidual -> GetYaxis() -> SetTitle("Residual [MeV]");
    GraphResidual -> SetMarkerStyle(20);
    GraphResidual -> SetMarkerSize(1.5);
    GraphResidual -> SetMarkerColor(kBlack);


    // Mev / mm 
    cout << "Energy Loss: "<< LinFit->GetParameter(1)*1e3 << " +/- "<< LinFit->GetParError(1)*1e3 <<" MeV / mm" << endl;

    





    return 0;
}