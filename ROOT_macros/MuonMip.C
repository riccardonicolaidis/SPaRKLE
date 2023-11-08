#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TFormula.h"


using namespace std;





void MuonMip()
{

    TString FileName = "MIP_14_08_2022_MUON.root";
    TFile *file = new TFile(FileName);
    TTree *Edep = (TTree*)file->Get("Edep");


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
    E_min_thin   = 0.04;
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
            
            ++j;
            ++CopyNumber;
        }
    }

    cout << "Branches Names : \n############################" << endl;
    for(int i = 0; i < 47; ++i)
    {
        cout << i << " " << BranchName[i].Data() << endl;
    }

    

    //  |-----|-----|-----|-----|
    //  |  0  |  1  |  2  |  3  |
    //  |-----|-----|-----|-----|
    //  |  4  |  5  |  6  |  7  |
    //  |-----|-----|-----|-----|
    //  |  8  |  9  | 10  | 11  |
    //  |-----|-----|-----|-----|
    //  | 12  | 13  | 14  | 15  |
    //  |-----|-----|-----|-----|

    int NAngle1 [4] = {0,3,12,15};
    int NAngle2 [4] = {1,7,14,8};
    int NAngle3 [4] = {5,6,9,10};

    TH1D *h1[3];

    h1[0] = new TH1D("h10","h10",300, 0.1,0.5);
    h1[1] = new TH1D("h11","h11",300, 0.1,0.5);
    h1[2] = new TH1D("h12","h12",300, 0.1,0.5);

    h1[0] -> SetLineColor(kRed);
    h1[1] -> SetLineColor(kBlue);
    h1[2] -> SetLineColor(kBlack);
    h1[0] -> SetLineWidth(4);
    h1[1] -> SetLineWidth(4);
    h1[2] -> SetLineWidth(4);
    h1[0] -> SetLineStyle(1);
    h1[1] -> SetLineStyle(1);
    h1[2] -> SetLineStyle(1);
    h1[0] -> SetMarkerColor(kRed);
    h1[1] -> SetMarkerColor(kBlue);
    h1[2] -> SetMarkerColor(kBlack);
    h1[0] -> SetMarkerStyle(7);
    h1[1] -> SetMarkerStyle(7);
    h1[2] -> SetMarkerStyle(7);
    h1[0] -> SetMarkerSize(1.);
    h1[1] -> SetMarkerSize(1.);
    h1[2] -> SetMarkerSize(1.);
    h1[0] -> SetTitle(";Energy [MeV];Counts");
    h1[1] -> SetTitle(";Energy [MeV];Counts");
    h1[2] -> SetTitle(";Energy [MeV];Counts");
    h1[0] -> SetStats(0);
    h1[1] -> SetStats(0);
    h1[2] -> SetStats(0);


    for(int i= 0; i < 4; ++i)
    {
        cout << "i = " << i << endl;
        if(i == 0)
        {
            Edep -> Draw(Form("%s>>h12",BranchName[11+NAngle3[i]+Ntot].Data()), "", "");
            Edep -> Draw(Form("%s>>h10",BranchName[11+NAngle1[i]+Ntot].Data()), "", "");
            Edep -> Draw(Form("%s>>h11",BranchName[11+NAngle2[i]+Ntot].Data()), "", "");
        }
        else
        {
            Edep -> Draw(Form("%s>>+h12",BranchName[11+NAngle3[i]+Ntot].Data()), "", "");
            Edep -> Draw(Form("%s>>+h10",BranchName[11+NAngle1[i]+Ntot].Data()), "", "");
            Edep -> Draw(Form("%s>>+h11",BranchName[11+NAngle2[i]+Ntot].Data()), "", "");
        }
    }

    //h1[0] -> Scale(1./((h1[0]->Integral())*(h1[0]->GetBinWidth(1))));
    //h1[1] -> Scale(1./((h1[1]->Integral())*(h1[1]->GetBinWidth(1))));
    //h1[2] -> Scale(1./((h1[2]->Integral())*(h1[2]->GetBinWidth(1))));

    
    //h1[0] -> GetYaxis() -> SetRangeUser(0,1.);
    //h1[1] -> GetYaxis() -> SetRangeUser(0,100000);
    //h1[2] -> GetYaxis() -> SetRangeUser(0,100000);


    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    h1[2] -> Draw("HISTsame");
    h1[0] -> Draw("HISTsame");
    h1[1] -> Draw("HISTsame");

    /*
    TF1 *f1[3];

    f1[0] = new TF1("f10","landau(0) + [3] + [4]*x + [5]*pow(x,2.)",0.001,0.35);
    f1[1] = new TF1("f11","landau(0) + [3] + [4]*x + [5]*pow(x,2.)",0.001,0.35);
    f1[2] = new TF1("f12","landau(0) + [3] + [4]*x + [5]*pow(x,2.)",0.001,0.35);

    f1[0] -> SetLineColor(kRed);
    f1[1] -> SetLineColor(kBlue);
    f1[2] -> SetLineColor(kBlack);
    f1[0] -> SetLineWidth(3);
    f1[1] -> SetLineWidth(3);
    f1[2] -> SetLineWidth(3);
    f1[0] -> SetLineStyle(1);
    f1[1] -> SetLineStyle(1);
    f1[2] -> SetLineStyle(1);
    f1[0] -> SetMarkerColor(kRed);
    f1[1] -> SetMarkerColor(kBlue);
    f1[2] -> SetMarkerColor(kBlack);
    f1[0] -> SetMarkerStyle(7);
    f1[1] -> SetMarkerStyle(7);
    f1[2] -> SetMarkerStyle(7);

    f1[0] -> FixParameter(3, 0.001);
    f1[0] -> FixParameter(4, 0.1);
    f1[0] -> FixParameter(5, 0.001);
    h1[0] -> Fit("f10","","L",0.001,0.35);
    f1[0] -> ReleaseParameter(3);
    f1[0] -> ReleaseParameter(4);
    f1[0] -> ReleaseParameter(5);
    h1[0] -> Fit("f10","","L",0.001,0.35);

    f1[1] -> FixParameter(3, 0.001);
    f1[1] -> FixParameter(4, 0.1);
    f1[1] -> FixParameter(5, 0.001);
    h1[1] -> Fit("f11","","L",0.001,0.35);
    f1[1] -> ReleaseParameter(3);
    f1[1] -> ReleaseParameter(4);
    f1[1] -> ReleaseParameter(5);
    h1[1] -> Fit("f11","","L",0.001,0.35);

    f1[2] -> FixParameter(3, 0.001);
    f1[2] -> FixParameter(4, 0.1);
    f1[2] -> FixParameter(5, 0.001);
    h1[2] -> Fit("f12","","L",0.001,0.35);
    f1[2] -> ReleaseParameter(3);
    f1[2] -> ReleaseParameter(4);
    f1[2] -> ReleaseParameter(5);
    h1[2] -> Fit("f12","","L",0.001,0.35);


    cout << "FIT PARAMETERS TABLE" << endl;
    cout << "====================" << endl;
    cout << "Parameter 0" << endl;
    cout << f1[0] -> GetParameter(0) << f1[0] -> GetParError(0) <<  endl;
    cout << f1[1] -> GetParameter(0) << f1[1] -> GetParError(0) <<  endl;
    cout << f1[2] -> GetParameter(0) << f1[2] -> GetParError(0) <<  endl;
    cout << "Parameter 1" << endl;
    cout << f1[0] -> GetParameter(1) << f1[0] -> GetParError(1) <<  endl;
    cout << f1[1] -> GetParameter(1) << f1[1] -> GetParError(1) <<  endl;
    cout << f1[2] -> GetParameter(1) << f1[2] -> GetParError(1) <<  endl;
    cout << "Parameter 2" << endl;
    cout << f1[0] -> GetParameter(2) << f1[0] -> GetParError(2) <<  endl;
    cout << f1[1] -> GetParameter(2) << f1[1] -> GetParError(2) <<  endl;
    cout << f1[2] -> GetParameter(2) << f1[2] -> GetParError(2) <<  endl;

    */
    

    gPad -> SetGrid();

    return;
}