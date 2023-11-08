#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"


#include "../ROOT_macros/Custom/TH2DLog.h"

#define Nx 3
#define Ny 3
#define N_PLASTIC_NO_VETO 5
#define N_PASTIC_LATERAL 4
#define Ntot Nx*Ny
#define N_BRANCHES 9+2*Ntot+N_PLASTIC_NO_VETO+N_PASTIC_LATERAL+2



using namespace std;

int ExtendedGeometry()
{
    int Nfiles = 3;
    int j;
    int CopyNumber;

    double ResSilicon;
    double ResPlastic;
    double E_min_thin;   // Thick layer
    double E_min_thick;  // Thin layer
    double E_th_Vetoed;  // Energy dispersion (Veto threshold)
    double E_th_plastic;
    
    double ERangeGen[3][2];

    ERangeGen[0][0] = 0.09;
    ERangeGen[0][1] = 50.0;
    ERangeGen[1][0] = 3.;
    ERangeGen[1][1] = 140.0;
    ERangeGen[2][0] = 10.;
    ERangeGen[2][1] = 500.0;
    
    // ********************************************
    // ********************************************
    //                 SETTINGS
    // ********************************************
    // ********************************************
    Nfiles       = 3;
    // Geometric parameter
    // Smearing parameters
    ResSilicon   = 0.01;
    ResPlastic   = 0.05;
    // Trigger and veto thresholds
    E_min_thin   = 0.02;
    E_min_thick  = 0.04;
    E_th_Vetoed  = 0.01;
    E_th_plastic = 0.2;
    
    // ********************************************
    // ********************************************

    ofstream AliasFile; // overwriting the file if it exists
    AliasFile.open("../ROOT_macro_diagnostics/AliasFile.txt");

    ofstream BranchFile;
    BranchFile.open("../ROOT_macro_diagnostics/BranchFile.txt");

    ofstream ConditionFile;
    ConditionFile.open("../ROOT_macro_diagnostics/ConditionFile.txt");



    // Definitions of the files 
    TString FileName[Nfiles];

    // Read the FileNames.txt and for each line of the file fill the TString FileName[]
    ifstream FileNames("FileNames.txt");
    if(!FileNames.is_open())
    {
        std::cout << "FileNames.txt not found" << std::endl;
        return 0;
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
    }

    // Definitions of the initial quantities from the data of the simulation


    TString BranchName[N_BRANCHES];
    TString TotEnergyName[Ntot];
    TString TotEnergyCondition[Ntot];
    TString PIDName[Ntot];
    TString DirName[3];
    TString PolarAngle[2];
    TString NewPolarAngle[2];
    TString AliasETot;
    TString ThinTot;
    TString ThickTot;
    TString gThinTot;
    TString gThickTot;

    TString EnergyPlasticScint[N_PLASTIC_NO_VETO+1];
    TString gEnergyPlasticScint[N_PLASTIC_NO_VETO+1];
    TString TotalEnergyPlasticScint;
    TString gTotalEnergyPlasticScint;

    TString EnergyPlastic2Veto[N_PLASTIC_NO_VETO+1];
    TString gEnergyPlastic2Veto[N_PLASTIC_NO_VETO+1];

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

    const int jPlastic = 9;
    j = jPlastic;

    for(int i = 0; i < (N_PLASTIC_NO_VETO+1); ++i)
    {
        BranchName[j] = Form("Ed_Veto%d",i);
        j++;
    }

    const int jDrilledVeto = j;

    BranchName[j++] = "Ed_DrilledVeto";

    const int jLateral = j;

    BranchName[j++] = "Ed_Lateral0";
    BranchName[j++] = "Ed_Lateral1";
    BranchName[j++] = "Ed_Lateral2";
    BranchName[j++] = "Ed_Lateral3";

    TString LateralDrilled = "(Ed_Lateral0 + Ed_Lateral1 + Ed_Lateral2 + Ed_Lateral3 + Ed_DrilledVeto)";

    const int jSilicon = j;
    
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
            ++j;
            ++CopyNumber;
        }
    }

    std::cout << "Branches Names : \n############################" << std::endl;
    BranchFile << "Branches Names : \n############################" << std::endl;
    for(int i = 0; i < N_BRANCHES; ++i)
    {
        std::cout << i << " " << BranchName[i].Data() << std::endl;
        BranchFile << i << " " << BranchName[i].Data() << std::endl;
    }


    /* -------------------------------------------------------------------------- */
    /*                              Alias definitions                             */
    /* -------------------------------------------------------------------------- */

    /* -------------------------------------------------------------------------- */
    /*                              ALIAS DEFINITIONS                             */
    /* -------------------------------------------------------------------------- */

    for(int i = 0; i < Nfiles ; ++i)
    {
        std::cout << "Setting alias in File " << i << " : " << FileName[i].Data() << std::endl;
        AliasFile << "Setting alias in File " << i << " : " << FileName[i].Data() << std::endl;

        Edep[i] -> SetAlias("ThinTot",ThinTot.Data());
        Edep[i] -> SetAlias("ThickTot",ThickTot.Data());
        Edep[i] -> SetAlias("gThinTot",gThinTot.Data());
        Edep[i] -> SetAlias("gThickTot",gThickTot.Data());
        Edep[i] -> SetAlias("PID", "TMath::Log10(ThinTot*(ThinTot + ThickTot))");
        Edep[i] -> SetAlias("gPID", "TMath::Log10(gThinTot*(gThinTot + gThickTot))");
        Edep[i] -> SetAlias("PID_NoLog", "(ThinTot*(ThinTot + ThickTot))");
        Edep[i] -> SetAlias("gPID_NoLog", "(gThinTot*(gThinTot + gThickTot))");
        
        Edep[i] -> SetAlias("ETot", "(ThinTot + ThickTot)");
        Edep[i] -> SetAlias("gETot", "(gThinTot + gThickTot)");

        AliasFile << "ThinTot = " << ThinTot.Data() << std::endl;
        AliasFile << "ThickTot = " << ThickTot.Data() << std::endl;
        AliasFile << "gThinTot = " << gThinTot.Data() << std::endl;
        AliasFile << "gThickTot = " << gThickTot.Data() << std::endl;
        AliasFile << "PID = " << "TMath::Log10(ThinTot*(ThinTot + ThickTot))" << std::endl;
        AliasFile << "gPID = " << "TMath::Log10(gThinTot*(gThinTot + gThickTot))" << std::endl;
        AliasFile << "ETot = " << "(ThinTot + ThickTot)" << std::endl;
        AliasFile << "gETot = " << "(gThinTot + gThickTot)" << std::endl;


        TotalEnergyPlasticScint = "(";
        gTotalEnergyPlasticScint = "(";
        for(int ii = 0; ii < (N_PLASTIC_NO_VETO+1); ++ii)
        {
            EnergyPlasticScint[ii] = "(";
            gEnergyPlasticScint[ii] = "(";
            
            if(ii==0)
            {
                EnergyPlasticScint[ii]+="0";
                gEnergyPlasticScint[ii]+="0";
            }

            EnergyPlastic2Veto[ii] = "(";
            gEnergyPlastic2Veto[ii] = "(";

            if(ii == 0)
            {
                TotalEnergyPlasticScint += Form("%s",BranchName[jPlastic + ii].Data());
                gTotalEnergyPlasticScint += Form("g%s",BranchName[jPlastic + ii].Data());
            }
            else
            {
                TotalEnergyPlasticScint += Form(" + %s",BranchName[jPlastic + ii].Data());
                gTotalEnergyPlasticScint += Form(" + g%s",BranchName[jPlastic + ii].Data());
            }


            for(int k = 0; k < ii; ++k)
            {
                if(k == 0)
                {
                    EnergyPlasticScint[ii] += Form("%s",BranchName[jPlastic + k].Data());
                    gEnergyPlasticScint[ii] += Form("g%s",BranchName[jPlastic + k].Data());
                }
                else
                {
                    EnergyPlasticScint[ii] += Form(" + %s",BranchName[jPlastic + k].Data());
                    gEnergyPlasticScint[ii] += Form(" + g%s",BranchName[jPlastic + k].Data());
                }
            }

            for(int k = ii; k < (N_PLASTIC_NO_VETO+1); ++k)
            {
                if(k == ii)
                {
                    EnergyPlastic2Veto[ii] += Form("%s",BranchName[jPlastic + k].Data());
                    gEnergyPlastic2Veto[ii] += Form("g%s",BranchName[jPlastic + k].Data());
                }
                else
                {
                    EnergyPlastic2Veto[ii] += Form(" + %s",BranchName[jPlastic + k].Data());
                    gEnergyPlastic2Veto[ii] += Form(" + g%s",BranchName[jPlastic + k].Data());
                }   
            }

            


            EnergyPlasticScint[ii] += ")";
            gEnergyPlasticScint[ii] += ")";
            EnergyPlastic2Veto[ii] += ")";
            gEnergyPlastic2Veto[ii] += ")";

            Edep[i] -> SetAlias(Form("EnergyPlasticScint_until_%d",ii),EnergyPlasticScint[ii].Data());
            Edep[i] -> SetAlias(Form("gEnergyPlasticScint_until_%d",ii),gEnergyPlasticScint[ii].Data());
            Edep[i] -> SetAlias(Form("EnergyPlastic2Veto_%d",ii),EnergyPlastic2Veto[ii].Data());
            Edep[i] -> SetAlias(Form("gEnergyPlastic2Veto_%d",ii),gEnergyPlastic2Veto[ii].Data());

            
            //AliasFile << Form("gEnergyPlasticScint_until_%d",ii) << " = " << gEnergyPlasticScint[ii].Data() << std::endl;
            //AliasFile << Form("EnergyPlastic2Veto_until_%d",ii) << " = " << EnergyPlastic2Veto[ii].Data() << std::endl;
            //AliasFile << Form("gEnergyPlastic2Veto_until_%d",ii) << " = " << gEnergyPlastic2Veto[ii].Data() << std::endl;
        }

        for(int ii = 0; ii < (N_PLASTIC_NO_VETO+1); ++ii)
        {
            AliasFile << Form("EnergyPlasticScint_until_%d",ii) << " = " << EnergyPlasticScint[ii].Data() << std::endl;
            
        }

        for(int ii = 0; ii < (N_PLASTIC_NO_VETO+1); ++ii)
        {
            AliasFile << Form("EnergyPlastic2Veto_%d",ii) << " = " << EnergyPlastic2Veto[ii].Data() << std::endl;
        }


        TotalEnergyPlasticScint += ")";
        gTotalEnergyPlasticScint += ")";

        Edep[i] -> SetAlias("TotalEnergyPlasticScint",TotalEnergyPlasticScint.Data());
        Edep[i] -> SetAlias("gTotalEnergyPlasticScint",gTotalEnergyPlasticScint.Data());
        AliasFile << "TotalEnergyPlasticScint = " << TotalEnergyPlasticScint.Data() << std::endl;
        AliasFile << "gTotalEnergyPlasticScint = " << gTotalEnergyPlasticScint.Data() << std::endl;



        /* -------------------------------- SMEARING -------------------------------- */
        /* ------------------------ YOU ONLY NEED TO ADD A g ------------------------ */
        for(int k = jPlastic; k < N_BRANCHES; ++k)
        {
            Edep[i] -> SetAlias(Form("wnorm%d",k),"(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))");
            AliasFile << Form("wnorm%d",k) << " = " << "(sin(2 *pi*rndm)*sqrt(-2*log(rndm)))" << std::endl;
            if((k >= jPlastic && k < jSilicon))
            {
                Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResPlastic));
                AliasFile << Form("g%s",BranchName[k].Data()) << " = " << Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResPlastic) << std::endl;
            }
            else 
            {
                Edep[i] -> SetAlias(Form("g%s",BranchName[k].Data()), Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon));
                AliasFile << Form("g%s",BranchName[k].Data()) << " = " << Form("((%s)*(1 + (wnorm%d * %f)))",BranchName[k].Data(),k,ResSilicon) << std::endl;
            }
        }

        /* ------------------- Incident direction of the particle ------------------- */

        Edep[i] -> SetAlias("NormP", Form("(TMath::Sqrt(TMath::Power(%s,2) + TMath::Power(%s,2) + TMath::Power(%s,2)))", BranchName[6].Data(), BranchName[7].Data(), BranchName[8].Data()));
        AliasFile << "NormP = " << Form("(TMath::Sqrt(TMath::Power(%s,2) + TMath::Power(%s,2) + TMath::Power(%s,2)))", BranchName[6].Data(), BranchName[7].Data(), BranchName[8].Data()) << std::endl;
        for(int k = 6; k <= 8; ++k)
        {
            Edep[i] -> SetAlias(DirName[k-6].Data(), Form("((%s)/(NormP))", BranchName[k].Data()));
            AliasFile << DirName[k-6].Data() << " = " << Form("((%s)/(NormP))", BranchName[k].Data()) << std::endl;
        }

        // Legend
        // [0] : Theta
        // [1] : Phi
        Edep[i] -> SetAlias(PolarAngle[0].Data()   , Form("((TMath::ACos(%s)))", DirName[2].Data())); // *(180/3.415927))
        Edep[i] -> SetAlias(PolarAngle[1].Data()   , Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data())); // *(180/3.415927))
        Edep[i] -> SetAlias(NewPolarAngle[1].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(TMath::Cos(%s))))",              PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() )); // *(180/3.415927))
        Edep[i] -> SetAlias(NewPolarAngle[0].Data(), Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)),TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data())); // *(180/3.415927))
        
        AliasFile << PolarAngle[0].Data() << " = " << Form("((TMath::ACos(%s)))", DirName[2].Data()) << std::endl;
        AliasFile << PolarAngle[1].Data() << " = " << Form("(TMath::ATan2(%s,%s))", DirName[1].Data(), DirName[0].Data()) << std::endl;
        AliasFile << NewPolarAngle[1].Data() << " = " << Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Cos(%s)),(TMath::Cos(%s))))",              PolarAngle[0].Data(), PolarAngle[0].Data(), PolarAngle[1].Data() ) << std::endl;
        AliasFile << NewPolarAngle[0].Data() << " = " << Form("(TMath::ATan2((TMath::Sin(%s)*TMath::Sin(%s)*TMath::Sin(%s)),TMath::Cos(%s)))", PolarAngle[0].Data(), PolarAngle[1].Data(), NewPolarAngle[1].Data(), PolarAngle[0].Data()) << std::endl;



        /* -------------------- Particle identification parameter ------------------- */

        j = jSilicon;
        CopyNumber = 0;
        std::cout << "Defining total energy in File " << i << " : " << FileName[i].Data() << std::endl;
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
                AliasFile << TotEnergyName[CopyNumber] << " = " << TotEnergyCondition[CopyNumber] << std::endl;
                AliasFile << PIDName[CopyNumber] << " = " << Form("(TMath::Log10(g%s * %s))",BranchName[j].Data(), TotEnergyName[CopyNumber].Data()) << std::endl;
                ++j;
                ++CopyNumber;
                std::cout << Form("ix= %d iy= %d CopyNumber= %d",ix, iy, CopyNumber) << std::endl;
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
    TString ConditionVetoBottom[N_PLASTIC_NO_VETO+1];
    TString ConditionLateralDrilled = Form("(%s < %g)", LateralDrilled.Data(), E_th_Vetoed);

    ConditionFile << "Lateral Drilled condition: \n" << ConditionLateralDrilled.Data() << std::endl;


    for(int i = 0; i< Ntot; ++i)
    {
        std::cout << "Defining good events for pair" << BranchName[jSilicon+i].Data() << " & " <<BranchName[jSilicon+ Ntot+i].Data() << endl;
        ConditionPairSilicon[i] = Form("((%s > %g) && (%s > %g))", BranchName[jSilicon+i].Data(), E_min_thin, BranchName[jSilicon + Ntot + i].Data(), E_min_thick);
        ConditionFile << "Good Events for Pair  "<< i << "\n" << ConditionPairSilicon[i].Data() << std::endl;
        
        std::cout << "ConditionPairSilicon[" << i << "] = " << ConditionPairSilicon[i].Data() << std::endl;
        if(i == 0)
        {
            ConditionPairSiliconAll = Form("((%s > %g) && (%s > %g))", BranchName[jSilicon+i].Data(), E_min_thin, BranchName[jSilicon + Ntot + i].Data(), E_min_thick);
        }
        else 
        {
            ConditionPairSiliconAll += Form("|| ((%s > %g) && (%s > %g))", BranchName[jSilicon+i].Data(), E_min_thin, BranchName[jSilicon + Ntot + i].Data(), E_min_thick);
        }
    }

    for(int i = 0; i < (N_PLASTIC_NO_VETO+1); ++i)
    {   // Form("EnergyPlasticScint_until_%d",ii)
        ConditionVetoBottom[i] = Form("(EnergyPlastic2Veto_%d) < %g", i, E_th_Vetoed);
        ConditionFile << "Condition for Veto " << i << "\n" << ConditionVetoBottom[i].Data() << std::endl;
    }
    

    std::cout << "ConditionPairSiliconAll = " << ConditionPairSiliconAll.Data() << std::endl;
    ConditionFile << "Good Events for all pairs \n" << ConditionPairSiliconAll.Data() << std::endl;

    ConditionDrilledVeto       = Form("(%s < %g)", BranchName[10].Data(), E_th_Vetoed);
    ConditionFile << "Good Events for drilled Veto \n" << ConditionDrilledVeto.Data() << std::endl;
    //ConditionEnergyDispersion  = Form("((%s) - (%s) - (%s) - (%s) - (%s) - (%s) - (%s)) < %g", BranchName[1].Data(), ThinTot.Data(), ThickTot.Data(), BranchName[9].Data(),BranchName[43].Data(), BranchName[44].Data(), BranchName[45].Data(), E_th_Vetoed);
    ConditionEnergyDispersion  = Form("((%s)" , BranchName[1].Data());
    ConditionEnergyDispersion += Form("- (%s)", ThinTot.Data());
    ConditionEnergyDispersion += Form("- (%s)", ThickTot.Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[9].Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[43].Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[44].Data());
    //ConditionEnergyDispersion += Form("- (%s)", BranchName[45].Data());
    ConditionEnergyDispersion += Form(") < %g", E_th_Vetoed);
    ConditionFile << "Condition: Energy dispersion \n" << ConditionEnergyDispersion.Data() << std::endl;
    
    ConditionGoodEvents        = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
    ConditionFile << "Good Events \n" << ConditionGoodEvents.Data() << std::endl;
    
    
    for(int i = 0; i < Ntot; ++i)
    {
        std::cout << "Defining good events for pair" << BranchName[11+i].Data() << " & " <<BranchName[11+ Ntot+i].Data() << std::endl;
        ConditionGoodEventsSinglePair[i] = Form("(%s) && (%s) && (%s)", ConditionPairSilicon[i].Data(), ConditionDrilledVeto.Data(), ConditionEnergyDispersion.Data());
        ConditionFile << "Good Events for Pair  "<< i << "\n" << ConditionGoodEventsSinglePair[i].Data() << std::endl;
    }


    TString ConditionEnergyDispUntilPlastic[N_PLASTIC_NO_VETO+1];

    cout << "Building the conditions for the energy dispersion until plastic\n------------\n" << endl;

    for(int i = 0; i < (N_PLASTIC_NO_VETO+1); ++i)
    {
        ConditionEnergyDispUntilPlastic[i] = Form("((%s)-(ThinTot) - (ThickTot) - (EnergyPlasticScint_until_%d)) < %g", BranchName[1].Data(), i, E_th_Vetoed);
        cout << "ConditionEnergyDispUntilPlastic[" << i << "] = " << ConditionEnergyDispUntilPlastic[i] << endl;
        ConditionFile << "ConditionEnergyDispUntilPlastic[" << i << "] = " << ConditionEnergyDispUntilPlastic[i] << endl;
    }


    TH1D *hEfficiency[Nfiles][N_PLASTIC_NO_VETO+2];
    THStack *hsEfficiency[Nfiles];

    TH1D *hEfficiencyRaw[Nfiles];

    TCanvas *cEfficiency[Nfiles];


    for(int i = 0; i < Nfiles; ++i)
    {
        cout << "Starting the loop over the files" << endl;

            cout << "Opening file " << i << endl;
            cEfficiency[i] = new TCanvas(Form("cEfficiency_%d",i), Form("cEfficiency_%d",i), 1200, 600);
            hsEfficiency[i] = new THStack(Form("hsEfficiency_%d",i), Form("hsEfficiency_%d",i));
            TLegend *leg;
            switch (i)
            {
            case 0:
                leg = new TLegend(0.6,0.6,0.9,0.9);
                break;
            case 1:
                leg = new TLegend(0.1,0.1,0.3,0.3);
                break;
            case 2:
                leg = new TLegend(0.1,0.1,0.3,0.3);
                break;
            default:
                break;
            }


            int NBin = 140;

            //hEfficiency[i][0] = new TH1D(Form("hEfficiency_%d_%d",i,j), Form("hEfficiency_%d_%d",i,j), NBin, ERangeGen[i][0], ERangeGen[i][1]);
            //Edep[i] -> Draw(Form("%s >> hEfficiency_%d_%d", BranchName[1].Data(), i, j), Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispersion.Data()), "goff");
            //hEfficiency[i][0] -> SetLineColor(1);
            //hEfficiency[i][0] -> SetMarkerColor(1);
            //hEfficiency[i][0] -> SetMarkerStyle(20);
            //hEfficiency[i][0] -> SetLineWidth(3);


            //hEfficiency[i][0] -> Scale(1./NPerBin);
            //hsEfficiency[i] -> Add(hEfficiency[i][0], "histolp");
            //leg -> AddEntry(hEfficiency[i][0], "Si + CZT", "lp");

            hEfficiencyRaw[i] = new TH1D(Form("hEfficiencyRaw_%d",i), Form("hEfficiencyRaw_%d",i), NBin, ERangeGen[i][0], ERangeGen[i][1]);
            Edep[i] -> Draw(Form("%s >> hEfficiencyRaw_%d", BranchName[1].Data(), i), "", "goff");


            for(j = 0; j < (N_PLASTIC_NO_VETO+1); ++j)
            {
                //j = 4;

                cout << "Opening file " << i << " and plastic " << j << endl;
                hEfficiency[i][j+1] = new TH1D(Form("hEfficiency_%d_%d",i,j), Form("hEfficiency_%d_%d",i,j), NBin, ERangeGen[i][0], ERangeGen[i][1]);
                Edep[i] -> Draw(Form("%s >> hEfficiency_%d_%d", BranchName[1].Data(), i, j), Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[j].Data(), ConditionLateralDrilled.Data()), "goff");
                hEfficiency[i][j+1] -> SetLineColor(j+2);
                hEfficiency[i][j+1] -> SetMarkerColor(j+2);
                hEfficiency[i][j+1] -> SetMarkerStyle(20);


                if((j+2) == 10)
                {
                    hEfficiency[i][j+1] -> SetLineColor(30);
                    hEfficiency[i][j+1] -> SetMarkerColor(30);
                }

                hEfficiency[i][j+1] -> SetLineWidth(3);
                hEfficiency[i][j+1] -> Divide(hEfficiencyRaw[i]);
                
                hsEfficiency[i] -> Add(hEfficiency[i][j+1], "histolp");
                


                int NPlastic = j;
                if(NPlastic == 0)
                {
                    leg -> AddEntry(hEfficiency[i][j+1], "Si+CZT", "lp");
                }
                else
                {
                    leg -> AddEntry(hEfficiency[i][j+1], Form("Si+CZT+PS %d mm", 10 *NPlastic), "lp");
                }
                
            }

            hsEfficiency[i] -> Draw("nostack");
            hsEfficiency[i] -> GetXaxis() -> SetTitle("Energy [MeV]");
            hsEfficiency[i] -> GetYaxis() -> SetTitle("Efficiency");
            leg -> Draw();

            //gPad -> SetLogy();
            gPad -> SetGrid();

            cEfficiency[i] -> SaveAs(Form("../docs/LEM_100MeV/Efficiency_%d.pdf", i));
            cEfficiency[i] -> SaveAs(Form("../docs/LEM_100MeV/Efficiency_%d.png", i));
            cEfficiency[i] -> SaveAs(Form("../docs/LEM_100MeV/Efficiency_%d.root", i));
            

    }

 

    TCanvas *cPID = new TCanvas("cPID", "cPID", 1200, 600);
    cPID -> Divide(2,1);

    TH2D *hPID[2];

    int NbinsX = 300;
    int NbinsY = 300;
    double Xmin = log10(0.05);
    double Xmax = log10(600);
    double Ymin = -3;
    double Ymax = 3;

    hPID[0] = new TH2D("hPID_0", "hPID_0", NbinsX, Xmin, Xmax, NbinsY, Ymin, Ymax);
    hPID[1] = new TH2D("hPID_1", "hPID_1", NbinsX, Xmin, Xmax, NbinsY, Ymin, Ymax);


    cPID -> cd(1);

    for(int i = 0; i < Nfiles; ++i)
    {
        if(i == 0)
        {
            // Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[j].Data(), ConditionLateralDrilled.Data())
            //Edep[i] -> Draw("gPID:log10(gETot)>>hPID_0", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispersion.Data()), "colz");
            Edep[i] -> Draw("gPID:log10(gETot)>>hPID_0", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[0].Data(), ConditionLateralDrilled.Data()) , "colz");

        }
        else
        {
            //Edep[i] -> Draw("gPID:log10(gETot)>>+hPID_0", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispersion.Data()), "colz");
            Edep[i] -> Draw("gPID:log10(gETot)>>+hPID_0", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[0].Data(), ConditionLateralDrilled.Data()) , "colz");
        }
    }

    hPID[0] -> GetXaxis() -> SetTitle("log_{10}(E_{tot} / 1 MeV)");
    hPID[0] -> GetYaxis() -> SetTitle("PID");
    hPID[0] -> SetTitle("Si 100 um + CZT 5 mm");
    hPID[0] -> SetStats(0);

    gPad -> SetLogz();
    gPad -> SetGrid();

    cPID -> cd(2);

    for(int i = 0; i < Nfiles; ++i)
    {
        if(i == 0)
        {
            //Edep[i] -> Draw("gPID:log10(gETot)>>hPID_1", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispUntilPlastic[8].Data()), "colz");
            Edep[i] -> Draw("log10((gThinTot)*(gThinTot+gThickTot+TotalEnergyPlasticScint)):log10((gThinTot+gThickTot+TotalEnergyPlasticScint))>>hPID_1", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[N_PLASTIC_NO_VETO].Data(), ConditionLateralDrilled.Data()) , "colz");
            ConditionFile << Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[9].Data(), ConditionLateralDrilled.Data()) << endl;
        }
        else
        {// TotalEnergyPlasticScint
            //Edep[i] -> Draw("log10((gThinTot)*(gThinTot+gThickTot+gEnergyPlasticScint_until_7)):log10((gThinTot+gThickTot+gEnergyPlasticScint_until_7))>>+hPID_1", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispUntilPlastic[7].Data()), "colz");
            Edep[i] -> Draw("log10((gThinTot)*(gThinTot+gThickTot+TotalEnergyPlasticScint)):log10((gThinTot+gThickTot+TotalEnergyPlasticScint))>>+hPID_1", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[N_PLASTIC_NO_VETO].Data(), ConditionLateralDrilled.Data()) , "colz");
        }
    }

    hPID[1] -> GetXaxis() -> SetTitle("log_{10}(E_{tot} / 1 MeV)");
    hPID[1] -> GetYaxis() -> SetTitle("PID");
    hPID[1] -> SetTitle("Si 100 um + CZT 5 mm + PS 50 mm");
    hPID[1] -> SetStats(0);

    gPad -> SetLogz();
    gPad -> SetGrid();

    cPID -> SaveAs("../docs/LEM_100MeV/PID.pdf");
    cPID -> SaveAs("../docs/LEM_100MeV/PID.png");
    cPID -> SaveAs("../docs/LEM_100MeV/PID.root");



    TCanvas* cPIDLog = new TCanvas("cPIDLog", "cPIDLog", 1200, 600);

    cPIDLog -> Divide(2,1);

    cPIDLog -> cd(1);

    NbinsX = 300;
    NbinsY = 300;
    Xmin = 0.05;
    Xmax = 600.;
    Ymin = 0.001;
    Ymax = 1000.;

    
    
    
    TH2DLog *hPIDLog[2];
    TH2D *hPID2[2];

    hPIDLog[0] = new TH2DLog();
    hPIDLog[0] -> SetName("hPIDLog_0");
    hPIDLog[0] -> SetTitle("Si 100 um + CZT 5 mm");
    hPIDLog[0] -> SetXTitle("E_{tot} [MeV]");
    hPIDLog[0] -> SetYTitle("PID");
    hPIDLog[0] -> SetXAxis(Xmin, Xmax, NbinsX);
    hPIDLog[0] -> SetYAxis(Ymin, Ymax, NbinsY);
    hPIDLog[0] -> GenerateHistogram();
    hPID2[0] = (TH2D*) hPIDLog[0] -> GetHistogram();
    hPID2[0] -> SetStats(0);


    for(int i = 0; i < Nfiles; ++i)
    {
        if(i == 0)
        {
            //Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot)):gETot>>hPIDLog_0", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispersion.Data()), "colz");
            Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot)):gETot>>hPIDLog_0", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[0].Data(), ConditionLateralDrilled.Data()) , "colz");
            ConditionFile << Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[0].Data(), ConditionLateralDrilled.Data()) << endl;
        
        }
        else
        {
            //Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot)):gETot>>+hPIDLog_0", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispersion.Data()), "colz");
            Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot)):gETot>>+hPIDLog_0", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[0].Data(), ConditionLateralDrilled.Data()) , "colz");
        }
    }
    hPID2[0] -> Draw("colz");
    gPad -> SetLogz();
    gPad -> SetGrid();
    gPad -> SetLogx();
    gPad -> SetLogy();

    cPIDLog -> cd(2);

    hPIDLog[1] = new TH2DLog();
    hPIDLog[1] -> SetName("hPIDLog_1");
    hPIDLog[1] -> SetTitle("Si 100 um + CZT 5 mm + PS 40 mm");
    hPIDLog[1] -> SetXTitle("E_{tot} [MeV]");
    hPIDLog[1] -> SetYTitle("PID");
    hPIDLog[1] -> SetXAxis(Xmin, Xmax, NbinsX);
    hPIDLog[1] -> SetYAxis(Ymin, Ymax, NbinsY);
    hPIDLog[1] -> GenerateHistogram();
    hPID2[1] = (TH2D*) hPIDLog[1] -> GetHistogram();
    hPID2[1] -> SetStats(0);

    for(int i = 0; i < Nfiles; ++i)
    {
        if(i == 0)
        {
            //Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot+gEnergyPlasticScint_until_9)):(gThinTot+gThickTot+gEnergyPlasticScint_until_9)>>hPIDLog_1", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispUntilPlastic[8].Data()), "colz");
            Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot+TotalEnergyPlasticScint)):(gThinTot+gThickTot+TotalEnergyPlasticScint)>>hPIDLog_1", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[N_PLASTIC_NO_VETO].Data(), ConditionLateralDrilled.Data()) , "colz");
            ConditionFile << Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[8].Data(), ConditionLateralDrilled.Data()) << endl;
        }
        else
        {
            //Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot+gEnergyPlasticScint_until_7)):(gThinTot+gThickTot+gEnergyPlasticScint_until_7)>>+hPIDLog_1", Form("(%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionEnergyDispUntilPlastic[8].Data()), "colz");
            Edep[i] -> Draw("((gThinTot)*(gThinTot+gThickTot+TotalEnergyPlasticScint)):(gThinTot+gThickTot+TotalEnergyPlasticScint)>>+hPIDLog_1", Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[N_PLASTIC_NO_VETO].Data(), ConditionLateralDrilled.Data()) , "colz");
        }
    }

    hPID2[1] -> Draw("colz");
    gPad -> SetLogz();
    gPad -> SetGrid();
    gPad -> SetLogx();
    gPad -> SetLogy();






    return 0;
}