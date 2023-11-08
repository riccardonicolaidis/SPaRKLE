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
#include "../ROOT_macros/Custom/TH1DLog.h"


#define Nx 3
#define Ny 3
#define N_PLASTIC_NO_VETO 5
#define N_PASTIC_LATERAL 4
#define Ntot Nx*Ny
#define N_BRANCHES 9+2*Ntot+N_PLASTIC_NO_VETO+N_PASTIC_LATERAL+2



using namespace std;




int Acceptance_LEM100()
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

    TH1DLog *hGenLog[Nfiles];
    TH1DLog *hAccLog[Nfiles];

    TFile *fOut = new TFile("out.root", "RECREATE");

    TH1D *hGen[Nfiles];
    TH1D *hAcc[Nfiles];

    TCanvas *cAcceptance[Nfiles];

    for(int i = 0; i < Nfiles; ++i)
    {
        hGenLog[i] = new TH1DLog();
        hAccLog[i] = new TH1DLog();

        hGenLog[i] -> SetName(Form("hGenLog_%d", i));
        hGenLog[i] -> SetTitle("Generated");
        hGenLog[i] -> SetXTitle("Energy [MeV]");
        hGenLog[i] -> SetYTitle("Counts");
        hGenLog[i] -> SetXAxis(ERangeGen[i][0], ERangeGen[i][1], 100);


        hAccLog[i] -> SetName(Form("hAccLog_%d", i));
        hAccLog[i] -> SetTitle("Accepted");
        hAccLog[i] -> SetXTitle("Energy [MeV]");
        hAccLog[i] -> SetYTitle("Geometric Factor [cm^{2} sr]");
        hAccLog[i] -> SetXAxis(ERangeGen[i][0], ERangeGen[i][1], 100);

        hGenLog[i] -> GenerateHistogram();
        hAccLog[i] -> GenerateHistogram();

        hGen[i] = hGenLog[i]->GetHistogram();
        hAcc[i] = hAccLog[i]->GetHistogram();

        cout << "Drawing Generated for file " << i << endl;
        Edep[i] -> Draw(Form("%s>>hGenLog_%d", BranchName[1].Data(), i), "", "goff");
        TString condition = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[N_PLASTIC_NO_VETO].Data(), ConditionLateralDrilled.Data());
        cout << "Condition for file " << i << " = " << condition.Data() << endl;
        cout << "Drawing Accepted for file " << i << endl;
        Edep[i] -> Draw(Form("%s>>hAccLog_%d", BranchName[1].Data(), i), condition.Data() , "goff");

        hAcc[i] -> Sumw2();
        hGen[i] -> Sumw2();

        hAcc[i] -> Divide(hGen[i]);

        double LGeneration = 4.; // cm

        hAcc[i] -> Scale(3.*LGeneration*LGeneration*TMath::Pi());
        hAcc[i] -> Write();

        cAcceptance[i] = new TCanvas(Form("cAcceptance_%d", i), Form("cAcceptance_%d", i));

        hAcc[i] -> SetMarkerStyle(20);
        hAcc[i] -> SetMarkerSize(0.5);
        hAcc[i] -> SetMarkerColor(kBlue);
        hAcc[i] -> SetLineColor(kBlue);
        hAcc[i] -> SetLineWidth(2);
        hAcc[i] -> SetStats(0);

        hAcc[i] -> Draw("lp");
        gPad -> SetLogx();
        gPad -> SetLogy();
        gPad -> SetGrid();

        cAcceptance[i] -> SaveAs(Form("Acceptance_%d.png", i));
        cAcceptance[i] -> Write();

    }


    int Colors[3] = {kAzure+9, kRed-4, kViolet-5};
    TString ParticleNames[3] = {"Electrons", "Protons", "Alpha"};

    TCanvas *cAcceptanceAll = new TCanvas("cAcceptanceAll", "cAcceptanceAll");
    TLegend *lAcceptanceAll = new TLegend(0.1, 0.1, 0.3, 0.4);


    TH1DLog *hAccLogAll[Nfiles];
    TH1DLog *hGenLogAll[Nfiles];
    TH1D *hAccAll[Nfiles];
    TH1D *hGenAll[Nfiles];

    THStack *hsAcceptance = new THStack("hsAcceptance", "hsAcceptance");
    

    for(int i = 0; i < Nfiles; ++i)
    {




        hGenLogAll[i] = new TH1DLog();
        hGenLogAll[i] -> SetName(Form("hGenLog2_%d", i));
        hGenLogAll[i] -> SetTitle("Generated");
        hGenLogAll[i] -> SetXTitle("Energy [MeV]");
        hGenLogAll[i] -> SetYTitle("Counts");
        hGenLogAll[i] -> SetXAxis(0.09, 500, 110);
        hGenLogAll[i] -> GenerateHistogram();
        hGenAll[i] = hGenLogAll[i]->GetHistogram();

        hAccLogAll[i] = new TH1DLog();
        hAccLogAll[i] -> SetName(Form("hAccLog2_%d", i));
        hAccLogAll[i] -> SetTitle("Accepted");
        hAccLogAll[i] -> SetXTitle("Energy [MeV]");
        hAccLogAll[i] -> SetYTitle("Geometric Factor [cm^{2} sr]");
        hAccLogAll[i] -> SetXAxis(0.09, 500, 110);
        hAccLogAll[i] -> GenerateHistogram();
        hAccAll[i] = hAccLogAll[i]->GetHistogram();


        cout << "Drawing Generated for file " << i << endl;
        Edep[i] -> Draw(Form("%s>>hGenLog2_%d", BranchName[1].Data(), i), "", "goff");
        
        TString condition = Form("(%s) && (%s) && (%s)", ConditionPairSiliconAll.Data(), ConditionVetoBottom[N_PLASTIC_NO_VETO].Data(), ConditionLateralDrilled.Data());
        cout << "Condition for file " << i << " = " << condition.Data() << endl;
        cout << "Drawing Accepted for file " << i << endl;

        Edep[i] -> Draw(Form("%s>>hAccLog2_%d", BranchName[1].Data(), i), condition.Data() , "goff");

        hAccAll[i] -> Sumw2();
        hGenAll[i] -> Sumw2();

        hAccAll[i] -> Divide(hGenAll[i]);

        double LGeneration = 4.; // cm

        hAccAll[i] -> Scale(3.*LGeneration*LGeneration*TMath::Pi());
        hAccAll[i] -> Write();


        hAccAll[i] -> SetMarkerStyle(20);
        hAccAll[i] -> SetMarkerSize(0.7);
        hAccAll[i] -> SetMarkerColor(Colors[i]);
        hAccAll[i] -> SetLineColor(Colors[i]);
        hAccAll[i] -> SetLineWidth(3);
        hAccAll[i] -> SetStats(0);

        hsAcceptance -> Add(hAccAll[i]);
        lAcceptanceAll -> AddEntry(hAccAll[i], ParticleNames[i].Data(), "lp");

    }

    //hsAcceptance -> GetXaxis() -> SetTitle("Energy [MeV]");
    //hsAcceptance -> GetYaxis() -> SetTitle("Geometric Factor [cm^{2} sr]");

    hsAcceptance -> SetTitle("Geometric Factor; Energy [MeV]; Geometric Factor [cm^{2} sr]");

    hsAcceptance -> Draw("nostack");
    gPad -> SetLogx();
    gPad -> SetLogy();
    gPad -> SetGrid();
    lAcceptanceAll -> Draw("same");

    hsAcceptance -> Write();
    cAcceptanceAll -> Write();

    //fOut -> Close();


    return 0;
}