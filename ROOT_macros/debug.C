void debug()
{

    TString FileName = "100um_500um_mu.root";
    TFile *file_debug = TFile::Open(FileName);
    TTree *Edep = (TTree*) file_debug -> Get("Edep");
    TString DrawInstruction;

    Edep -> Print();

    // ##################################
    //  STRING DEFINITIONS FROM THE TREE
    // ##################################
    
    TString EnergyEntries[4];
    EnergyEntries[0] = "RandNumber";
    EnergyEntries[1] = "RandEnergy";
    EnergyEntries[2] = "Ed_Veto";
    EnergyEntries[3] = "Ed_DrilledVeto";
    

    TString PositionEntries[3];
    PositionEntries[0] = "Xgen";
    PositionEntries[1] = "YGen";
    PositionEntries[2] = "ZGen";
    
    TString MomentumEntries[3];
    MomentumEntries[0] = "pDirX";
    MomentumEntries[1] = "pDirY";
    MomentumEntries[2] = "pDirZ";

    int Nx   = 4;
    int Ny   = 4;
    int Ntot = Nx*Ny;

    TString ThinEntries[Ntot];
    TString ThickEntries[Ntot];
    
    int nID = 0;

    for(int ix = 0; ix < Nx ; ++ix)
    {
        for(int iy = 0; iy < Ny ; ++iy)
        {
            ThinEntries[nID] = "Thin_x" + std::to_string(ix) + "_y" + std::to_string(iy) + "_ID" + std::to_string(nID);
            ThickEntries[nID] = "Thick_x" + std::to_string(ix) + "_y" + std::to_string(iy) + "_ID" + std::to_string(nID);
            cout << ThinEntries[nID] << endl;
            cout << ThickEntries[nID] << endl;
            ++nID;
        }
    }

    
    // ##################################
    // ##################################
    

    TCanvas *c1 = new TCanvas();
    TH1D *h1 = new TH1D("h1", "", 100, 1e-5, 10);
    DrawInstruction = EnergyEntries[2] + ">>h1";
    Edep -> Draw(DrawInstruction, "");

    TF1 *f1 = new TF1("f1", "landau +[3]*TMath::ATan([4]*x - [5]) + [6]", 0.001, 10);
    f1 -> SetParameter(0,100);
    f1 -> SetParameter(1,2.);
    f1 -> SetParameter(3,50);
    h1 -> Fit("f1");



    return;
}