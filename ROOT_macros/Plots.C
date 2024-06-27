#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TMath.h"
#include "TLegend.h"
#include "TH1DLog.h"
#include "TChain.h"
#include "TLine.h"
#include "TText.h"


int Plots()
{
    TString path_geomfact = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/ROOT_macros/GeomFactor.root";
    TString path_EffArea = "/home/riccardo/Documenti/GeantProjects/SPaRKLE/ROOT_macros/EffArea.root";

    TFile *file_geomfact = new TFile(path_geomfact);
    TFile *file_EffArea = new TFile(path_EffArea);

    TCanvas *cgeom = (TCanvas*)file_geomfact->Get("c1");
    // Change the name of the canvas
    cgeom->SetName("cgeom");
    cgeom->Draw();

    TCanvas *cEffArea = (TCanvas*)file_EffArea->Get("c");
    // Change the name of the canvas

    cEffArea->SetName("cEffArea");
    cEffArea->Draw();


    TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
    c2->Divide(2,1);
    c2->cd(1);
    cgeom->DrawClonePad();
    c2->cd(2);
    cEffArea->DrawClonePad();
    



    return 0;
}