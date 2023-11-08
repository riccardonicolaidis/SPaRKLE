#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"

#include "TRandom2.h"


using namespace std;


int TestLogUniform()
{
    int N = 1000000;

    double x0 = 2.;
    double x1 = 1000.;

    double logx0 = TMath::Log10(x0);
    double logx1 = TMath::Log10(x1);

    TRandom2 *rnd = new TRandom2(0);
    TTree *tree = new TTree("tree", "tree");
    double rndVal;

    tree->Branch("rndVal", &rndVal, "rndVal/D");




    for(int i = 0; i < N; ++i)
    {
        double logRnd = logx0 + rnd->Rndm()*(logx1 - logx0);
        rndVal = pow(10.,logRnd);
        tree->Fill();
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TH1D *h1 = new TH1D("h1", "h1", 2000, x0, x1);
    tree->Draw("rndVal>>h1");
    gPad -> SetLogx();
    gPad -> SetLogy();

    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    TH1D *h2 = new TH1D("h2", "h2", 1000, logx0, logx1);
    tree->Draw("TMath::Log10(rndVal)>>h2", "");





    return 0;
}