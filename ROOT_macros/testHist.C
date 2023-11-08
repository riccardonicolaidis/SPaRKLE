#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"



void testHist()
{
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TH1F *h1 = new TH1F("h1","h1",100,-5.,5.);
  h1->FillRandom("gaus",100000);
  h1->Draw();
  gPad -> SetGrid();

    // Resize font
    h1->GetXaxis()->SetTitleSize(0.1);
    h1->GetYaxis()->SetTitleSize(10.3);
    h1->GetXaxis()->SetLabelSize(0.06);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetTitleOffset(1.2);
    h1->GetYaxis()->SetTitleOffset(1.2);
    


}