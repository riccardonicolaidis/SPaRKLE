#include "TH1DLog.h"

#include <iostream>
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TF1.h"




int TestTH1DLog()
{
    // This code generates a histogram with a log-uniform x-axis.
// It then fills the histogram with log-uniform random numbers.
// Finally, it fits a function to the histogram.

// Create the histogram
TH1DLog *hist = new TH1DLog();
double Xmin = 0.1;
double Xmax = 1000;
int NbinsX = 100;
hist->SetXAxis(Xmin, Xmax, NbinsX);
hist->SetName("hist");
hist->SetTitle("hist Test");
hist->SetXTitle("Random Number Log Uniform");
hist->SetYTitle("Counts");
hist->GenerateHistogram();

// Fill the histogram
TH1D *histogram = hist->GetHistogram();
int TotNumber = 1000000;
double XminGen = 0.2;
double XmaxGen = 700;
for (int i = 0; i < TotNumber; i++)
{
    double x = pow(10., log10(XminGen) + (log10(XmaxGen) - log10(XminGen)) * gRandom->Rndm());
    histogram->Fill(x);
}

// Plot the histogram
TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
histogram->Draw("lep");
gPad -> SetGrid();
gPad -> SetLogx();

// Fit a linear function to the histogram
TF1 *f1 = new TF1("f1", "[0]+[1]*x", XminGen+10, XmaxGen-10);
histogram->Fit("f1", "R", "", XminGen+10, XmaxGen-10);
f1->Draw("same");

return 0;
}
