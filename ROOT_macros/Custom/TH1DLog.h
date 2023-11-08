#ifndef TH1DLOG
#define TH1DLOG

#include <iostream>
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"

class TH1DLog
{
    public:

        TH1DLog();
        ~TH1DLog();

        void  SetXAxis(double, double, int);
        void  SetName(TString);
        void  SetTitle(TString);
        void  SetXTitle(TString);
        void SetYTitle(TString);
        
        void  GenerateHistogram();
        TH1D* GetHistogram();


    private:
        TH1D *histogram;

        TString Name;
        TString Title;
        TString XTitle;
        TString YTitle;
        
        double xMin;
        double xMax;
        int NbinsX;
        double xMinLog;
        double xMaxLog;
};

TH1DLog::TH1DLog()
{}

TH1DLog::~TH1DLog()
{}

void TH1DLog::SetXAxis(double xMinV, double xMaxV, int NbinsXV)
{
    xMin = xMinV;
    xMax = xMaxV;
    NbinsX = NbinsXV;

    xMinLog = TMath::Log10(xMin);
    xMaxLog = TMath::Log10(xMax);
}


void TH1DLog::SetName(TString NameV)
{
    Name = NameV;
}

void TH1DLog::SetTitle(TString TitleV)
{
    Title = TitleV;
}

void TH1DLog::SetXTitle(TString XTitleV)
{
    XTitle = XTitleV;
}

void TH1DLog::SetYTitle(TString YTitleV)
{
    YTitle = YTitleV;
}

void TH1DLog::GenerateHistogram()
{
    std::vector<double> xBins;
    std::vector<double> xBinsLog;
    
    int NEdgesX = NbinsX + 1;
    
    double xBinWidthLog = (xMaxLog - xMinLog) / NbinsX;
    
    for (int i = 0; i < NEdgesX; i++)
    {
        xBinsLog.push_back(xMinLog + i * xBinWidthLog);
        xBins.push_back(TMath::Power(10, xBinsLog[i]));
    }

    
    histogram = new TH1D(Name, Title, NbinsX, &xBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);
}


TH1D* TH1DLog::GetHistogram()
{
    return histogram;
}





#endif