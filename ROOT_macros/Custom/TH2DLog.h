#ifndef TH2DLOG
#define TH2DLOG

#include <iostream>
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

class TH2DLog
{
    public:

        TH2DLog();
        ~TH2DLog();

        void  SetXAxis(double, double, int);
        void  SetYAxis(double, double, int);
        void  SetName(TString);
        void  SetTitle(TString);
        void  SetXTitle(TString);
        void  SetYTitle(TString);
        
        void  GenerateHistogram();
        TH2D* GetHistogram();


    private:
        TH2D *histogram;

        TString Name;
        TString Title;
        TString XTitle;
        TString YTitle;

        double xMin;
        double xMax;
        double yMin;
        double yMax;
        int NbinsX;
        int NbinsY;
        double xMinLog;
        double xMaxLog;
        double yMinLog;
        double yMaxLog;

};

TH2DLog::TH2DLog()
{}

TH2DLog::~TH2DLog()
{}

void TH2DLog::SetXAxis(double xMinV, double xMaxV, int NbinsXV)
{
    xMin = xMinV;
    xMax = xMaxV;
    NbinsX = NbinsXV;

    xMinLog = TMath::Log10(xMin);
    xMaxLog = TMath::Log10(xMax);
}

void TH2DLog::SetYAxis(double yMinV, double yMaxV, int NbinsYV)
{
    yMin = yMinV;
    yMax = yMaxV;
    NbinsY = NbinsYV;

    yMinLog = TMath::Log10(yMin);
    yMaxLog = TMath::Log10(yMax);
}


void TH2DLog::SetName(TString NameV)
{
    Name = NameV;
}

void TH2DLog::SetTitle(TString TitleV)
{
    Title = TitleV;
}

void TH2DLog::SetXTitle(TString XTitleV)
{
    XTitle = XTitleV;
}

void TH2DLog::SetYTitle(TString YTitleV)
{
    YTitle = YTitleV;
}

void TH2DLog::GenerateHistogram()
{
    std::vector<double> xBins;
    std::vector<double> yBins;

    std::vector<double> xBinsLog;
    std::vector<double> yBinsLog;

    int NEdgesX = NbinsX + 1;
    int NEdgesY = NbinsY + 1;

    double xBinWidthLog = (xMaxLog - xMinLog) / NbinsX;
    double yBinWidthLog = (yMaxLog - yMinLog) / NbinsY;

    for (int i = 0; i < NEdgesX; i++)
    {
        xBinsLog.push_back(xMinLog + i * xBinWidthLog);
        xBins.push_back(TMath::Power(10, xBinsLog[i]));
    }

    for (int i = 0; i < NEdgesY; i++)
    {
        yBinsLog.push_back(yMinLog + i * yBinWidthLog);
        yBins.push_back(TMath::Power(10, yBinsLog[i]));
    }

    histogram = new TH2D(Name, Title, NbinsX, &xBins[0], NbinsY, &yBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);
}


TH2D* TH2DLog::GetHistogram()
{
    return histogram;
}





#endif