#ifndef TH2DLOGX
#define TH2DLOGX

#include <iostream>
#include "TH2D.h"
#include "TMath.h"
#include "TString.h"

class TH2DLogX
{
    public:

        TH2DLogX();
        ~TH2DLogX();

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
        

};

TH2DLogX::TH2DLogX()
{}

TH2DLogX::~TH2DLogX()
{}

void TH2DLogX::SetXAxis(double xMinV, double xMaxV, int NbinsXV)
{
    xMin = xMinV;
    xMax = xMaxV;
    NbinsX = NbinsXV;

    xMinLog = TMath::Log10(xMin);
    xMaxLog = TMath::Log10(xMax);
}

void TH2DLogX::SetYAxis(double yMinV, double yMaxV, int NbinsYV)
{
    yMin = yMinV;
    yMax = yMaxV;
    NbinsY = NbinsYV;
}


void TH2DLogX::SetName(TString NameV)
{
    Name = NameV;
}

void TH2DLogX::SetTitle(TString TitleV)
{
    Title = TitleV;
}

void TH2DLogX::SetXTitle(TString XTitleV)
{
    XTitle = XTitleV;
}

void TH2DLogX::SetYTitle(TString YTitleV)
{
    YTitle = YTitleV;
}

void TH2DLogX::GenerateHistogram()
{
    std::vector<double> xBins;
    std::vector<double> yBins;

    std::vector<double> xBinsLog;

    int NEdgesX = NbinsX + 1;
    int NEdgesY = NbinsY + 1;

    double xBinWidthLog = (xMaxLog - xMinLog) / NbinsX;

    for (int i = 0; i < NEdgesX; i++)
    {
        xBinsLog.push_back(xMinLog + i * xBinWidthLog);
        xBins.push_back(TMath::Power(10, xBinsLog[i]));
    }

    for (int i = 0; i < NEdgesY; i++)
    {
        yBins.push_back(yMin + i * (yMax - yMin) / NbinsY);
    }   

    histogram = new TH2D(Name, Title, NbinsX, &xBins[0], NbinsY, &yBins[0]);
    histogram->GetXaxis()->SetTitle(XTitle);
    histogram->GetYaxis()->SetTitle(YTitle);
}


TH2D* TH2DLogX::GetHistogram()
{
    return histogram;
}





#endif