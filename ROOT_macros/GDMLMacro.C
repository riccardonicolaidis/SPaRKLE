#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

#include "TGeoManager.h"


int GDMLMacro()
{

    TGeoManager::Import("LEM_geometry_08-12-2022_16-27-45.gdml");


    gGeoManager->GetTopVolume()->Draw("ogl");
    return 0;
}