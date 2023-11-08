void Resolutions()
{

    TString FileNames[3];
    FileNames[0] = "Res0.csv";
    FileNames[1] = "Res1.csv";
    FileNames[2] = "Res2.csv";

    TTree *tree[3];
    for (int i = 0; i < 3; i++)
    {
        tree[i] = new TTree("tree", "tree");
        tree[i] -> ReadFile(FileNames[i].Data(), "CopyNO/I:StdX/D:StdY/D", ',');
    }

    TGraph *grX[3];
    TGraph *grY[3];

    TMultiGraph *mgX = new TMultiGraph();
    TMultiGraph *mgY = new TMultiGraph();

    TLegend *legX = new TLegend(0.1, 0.1, 0.4, 0.3);

    for (int i = 0; i < 3; i++)
    {
        tree[i] -> Draw("StdX:CopyNO", "", "");
        grX[i] = new TGraph(tree[i]-> GetSelectedRows(), tree[i]->GetV2(), tree[i]->GetV1());
        grX[i] -> SetMarkerStyle(20);
        grX[i] -> SetMarkerColor(i+1);
        grX[i] -> SetLineColor(i+1);
        grX[i] -> SetLineWidth(2);
        grX[i] -> SetTitle(" ;Copy Number;Standard Deviation of X [deg]");
        grX[i] -> GetYaxis() -> SetRangeUser(4.5, 9.);
        mgX -> Add(grX[i]);


        tree[i] -> Draw("StdY:CopyNO", "", "");
        grY[i] = new TGraph(tree[i]-> GetSelectedRows(), tree[i]->GetV2(), tree[i]->GetV1());
        grY[i] -> SetMarkerStyle(20);
        grY[i] -> SetMarkerColor(i+1);
        grY[i] -> SetLineColor(i+1);
        grY[i] -> SetLineWidth(2);
        grY[i] -> SetTitle(" ;Copy Number;Standard Deviation of Y [deg]");
        grY[i] -> GetYaxis() -> SetRangeUser(4.5, 9.);
        mgY -> Add(grY[i]);

    }

    legX -> AddEntry(grX[0], "Electrons", "lp");
    legX -> AddEntry(grX[1], "Protons", "lp");
    legX -> AddEntry(grX[2], "Alpha", "lp");


    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
    c1 -> Divide(2,1);
    c1 -> cd(1);
    mgX -> SetTitle(" ;Copy Number;Standard Deviation of X [deg]");
    mgX -> Draw("APL");
    mgX -> GetYaxis() -> SetRangeUser(4.2, 8.6);
    legX -> Draw();
    gPad -> SetGrid();
    c1 -> cd(2);
    mgY -> SetTitle(" ;Copy Number;Standard Deviation of Y [deg]");
    mgY -> Draw("APL");
    mgY -> GetYaxis() -> SetRangeUser(4.2, 8.6);
    legX -> Draw();
    gPad -> SetGrid();



}