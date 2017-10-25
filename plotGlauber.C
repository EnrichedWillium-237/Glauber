# include "TAttAxis.h"
# include "TAttFill.h"
# include "TCanvas.h"
# include "TDirectory.h"
# include "TFile.h"
# include "TGraph.h"
# include "TGraphErrors.h"
# include "TH1.h"
# include "TH2.h"
# include "TH3.h"
# include "TLegend.h"
# include "TLine.h"
# include "TPaveText.h"
# include "TString.h"
# include "TStyle.h"
# include <fstream>
# include <iostream>

# include "../2015_PbPb/published_results/PhysRevC90_024902.h"
# include "../2015_PbPb/published_results/EurPhysJ_C73_2510.h"

using namespace std;

Bool_t print_plots = kTRUE;
Bool_t close_plots = kFALSE;
Bool_t minorAxes = kFALSE;

// static const int ncentbins = 14;
// static const double centbins[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80,  90, 100};

static const int ncentbins = 20;
static const double centbins[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};

static const int colors[] = {kTeal+4, kBlack, kMagenta, kRed, kAzure+10, kBlue};
static const int Gstyle[] = {29, 20, 31, 22, 34, 33};
static const float Gsize[] = {1.9, 1.2, 1.9, 1.5, 1.6, 1.9};

TH2D * xy1[ncentbins];
TH2D * xy2[ncentbins];
TH2D * xy3[ncentbins];
TH2D * xy4[ncentbins];
TH2D * xy5[ncentbins];
TH2D * xy6[ncentbins];

TH1D * hpsi[7][ncentbins];

TH2D * Psi2Psi1[ncentbins];
TH2D * Psi3Psi1[ncentbins];
TH1D * CosPsi2Psi1[ncentbins];
TH1D * CosPsi3Psi1[ncentbins];

TH2D * xyC1[ncentbins];
TH2D * xyC2[ncentbins];
TH2D * xyC3[ncentbins];
TH2D * xyC4[ncentbins];
TH2D * xyC5[ncentbins];
TH2D * xyC6[ncentbins];

TH1D * hpsiC[7][ncentbins];

TH2D * Psi2Psi1C[ncentbins];
TH2D * Psi3Psi1C[ncentbins];
TH1D * CosPsi2Psi1C[ncentbins];
TH1D * CosPsi3Psi1C[ncentbins];

TGraphErrors * g[7];
TGraphErrors * g_2_1;
TGraphErrors * g_3_1;
TGraphErrors * g_3_2;
TGraphErrors * g_4_2;
TGraphErrors * sg[7];
TGraphErrors * sg_2_1;
TGraphErrors * sg_3_1;
TGraphErrors * sg_3_2;
TGraphErrors * sg_4_2;

TGraphErrors * cg_4_p4_p2;
TGraphErrors * cg_8_p4_p2;
TGraphErrors * cg_12_p4_p2;
TGraphErrors * cg_6_p3_p2;
TGraphErrors * cg_6_p2_p6;
TGraphErrors * cg_6_p3_p6;
TGraphErrors * cg_12_p3_p4;
TGraphErrors * cg_10_p2_p5;

TGraphErrors * cg_2p2_3p3_5p5;
TGraphErrors * cg_2p2_4p4_6p6;
TGraphErrors * cg_2p2_6p3_4p4;
TGraphErrors * cg_8p2_3p3_5p5;
TGraphErrors * cg_10p2_4p4_6p6;
TGraphErrors * cg_10p2_6p3_4p4;

TGraphErrors * cg_6p2_p3;
TGraphErrors * cg_4p2_p4;
TGraphErrors * cg_6p2_p6;
TGraphErrors * cg_6p3_p6;
TGraphErrors * cg_10p2_p5;
TGraphErrors * cg_15p3_p5;

TGraphErrors * cg_2p1_p2;
TGraphErrors * cg_3p1_p3;
TGraphErrors * cg_4p1_p4;
TGraphErrors * cg_5p1_p5;
TGraphErrors * cg_6p1_p6;

TGraphErrors * cg_1p1_2p2_3p3;
TGraphErrors * cg_4p1_2p2_6p3;
TGraphErrors * cg_2p1_4p2_6p3;
TGraphErrors * cg_5p1_2p2_3p3;
TGraphErrors * cg_2p1_2p2_4p4;
TGraphErrors * cg_2p1_6p2_4p4;
TGraphErrors * cg_2p1_6p2_8p4;

TGraphErrors * cg_1p1_3p3_4p4;
TGraphErrors * cg_2p1_6p3_4p4;
TGraphErrors * cg_2p1_4p2_5p5;
TGraphErrors * cg_3p1_2p2_5p5;
TGraphErrors * cg_1p1_6p3_5p5;
TGraphErrors * cg_2p1_3p3_5p5;
TGraphErrors * cg_4p1_9p3_5p5;

TGraphErrors * gC[7];
TGraphErrors * gC_2_1;
TGraphErrors * gC_3_1;
TGraphErrors * gC_3_2;
TGraphErrors * gC_4_2;

TGraphErrors * cgC_4_p4_p2;
TGraphErrors * cgC_8_p4_p2;
TGraphErrors * cgC_12_p4_p2;
TGraphErrors * cgC_6_p3_p2;
TGraphErrors * cgC_6_p2_p6;
TGraphErrors * cgC_6_p3_p6;
TGraphErrors * cgC_12_p3_p4;
TGraphErrors * cgC_10_p2_p5;

TGraphErrors * cgC_2p2_3p3_5p5;
TGraphErrors * cgC_2p2_4p4_6p6;
TGraphErrors * cgC_2p2_6p3_4p4;
TGraphErrors * cgC_8p2_3p3_5p5;
TGraphErrors * cgC_10p2_4p4_6p6;
TGraphErrors * cgC_10p2_6p3_4p4;

// Theory points (0 = moment, 1 = cumulant)
TGraph * cgthy_4_p4_p2[2];
TGraph * cgthy_8_p4_p2[2];
TGraph * cgthy_12_p4_p2[2];
TGraph * cgthy_6_p3_p2[2];
TGraph * cgthy_6_p2_p6[2];
TGraph * cgthy_6_p3_p6[2];
TGraph * cgthy_12_p3_p4[2];
TGraph * cgthy_10_p2_p5[2];

TGraph * cgthy_2p2_3p3_5p5[2];
TGraph * cgthy_2p2_4p4_6p6[2];
TGraph * cgthy_2p2_6p3_4p4[2];
TGraph * cgthy_8p2_3p3_5p5[2];
TGraph * cgthy_10p2_4p4_6p6[2];
TGraph * cgthy_10p2_6p3_4p4[2];

TGraphErrors * cg_EPJ_Fig2_cos6Phi2_Phi3;
TGraphErrors * cg_EPJ_Fig2_cos4Phi2_Phi4;
TGraphErrors * cg_EPJ_Fig2_cos6Phi2_Phi6;
TGraphErrors * cg_EPJ_Fig2_cos6Phi3_Phi6;
TGraphErrors * cg_EPJ_Fig2_cos10Phi2_Phi5;
TGraphErrors * cg_EPJ_Fig2_cos15Phi3_Phi5;

TGraphErrors * cg_EPJ_Fig2_cos2Phi1_Phi2;
TGraphErrors * cg_EPJ_Fig2_cos3Phi1_Phi3;
TGraphErrors * cg_EPJ_Fig2_cos4Phi1_Phi4;
TGraphErrors * cg_EPJ_Fig2_cos5Phi1_Phi5;
TGraphErrors * cg_EPJ_Fig2_cos6Phi1_Phi6;

TGraphErrors * cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3;
TGraphErrors * cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3;
TGraphErrors * cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3;
TGraphErrors * cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3;

TGraphErrors * cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4;
TGraphErrors * cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4;
TGraphErrors * cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5;
TGraphErrors * cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5;

TGraphErrors * cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4;
TGraphErrors * cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4;
TGraphErrors * cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4;

TGraphErrors * cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5;
TGraphErrors * cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5;
TGraphErrors * cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5;


TFile * tfin;
TDirectory * tdmom;
TDirectory * tdcumu;

void plotGlauber()
{
    // if (minorAxes) tfin = new TFile("results/glauber_minorAxes.root");
    // else tfin = new TFile("results/glauber.root");
    // if (minorAxes) tfin = new TFile("results/results_AuAu/results_r2_weighting/glauber_minorAxes.root");
    // else tfin = new TFile("results/results_AuAu/results_r2_weighting/glauber.root");
    if (minorAxes) tfin = new TFile("results/results_AuAu/results_rn_weighting/glauber_minorAxes.root");
    else tfin = new TFile("results/results_AuAu/results_rn_weighting/glauber.root");

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        int lcmin = (int)centbins[cbin];
        int lcmax = (int)centbins[cbin+1];
        tdmom = (TDirectory *) tfin->Get(Form("moments/%d-%d",lcmin,lcmax));
        tdcumu = (TDirectory *) tfin->Get(Form("cumulants/%d-%d",lcmin,lcmax));

        xy1[cbin] = (TH2D *) tdmom->Get(Form("x_vs_y/xy1_%d-%d",lcmin,lcmax));
        xy2[cbin] = (TH2D *) tdmom->Get(Form("x_vs_y/xy2_%d-%d",lcmin,lcmax));
        xy3[cbin] = (TH2D *) tdmom->Get(Form("x_vs_y/xy3_%d-%d",lcmin,lcmax));
        xy4[cbin] = (TH2D *) tdmom->Get(Form("x_vs_y/xy4_%d-%d",lcmin,lcmax));
        xy5[cbin] = (TH2D *) tdmom->Get(Form("x_vs_y/xy5_%d-%d",lcmin,lcmax));
        xy6[cbin] = (TH2D *) tdmom->Get(Form("x_vs_y/xy6_%d-%d",lcmin,lcmax));

        xyC1[cbin] = (TH2D *) tdcumu->Get(Form("x_vs_y/xyC1_%d-%d",lcmin,lcmax));
        xyC2[cbin] = (TH2D *) tdcumu->Get(Form("x_vs_y/xyC2_%d-%d",lcmin,lcmax));
        xyC3[cbin] = (TH2D *) tdcumu->Get(Form("x_vs_y/xyC3_%d-%d",lcmin,lcmax));
        xyC4[cbin] = (TH2D *) tdcumu->Get(Form("x_vs_y/xyC4_%d-%d",lcmin,lcmax));
        xyC5[cbin] = (TH2D *) tdcumu->Get(Form("x_vs_y/xyC5_%d-%d",lcmin,lcmax));
        xyC6[cbin] = (TH2D *) tdcumu->Get(Form("x_vs_y/xyC6_%d-%d",lcmin,lcmax));

        xy1[cbin]->SetTitle("");
        xy2[cbin]->SetTitle("");
        xy3[cbin]->SetTitle("");
        xy4[cbin]->SetTitle("");
        xy5[cbin]->SetTitle("");
        xy6[cbin]->SetTitle("");

        xyC1[cbin]->SetTitle("");
        xyC2[cbin]->SetTitle("");
        xyC3[cbin]->SetTitle("");
        xyC4[cbin]->SetTitle("");
        xyC5[cbin]->SetTitle("");
        xyC6[cbin]->SetTitle("");

        xy1[cbin]->SetStats(0);
        xy2[cbin]->SetStats(0);
        xy3[cbin]->SetStats(0);
        xy4[cbin]->SetStats(0);
        xy5[cbin]->SetStats(0);
        xy6[cbin]->SetStats(0);

        xyC1[cbin]->SetStats(0);
        xyC2[cbin]->SetStats(0);
        xyC3[cbin]->SetStats(0);
        xyC4[cbin]->SetStats(0);
        xyC5[cbin]->SetStats(0);
        xyC6[cbin]->SetStats(0);

        xy1[cbin]->GetXaxis()->SetLabelSize(0.);
        xy2[cbin]->GetXaxis()->SetLabelSize(0.);
        xy3[cbin]->GetXaxis()->SetLabelSize(0.);
        xy4[cbin]->GetXaxis()->SetLabelSize(0.);
        xy5[cbin]->GetXaxis()->SetLabelSize(0.);
        xy6[cbin]->GetXaxis()->SetLabelSize(0.);

        xy1[cbin]->GetYaxis()->SetLabelSize(0.);
        xy2[cbin]->GetYaxis()->SetLabelSize(0.);
        xy3[cbin]->GetYaxis()->SetLabelSize(0.);
        xy4[cbin]->GetYaxis()->SetLabelSize(0.);
        xy5[cbin]->GetYaxis()->SetLabelSize(0.);
        xy6[cbin]->GetYaxis()->SetLabelSize(0.);

        xyC1[cbin]->GetXaxis()->SetLabelSize(0.);
        xyC2[cbin]->GetXaxis()->SetLabelSize(0.);
        xyC3[cbin]->GetXaxis()->SetLabelSize(0.);
        xyC4[cbin]->GetXaxis()->SetLabelSize(0.);
        xyC5[cbin]->GetXaxis()->SetLabelSize(0.);
        xyC6[cbin]->GetXaxis()->SetLabelSize(0.);

        xyC1[cbin]->GetYaxis()->SetLabelSize(0.);
        xyC2[cbin]->GetYaxis()->SetLabelSize(0.);
        xyC3[cbin]->GetYaxis()->SetLabelSize(0.);
        xyC4[cbin]->GetYaxis()->SetLabelSize(0.);
        xyC5[cbin]->GetYaxis()->SetLabelSize(0.);
        xyC6[cbin]->GetYaxis()->SetLabelSize(0.);

        xy1[cbin]->SetOption("CONT1");
        xy2[cbin]->SetOption("CONT1");
        xy3[cbin]->SetOption("CONT1");
        xy4[cbin]->SetOption("CONT1");
        xy5[cbin]->SetOption("CONT1");
        xy6[cbin]->SetOption("CONT1");

        xyC1[cbin]->SetOption("CONT1");
        xyC2[cbin]->SetOption("CONT1");
        xyC3[cbin]->SetOption("CONT1");
        xyC4[cbin]->SetOption("CONT1");
        xyC5[cbin]->SetOption("CONT1");
        xyC6[cbin]->SetOption("CONT1");

        for (int iorder = 1; iorder<=6; iorder++) {
            hpsi[iorder][cbin] = (TH1D *) tdmom->Get(Form("Psi/psi%d_%d-%d",iorder,lcmin,lcmax));
            hpsi[iorder][cbin]->SetTitle("");
            hpsi[iorder][cbin]->SetStats(0);
            hpsi[iorder][cbin]->SetXTitle(Form("#psi_{%d} Centrality %d-%d%%",iorder,lcmin,lcmax));
            hpsi[iorder][cbin]->GetXaxis()->SetLabelSize(0.06);
            hpsi[iorder][cbin]->GetYaxis()->SetLabelSize(0.06);
            hpsi[iorder][cbin]->GetXaxis()->SetTitleSize(0.06);
            hpsi[iorder][cbin]->GetXaxis()->SetTitleOffset(1.2);

            hpsiC[iorder][cbin] = (TH1D *) tdcumu->Get(Form("Psi/psiC%d_%d-%d",iorder,lcmin,lcmax));
            hpsiC[iorder][cbin]->SetTitle("");
            hpsiC[iorder][cbin]->SetStats(0);
            hpsiC[iorder][cbin]->SetXTitle(Form("#psi_{%d} Centrality %d-%d%%",iorder,lcmin,lcmax));
            hpsiC[iorder][cbin]->GetXaxis()->SetLabelSize(0.06);
            hpsiC[iorder][cbin]->GetYaxis()->SetLabelSize(0.06);
            hpsiC[iorder][cbin]->GetXaxis()->SetTitleSize(0.06);
            hpsiC[iorder][cbin]->GetXaxis()->SetTitleOffset(1.2);
        }

        Psi2Psi1[cbin] = (TH2D *) tdmom->Get(Form("Psi/Psi2Psi1_%d-%d",lcmin,lcmax));
        Psi3Psi1[cbin] = (TH2D *) tdmom->Get(Form("Psi/Psi3Psi1_%d-%d",lcmin,lcmax));
        CosPsi2Psi1[cbin] = (TH1D *) tdmom->Get(Form("Psi/CosPsi2Psi1_%d-%d",lcmin,lcmax));
        CosPsi3Psi1[cbin] = (TH1D *) tdmom->Get(Form("Psi/CosPsi3Psi1_%d-%d",lcmin,lcmax));

        Psi2Psi1C[cbin] = (TH2D *) tdcumu->Get(Form("Psi/Psi2Psi1C_%d-%d",lcmin,lcmax));
        Psi3Psi1C[cbin] = (TH2D *) tdcumu->Get(Form("Psi/Psi3Psi1C_%d-%d",lcmin,lcmax));
        CosPsi2Psi1C[cbin] = (TH1D *) tdcumu->Get(Form("Psi/CosPsi2Psi1C_%d-%d",lcmin,lcmax));
        CosPsi3Psi1C[cbin] = (TH1D *) tdcumu->Get(Form("Psi/CosPsi3Psi1C_%d-%d",lcmin,lcmax));
    }


    for (int iorder = 1; iorder<=6; iorder++) {
        g[iorder] = (TGraphErrors *) tfin->Get(Form("moments/n_%d/g%d",iorder,iorder));
        sg[iorder] = (TGraphErrors *) tfin->Get(Form("moments/n_%d/sg%d",iorder,iorder));
        gC[iorder] = (TGraphErrors *) tfin->Get(Form("cumulants/n_%d/g%d",iorder,iorder));
    }

    g_2_1  = (TGraphErrors *) tfin->Get("moments/g_2_1");
    g_3_1  = (TGraphErrors *) tfin->Get("moments/g_3_1");
    g_3_2  = (TGraphErrors *) tfin->Get("moments/g_3_2");
    g_4_2  = (TGraphErrors *) tfin->Get("moments/g_4_2");

    sg_2_1  = (TGraphErrors *) tfin->Get("moments/sg_2_1");
    sg_3_1  = (TGraphErrors *) tfin->Get("moments/sg_3_1");
    sg_3_2  = (TGraphErrors *) tfin->Get("moments/sg_3_2");
    sg_4_2  = (TGraphErrors *) tfin->Get("moments/sg_4_2");

    cg_4_p4_p2  = (TGraphErrors *) tfin->Get("moments/correlations/<cos4(psi4-psi2)>");
    cg_8_p4_p2  = (TGraphErrors *) tfin->Get("moments/correlations/<cos8(psi4-psi2)>");
    cg_12_p4_p2 = (TGraphErrors *) tfin->Get("moments/correlations/<cos12(psi4-psi2)>");
    cg_6_p3_p2  = (TGraphErrors *) tfin->Get("moments/correlations/<cos6(psi3-psi2)>");
    cg_6_p2_p6  = (TGraphErrors *) tfin->Get("moments/correlations/<cos6(psi2-psi6)>");
    cg_6_p3_p6  = (TGraphErrors *) tfin->Get("moments/correlations/<cos6(psi3-psi6)>");
    cg_12_p3_p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos12(psi3-psi4)>");
    cg_10_p2_p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos10(psi2-psi5)>");

    cg_2p2_3p3_5p5  = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi2+3psi3-5psi5)>");
    cg_2p2_4p4_6p6  = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi2+4psi4-6psi6)>");
    cg_2p2_6p3_4p4  = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi2-6psi3+4psi4)>");
    cg_8p2_3p3_5p5  = (TGraphErrors *) tfin->Get("moments/correlations/<cos(-8psi2+3psi3+5psi5)>");
    cg_10p2_4p4_6p6 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(-10psi2+4psi4+6psi6)>");
    cg_10p2_6p3_4p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(-10psi2+6psi3+4psi4)>");

    cgC_4_p4_p2  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos4(psi4-psi2)>");
    cgC_8_p4_p2  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos8(psi4-psi2)>");
    cgC_12_p4_p2 = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos12(psi4-psi2)>");
    cgC_6_p3_p2  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos6(psi3-psi2)>");
    cgC_6_p2_p6  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos6(psi2-psi6)>");
    cgC_6_p3_p6  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos6(psi3-psi6)>");
    cgC_12_p3_p4 = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos12(psi3-psi4)>");
    cgC_10_p2_p5 = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos10(psi2-psi5)>");

    cgC_2p2_3p3_5p5  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos(2psi2+3psi3-5psi5)>");
    cgC_2p2_4p4_6p6  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos(2psi2+4psi4-6psi6)>");
    cgC_2p2_6p3_4p4  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos(2psi2-6psi3+4psi4)>");
    cgC_8p2_3p3_5p5  = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos(-8psi2+3psi3+5psi5)>");
    cgC_10p2_4p4_6p6 = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos(-10psi2+4psi4+6psi6)>");
    cgC_10p2_6p3_4p4 = (TGraphErrors *) tfin->Get("cumulants/correlations/<cos(-10psi2+6psi3+4psi4)>");

    cg_6p2_p3 = (TGraphErrors *) tfin->Get("moments/correlations/<cos6(psi2-psi3)>");
    cg_4p2_p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos4(psi2-psi4)>");
    cg_6p2_p6 = (TGraphErrors *) tfin->Get("moments/correlations/<cos6(psi2-psi6)>");
    cg_6p3_p6 = (TGraphErrors *) tfin->Get("moments/correlations/<cos6(psi3-psi6)>");
    cg_10p2_p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos10(psi2-psi5)>");
    cg_15p3_p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos15(psi3-psi5)>");

    cg_2p1_p2 = (TGraphErrors *) tfin->Get("moments/correlations/<cos2(psi1-psi2)>");
    cg_3p1_p3 = (TGraphErrors *) tfin->Get("moments/correlations/<cos3(psi1-psi3)>");
    cg_4p1_p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos4(psi1-psi4)>");
    cg_5p1_p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos5(psi1-psi5)>");
    cg_6p1_p6 = (TGraphErrors *) tfin->Get("moments/correlations/<cos6(psi1-psi6)>");

    cg_1p1_2p2_3p3 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(1psi1+2psi2-3psi3)>");
    cg_4p1_2p2_6p3 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(4psi1+2psi2-6psi3)>");
    cg_2p1_4p2_6p3 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi1+4psi2-6psi3)>");
    cg_5p1_2p2_3p3 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(5psi1-2psi2-3psi3)>");
    cg_2p1_2p2_4p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi1+2psi2-4psi4)>");
    cg_2p1_6p2_4p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi1-6psi2+4psi4)>");
    cg_2p1_6p2_8p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi1+6psi2-8psi4)>");

    cg_1p1_3p3_4p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(1psi1+3psi3-4psi4)>");
    cg_2p1_6p3_4p4 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi1-6psi3+4psi4)>");
    cg_2p1_4p2_5p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(1psi1+4psi2-5psi5)>");
    cg_3p1_2p2_5p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(3psi1+2psi2-5psi5)>");
    cg_1p1_6p3_5p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(1psi1-6psi3+5psi5)>");
    cg_2p1_3p3_5p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(2psi1+3psi3-5psi5)>");
    cg_4p1_9p3_5p5 = (TGraphErrors *) tfin->Get("moments/correlations/<cos(4psi1-9psi3+5psi5)>");

    for (int iorder = 1; iorder<=6; iorder++) {
        g[iorder]->SetMarkerColor(colors[iorder-1]);
        g[iorder]->SetLineColor(colors[iorder-1]);
        g[iorder]->SetMarkerStyle(Gstyle[iorder-1]);
        g[iorder]->SetMarkerSize(Gsize[iorder-1]);

        sg[iorder]->SetMarkerColor(colors[iorder-1]);
        sg[iorder]->SetLineColor(colors[iorder-1]);
        sg[iorder]->SetMarkerStyle(Gstyle[iorder-1]);
        sg[iorder]->SetMarkerSize(Gsize[iorder-1]);

        gC[iorder]->SetMarkerColor(colors[iorder-1]);
        gC[iorder]->SetLineColor(colors[iorder-1]);
        gC[iorder]->SetMarkerStyle(Gstyle[iorder-1]);
        gC[iorder]->SetMarkerSize(Gsize[iorder-1]);
    }

    cg_4_p4_p2->SetMarkerColor(kBlue);
    cg_4_p4_p2->SetLineColor(kBlue);
    cg_4_p4_p2->SetMarkerStyle(20);
    cg_4_p4_p2->SetMarkerSize(1.1);
    cg_8_p4_p2->SetMarkerColor(kBlue);
    cg_8_p4_p2->SetLineColor(kBlue);
    cg_8_p4_p2->SetMarkerStyle(20);
    cg_8_p4_p2->SetMarkerSize(1.1);
    cg_12_p4_p2->SetMarkerColor(kBlue);
    cg_12_p4_p2->SetLineColor(kBlue);
    cg_12_p4_p2->SetMarkerStyle(20);
    cg_12_p4_p2->SetMarkerSize(1.1);
    cg_6_p3_p2->SetMarkerColor(kBlue);
    cg_6_p3_p2->SetLineColor(kBlue);
    cg_6_p3_p2->SetMarkerStyle(20);
    cg_6_p3_p2->SetMarkerSize(1.1);
    cg_6_p2_p6->SetMarkerColor(kBlue);
    cg_6_p2_p6->SetLineColor(kBlue);
    cg_6_p2_p6->SetMarkerStyle(20);
    cg_6_p2_p6->SetMarkerSize(1.1);
    cg_6_p3_p6->SetMarkerColor(kBlue);
    cg_6_p3_p6->SetLineColor(kBlue);
    cg_6_p3_p6->SetMarkerStyle(20);
    cg_6_p3_p6->SetMarkerSize(1.1);
    cg_12_p3_p4->SetMarkerColor(kBlue);
    cg_12_p3_p4->SetLineColor(kBlue);
    cg_12_p3_p4->SetMarkerStyle(20);
    cg_12_p3_p4->SetMarkerSize(1.1);
    cg_10_p2_p5->SetMarkerColor(kBlue);
    cg_10_p2_p5->SetLineColor(kBlue);
    cg_10_p2_p5->SetMarkerStyle(20);
    cg_10_p2_p5->SetMarkerSize(1.1);

    cg_2p2_3p3_5p5->SetMarkerColor(kBlue);
    cg_2p2_3p3_5p5->SetLineColor(kBlue);
    cg_2p2_3p3_5p5->SetMarkerStyle(20);
    cg_2p2_3p3_5p5->SetMarkerSize(1.1);
    cg_2p2_4p4_6p6->SetMarkerColor(kBlue);
    cg_2p2_4p4_6p6->SetLineColor(kBlue);
    cg_2p2_4p4_6p6->SetMarkerStyle(20);
    cg_2p2_4p4_6p6->SetMarkerSize(1.1);
    cg_2p2_6p3_4p4->SetMarkerColor(kBlue);
    cg_2p2_6p3_4p4->SetLineColor(kBlue);
    cg_2p2_6p3_4p4->SetMarkerStyle(20);
    cg_2p2_6p3_4p4->SetMarkerSize(1.1);
    cg_8p2_3p3_5p5->SetMarkerColor(kBlue);
    cg_8p2_3p3_5p5->SetLineColor(kBlue);
    cg_8p2_3p3_5p5->SetMarkerStyle(20);
    cg_8p2_3p3_5p5->SetMarkerSize(1.1);
    cg_10p2_4p4_6p6->SetMarkerColor(kBlue);
    cg_10p2_4p4_6p6->SetLineColor(kBlue);
    cg_10p2_4p4_6p6->SetMarkerStyle(20);
    cg_10p2_4p4_6p6->SetMarkerSize(1.1);
    cg_10p2_6p3_4p4->SetMarkerColor(kBlue);
    cg_10p2_6p3_4p4->SetLineColor(kBlue);
    cg_10p2_6p3_4p4->SetMarkerStyle(20);
    cg_10p2_6p3_4p4->SetMarkerSize(1.1);

    cg_6p2_p3->SetMarkerColor(kBlack);
    cg_6p2_p3->SetLineColor(kBlack);
    cg_6p2_p3->SetMarkerStyle(20);
    cg_6p2_p3->SetMarkerSize(1.0);
    cg_4p2_p4->SetMarkerColor(kRed);
    cg_4p2_p4->SetLineColor(kRed);
    cg_4p2_p4->SetMarkerStyle(24);
    cg_4p2_p4->SetMarkerSize(1.0);
    cg_6p2_p6->SetMarkerColor(kCyan+2);
    cg_6p2_p6->SetLineColor(kCyan+2);
    cg_6p2_p6->SetMarkerStyle(20);
    cg_6p2_p6->SetMarkerSize(1.0);
    cg_6p3_p6->SetMarkerColor(kBlue);
    cg_6p3_p6->SetLineColor(kBlue);
    cg_6p3_p6->SetMarkerStyle(24);
    cg_6p3_p6->SetMarkerSize(1.0);
    cg_10p2_p5->SetMarkerColor(kMagenta);
    cg_10p2_p5->SetLineColor(kMagenta);
    cg_10p2_p5->SetMarkerStyle(20);
    cg_10p2_p5->SetMarkerSize(1.0);
    cg_15p3_p5->SetMarkerColor(kTeal+3);
    cg_15p3_p5->SetLineColor(kTeal+3);
    cg_15p3_p5->SetMarkerStyle(24);
    cg_15p3_p5->SetMarkerSize(1.0);

    cg_2p1_p2->SetMarkerColor(kBlack);
    cg_2p1_p2->SetLineColor(kBlack);
    cg_2p1_p2->SetMarkerStyle(20);
    cg_2p1_p2->SetMarkerSize(1.0);
    cg_3p1_p3->SetMarkerColor(kRed);
    cg_3p1_p3->SetLineColor(kRed);
    cg_3p1_p3->SetMarkerStyle(24);
    cg_3p1_p3->SetMarkerSize(1.0);
    cg_4p1_p4->SetMarkerColor(kCyan+2);
    cg_4p1_p4->SetLineColor(kCyan+2);
    cg_4p1_p4->SetMarkerStyle(20);
    cg_4p1_p4->SetMarkerSize(1.0);
    cg_5p1_p5->SetMarkerColor(kBlue);
    cg_5p1_p5->SetLineColor(kBlue);
    cg_5p1_p5->SetMarkerStyle(24);
    cg_5p1_p5->SetMarkerSize(1.0);
    cg_6p1_p6->SetMarkerColor(kMagenta);
    cg_6p1_p6->SetLineColor(kMagenta);
    cg_6p1_p6->SetMarkerStyle(20);
    cg_6p1_p6->SetMarkerSize(1.0);

    cg_1p1_2p2_3p3->SetMarkerColor(kBlack);
    cg_1p1_2p2_3p3->SetLineColor(kBlack);
    cg_1p1_2p2_3p3->SetMarkerStyle(24);
    cg_1p1_2p2_3p3->SetMarkerSize(1.1);
    cg_4p1_2p2_6p3->SetMarkerColor(kRed);
    cg_4p1_2p2_6p3->SetLineColor(kRed);
    cg_4p1_2p2_6p3->SetMarkerStyle(30);
    cg_4p1_2p2_6p3->SetMarkerSize(1.7);
    cg_2p1_4p2_6p3->SetMarkerColor(kBlue);
    cg_2p1_4p2_6p3->SetLineColor(kBlue);
    cg_2p1_4p2_6p3->SetMarkerStyle(28);
    cg_2p1_4p2_6p3->SetMarkerSize(1.6);
    cg_5p1_2p2_3p3->SetMarkerColor(kMagenta);
    cg_5p1_2p2_3p3->SetLineColor(kMagenta);
    cg_5p1_2p2_3p3->SetMarkerStyle(25);
    cg_5p1_2p2_3p3->SetMarkerSize(1.1);
    cg_2p1_2p2_4p4->SetMarkerColor(kBlack);
    cg_2p1_2p2_4p4->SetLineColor(kBlack);
    cg_2p1_2p2_4p4->SetMarkerStyle(24);
    cg_2p1_2p2_4p4->SetMarkerSize(1.1);
    cg_2p1_6p2_4p4->SetMarkerColor(kRed);
    cg_2p1_6p2_4p4->SetLineColor(kRed);
    cg_2p1_6p2_4p4->SetMarkerStyle(30);
    cg_2p1_6p2_4p4->SetMarkerSize(1.7);
    cg_2p1_6p2_8p4->SetMarkerColor(kBlue);
    cg_2p1_6p2_8p4->SetLineColor(kBlue);
    cg_2p1_6p2_8p4->SetMarkerStyle(28);
    cg_2p1_6p2_8p4->SetMarkerSize(1.6);

    cg_1p1_3p3_4p4->SetMarkerColor(kBlack);
    cg_1p1_3p3_4p4->SetLineColor(kBlack);
    cg_1p1_3p3_4p4->SetMarkerStyle(24);
    cg_1p1_3p3_4p4->SetMarkerSize(1.1);
    cg_2p1_6p3_4p4->SetMarkerColor(kRed);
    cg_2p1_6p3_4p4->SetLineColor(kRed);
    cg_2p1_6p3_4p4->SetMarkerStyle(30);
    cg_2p1_6p3_4p4->SetMarkerSize(1.7);
    cg_2p1_4p2_5p5->SetMarkerColor(kBlue);
    cg_2p1_4p2_5p5->SetLineColor(kBlue);
    cg_2p1_4p2_5p5->SetMarkerStyle(28);
    cg_2p1_4p2_5p5->SetMarkerSize(1.6);
    cg_3p1_2p2_5p5->SetMarkerColor(kMagenta);
    cg_3p1_2p2_5p5->SetLineColor(kMagenta);
    cg_3p1_2p2_5p5->SetMarkerStyle(25);
    cg_3p1_2p2_5p5->SetMarkerSize(1.1);
    cg_1p1_6p3_5p5->SetMarkerColor(kBlack);
    cg_1p1_6p3_5p5->SetLineColor(kBlack);
    cg_1p1_6p3_5p5->SetMarkerStyle(24);
    cg_1p1_6p3_5p5->SetMarkerSize(1.1);
    cg_2p1_3p3_5p5->SetMarkerColor(kRed);
    cg_2p1_3p3_5p5->SetLineColor(kRed);
    cg_2p1_3p3_5p5->SetMarkerStyle(30);
    cg_2p1_3p3_5p5->SetMarkerSize(1.7);
    cg_4p1_9p3_5p5->SetMarkerColor(kBlue);
    cg_4p1_9p3_5p5->SetLineColor(kBlue);
    cg_4p1_9p3_5p5->SetMarkerStyle(28);
    cg_4p1_9p3_5p5->SetMarkerSize(1.6);

    cgC_4_p4_p2->SetMarkerColor(kRed);
    cgC_4_p4_p2->SetLineColor(kRed);
    cgC_4_p4_p2->SetMarkerStyle(20);
    cgC_4_p4_p2->SetMarkerSize(1.1);
    cgC_8_p4_p2->SetMarkerColor(kRed);
    cgC_8_p4_p2->SetLineColor(kRed);
    cgC_8_p4_p2->SetMarkerStyle(20);
    cgC_8_p4_p2->SetMarkerSize(1.1);
    cgC_12_p4_p2->SetMarkerColor(kRed);
    cgC_12_p4_p2->SetLineColor(kRed);
    cgC_12_p4_p2->SetMarkerStyle(20);
    cgC_12_p4_p2->SetMarkerSize(1.1);
    cgC_6_p3_p2->SetMarkerColor(kRed);
    cgC_6_p3_p2->SetLineColor(kRed);
    cgC_6_p3_p2->SetMarkerStyle(20);
    cgC_6_p3_p2->SetMarkerSize(1.1);
    cgC_6_p2_p6->SetMarkerColor(kRed);
    cgC_6_p2_p6->SetLineColor(kRed);
    cgC_6_p2_p6->SetMarkerStyle(20);
    cgC_6_p2_p6->SetMarkerSize(1.1);
    cgC_6_p3_p6->SetMarkerColor(kRed);
    cgC_6_p3_p6->SetLineColor(kRed);
    cgC_6_p3_p6->SetMarkerStyle(20);
    cgC_6_p3_p6->SetMarkerSize(1.1);
    cgC_12_p3_p4->SetMarkerColor(kRed);
    cgC_12_p3_p4->SetLineColor(kRed);
    cgC_12_p3_p4->SetMarkerStyle(20);
    cgC_12_p3_p4->SetMarkerSize(1.1);
    cgC_10_p2_p5->SetMarkerColor(kRed);
    cgC_10_p2_p5->SetLineColor(kRed);
    cgC_10_p2_p5->SetMarkerStyle(20);
    cgC_10_p2_p5->SetMarkerSize(1.1);

    cgC_2p2_3p3_5p5->SetMarkerColor(kRed);
    cgC_2p2_3p3_5p5->SetLineColor(kRed);
    cgC_2p2_3p3_5p5->SetMarkerStyle(20);
    cgC_2p2_3p3_5p5->SetMarkerSize(1.1);
    cgC_2p2_4p4_6p6->SetMarkerColor(kRed);
    cgC_2p2_4p4_6p6->SetLineColor(kRed);
    cgC_2p2_4p4_6p6->SetMarkerStyle(20);
    cgC_2p2_4p4_6p6->SetMarkerSize(1.1);
    cgC_2p2_6p3_4p4->SetMarkerColor(kRed);
    cgC_2p2_6p3_4p4->SetLineColor(kRed);
    cgC_2p2_6p3_4p4->SetMarkerStyle(20);
    cgC_2p2_6p3_4p4->SetMarkerSize(1.1);
    cgC_8p2_3p3_5p5->SetMarkerColor(kRed);
    cgC_8p2_3p3_5p5->SetLineColor(kRed);
    cgC_8p2_3p3_5p5->SetMarkerStyle(20);
    cgC_8p2_3p3_5p5->SetMarkerSize(1.1);
    cgC_10p2_4p4_6p6->SetMarkerColor(kRed);
    cgC_10p2_4p4_6p6->SetLineColor(kRed);
    cgC_10p2_4p4_6p6->SetMarkerStyle(20);
    cgC_10p2_4p4_6p6->SetMarkerSize(1.1);
    cgC_10p2_6p3_4p4->SetMarkerColor(kRed);
    cgC_10p2_6p3_4p4->SetLineColor(kRed);
    cgC_10p2_6p3_4p4->SetMarkerStyle(20);
    cgC_10p2_6p3_4p4->SetMarkerSize(1.1);

    cgthy_4_p4_p2[0] = new TGraph(13, cdat_4_p4_p2_moment_x, cdat_4_p4_p2_moment_y);
    cgthy_8_p4_p2[0] = new TGraph(13, cdat_8_p4_p2_moment_x, cdat_8_p4_p2_moment_y);
    cgthy_12_p4_p2[0] = new TGraph(13, cdat_12_p4_p2_moment_x, cdat_12_p4_p2_moment_y);
    cgthy_6_p3_p2[0] = new TGraph(13, cdat_6_p3_p2_moment_x, cdat_6_p3_p2_moment_y);
    cgthy_6_p2_p6[0] = new TGraph(13, cdat_6_p2_p6_moment_x, cdat_6_p2_p6_moment_y);
    cgthy_6_p3_p6[0] = new TGraph(13, cdat_6_p3_p6_moment_x, cdat_6_p3_p6_moment_y);
    cgthy_12_p3_p4[0] = new TGraph(13, cdat_12_p3_p4_moment_x, cdat_12_p3_p4_moment_y);
    cgthy_10_p2_p5[0] = new TGraph(13, cdat_10_p2_p5_moment_x, cdat_10_p2_p5_moment_y);

    cgthy_2p2_3p3_5p5[0] = new TGraph(13, cdat_2p2_3p3_5p5_moment_x, cdat_2p2_3p3_5p5_moment_y);
    cgthy_2p2_4p4_6p6[0] = new TGraph(13, cdat_2p2_4p4_6p6_moment_x, cdat_2p2_4p4_6p6_moment_y);
    cgthy_2p2_6p3_4p4[0] = new TGraph(13, cdat_2p2_6p3_4p4_moment_x, cdat_2p2_6p3_4p4_moment_y);
    cgthy_8p2_3p3_5p5[0] = new TGraph(13, cdat_8p2_3p3_5p5_moment_x, cdat_8p2_3p3_5p5_moment_y);
    cgthy_10p2_4p4_6p6[0] = new TGraph(13, cdat_10p2_4p4_6p6_moment_x, cdat_10p2_4p4_6p6_moment_y);
    cgthy_10p2_6p3_4p4[0] = new TGraph(13, cdat_10p2_6p3_4p4_moment_x, cdat_10p2_6p3_4p4_moment_y);

    cgthy_4_p4_p2[1] = new TGraph(13, cdat_4_p4_p2_cumu_x, cdat_4_p4_p2_cumu_y);
    cgthy_8_p4_p2[1] = new TGraph(13, cdat_8_p4_p2_cumu_x, cdat_8_p4_p2_cumu_y);
    cgthy_12_p4_p2[1] = new TGraph(13, cdat_12_p4_p2_cumu_x, cdat_12_p4_p2_cumu_y);
    cgthy_6_p2_p6[1] = new TGraph(13, cdat_6_p2_p6_cumu_x, cdat_6_p2_p6_cumu_y);
    cgthy_6_p3_p6[1] = new TGraph(13, cdat_6_p3_p6_cumu_x, cdat_6_p3_p6_cumu_y);
    cgthy_12_p3_p4[1] = new TGraph(13, cdat_12_p3_p4_cumu_x, cdat_12_p3_p4_cumu_y);
    cgthy_10_p2_p5[1] = new TGraph(13, cdat_10_p2_p5_cumu_x, cdat_10_p2_p5_cumu_y);

    cgthy_2p2_3p3_5p5[1] = new TGraph(13, cdat_2p2_3p3_5p5_cumu_x, cdat_2p2_3p3_5p5_cumu_y);
    cgthy_2p2_4p4_6p6[1] = new TGraph(13, cdat_2p2_4p4_6p6_cumu_x, cdat_2p2_4p4_6p6_cumu_y);
    cgthy_2p2_6p3_4p4[1] = new TGraph(13, cdat_2p2_6p3_4p4_cumu_x, cdat_2p2_6p3_4p4_cumu_y);
    cgthy_8p2_3p3_5p5[1] = new TGraph(13, cdat_8p2_3p3_5p5_cumu_x, cdat_8p2_3p3_5p5_cumu_y);
    cgthy_10p2_4p4_6p6[1] = new TGraph(13, cdat_10p2_4p4_6p6_cumu_x, cdat_10p2_4p4_6p6_cumu_y);
    cgthy_10p2_6p3_4p4[1] = new TGraph(13, cdat_10p2_6p3_4p4_cumu_x, cdat_10p2_6p3_4p4_cumu_y);


    cg_EPJ_Fig2_cos6Phi2_Phi3 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos6Phi2_Phi3, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos4Phi2_Phi4 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos4Phi2_Phi4, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos6Phi2_Phi6 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos6Phi2_Phi6, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos6Phi3_Phi6 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos6Phi3_Phi6, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos10Phi2_Phi5 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos10Phi2_Phi5, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos15Phi3_Phi5 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos15Phi3_Phi5, 0, EPJerr_Fig2);

    cg_EPJ_Fig2_cos2Phi1_Phi2 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos2Phi1_Phi2, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos3Phi1_Phi3 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos3Phi1_Phi3, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos4Phi1_Phi4 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos4Phi1_Phi4, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos5Phi1_Phi5 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos5Phi1_Phi5, 0, EPJerr_Fig2);
    cg_EPJ_Fig2_cos6Phi1_Phi6 = new TGraphErrors(EPJ_Fig2_NumNpart, EPJ_Fig2_Npart, EPJ_Fig2_cos6Phi1_Phi6, 0, EPJerr_Fig2);

    cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos1Phi1_2Phi2_3Phi3, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos4Phi1_2Phi2_6Phi3, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos2Phi1_4Phi2_6Phi3, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos5Phi1_2Phi2_3Phi3, 0, EPJerr_Fig12);

    cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos1Phi1_3Phi3_4Phi4, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos2Phi1_6Phi3_4Phi4, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos1Phi1_4Phi2_5Phi5, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos3Phi1_2Phi2_5Phi5, 0, EPJerr_Fig12);

    cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos2Phi1_2Phi2_4Phi4, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos2Phi1_6Phi2_4Phi4, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos2Phi1_6Phi2_8Phi4, 0, EPJerr_Fig12);

    cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos1Phi1_6Phi3_5Phi5, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos2Phi1_3Phi3_5Phi5, 0, EPJerr_Fig12);
    cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5 = new TGraphErrors(EPJ_Fig12_NumNpart, EPJ_Fig12_Npart, EPJ_Fig12_cos4Phi1_9Phi3_5Phi5, 0, EPJerr_Fig12);

    for (int i = 0; i<2; i++) {
        if (i == 0) cgthy_4_p4_p2[i]->SetMarkerStyle(20);
        else cgthy_4_p4_p2[i]->SetMarkerStyle(24);
        cgthy_4_p4_p2[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_8_p4_p2[i]->SetMarkerStyle(20);
        else cgthy_8_p4_p2[i]->SetMarkerStyle(24);
        cgthy_8_p4_p2[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_12_p4_p2[i]->SetMarkerStyle(20);
        else cgthy_12_p4_p2[i]->SetMarkerStyle(24);
        cgthy_12_p4_p2[i]->SetMarkerSize(1.1);

        cgthy_6_p3_p2[0]->SetMarkerStyle(20);
        cgthy_6_p3_p2[0]->SetMarkerSize(1.1);

        if (i == 0) cgthy_6_p2_p6[i]->SetMarkerStyle(20);
        else cgthy_6_p2_p6[i]->SetMarkerStyle(24);
        cgthy_6_p2_p6[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_6_p3_p6[i]->SetMarkerStyle(20);
        else cgthy_6_p3_p6[i]->SetMarkerStyle(24);
        cgthy_6_p3_p6[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_12_p3_p4[i]->SetMarkerStyle(20);
        else cgthy_12_p3_p4[i]->SetMarkerStyle(24);
        cgthy_12_p3_p4[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_10_p2_p5[i]->SetMarkerStyle(20);
        else cgthy_10_p2_p5[i]->SetMarkerStyle(24);
        cgthy_10_p2_p5[i]->SetMarkerSize(1.1);


        if (i == 0) cgthy_2p2_3p3_5p5[i]->SetMarkerStyle(20);
        else cgthy_2p2_3p3_5p5[i]->SetMarkerStyle(24);
        cgthy_2p2_3p3_5p5[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_2p2_4p4_6p6[i]->SetMarkerStyle(20);
        else cgthy_2p2_4p4_6p6[i]->SetMarkerStyle(24);
        cgthy_2p2_4p4_6p6[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_2p2_6p3_4p4[i]->SetMarkerStyle(20);
        else cgthy_2p2_6p3_4p4[i]->SetMarkerStyle(24);
        cgthy_2p2_6p3_4p4[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_8p2_3p3_5p5[i]->SetMarkerStyle(20);
        else cgthy_8p2_3p3_5p5[i]->SetMarkerStyle(24);
        cgthy_8p2_3p3_5p5[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_10p2_4p4_6p6[i]->SetMarkerStyle(20);
        else cgthy_10p2_4p4_6p6[i]->SetMarkerStyle(24);
        cgthy_10p2_4p4_6p6[i]->SetMarkerSize(1.1);

        if (i == 0) cgthy_10p2_6p3_4p4[i]->SetMarkerStyle(20);
        else cgthy_10p2_6p3_4p4[i]->SetMarkerStyle(24);
        cgthy_10p2_6p3_4p4[i]->SetMarkerSize(1.1);
    }

    float transfill = 0.85;

    cg_EPJ_Fig2_cos6Phi2_Phi3->SetLineColor(kBlack);
    cg_EPJ_Fig2_cos6Phi2_Phi3->SetFillColorAlpha(kBlack, transfill);
    cg_EPJ_Fig2_cos6Phi2_Phi3->SetFillStyle(3002);
    cg_EPJ_Fig2_cos6Phi2_Phi3->SetLineWidth(2);
    cg_EPJ_Fig2_cos6Phi2_Phi3->SetLineStyle(2);
    cg_EPJ_Fig2_cos4Phi2_Phi4->SetLineColor(kRed);
    cg_EPJ_Fig2_cos4Phi2_Phi4->SetFillColorAlpha(kRed, transfill);
    cg_EPJ_Fig2_cos4Phi2_Phi4->SetFillStyle(3002);
    cg_EPJ_Fig2_cos4Phi2_Phi4->SetLineWidth(2);
    cg_EPJ_Fig2_cos4Phi2_Phi4->SetLineStyle(2);
    cg_EPJ_Fig2_cos6Phi2_Phi6->SetLineColor(kCyan+2);
    cg_EPJ_Fig2_cos6Phi2_Phi6->SetFillColorAlpha(kCyan+2, transfill);
    cg_EPJ_Fig2_cos6Phi2_Phi6->SetFillStyle(3002);
    cg_EPJ_Fig2_cos6Phi2_Phi6->SetLineWidth(2);
    cg_EPJ_Fig2_cos6Phi2_Phi6->SetLineStyle(2);
    cg_EPJ_Fig2_cos6Phi3_Phi6->SetLineColor(kBlue);
    cg_EPJ_Fig2_cos6Phi3_Phi6->SetFillColorAlpha(kBlue, transfill);
    cg_EPJ_Fig2_cos6Phi3_Phi6->SetFillStyle(3002);
    cg_EPJ_Fig2_cos6Phi3_Phi6->SetLineWidth(2);
    cg_EPJ_Fig2_cos6Phi3_Phi6->SetLineStyle(2);
    cg_EPJ_Fig2_cos10Phi2_Phi5->SetLineColor(kMagenta);
    cg_EPJ_Fig2_cos10Phi2_Phi5->SetFillColorAlpha(kMagenta, transfill);
    cg_EPJ_Fig2_cos10Phi2_Phi5->SetFillStyle(3002);
    cg_EPJ_Fig2_cos10Phi2_Phi5->SetLineWidth(2);
    cg_EPJ_Fig2_cos10Phi2_Phi5->SetLineStyle(2);
    cg_EPJ_Fig2_cos15Phi3_Phi5->SetLineColor(kTeal+3);
    cg_EPJ_Fig2_cos15Phi3_Phi5->SetFillColorAlpha(kTeal+3, transfill);
    cg_EPJ_Fig2_cos15Phi3_Phi5->SetFillStyle(3002);
    cg_EPJ_Fig2_cos15Phi3_Phi5->SetLineWidth(2);
    cg_EPJ_Fig2_cos15Phi3_Phi5->SetLineStyle(2);

    cg_EPJ_Fig2_cos2Phi1_Phi2->SetLineColor(kBlack);
    cg_EPJ_Fig2_cos2Phi1_Phi2->SetFillColorAlpha(kBlack, transfill);
    cg_EPJ_Fig2_cos2Phi1_Phi2->SetFillStyle(3002);
    cg_EPJ_Fig2_cos2Phi1_Phi2->SetLineWidth(2);
    cg_EPJ_Fig2_cos2Phi1_Phi2->SetLineStyle(2);
    cg_EPJ_Fig2_cos3Phi1_Phi3->SetLineColor(kRed);
    cg_EPJ_Fig2_cos3Phi1_Phi3->SetFillColorAlpha(kRed, transfill);
    cg_EPJ_Fig2_cos3Phi1_Phi3->SetFillStyle(3002);
    cg_EPJ_Fig2_cos3Phi1_Phi3->SetLineWidth(2);
    cg_EPJ_Fig2_cos3Phi1_Phi3->SetLineStyle(2);
    cg_EPJ_Fig2_cos4Phi1_Phi4->SetLineColor(kCyan+2);
    cg_EPJ_Fig2_cos4Phi1_Phi4->SetFillColorAlpha(kCyan+2, transfill);
    cg_EPJ_Fig2_cos4Phi1_Phi4->SetFillStyle(3002);
    cg_EPJ_Fig2_cos4Phi1_Phi4->SetLineWidth(2);
    cg_EPJ_Fig2_cos4Phi1_Phi4->SetLineStyle(2);
    cg_EPJ_Fig2_cos5Phi1_Phi5->SetLineColor(kBlue);
    cg_EPJ_Fig2_cos5Phi1_Phi5->SetFillColorAlpha(kBlue, transfill);
    cg_EPJ_Fig2_cos5Phi1_Phi5->SetFillStyle(3002);
    cg_EPJ_Fig2_cos5Phi1_Phi5->SetLineWidth(2);
    cg_EPJ_Fig2_cos5Phi1_Phi5->SetLineStyle(2);
    cg_EPJ_Fig2_cos6Phi1_Phi6->SetLineColor(kMagenta);
    cg_EPJ_Fig2_cos6Phi1_Phi6->SetFillColorAlpha(kMagenta, transfill);
    cg_EPJ_Fig2_cos6Phi1_Phi6->SetFillStyle(3002);
    cg_EPJ_Fig2_cos6Phi1_Phi6->SetLineWidth(2);
    cg_EPJ_Fig2_cos6Phi1_Phi6->SetLineStyle(2);

    cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3->SetLineColor(kBlack);
    cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3->SetFillColorAlpha(kBlack, transfill);
    cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3->SetFillStyle(3002);
    cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3->SetLineWidth(2);
    cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3->SetLineStyle(2);
    cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3->SetLineColor(kRed);
    cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3->SetFillColorAlpha(kRed, transfill);
    cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3->SetFillStyle(3002);
    cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3->SetLineWidth(2);
    cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3->SetLineStyle(2);
    cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3->SetLineColor(kBlue);
    cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3->SetFillColorAlpha(kBlue, transfill);
    cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3->SetFillStyle(3002);
    cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3->SetLineWidth(2);
    cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3->SetLineStyle(2);
    cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3->SetLineColor(kMagenta);
    cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3->SetFillColorAlpha(kMagenta, transfill);
    cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3->SetFillStyle(3002);
    cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3->SetLineWidth(2);
    cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3->SetLineStyle(2);

    cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4->SetLineColor(kBlack);
    cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4->SetFillColorAlpha(kBlack, transfill);
    cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4->SetFillStyle(3002);
    cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4->SetLineWidth(2);
    cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4->SetLineStyle(2);
    cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4->SetLineColor(kRed);
    cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4->SetFillColorAlpha(kRed, transfill);
    cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4->SetFillStyle(3002);
    cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4->SetLineWidth(2);
    cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4->SetLineStyle(2);
    cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5->SetLineColor(kBlue);
    cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5->SetFillColorAlpha(kBlue, transfill);
    cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5->SetFillStyle(3002);
    cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5->SetLineWidth(2);
    cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5->SetLineStyle(2);
    cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5->SetLineColor(kMagenta);
    cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5->SetFillColorAlpha(kMagenta, transfill);
    cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5->SetFillStyle(3002);
    cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5->SetLineWidth(2);
    cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5->SetLineStyle(2);

    cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4->SetLineColor(kBlack);
    cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4->SetFillColorAlpha(kBlack, transfill);
    cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4->SetFillStyle(3002);
    cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4->SetLineWidth(2);
    cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4->SetLineStyle(2);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4->SetLineColor(kRed);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4->SetFillColorAlpha(kRed, transfill);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4->SetFillStyle(3002);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4->SetLineWidth(2);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4->SetLineStyle(2);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4->SetLineColor(kBlue);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4->SetFillColorAlpha(kBlue, transfill);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4->SetFillStyle(3002);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4->SetLineWidth(2);
    cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4->SetLineStyle(2);

    cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5->SetLineColor(kBlack);
    cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5->SetFillColorAlpha(kBlack, transfill);
    cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5->SetFillStyle(3002);
    cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5->SetLineWidth(2);
    cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5->SetLineStyle(2);
    cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5->SetLineColor(kRed);
    cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5->SetFillColorAlpha(kRed, transfill);
    cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5->SetFillStyle(3002);
    cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5->SetLineWidth(2);
    cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5->SetLineStyle(2);
    cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5->SetLineColor(kBlue);
    cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5->SetFillColorAlpha(kBlue, transfill);
    cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5->SetFillStyle(3002);
    cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5->SetLineWidth(2);
    cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5->SetLineStyle(2);


    if (!fopen("plots","r")) system("mkdir plots");
    if (!fopen("plots/minorAxes/","r")) system("mkdir plots/minorAxes");
    if (!fopen("plots/majorAxes/","r")) system("mkdir plots/majorAxes");


    //-- scan of n = 1 to 6 for different cent values

    TCanvas * cxySpread = new TCanvas("cxySpread","cxySpread",740,850);
    cxySpread->Divide(5,6,0,0);
    // 1st row
    cxySpread->cd(1);
    xy1[0]->Draw();
    cxySpread->cd(2);
    xy1[2]->Draw();
    cxySpread->cd(3);
    xy1[6]->Draw();
    cxySpread->cd(4);
    xy1[8]->Draw();
    cxySpread->cd(5);
    xy1[11]->Draw();
    // 2nd row
    cxySpread->cd(6);
    xy2[0]->Draw();
    cxySpread->cd(7);
    xy2[2]->Draw();
    cxySpread->cd(8);
    xy2[6]->Draw();
    cxySpread->cd(9);
    xy2[8]->Draw();
    cxySpread->cd(10);
    xy2[11]->Draw();
    // 3rd row
    cxySpread->cd(11);
    xy3[0]->Draw();
    cxySpread->cd(12);
    xy3[2]->Draw();
    cxySpread->cd(13);
    xy3[6]->Draw();
    cxySpread->cd(14);
    xy3[8]->Draw();
    cxySpread->cd(15);
    xy3[11]->Draw();
    // 4th row
    cxySpread->cd(16);
    xy4[0]->Draw();
    cxySpread->cd(17);
    xy4[2]->Draw();
    cxySpread->cd(18);
    xy4[6]->Draw();
    cxySpread->cd(19);
    xy4[8]->Draw();
    cxySpread->cd(20);
    xy4[11]->Draw();
    // 5th row
    cxySpread->cd(21);
    xy5[0]->Draw();
    cxySpread->cd(22);
    xy5[2]->Draw();
    cxySpread->cd(23);
    xy5[6]->Draw();
    cxySpread->cd(24);
    xy5[8]->Draw();
    cxySpread->cd(25);
    xy5[11]->Draw();
    // 6th row
    cxySpread->cd(26);
    xy6[0]->Draw();
    cxySpread->cd(27);
    xy6[2]->Draw();
    cxySpread->cd(28);
    xy6[6]->Draw();
    cxySpread->cd(29);
    xy6[8]->Draw();
    cxySpread->cd(30);
    xy6[11]->Draw();

    int xySpreadtxtcol = kRed;

    TPaveText * txxyspread0 = new TPaveText(0.24, 0.78, 0.55, 0.99,"NDC");
    txxyspread0->SetFillColorAlpha(0,0);
    txxyspread0->SetBorderSize(0);
    txxyspread0->SetTextFont(43);
    txxyspread0->SetTextSize(24);
    txxyspread0->SetTextColor(xySpreadtxtcol);

    txxyspread0->AddText(Form("#bf{%d-%d%%}",(int)centbins[0],(int)centbins[1]));
    cxySpread->cd(1);
    txxyspread0->Draw();

    TPaveText * txxyspread1 = new TPaveText(0.16, 0.78, 0.56, 0.99,"NDC");
    txxyspread1->SetFillColorAlpha(0,0);
    txxyspread1->SetBorderSize(0);
    txxyspread1->SetTextFont(43);
    txxyspread1->SetTextSize(20);
    txxyspread1->SetTextColor(xySpreadtxtcol);
    TPaveText * txxyspread2 = (TPaveText *) txxyspread1->Clone();
    TPaveText * txxyspread3 = (TPaveText *) txxyspread1->Clone();
    TPaveText * txxyspread4 = (TPaveText *) txxyspread1->Clone();
    TPaveText * txxyspread5 = (TPaveText *) txxyspread1->Clone();
    txxyspread1->AddText(Form("#bf{%d-%d%%}",(int)centbins[2],(int)centbins[3]));
    cxySpread->cd(2);
    txxyspread1->Draw();
    txxyspread2->AddText(Form("#bf{%d-%d%%}",(int)centbins[6],(int)centbins[7]));
    cxySpread->cd(3);
    txxyspread2->Draw();
    txxyspread3->AddText(Form("#bf{%d-%d%%}",(int)centbins[8],(int)centbins[9]));
    cxySpread->cd(4);
    txxyspread3->Draw();
    txxyspread4->AddText(Form("#bf{%d-%d%%}",(int)centbins[11],(int)centbins[12]));
    cxySpread->cd(5);
    txxyspread4->Draw();

    TPaveText * txyPsiN1 = new TPaveText(0.62, 0.05, 0.95, 0.23,"NDC");
    txyPsiN1->SetFillColorAlpha(0,0);
    txyPsiN1->SetBorderSize(0);
    txyPsiN1->SetTextFont(43);
    txyPsiN1->SetTextSize(24);
    txyPsiN1->SetTextColor(kBlack);
    TPaveText * txyPsiN2 = (TPaveText *) txyPsiN1->Clone();
    TPaveText * txyPsiN3 = (TPaveText *) txyPsiN1->Clone();
    TPaveText * txyPsiN4 = (TPaveText *) txyPsiN1->Clone();
    TPaveText * txyPsiN5 = (TPaveText *) txyPsiN1->Clone();
    cxySpread->cd(5);
    txyPsiN1->AddText("#bf{n = 1}");
    txyPsiN1->Draw();
    cxySpread->cd(10);
    txyPsiN2->AddText("#bf{n = 2}");
    txyPsiN2->Draw();
    cxySpread->cd(15);
    txyPsiN3->AddText("#bf{n = 3}");
    txyPsiN3->Draw();
    cxySpread->cd(20);
    txyPsiN4->AddText("#bf{n = 4}");
    txyPsiN4->Draw();
    cxySpread->cd(25);
    txyPsiN5->AddText("#bf{n = 5}");
    txyPsiN5->Draw();
    TPaveText * txyPsiN6 = new TPaveText(0.62, 0.20, 0.95, 0.38,"NDC");
    txyPsiN6->SetFillColorAlpha(0,0);
    txyPsiN6->SetBorderSize(0);
    txyPsiN6->SetTextFont(43);
    txyPsiN6->SetTextSize(24);
    txyPsiN6->SetTextColor(kBlack);
    cxySpread->cd(30);
    txyPsiN6->AddText("#bf{n = 6}");
    txyPsiN6->Draw();

    if (print_plots && minorAxes) cxySpread->Print("plots/minorAxes/xyScan.png","png");
    else if (print_plots && !minorAxes) cxySpread->Print("plots/majorAxes/xyScan.png","png");
    if (close_plots) cxySpread->Close();


    //-- cumulants
    TCanvas * cxySpreadCumu = new TCanvas("cxySpreadCumu","cxySpreadCumu",1100,650);
    cxySpreadCumu->Divide(3,2,0,0);
    int centcumubin = 8;
    TPaveText * txyPsiN4cumu = new TPaveText(0.27, 0.06, 0.44, 0.16,"NDC");
    txyPsiN4cumu->SetFillColorAlpha(0,0);
    txyPsiN4cumu->SetBorderSize(0);
    txyPsiN4cumu->SetTextFont(43);
    txyPsiN4cumu->SetTextSize(24);
    txyPsiN4cumu->SetTextColor(kBlack);
    TPaveText * txyPsiN5cumu = (TPaveText *) txyPsiN4cumu->Clone();
    TPaveText * txyPsiN6cumu = (TPaveText *) txyPsiN4cumu->Clone();
    cxySpreadCumu->cd(1);
    xy4[centcumubin]->Draw();
    txyPsiN4cumu->AddText("n = 4");
    txyPsiN4cumu->Draw();
    TPaveText * txymoment = new TPaveText(0.25, 0.82, 0.56, 0.95,"NDC");
    txymoment->SetFillColorAlpha(0,0);
    txymoment->SetBorderSize(0);
    txymoment->SetTextFont(43);
    txymoment->SetTextSize(24);
    txymoment->SetTextColor(kBlack);
    TPaveText * txycumulant = (TPaveText *) txymoment->Clone();
    txymoment->AddText("Moments");
    txymoment->Draw();
    cxySpreadCumu->cd(2);
    xy5[centcumubin]->Draw();
    txyPsiN5cumu->AddText("n = 5");
    txyPsiN5cumu->Draw();
    cxySpreadCumu->cd(3);
    xy6[centcumubin]->Draw();
    txyPsiN6cumu->AddText("n = 6");
    txyPsiN6cumu->Draw();
    cxySpreadCumu->cd(4);
    xyC4[centcumubin]->Draw();
    txycumulant->AddText("Cumulants");
    txycumulant->Draw();
    cxySpreadCumu->cd(5);
    xyC5[centcumubin]->Draw();
    cxySpreadCumu->cd(6);
    xyC6[centcumubin]->Draw();

    if (print_plots && minorAxes) cxySpreadCumu->Print("plots/minorAxes/xyScanCumu.png","png");
    else if (print_plots && !minorAxes) cxySpreadCumu->Print("plots/majorAxes/xyScanCumu.png","png");
    if (close_plots) cxySpreadCumu->Close();


    //-- individual xy cartoons of harmonics
    int cxybin = 9;

    TCanvas * cXY1 = new TCanvas("cXY1","cXY1",400,400);
    TPad * padXY1 = (TPad *) cXY1->cd();
    padXY1->SetLogz();
    TH1D * xy1tmp = (TH1D *) xy1[cxybin]->Clone("xy1tmp");
    xy1[cxybin]->GetXaxis()->SetTicks("-");
    xy1[cxybin]->GetYaxis()->SetTicks("-");
    xy1[cxybin]->Draw("CONT1");
    if (print_plots && minorAxes) cXY1->Print("plots/minorAxes/v1.png","png");
    else if (print_plots && !minorAxes) cXY1->Print("plots/majorAxes/v1.png","png");
    if (close_plots) cXY1->Close();

    TCanvas * cXY2 = new TCanvas("cXY2","cXY2",400,400);
    TPad * padXY2 = (TPad *) cXY2->cd();
    padXY2->SetLogz();
    TH1D * xy2tmp = (TH1D *) xy2[cxybin]->Clone("xy1tmp");
    xy2[cxybin]->GetXaxis()->SetTicks("-");
    xy2[cxybin]->GetYaxis()->SetTicks("-");
    xy2[cxybin]->Draw("CONT1");
    if (print_plots && minorAxes) cXY2->Print("plots/minorAxes/v2.png","png");
    else if (print_plots && !minorAxes) cXY2->Print("plots/majorAxes/v2.png","png");
    if (close_plots) cXY2->Close();

    TCanvas * cXY3 = new TCanvas("cXY3","cXY3",400,400);
    TPad * padXY3 = (TPad *) cXY3->cd();
    padXY3->SetLogz();
    TH2D * xy3tmp = (TH2D *) xy3[cxybin]->Clone("xy3tmp");
    xy3[cxybin]->GetXaxis()->SetTicks("-");
    xy3[cxybin]->GetYaxis()->SetTicks("-");
    xy3[cxybin]->Draw("CONT1");
    if (print_plots && minorAxes) cXY3->Print("plots/minorAxes/v3.png","png");
    else if (print_plots && !minorAxes) cXY3->Print("plots/majorAxes/v3.png","png");
    if (close_plots) cXY3->Close();

    TCanvas * cXY4 = new TCanvas("cXY4","cXY4",400,400);
    TPad * padXY4 = (TPad *) cXY4->cd();
    padXY4->SetLogz();
    TH2D * xy4tmp = (TH2D *) xy4[cxybin]->Clone("xy4tmp");
    xy4[cxybin]->GetXaxis()->SetTicks("-");
    xy4[cxybin]->GetYaxis()->SetTicks("-");
    xy4[cxybin]->Draw("CONT1");
    if (print_plots && minorAxes) cXY4->Print("plots/minorAxes/v4.png","png");
    else if (print_plots && !minorAxes) cXY4->Print("plots/majorAxes/v4.png","png");
    if (close_plots) cXY4->Close();

    TCanvas * cXY5 = new TCanvas("cXY5","cXY5",400,400);
    TPad * padXY5 = (TPad *) cXY5->cd();
    padXY5->SetLogz();
    TH2D * xy5tmp = (TH2D *) xy5[cxybin]->Clone("xy5tmp");
    xy5[cxybin]->GetXaxis()->SetTicks("-");
    xy5[cxybin]->GetYaxis()->SetTicks("-");
    xy5[cxybin]->Draw("CONT1");
    if (print_plots && minorAxes) cXY5->Print("plots/minorAxes/v5.png","png");
    else if (print_plots && !minorAxes) cXY5->Print("plots/majorAxes/v5.png","png");
    if (close_plots) cXY5->Close();

    TCanvas * cXY6 = new TCanvas("cXY6","cXY6",400,400);
    TPad * padXY6 = (TPad *) cXY6->cd();
    padXY6->SetLogz();
    TH2D * xy6tmp = (TH2D *) xy6[cxybin]->Clone("xy6tmp");
    xy6[cxybin]->GetXaxis()->SetTicks("-");
    xy6[cxybin]->GetYaxis()->SetTicks("-");
    xy6[cxybin]->Draw("CONT1");
    if (print_plots && minorAxes) cXY6->Print("plots/minorAxes/v6.png","png");
    else if (print_plots && !minorAxes) cXY6->Print("plots/majorAxes/v6.png","png");
    if (close_plots) cXY6->Close();


    //-- eccentricity as a function of npart

    TCanvas * ceccCent = new TCanvas("ceccCent","ceccCent",700,650);
    ceccCent->cd();
    g[2]->GetXaxis()->SetTitle("#LTN_{part}#GT");
    g[2]->GetYaxis()->SetTitle("#epsilon");
    g[2]->GetXaxis()->SetTitleSize(0.05);
    g[2]->GetYaxis()->SetTitleSize(0.06);
    g[2]->GetXaxis()->SetLabelSize(0.05);
    g[2]->GetYaxis()->SetLabelSize(0.05);
    g[2]->GetXaxis()->SetTitleOffset(1.12);
    g[2]->GetYaxis()->SetTitleOffset(1.08);
    g[2]->Draw("alp");
    for (int iorder = 3; iorder<=6; iorder++) {
        g[iorder]->Draw("lp");
    }
    g[1]->Draw("lp");
    TLegend * lgeccCent = new TLegend(0.82, 0.59, 0.90, 0.91);
    lgeccCent->SetFillColorAlpha(0,0);
    lgeccCent->SetBorderSize(0);
    lgeccCent->SetTextFont(43);
    lgeccCent->SetTextSize(25);
    lgeccCent->AddEntry(g[1]," #epsilon_{1}","p");
    lgeccCent->AddEntry(g[2]," #epsilon_{2}","p");
    lgeccCent->AddEntry(g[3]," #epsilon_{3}","p");
    lgeccCent->AddEntry(g[4]," #epsilon_{4}","p");
    lgeccCent->AddEntry(g[5]," #epsilon_{5}","p");
    lgeccCent->AddEntry(g[6]," #epsilon_{6}","p");
    lgeccCent->Draw();
    if (print_plots && minorAxes) ceccCent->Print("plots/minorAxes/eccCent.png","png");
    else if (print_plots && !minorAxes) ceccCent->Print("plots/majorAxes/eccCent.png","png");
    if (close_plots) ceccCent->Close();


    // overlap area as a function of Npart
    TCanvas * cSCent = new TCanvas("cSCent","cSCent",700,650);
    cSCent->cd();
    sg[2]->GetXaxis()->SetTitle("#LTN_{part}#GT");
    sg[2]->GetYaxis()->SetTitle("S fm^{2}");
    sg[2]->GetXaxis()->SetTitleSize(0.05);
    sg[2]->GetYaxis()->SetTitleSize(0.06);
    sg[2]->GetXaxis()->SetLabelSize(0.05);
    sg[2]->GetYaxis()->SetLabelSize(0.05);
    sg[2]->GetXaxis()->SetTitleOffset(1.12);
    sg[2]->GetYaxis()->SetTitleOffset(1.08);
    sg[2]->Draw("alp");
    for (int iorder = 3; iorder<=6; iorder++) {
        sg[iorder]->Draw("lp");
    }
    g[1]->Draw("lp");
    TLegend * lgSCent = new TLegend(0.79, 0.20, 0.91, 0.52);
    lgSCent->SetFillColorAlpha(0,0);
    lgSCent->SetBorderSize(0);
    lgSCent->SetTextFont(43);
    lgSCent->SetTextSize(25);
    lgSCent->AddEntry(sg[1]," n = 1","p");
    lgSCent->AddEntry(sg[2]," n = 2","p");
    lgSCent->AddEntry(sg[3]," n = 3","p");
    lgSCent->AddEntry(sg[4]," n = 4","p");
    lgSCent->AddEntry(sg[5]," n = 5","p");
    lgSCent->AddEntry(sg[6]," n = 6","p");
    lgSCent->Draw();
    if (print_plots && minorAxes) cSCent->Print("plots/minorAxes/SCent.png","png");
    else if (print_plots && !minorAxes) cSCent->Print("plots/majorAxes/SCent.png","png");
    if (close_plots) cSCent->Close();


    //-- Participant 2-plane correlations

    if (minorAxes) {
        TCanvas * c2pcorr = new TCanvas("c2pcorr","c2pcorr",1250,600);
        c2pcorr->Divide(4,2);

        TPaveText * tx2pcorr = new TPaveText(0.38, 0.81, 0.88, 0.90,"NDC");
        tx2pcorr->SetFillColorAlpha(0,0);
        tx2pcorr->SetBorderSize(0);
        tx2pcorr->SetTextFont(43);
        tx2pcorr->SetTextSize(18);

        TH1D * hdummy_2corr = new TH1D("hdummy_2corr","hdummy_2corr", 50, 0, 425);
        hdummy_2corr->SetTitle("");
        hdummy_2corr->SetStats(0);
        hdummy_2corr->GetYaxis()->SetRangeUser(-1, 1);
        hdummy_2corr->SetXTitle("#LTN_{part}#GT");
        hdummy_2corr->GetXaxis()->CenterTitle(kTRUE);
        hdummy_2corr->GetXaxis()->SetTitleSize(0.06);
        hdummy_2corr->GetXaxis()->SetTitleOffset(1.08);

        TLegend * leg2pcorr = new TLegend(0.41, 0.19, 0.85, 0.33);
        leg2pcorr->SetFillColorAlpha(0,0);
        leg2pcorr->SetBorderSize(0);
        leg2pcorr->SetTextFont(43);
        leg2pcorr->SetTextSize(16);
        leg2pcorr->AddEntry(cg_2p2_3p3_5p5,"Toy Glauber","p");
        leg2pcorr->AddEntry(cgthy_2p2_3p3_5p5[0],"PRC 90, 024902","p");

        c2pcorr->cd(1);
        TH1D * h2pcorr_1 = (TH1D *) hdummy_2corr->Clone("h2pcorr_1");
        h2pcorr_1->GetYaxis()->SetRangeUser(-1, 1);
        h2pcorr_1->Draw();
        cg_4_p4_p2->Draw("same p");
        cgC_4_p4_p2->Draw("same p");
        cgthy_4_p4_p2[0]->Draw("same p");
        cgthy_4_p4_p2[1]->Draw("same p");
        TPaveText * tx2pcorr_1 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_1");
        tx2pcorr_1->AddText("#LT#LTcos4(#Psi_{4} - #Psi_{2})#GT#GT");
        tx2pcorr_1->Draw();
        c2pcorr->cd(2);
        TH1D * h2pcorr_2 = (TH1D *) hdummy_2corr->Clone("h2pcorr_2");
        h2pcorr_2->GetYaxis()->SetRangeUser(-0.1, 0.6);
        h2pcorr_2->Draw();
        cg_8_p4_p2->Draw("same p");
        cgC_8_p4_p2->Draw("same p");
        cgthy_8_p4_p2[0]->Draw("same p");
        cgthy_8_p4_p2[1]->Draw("same p");
        TPaveText * tx2pcorr_2 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_2");
        tx2pcorr_2->AddText("#LT#LTcos8(#Psi_{4} - #Psi_{2})#GT#GT");
        tx2pcorr_2->Draw();
        c2pcorr->cd(3);
        TH1D * h2pcorr_3 = (TH1D *) hdummy_2corr->Clone("h2pcorr_3");
        h2pcorr_3->GetYaxis()->SetRangeUser(-0.5, 0.5);
        h2pcorr_3->Draw();
        cg_12_p4_p2->Draw("same p");
        cgC_12_p4_p2->Draw("same p");
        cgthy_12_p4_p2[0]->Draw("same p");
        cgthy_12_p4_p2[1]->Draw("same p");
        TPaveText * tx2pcorr_3 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_3");
        tx2pcorr_3->AddText("#LT#LTcos12(#Psi_{4} - #Psi_{2})#GT#GT");
        tx2pcorr_3->Draw();
        c2pcorr->cd(4);
        TH1D * h2pcorr_4 = (TH1D *) hdummy_2corr->Clone("h2pcorr_4");
        h2pcorr_4->GetYaxis()->SetRangeUser(-0.6, 0.6);
        h2pcorr_4->Draw();
        cg_6_p3_p2->Draw("same p");
        cgthy_6_p3_p2[0]->Draw("same p");
        TPaveText * tx2pcorr_4 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_4");
        tx2pcorr_4->AddText("#LT#LTcos6(#Psi_{3} - #Psi_{2})#GT#GT#times10");
        tx2pcorr_4->Draw();
        c2pcorr->cd(5);
        TH1D * h2pcorr_5 = (TH1D *) hdummy_2corr->Clone("h2pcorr_5");
        h2pcorr_5->GetYaxis()->SetRangeUser(-0.2, 0.6);
        h2pcorr_5->Draw();
        cg_6_p2_p6->Draw("same p");
        cgC_6_p2_p6->Draw("same p");
        cgthy_6_p2_p6[0]->Draw("same p");
        cgthy_6_p2_p6[1]->Draw("same p");
        TPaveText * tx2pcorr_5 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_5");
        tx2pcorr_5->AddText("#LT#LTcos6(#Psi_{2} - #Psi_{6})#GT#GT");
        tx2pcorr_5->Draw();
        c2pcorr->cd(6);
        TH1D * h2pcorr_6 = (TH1D *) hdummy_2corr->Clone("h2pcorr_6");
        h2pcorr_6->GetYaxis()->SetRangeUser(-0.5, 0.5);
        h2pcorr_6->Draw();
        cg_6_p3_p6->Draw("same p");
        cgC_6_p3_p6->Draw("same p");
        cgthy_6_p3_p6[0]->Draw("same p");
        cgthy_6_p3_p6[1]->Draw("same p");
        TPaveText * tx2pcorr_6 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_6");
        tx2pcorr_6->AddText("#LT#LTcos6(#Psi_{3} - #Psi_{2})#GT#GT");
        tx2pcorr_6->Draw();
        c2pcorr->cd(7);
        TH1D * h2pcorr_7 = (TH1D *) hdummy_2corr->Clone("h2pcorr_7");
        h2pcorr_7->GetYaxis()->SetRangeUser(-0.1, 0.1);
        h2pcorr_7->Draw();
        cg_12_p3_p4->Draw("same p");
        cgC_12_p3_p4->Draw("same p");
        cgthy_12_p3_p4[0]->Draw("same p");
        cgthy_12_p3_p4[1]->Draw("same p");
        TPaveText * tx2pcorr_7 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_7");
        tx2pcorr_7->AddText("#LT#LTcos12(#Psi_{3} - #Psi_{4})#GT#GT");
        tx2pcorr_7->Draw();
        c2pcorr->cd(8);
        TH1D * h2pcorr_8 = (TH1D *) hdummy_2corr->Clone("h2pcorr_8");
        h2pcorr_8->GetYaxis()->SetRangeUser(-1.5, 1.5);
        h2pcorr_8->Draw();
        cg_10_p2_p5->Draw("same p");
        cgC_10_p2_p5->Draw("same p");
        cgthy_10_p2_p5[0]->Draw("same p");
        cgthy_10_p2_p5[1]->Draw("same p");
        TPaveText * tx2pcorr_8 = (TPaveText *) tx2pcorr->Clone("tx2pcorr_8");
        tx2pcorr_8->AddText("#LT#LTcos10(#Psi_{2} - #Psi_{5})#GT#GT#times10");
        tx2pcorr_8->Draw();
        leg2pcorr->Draw();

        if (print_plots && minorAxes) c2pcorr->Print("plots/minorAxes/2PlaneCorr.png","png");
        else if (print_plots && !minorAxes) c2pcorr->Print("plots/majorAxes/2PlaneCorr.png","png");
        if (close_plots) c2pcorr->Close();


        //-- Participant 3-plane correlations

        TCanvas * c3pcorr = new TCanvas("c3pcorr","c3pcorr",1100,700);
        c3pcorr->Divide(3,2);

        TPaveText * tx3pcorr = new TPaveText(0.38, 0.81, 0.88, 0.90,"NDC");
        tx3pcorr->SetFillColorAlpha(0,0);
        tx3pcorr->SetBorderSize(0);
        tx3pcorr->SetTextFont(43);
        tx3pcorr->SetTextSize(18);

        TH1D * hdummy_3corr = new TH1D("hdummy_3corr","hdummy_3corr", 50, 0, 425);
        hdummy_3corr->SetTitle("");
        hdummy_3corr->SetStats(0);
        hdummy_3corr->GetYaxis()->SetRangeUser(-1, 1);
        hdummy_3corr->SetXTitle("#LTN_{part}#GT");
        hdummy_3corr->GetXaxis()->CenterTitle(kTRUE);
        hdummy_3corr->GetXaxis()->SetTitleSize(0.06);
        hdummy_3corr->GetXaxis()->SetTitleOffset(1.08);

        TLegend * leg3pcorr = new TLegend(0.48, 0.65, 0.92, 0.79);
        leg3pcorr->SetFillColorAlpha(0,0);
        leg3pcorr->SetBorderSize(0);
        leg3pcorr->SetTextFont(43);
        leg3pcorr->SetTextSize(16);
        leg3pcorr->AddEntry(cg_2p2_3p3_5p5,"Toy Glauber","p");
        leg3pcorr->AddEntry(cgthy_2p2_3p3_5p5[0],"PRC 90, 024902","p");

        c3pcorr->cd(1);
        TH1D * h3pcorr_1 = (TH1D *) hdummy_3corr->Clone("h3pcorr_1");
        h3pcorr_1->GetYaxis()->SetRangeUser(-1, 1);
        h3pcorr_1->Draw();
        cg_2p2_3p3_5p5->Draw("same p");
        cgC_2p2_3p3_5p5->Draw("same p");
        cgthy_2p2_3p3_5p5[0]->Draw("same p");
        cgthy_2p2_3p3_5p5[1]->Draw("same p");
        TPaveText * tx3pcorr_1 = (TPaveText *) tx2pcorr->Clone("tx3pcorr_1");
        tx3pcorr_1->AddText("#LT#LTcos(2#Psi_{2} + 3#Psi_{3} - 5#Psi_{5})#GT#GT");
        tx3pcorr_1->Draw();
        leg3pcorr->Draw();
        c3pcorr->cd(2);
        TH1D * h3pcorr_2 = (TH1D *) hdummy_3corr->Clone("h3pcorr_2");
        h3pcorr_2->GetYaxis()->SetRangeUser(-1, 1);
        h3pcorr_2->Draw();
        cg_2p2_4p4_6p6->Draw("same p");
        cgC_2p2_4p4_6p6->Draw("same p");
        cgthy_2p2_4p4_6p6[0]->Draw("same p");
        cgthy_2p2_4p4_6p6[1]->Draw("same p");
        TPaveText * tx3pcorr_2 = (TPaveText *) tx2pcorr->Clone("tx3pcorr_2");
        tx3pcorr_2->AddText("#LT#LTcos(2#Psi_{2} + 4#Psi_{4} - 6#Psi_{6})#GT#GT");
        tx3pcorr_2->Draw();
        c3pcorr->cd(3);
        TH1D * h3pcorr_3 = (TH1D *) hdummy_3corr->Clone("h3pcorr_3");
        h3pcorr_3->GetYaxis()->SetRangeUser(-0.3, 0.3);
        h3pcorr_3->Draw();
        cg_2p2_6p3_4p4->Draw("same p");
        cgC_2p2_6p3_4p4->Draw("same p");
        cgthy_2p2_6p3_4p4[0]->Draw("same p");
        cgthy_2p2_6p3_4p4[1]->Draw("same p");
        TPaveText * tx3pcorr_3 = (TPaveText *) tx2pcorr->Clone("tx3pcorr_3");
        tx3pcorr_3->AddText("#LT#LTcos(2#Psi_{2} - 6#Psi_{3} + 4#Psi_{4})#GT#GT");
        tx3pcorr_3->Draw();
        c3pcorr->cd(4);
        TH1D * h3pcorr_4 = (TH1D *) hdummy_3corr->Clone("h3pcorr_4");
        h3pcorr_4->GetYaxis()->SetRangeUser(-0.1, 0.2);
        h3pcorr_4->Draw();
        cg_8p2_3p3_5p5->Draw("same p");
        cgC_8p2_3p3_5p5->Draw("same p");
        cgthy_8p2_3p3_5p5[0]->Draw("same p");
        cgthy_8p2_3p3_5p5[1]->Draw("same p");
        TPaveText * tx3pcorr_4 = (TPaveText *) tx2pcorr->Clone("tx3pcorr_4");
        tx3pcorr_4->AddText("#LT#LTcos(-8#Psi_{2} + 3#Psi_{3} + 5#Psi_{5})#GT#GT");
        tx3pcorr_4->Draw();
        c3pcorr->cd(5);
        TH1D * h3pcorr_5 = (TH1D *) hdummy_3corr->Clone("h3pcorr_5");
        h3pcorr_5->GetYaxis()->SetRangeUser(-0.5, 0.5);
        h3pcorr_5->Draw();
        cg_10p2_4p4_6p6->Draw("same p");
        cgC_10p2_4p4_6p6->Draw("same p");
        cgthy_10p2_4p4_6p6[0]->Draw("same p");
        cgthy_10p2_4p4_6p6[1]->Draw("same p");
        TPaveText * tx3pcorr_5 = (TPaveText *) tx2pcorr->Clone("tx3pcorr_5");
        tx3pcorr_5->AddText("#LT#LTcos(-10#Psi_{2} + 4#Psi_{4} + 6#Psi_{6})#GT#GT");
        tx3pcorr_5->Draw();
        c3pcorr->cd(6);
        TH1D * h3pcorr_6 = (TH1D *) hdummy_3corr->Clone("h3pcorr_6");
        h3pcorr_6->GetYaxis()->SetRangeUser(-0.1, 0.25);
        h3pcorr_6->Draw();
        cg_10p2_6p3_4p4->Draw("same p");
        cgC_10p2_6p3_4p4->Draw("same p");
        cgthy_10p2_6p3_4p4[0]->Draw("same p");
        cgthy_10p2_6p3_4p4[1]->Draw("same p");
        TPaveText * tx3pcorr_6 = (TPaveText *) tx2pcorr->Clone("tx3pcorr_6");
        tx3pcorr_6->AddText("#LT#LTcos(-10#Psi_{2} + 6#Psi_{3} + 4#Psi_{4})#GT#GT");
        tx3pcorr_6->Draw();

        if (print_plots && minorAxes) c3pcorr->Print("plots/minorAxes/3PlaneCorr.png","png");
        else if (print_plots && !minorAxes) c3pcorr->Print("plots/majorAxes/3PlaneCorr.png","png");
        if (close_plots) c3pcorr->Close();
    }


    //-- psi angle distribution

    int cangmin = 10;
    int cangmax = 12;
    TH1D * hpsi_new[7];
    for (int iorder = 1; iorder<=6; iorder++) {
        hpsi_new[iorder] = (TH1D *) hpsi[iorder][cangmin]->Clone(Form("psi%d_new",iorder));
        for (int cbin = cangmin+1; cbin<cangmax; cbin++) {
            hpsi_new[iorder]->Add(hpsi[iorder][cbin]);
            hpsi_new[iorder]->SetXTitle(Form("#psi_{%d} Centrality %d-%d%%",iorder,(int)centbins[cangmin],(int)centbins[cangmax]));
        }
    }
    for (int iorder = 1; iorder<=6; iorder++) {
        hpsi_new[iorder]->Scale(1/hpsi_new[iorder]->Integral("width"));
    }
    if (minorAxes) {
        TCanvas * cpsiAngDist = new TCanvas("cpsiAngDist","cpsiAngDist",1100,700);
        cpsiAngDist->Divide(3,2);
        for (int iorder = 1; iorder<=6; iorder++) {
            cpsiAngDist->cd(iorder);
            hpsi_new[iorder]->Draw("same");
        }
        if (print_plots && minorAxes) cpsiAngDist->Print(Form("plots/minorAxes/psiAngDist_%d-%d.png",(int)centbins[cangmin],(int)centbins[cangmax]),"png");
        else if (print_plots && !minorAxes) cpsiAngDist->Print(Form("plots/majorAxes/psiAngDist_%d-%d.png",(int)centbins[cangmin],(int)centbins[cangmax]),"png");
        if (close_plots) cpsiAngDist->Close();
    }


    //-- Psi2-Psi1 and Psi3-Psi1 correlations

    int ccorrmin = 4;
    int ccorrmax = 8;
    TH2D * Psi2Psi1_new = (TH2D *) Psi2Psi1[ccorrmin]->Clone("Psi2Psi1_new");
    TH2D * Psi3Psi1_new = (TH2D *) Psi3Psi1[ccorrmin]->Clone("Psi3Psi1_new");
    TH1D * CosPsi2Psi1_new = (TH1D *) CosPsi2Psi1[ccorrmin]->Clone("CosPsi2Psi1_new");
    TH1D * CosPsi3Psi1_new = (TH1D *) CosPsi3Psi1[ccorrmin]->Clone("CosPsi3Psi1_new");
    for (int cbin = ccorrmin+1; cbin<ccorrmax; cbin++) {
        Psi2Psi1_new->Add(Psi2Psi1[cbin]);
        Psi3Psi1_new->Add(Psi3Psi1[cbin]);
        CosPsi2Psi1_new->Add(CosPsi2Psi1[cbin]);
        CosPsi3Psi1_new->Add(CosPsi3Psi1[cbin]);
    }

    TCanvas * cPsi2Psi1 = new TCanvas("cPsi2Psi1","cPsi2Psi1",700,650);
    cPsi2Psi1->cd();
    Psi2Psi1_new->SetXTitle("#Psi_{1}");
    Psi2Psi1_new->GetXaxis()->SetTitleSize(0.06);
    Psi2Psi1_new->GetYaxis()->SetTitleSize(0.06);
    Psi2Psi1_new->GetYaxis()->SetTitleOffset(1.10);
    Psi2Psi1_new->Draw("col");
    if (print_plots && minorAxes) cPsi2Psi1->Print("plots/minorAxes/Psi2Psi1.png","png");
    else if (print_plots && !minorAxes) cPsi2Psi1->Print("plots/majorAxes/Psi2Psi1.png","png");
    if (close_plots) cPsi2Psi1->Close();

    TCanvas * cPsi3Psi1 = new TCanvas("cPsi3Psi1","cPsi3Psi1",700,650);
    cPsi3Psi1->cd();
    Psi3Psi1_new->SetXTitle("#Psi_{1}");
    Psi3Psi1_new->GetXaxis()->SetTitleSize(0.06);
    Psi3Psi1_new->GetYaxis()->SetTitleSize(0.06);
    Psi3Psi1_new->GetYaxis()->SetTitleOffset(1.10);
    Psi3Psi1_new->Draw("col");
    if (print_plots && minorAxes) cPsi3Psi1->Print("plots/minorAxes/Psi3Psi1.png","png");
    else if (print_plots && !minorAxes) cPsi3Psi1->Print("plots/majorAxes/Psi3Psi1.png","png");
    if (close_plots) cPsi3Psi1->Close();

    TCanvas * cCosPsi2Psi1 = new TCanvas("cCosPsi2Psi1","cCosPsi2Psi1",700,650);
    cCosPsi2Psi1->cd();
    CosPsi2Psi1_new->SetXTitle("cos(2#Psi_{2} - #Psi_{1})");
    CosPsi2Psi1_new->GetXaxis()->SetTitleSize(0.05);
    CosPsi2Psi1_new->GetXaxis()->SetTitleOffset(1.25);
    CosPsi2Psi1_new->Draw();
    TPaveText * txCosPsi2Psi1_new = new TPaveText(0.62, 0.83, 0.89, 0.89,"NDC");
    txCosPsi2Psi1_new->SetFillColorAlpha(0,0);
    txCosPsi2Psi1_new->SetBorderSize(0);
    txCosPsi2Psi1_new->SetTextFont(43);
    txCosPsi2Psi1_new->SetTextSize(24);
    txCosPsi2Psi1_new->AddText(Form("%d-%d%% Centrality",(int)centbins[ccorrmin],(int)centbins[ccorrmax]));
    txCosPsi2Psi1_new->Draw();
    if (print_plots && minorAxes) cCosPsi2Psi1->Print("plots/minorAxes/CosPsi2Psi1.png","png");
    else if (print_plots && !minorAxes) cCosPsi2Psi1->Print("plots/majorAxes/CosPsi2Psi1.png","png");
    if (close_plots) cCosPsi2Psi1->Close();

    TCanvas * cCosPsi3Psi1 = new TCanvas("cCosPsi3Psi1","cCosPsi3Psi1",700,650);
    cCosPsi3Psi1->cd();
    CosPsi3Psi1_new->SetXTitle("cos(3#Psi_{2} - #Psi_{1})");
    CosPsi3Psi1_new->GetXaxis()->SetTitleSize(0.05);
    CosPsi3Psi1_new->GetXaxis()->SetTitleOffset(1.25);
    CosPsi3Psi1_new->Draw();
    TPaveText * txCosPsi3Psi1_new = new TPaveText(0.62, 0.83, 0.89, 0.89,"NDC");
    txCosPsi3Psi1_new->SetFillColorAlpha(0,0);
    txCosPsi3Psi1_new->SetBorderSize(0);
    txCosPsi3Psi1_new->SetTextFont(43);
    txCosPsi3Psi1_new->SetTextSize(24);
    txCosPsi3Psi1_new->AddText(Form("%d-%d%% Centrality",(int)centbins[ccorrmin],(int)centbins[ccorrmax]));
    txCosPsi3Psi1_new->Draw();
    if (print_plots && minorAxes) cCosPsi3Psi1->Print("plots/minorAxes/CosPsi3Psi1.png","png");
    else if (print_plots && !minorAxes) cCosPsi3Psi1->Print("plots/majorAxes/CosPsi3Psi1.png","png");
    if (close_plots) cCosPsi3Psi1->Close();


    //-- more correlation plots
    if (!minorAxes) {
        TH1D * hdummy_corr_00 = new TH1D("hdummy_corr_00","hdummy_corr_00", 50, 0, 425);
        hdummy_corr_00->SetTitle("");
        hdummy_corr_00->SetStats(0);
        hdummy_corr_00->GetYaxis()->SetRangeUser(-0.1, 0.89);
        hdummy_corr_00->SetXTitle("#LTN_{part}#GT");
        hdummy_corr_00->GetXaxis()->CenterTitle(kTRUE);
        hdummy_corr_00->GetXaxis()->SetTitleSize(0.06);
        hdummy_corr_00->GetXaxis()->SetTitleOffset(1.08);

        TCanvas * cCorr_00 = new TCanvas("cCorr_00","cCorr_00",550,750);
        cCorr_00->Divide(1,2,0,0.01);
        cCorr_00->cd(1);
        TH1D * hCorr_00_0 = (TH1D *) hdummy_corr_00->Clone("hCorr_00_0");
        hCorr_00_0->Draw();
        cg_6p2_p3->Draw("same p");
        cg_4p2_p4->Draw("same p");
        cg_6p2_p6->Draw("same p");
        cg_6p3_p6->Draw("same p");
        cg_10p2_p5->Draw("same p");
        cg_15p3_p5->Draw("same p");
        // cg_EPJ_Fig2_cos6Phi2_Phi3->Draw("same 3");
        // cg_EPJ_Fig2_cos4Phi2_Phi4->Draw("same 3");
        // cg_EPJ_Fig2_cos6Phi2_Phi6->Draw("same 3");
        // cg_EPJ_Fig2_cos6Phi3_Phi6->Draw("same 3");
        // cg_EPJ_Fig2_cos10Phi2_Phi5->Draw("same 3");
        // cg_EPJ_Fig2_cos15Phi3_Phi5->Draw("same 3");
        TLegend * legCorr_00_0 = new TLegend(0.62, 0.56, 0.82, 0.96);
        legCorr_00_0->SetFillColorAlpha(0,0);
        legCorr_00_0->SetBorderSize(0);
        legCorr_00_0->SetTextFont(43);
        legCorr_00_0->SetTextSize(18);
        legCorr_00_0->AddEntry(cg_6p2_p3,"#LTcos6(#Phi_{2}* - #Phi_{3}*)#GT","lp");
        legCorr_00_0->AddEntry(cg_4p2_p4,"#LTcos4(#Phi_{2}* - #Phi_{4}*)#GT","lp");
        legCorr_00_0->AddEntry(cg_6p2_p6,"#LTcos6(#Phi_{2}* - #Phi_{6}*)#GT","lp");
        legCorr_00_0->AddEntry(cg_6p3_p6,"#LTcos6(#Phi_{3}* - #Phi_{6}*)#GT","lp");
        legCorr_00_0->AddEntry(cg_10p2_p5,"#LTcos10(#Phi_{2}* - #Phi_{5}*)#GT","lp");
        legCorr_00_0->AddEntry(cg_15p3_p5,"#LTcos15(#Phi_{3}* - #Phi_{5}*)#GT","lp");
        legCorr_00_0->Draw();

        cCorr_00->cd(2);
        TH1D * hCorr_00_1 = (TH1D *) hdummy_corr_00->Clone("hCorr_00_1");
        hCorr_00_1->GetYaxis()->SetRangeUser(-0.02, 0.62);
        hCorr_00_1->Draw();
        cg_2p1_p2->Draw("same p");
        cg_3p1_p3->Draw("same p");
        cg_4p1_p4->Draw("same p");
        cg_5p1_p5->Draw("same p");
        cg_6p1_p6->Draw("same p");
        // cg_EPJ_Fig2_cos2Phi1_Phi2->Draw("same 3");
        // cg_EPJ_Fig2_cos3Phi1_Phi3->Draw("same 3");
        // cg_EPJ_Fig2_cos4Phi1_Phi4->Draw("same 3");
        // cg_EPJ_Fig2_cos5Phi1_Phi5->Draw("same 3");
        // cg_EPJ_Fig2_cos6Phi1_Phi6->Draw("same 3");
        TLegend * legCorr_00_1 = new TLegend(0.62, 0.67, 0.82, 0.96);
        legCorr_00_1->SetFillColorAlpha(0,0);
        legCorr_00_1->SetBorderSize(0);
        legCorr_00_1->SetTextFont(43);
        legCorr_00_1->SetTextSize(18);
        legCorr_00_1->AddEntry(cg_2p1_p2,"#LTcos2(#Phi_{1}* - #Phi_{2}*)#GT","lp");
        legCorr_00_1->AddEntry(cg_3p1_p3,"#LTcos3(#Phi_{1}* - #Phi_{3}*)#GT","lp");
        legCorr_00_1->AddEntry(cg_4p1_p4,"#LTcos4(#Phi_{1}* - #Phi_{4}*)#GT","lp");
        legCorr_00_1->AddEntry(cg_5p1_p5,"#LTcos5(#Phi_{1}* - #Phi_{5}*)#GT","lp");
        legCorr_00_1->AddEntry(cg_6p1_p6,"#LTcos6(#Phi_{1}* - #Phi_{6}*)#GT","lp");
        legCorr_00_1->Draw();

        if (print_plots && !minorAxes) cCorr_00->Print("plots/majorAxes/Corr_00.png","png");
        if (close_plots) cCorr_00->Close();


        TH1D * hdummy_corr_01 = new TH1D("hdummy_corr_01","hdummy_corr_01", 50, 0, 425);
        hdummy_corr_01->SetTitle("");
        hdummy_corr_01->SetStats(0);
        hdummy_corr_01->GetYaxis()->SetRangeUser(-0.1, 0.89);
        hdummy_corr_01->SetXTitle("#LTN_{part}#GT");
        hdummy_corr_01->GetXaxis()->CenterTitle(kTRUE);
        hdummy_corr_01->GetXaxis()->SetTitleSize(0.06);
        hdummy_corr_01->GetXaxis()->SetTitleOffset(1.08);

        TCanvas * cCorr_01 = new TCanvas("cCorr_01","cCorr_01",800,800);
        cCorr_01->Divide(2,2,0,0);
        cCorr_01->cd(1);
        TH1D * hCorr_01_1 = (TH1D *) hdummy_corr_01->Clone("hCorr_01_1");
        hCorr_01_1->Draw();
        cg_1p1_2p2_3p3->Draw("same p");
        cg_4p1_2p2_6p3->Draw("same p");
        cg_2p1_4p2_6p3->Draw("same p");
        cg_5p1_2p2_3p3->Draw("same p");
        // cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3->Draw("same 3");
        // cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3->Draw("same 3");
        // cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3->Draw("same 3");
        // cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3->Draw("same 3");
        TLegend * legCorr_01_1 = new TLegend(0.43, 0.67, 0.92, 0.96);
        legCorr_01_1->SetFillColorAlpha(0,0);
        legCorr_01_1->SetBorderSize(0);
        legCorr_01_1->SetTextFont(43);
        legCorr_01_1->SetTextSize(17);
        legCorr_01_1->AddEntry(cg_1p1_2p2_3p3,"#LTcos(#Phi_{1}* + 2#Phi_{2}* - 3#Phi_{3}*)#GT","lp");
        legCorr_01_1->AddEntry(cg_4p1_2p2_6p3,"#LTcos(4#Phi_{1}* + 2#Phi_{2}* - 6#Phi_{3}*)#GT","lp");
        legCorr_01_1->AddEntry(cg_2p1_4p2_6p3,"#LTcos(2#Phi_{1}* + 4#Phi_{2}* - 6#Phi_{3}*)#GT","lp");
        legCorr_01_1->AddEntry(cg_5p1_2p2_3p3,"#LTcos(5#Phi_{1}* - 2#Phi_{2}* - 3#Phi_{3}*)#GT","lp");
        legCorr_01_1->Draw();
        cCorr_01->cd(2);
        TH1D * hCorr_01_2 = (TH1D *) hdummy_corr_01->Clone("hCorr_01_2");
        hCorr_01_2->GetYaxis()->SetLabelSize(0);
        hCorr_01_2->Draw();
        cg_2p1_2p2_4p4->Draw("same p");
        cg_2p1_6p2_4p4->Draw("same p");
        cg_2p1_6p2_8p4->Draw("same p");
        // cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4->Draw("same 3");
        // cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4->Draw("same 3");
        // cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4->Draw("same 3");
        TLegend * legCorr_01_2 = new TLegend(0.31, 0.74, 0.84, 0.96);
        legCorr_01_2->SetFillColorAlpha(0,0);
        legCorr_01_2->SetBorderSize(0);
        legCorr_01_2->SetTextFont(43);
        legCorr_01_2->SetTextSize(17);
        legCorr_01_2->AddEntry(cg_2p1_2p2_4p4,"#LTcos(2#Phi_{1}* + 2#Phi_{2}* - 4#Phi_{4}*)#GT","lp");
        legCorr_01_2->AddEntry(cg_2p1_6p2_4p4,"#LTcos(2#Phi_{1}* + 6#Phi_{2}* - 4#Phi_{4}*)#GT","lp");
        legCorr_01_2->AddEntry(cg_2p1_6p2_8p4,"#LTcos(2#Phi_{1}* + 6#Phi_{2}* - 8#Phi_{4}*)#GT","lp");
        legCorr_01_2->Draw();
        cCorr_01->cd(3);
        TH1D * hCorr_01_3 = (TH1D *) hdummy_corr_01->Clone("hCorr_01_3");
        hCorr_01_3->Draw();
        cg_1p1_3p3_4p4->Draw("same p");
        cg_2p1_6p3_4p4->Draw("same p");
        cg_2p1_4p2_5p5->Draw("same p");
        cg_3p1_2p2_5p5->Draw("same p");
        // cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4->Draw("same 3");
        // cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4->Draw("same 3");
        // cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5->Draw("same 3");
        // cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5->Draw("same 3");
        TLegend * legCorr_01_3 = new TLegend(0.43, 0.67, 0.92, 0.96);
        legCorr_01_3->SetFillColorAlpha(0,0);
        legCorr_01_3->SetBorderSize(0);
        legCorr_01_3->SetTextFont(43);
        legCorr_01_3->SetTextSize(17);
        legCorr_01_3->AddEntry(cg_1p1_3p3_4p4,"#LTcos(#Phi_{1}* + 3#Phi_{3}* - 4#Phi_{4}*)#GT","lp");
        legCorr_01_3->AddEntry(cg_2p1_6p3_4p4,"#LTcos(2#Phi_{1}* - 6#Phi_{3}* + 4#Phi_{4}*)#GT","lp");
        legCorr_01_3->AddEntry(cg_2p1_4p2_5p5,"#LTcos(#Phi_{1}* + 4#Phi_{2}* - 5#Phi_{5}*)#GT","lp");
        legCorr_01_3->AddEntry(cg_3p1_2p2_5p5,"#LTcos(3#Phi_{1}* + 2#Phi_{2}* - 5#Phi_{5}*)#GT","lp");
        legCorr_01_3->Draw();
        cCorr_01->cd(4);
        TH1D * hCorr_01_4 = (TH1D *) hdummy_corr_01->Clone("hCorr_01_4");
        hCorr_01_4->GetYaxis()->SetLabelSize(0);
        hCorr_01_4->Draw();
        cg_1p1_6p3_5p5->Draw("same p");
        cg_2p1_3p3_5p5->Draw("same p");
        cg_4p1_9p3_5p5->Draw("same p");
        // cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5->Draw("same 3");
        // cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5->Draw("same 3");
        // cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5->Draw("same 3");
        TLegend * legCorr_01_4 = new TLegend(0.31, 0.74, 0.84, 0.96);
        legCorr_01_4->SetFillColorAlpha(0,0);
        legCorr_01_4->SetBorderSize(0);
        legCorr_01_4->SetTextFont(43);
        legCorr_01_4->SetTextSize(17);
        legCorr_01_4->AddEntry(cg_1p1_6p3_5p5,"#LTcos(#Phi_{1}* - 6#Phi_{3}* - 5#Phi_{5}*)#GT","lp");
        legCorr_01_4->AddEntry(cg_2p1_3p3_5p5,"#LTcos(2#Phi_{1}* + 3#Phi_{2}* - 5#Phi_{5}*)#GT","lp");
        legCorr_01_4->AddEntry(cg_4p1_9p3_5p5,"#LTcos(4#Phi_{1}* - 9#Phi_{2}* - 5#Phi_{5}*)#GT","lp");
        legCorr_01_4->Draw();

        if (print_plots && !minorAxes) cCorr_01->Print("plots/majorAxes/Corr_01.png","png");
        if (close_plots) cCorr_01->Close();


        //--comparison between MC and theory predictions

        TPaveText * txCosComp_00_0 = new TPaveText(0.58, 0.80, 0.92, 0.94,"NDC");
        txCosComp_00_0->SetFillColorAlpha(0,0);
        txCosComp_00_0->SetBorderSize(0);
        txCosComp_00_0->SetTextFont(43);
        txCosComp_00_0->SetTextSize(20);

        TPaveText * txCosComp_00_1 = new TPaveText(0.52, 0.80, 0.86, 0.94,"NDC");
        txCosComp_00_1->SetFillColorAlpha(0,0);
        txCosComp_00_1->SetBorderSize(0);
        txCosComp_00_1->SetTextFont(43);
        txCosComp_00_1->SetTextSize(20);

        TCanvas * cCosComp_00 = new TCanvas("cCosComp_00","cCosComp_00",1100,650);
        cCosComp_00->Divide(3,2,0,0);
        cCosComp_00->cd(1);
        TH1D * hCosComp_00_0 = (TH1D *) hdummy_corr_00->Clone("hCosComp_00_0");
        hCosComp_00_0->Draw();
        cg_6p2_p3->Draw("same p");
        cg_EPJ_Fig2_cos6Phi2_Phi3->Draw("same 3");
        TPaveText * txCosComp_00_11 = (TPaveText *) txCosComp_00_1->Clone();
        txCosComp_00_11->AddText("#LTcos6(#Phi_{2}* - #Phi_{3}*)#GT");
        txCosComp_00_11->Draw();

        cCosComp_00->cd(2);
        TH1D * hCosComp_00_1 = (TH1D *) hdummy_corr_00->Clone("hCosComp_00_1");
        hCosComp_00_1->GetXaxis()->SetLabelSize(0);
        hCosComp_00_1->GetYaxis()->SetLabelSize(0);
        hCosComp_00_1->Draw();
        cg_4p2_p4->Draw("same p");
        cg_EPJ_Fig2_cos4Phi2_Phi4->Draw("same 3");
        TPaveText * txCosComp_00_2 = (TPaveText *) txCosComp_00_1->Clone();
        txCosComp_00_2->AddText("#LTcos4(#Phi_{2}* - #Phi_{4}*)#GT");
        txCosComp_00_2->Draw();

        cCosComp_00->cd(3);
        TH1D * hCosComp_00_2 = (TH1D *) hdummy_corr_00->Clone("hCosComp_00_2");
        hCosComp_00_2->GetXaxis()->SetLabelSize(0);
        hCosComp_00_2->GetYaxis()->SetLabelSize(0);
        hCosComp_00_2->Draw();
        cg_6p2_p6->Draw("same p");
        cg_EPJ_Fig2_cos6Phi2_Phi6->Draw("same 3");
        TPaveText * txCosComp_00_3 = (TPaveText *) txCosComp_00_1->Clone();
        txCosComp_00_3->AddText("#LTcos6(#Phi_{2}* - #Phi_{6}*)#GT");
        txCosComp_00_3->Draw();

        cCosComp_00->cd(4);
        TH1D * hCosComp_00_3 = (TH1D *) hdummy_corr_00->Clone("hCosComp_00_3");
        hCosComp_00_3->Draw();
        cg_6p3_p6->Draw("same p");
        cg_EPJ_Fig2_cos6Phi3_Phi6->Draw("same 3");
        TPaveText * txCosComp_00_4 = (TPaveText *) txCosComp_00_0->Clone();
        txCosComp_00_4->AddText("#LTcos6(#Phi_{3}* - #Phi_{6}*)#GT");
        txCosComp_00_4->Draw();

        cCosComp_00->cd(5);
        TH1D * hCosComp_00_4 = (TH1D *) hdummy_corr_00->Clone("hCosComp_00_4");
        hCosComp_00_4->GetYaxis()->SetLabelSize(0);
        hCosComp_00_4->Draw();
        cg_10p2_p5->Draw("same p");
        cg_EPJ_Fig2_cos10Phi2_Phi5->Draw("same 3");
        TPaveText * txCosComp_00_5 = (TPaveText *) txCosComp_00_1->Clone();
        txCosComp_00_5->AddText("#LTcos10(#Phi_{2}* - #Phi_{5}*)#GT");
        txCosComp_00_5->Draw();

        cCosComp_00->cd(6);
        TH1D * hCosComp_00_5 = (TH1D *) hdummy_corr_00->Clone("hCosComp_00_5");
        hCosComp_00_5->GetYaxis()->SetLabelSize(0);
        hCosComp_00_5->Draw();
        cg_15p3_p5->Draw("same p");
        cg_EPJ_Fig2_cos15Phi3_Phi5->Draw("same 3");
        TPaveText * txCosComp_00_6 = (TPaveText *) txCosComp_00_1->Clone();
        txCosComp_00_6->AddText("#LTcos15(#Phi_{3}* - #Phi_{5}*)#GT");
        txCosComp_00_6->Draw();

        cCosComp_00->cd(1);
        TLegend * legCosComp_00 = new TLegend(0.43, 0.49, 0.77, 0.80);
        legCosComp_00->SetFillColorAlpha(0,0);
        legCosComp_00->SetBorderSize(0);
        legCosComp_00->SetTextFont(43);
        legCosComp_00->SetTextSize(18);
        legCosComp_00->SetHeader("Au+Au");
        legCosComp_00->AddEntry(cg_6p2_p3," Toy Glauber","lp");
        legCosComp_00->AddEntry(cg_EPJ_Fig2_cos6Phi2_Phi3," Eur.Phys.J. C73, 2510","f");
        legCosComp_00->Draw();

        if (print_plots && !minorAxes) cCosComp_00->Print("plots/majorAxes/cCosComp_00.png","png");
        if (close_plots) cCosComp_00->Close();


        //-- more comparisons between Glauber and theory predictions

        TPaveText * txCosComp_01_0 = new TPaveText(0.44, 0.86, 0.90, 0.95,"NDC");
        txCosComp_01_0->SetFillColorAlpha(0,0);
        txCosComp_01_0->SetBorderSize(0);
        txCosComp_01_0->SetTextFont(43);
        txCosComp_01_0->SetTextSize(18);

        TPaveText * txCosComp_01_1 = new TPaveText(0.35, 0.86,  0.88, 0.95,"NDC");
        txCosComp_01_1->SetFillColorAlpha(0,0);
        txCosComp_01_1->SetBorderSize(0);
        txCosComp_01_1->SetTextFont(43);
        txCosComp_01_1->SetTextSize(18);

        TCanvas * cCosComp_01 = new TCanvas("cCosComp_01","cCosComp_01",1250,450);
        cCosComp_01->Divide(5,1,0,0);
        cCosComp_01->cd(1);
        TH1D * hCosComp_01_1 = (TH1D *) hdummy_corr_00->Clone("hCosComp_01_1");
        hCosComp_01_1->GetYaxis()->SetRangeUser(-0.1, 0.65);
        hCosComp_01_1->GetYaxis()->SetLabelSize(0.06);
        hCosComp_01_1->Draw();
        cg_2p1_p2->Draw("same p");
        cg_EPJ_Fig2_cos2Phi1_Phi2->Draw("same 3");
        TPaveText * txCosComp_01_11 = (TPaveText *) txCosComp_01_0->Clone();
        txCosComp_01_11->AddText("#LTcos2(#Phi_{1}* - #Phi_{2}*)#GT");
        txCosComp_01_11->Draw();

        cCosComp_01->cd(2);
        TH1D * hCosComp_01_2 = (TH1D *) hdummy_corr_00->Clone("hCosComp_01_2");
        hCosComp_01_2->GetYaxis()->SetRangeUser(-0.1, 0.65);
        hCosComp_01_2->GetXaxis()->SetLabelSize(0.06);
        hCosComp_01_2->GetXaxis()->SetLabelOffset(0.000);
        hCosComp_01_2->GetXaxis()->SetTitleSize(0.07);
        hCosComp_01_2->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_01_2->GetYaxis()->SetLabelSize(0);
        hCosComp_01_2->Draw();
        cg_3p1_p3->Draw("same p");
        cg_EPJ_Fig2_cos3Phi1_Phi3->Draw("same 3");
        TPaveText * txCosComp_01_2 = (TPaveText *) txCosComp_01_1->Clone();
        txCosComp_01_2->AddText("#LTcos3(#Phi_{1}* - #Phi_{3}*)#GT");
        txCosComp_01_2->Draw();

        cCosComp_01->cd(3);
        TH1D * hCosComp_01_3 = (TH1D *) hdummy_corr_00->Clone("hCosComp_01_3");
        hCosComp_01_3->GetYaxis()->SetRangeUser(-0.1, 0.65);
        hCosComp_01_3->GetXaxis()->SetLabelSize(0.06);
        hCosComp_01_3->GetXaxis()->SetLabelOffset(0.000);
        hCosComp_01_3->GetXaxis()->SetTitleSize(0.07);
        hCosComp_01_3->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_01_3->GetYaxis()->SetLabelSize(0);
        hCosComp_01_3->Draw();
        cg_4p1_p4->Draw("same p");
        cg_EPJ_Fig2_cos4Phi1_Phi4->Draw("same 3");
        TPaveText * txCosComp_01_3 = (TPaveText *) txCosComp_01_1->Clone();
        txCosComp_01_3->AddText("#LTcos4(#Phi_{1}* - #Phi_{4}*)#GT");
        txCosComp_01_3->Draw();

        cCosComp_01->cd(4);
        TH1D * hCosComp_01_4 = (TH1D *) hdummy_corr_00->Clone("hCosComp_01_4");
        hCosComp_01_4->GetYaxis()->SetRangeUser(-0.1, 0.65);
        hCosComp_01_4->GetYaxis()->SetLabelSize(0);
        hCosComp_01_4->GetXaxis()->SetLabelSize(0.06);
        hCosComp_01_4->GetXaxis()->SetLabelOffset(0.000);
        hCosComp_01_4->GetXaxis()->SetTitleSize(0.07);
        hCosComp_01_4->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_01_4->Draw();
        cg_5p1_p5->Draw("same p");
        cg_EPJ_Fig2_cos5Phi1_Phi5->Draw("same 3");
        TPaveText * txCosComp_01_4 = (TPaveText *) txCosComp_01_1->Clone();
        txCosComp_01_4->AddText("#LTcos5(#Phi_{1}* - #Phi_{5}*)#GT");
        txCosComp_01_4->Draw();

        cCosComp_01->cd(5);
        TH1D * hCosComp_01_5 = (TH1D *) hdummy_corr_00->Clone("hCosComp_01_5");
        hCosComp_01_5->GetYaxis()->SetRangeUser(-0.1, 0.65);
        hCosComp_01_5->GetYaxis()->SetLabelSize(0);
        hCosComp_01_5->GetXaxis()->SetLabelSize(0.06);
        hCosComp_01_5->GetXaxis()->SetLabelOffset(0.000);
        hCosComp_01_5->GetXaxis()->SetTitleSize(0.07);
        hCosComp_01_5->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_01_5->Draw();
        cg_6p1_p6->Draw("same p");
        cg_EPJ_Fig2_cos6Phi1_Phi6->Draw("same 3");
        TPaveText * txCosComp_01_5 = (TPaveText *) txCosComp_01_1->Clone();
        txCosComp_01_5->AddText("#LTcos6(#Phi_{1}* - #Phi_{6}*)#GT");
        txCosComp_01_5->Draw();

        cCosComp_01->cd(1);
        TLegend * legCosComp_01 = new TLegend(0.33, 0.69, 0.67, 0.87);
        legCosComp_01->SetFillColorAlpha(0,0);
        legCosComp_01->SetBorderSize(0);
        legCosComp_01->SetTextFont(43);
        legCosComp_01->SetTextSize(16);
        legCosComp_01->SetHeader("Au+Au");
        legCosComp_01->AddEntry(cg_2p1_p2," Toy Glauber","lp");
        legCosComp_01->AddEntry(cg_EPJ_Fig2_cos2Phi1_Phi2," Eur.Phys.J. C73, 2510","f");
        legCosComp_01->Draw();

        if (print_plots && !minorAxes) cCosComp_01->Print("plots/majorAxes/cCosComp_01.png","png");
        if (close_plots) cCosComp_01->Close();


        //-- even more comparisons between Glauber and theory predictions

        TPaveText * txCosComp_02_dum0 = new TPaveText(0.43, 0.85, 0.90, 0.94,"NDC");
        txCosComp_02_dum0->SetFillColorAlpha(0,0);
        txCosComp_02_dum0->SetBorderSize(0);
        txCosComp_02_dum0->SetTextFont(43);
        txCosComp_02_dum0->SetTextSize(18);

        TPaveText * txCosComp_02_dum1 = new TPaveText(0.38, 0.85, 0.84, 0.94,"NDC");
        txCosComp_02_dum1->SetFillColorAlpha(0,0);
        txCosComp_02_dum1->SetBorderSize(0);
        txCosComp_02_dum1->SetTextFont(43);
        txCosComp_02_dum1->SetTextSize(18);

        TCanvas * cCosComp_02 = new TCanvas("cCosComp_02","cCosComp_02",1200,650);
        cCosComp_02->Divide(4,2,0,0);
        cCosComp_02->cd(1);
        TH1D * hCosComp_02_1 = (TH1D *) hdummy_corr_00->Clone("hCosComp_02_1");
        hCosComp_02_1->GetYaxis()->SetLabelSize(0.06);
        hCosComp_02_1->Draw();
        cg_1p1_2p2_3p3->Draw("same p");
        cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3->Draw("same 3");
        TPaveText * txCosComp_02_01 = (TPaveText *) txCosComp_02_dum0->Clone();
        txCosComp_02_01->AddText("#LTcos(#Phi_{1}* + 2#Phi_{2}* - 3#Phi_{3}*)#GT");
        txCosComp_02_01->Draw();

        cCosComp_02->cd(2);
        TH1D * hCosComp_02_2 = (TH1D *) hdummy_corr_00->Clone("hCosComp_02_2");
        hCosComp_02_2->Draw();
        cg_4p1_2p2_6p3->Draw("same p");
        cg_EPJ_Fig12_cos4Phi1_2Phi2_6Phi3->Draw("same 3");
        TPaveText * txCosComp_02_02 = (TPaveText *) txCosComp_02_dum1->Clone();
        txCosComp_02_02->AddText("#LTcos(4#Phi_{1}* + 2#Phi_{2}* - 6#Phi_{3}*)#GT");
        txCosComp_02_02->Draw();

        cCosComp_02->cd(3);
        TH1D * hCosComp_02_3 = (TH1D *) hdummy_corr_00->Clone("hCosComp_02_3");
        hCosComp_02_3->Draw();
        cg_2p1_4p2_6p3->Draw("same p");
        cg_EPJ_Fig12_cos2Phi1_4Phi2_6Phi3->Draw("same 3");
        TPaveText * txCosComp_02_03 = (TPaveText *) txCosComp_02_dum1->Clone();
        txCosComp_02_03->AddText("#LTcos(2#Phi_{1}* + 4#Phi_{2}* - 6#Phi_{3}*)#GT");
        txCosComp_02_03->Draw();

        cCosComp_02->cd(4);
        TH1D * hCosComp_02_4 = (TH1D *) hdummy_corr_00->Clone("hCosComp_02_4");
        hCosComp_02_4->Draw();
        cg_5p1_2p2_3p3->Draw("same p");
        cg_EPJ_Fig12_cos5Phi1_2Phi2_3Phi3->Draw("same 3");
        TPaveText * txCosComp_02_04 = (TPaveText *) txCosComp_02_dum1->Clone();
        txCosComp_02_04->AddText("#LTcos(4#Phi_{1}* - 2#Phi_{2}* - 3#Phi_{3}*)#GT");
        txCosComp_02_04->Draw();

        cCosComp_02->cd(5);
        TH1D * hCosComp_02_5 = (TH1D *) hdummy_corr_00->Clone("hCosComp_02_5");
        hCosComp_02_5->Draw();
        cg_2p1_2p2_4p4->SetMarkerStyle(20);
        cg_2p1_2p2_4p4->Draw("same p");
        cg_EPJ_Fig12_cos2Phi1_2Phi2_4Phi4->Draw("same 3");
        TPaveText * txCosComp_02_05 = (TPaveText *) txCosComp_02_dum0->Clone();
        txCosComp_02_05->AddText("#LTcos(2#Phi_{1}* + 2#Phi_{2}* - 4#Phi_{4}*)#GT");
        txCosComp_02_05->Draw();

        cCosComp_02->cd(6);
        TH1D * hCosComp_02_6 = (TH1D *) hdummy_corr_00->Clone("hCosComp_02_6");
        hCosComp_02_6->GetXaxis()->SetTitleSize(0.07);
        hCosComp_02_6->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_02_6->Draw();
        cg_2p1_6p2_4p4->SetMarkerStyle(29);
        cg_2p1_6p2_4p4->Draw("same p");
        cg_EPJ_Fig12_cos2Phi1_6Phi2_4Phi4->Draw("same 3");
        TPaveText * txCosComp_02_06 = (TPaveText *) txCosComp_02_dum1->Clone();
        txCosComp_02_06->AddText("#LTcos(2#Phi_{1}* - 6#Phi_{2}* + 4#Phi_{4}*)#GT");
        txCosComp_02_06->Draw();

        cCosComp_02->cd(7);
        TH1D * hCosComp_02_7 = (TH1D *) hdummy_corr_00->Clone("hCosComp_02_7");
        hCosComp_02_7->GetXaxis()->SetTitleSize(0.07);
        hCosComp_02_7->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_02_7->Draw();
        cg_2p1_6p2_8p4->SetMarkerStyle(34);
        cg_2p1_6p2_8p4->Draw("same p");
        cg_EPJ_Fig12_cos2Phi1_6Phi2_8Phi4->Draw("same 3");
        TPaveText * txCosComp_02_07 = (TPaveText *) txCosComp_02_dum1->Clone();
        txCosComp_02_07->AddText("#LTcos(2#Phi_{1}* + 6#Phi_{2}* - 8#Phi_{4}*)#GT");
        txCosComp_02_07->Draw();

        cCosComp_02->cd(8);
        TLegend * legCosComp_02 = new TLegend(0.12, 0.44, 0.57, 0.71);
        legCosComp_02->SetFillColorAlpha(0,0);
        legCosComp_02->SetBorderSize(0);
        legCosComp_02->SetTextFont(43);
        legCosComp_02->SetTextSize(18);
        legCosComp_02->SetHeader("Au+Au");
        legCosComp_02->AddEntry(cg_1p1_2p2_3p3," Toy Glauber","lp");
        legCosComp_02->AddEntry(cg_EPJ_Fig12_cos1Phi1_2Phi2_3Phi3," Eur.Phys.J. C73, 2510","f");
        legCosComp_02->Draw();

        if (print_plots && !minorAxes) cCosComp_02->Print("plots/majorAxes/cCosComp_02.png","png");
        if (close_plots) cCosComp_02->Close();



        //-- still more comparisons between Glauber and theory predictions

        TPaveText * txCosComp_03_dum0 = new TPaveText(0.43, 0.85, 0.90, 0.94,"NDC");
        txCosComp_03_dum0->SetFillColorAlpha(0,0);
        txCosComp_03_dum0->SetBorderSize(0);
        txCosComp_03_dum0->SetTextFont(43);
        txCosComp_03_dum0->SetTextSize(18);

        TPaveText * txCosComp_03_dum1 = new TPaveText(0.38, 0.85, 0.84, 0.94,"NDC");
        txCosComp_03_dum1->SetFillColorAlpha(0,0);
        txCosComp_03_dum1->SetBorderSize(0);
        txCosComp_03_dum1->SetTextFont(43);
        txCosComp_03_dum1->SetTextSize(18);

        TCanvas * cCosComp_03 = new TCanvas("cCosComp_03","cCosComp_03",1200,650);
        cCosComp_03->Divide(4,2,0,0);
        cCosComp_03->cd(1);
        TH1D * hCosComp_03_1 = (TH1D *) hdummy_corr_00->Clone("hCosComp_03_1");
        hCosComp_03_1->GetYaxis()->SetLabelSize(0.06);
        hCosComp_03_1->Draw();
        cg_1p1_3p3_4p4->Draw("same p");
        cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4->Draw("same 3");
        TPaveText * txCosComp_03_01 = (TPaveText *) txCosComp_03_dum0->Clone();
        txCosComp_03_01->AddText("#LTcos(#Phi_{1}* + 3#Phi_{3}* - 4#Phi_{4}*)#GT");
        txCosComp_03_01->Draw();

        cCosComp_03->cd(2);
        TH1D * hCosComp_03_2 = (TH1D *) hdummy_corr_00->Clone("hCosComp_03_2");
        hCosComp_03_2->Draw();
        cg_2p1_6p3_4p4->Draw("same p");
        cg_EPJ_Fig12_cos2Phi1_6Phi3_4Phi4->Draw("same 3");
        TPaveText * txCosComp_03_02 = (TPaveText *) txCosComp_03_dum1->Clone();
        txCosComp_03_02->AddText("#LTcos(2#Phi_{1}* - 6#Phi_{3}* + 4#Phi_{4}*)#GT");
        txCosComp_03_02->Draw();

        cCosComp_03->cd(3);
        TH1D * hCosComp_03_3 = (TH1D *) hdummy_corr_00->Clone("hCosComp_03_3");
        hCosComp_03_3->Draw();
        cg_2p1_4p2_5p5->Draw("same p");
        cg_EPJ_Fig12_cos1Phi1_4Phi2_5Phi5->Draw("same 3");
        TPaveText * txCosComp_03_03 = (TPaveText *) txCosComp_03_dum1->Clone();
        txCosComp_03_03->AddText("#LTcos(#Phi_{1}* + 4#Phi_{2}* - 5#Phi_{5}*)#GT");
        txCosComp_03_03->Draw();

        cCosComp_03->cd(4);
        TH1D * hCosComp_03_4 = (TH1D *) hdummy_corr_00->Clone("hCosComp_03_4");
        hCosComp_03_4->Draw();
        cg_3p1_2p2_5p5->Draw("same p");
        cg_EPJ_Fig12_cos3Phi1_2Phi2_5Phi5->Draw("same 3");
        TPaveText * txCosComp_03_04 = (TPaveText *) txCosComp_03_dum1->Clone();
        txCosComp_03_04->AddText("#LTcos(3#Phi_{1}* + 2#Phi_{2}* - 5#Phi_{5}*)#GT");
        txCosComp_03_04->Draw();

        cCosComp_03->cd(5);
        TH1D * hCosComp_03_5 = (TH1D *) hdummy_corr_00->Clone("hCosComp_03_5");
        hCosComp_03_5->Draw();
        cg_1p1_6p3_5p5->SetMarkerStyle(20);
        cg_1p1_6p3_5p5->Draw("same p");
        cg_EPJ_Fig12_cos1Phi1_6Phi3_5Phi5->Draw("same 3");
        TPaveText * txCosComp_03_05 = (TPaveText *) txCosComp_03_dum0->Clone();
        txCosComp_03_05->AddText("#LTcos(#Phi_{1}* - 6#Phi_{3}* + 5#Phi_{5}*)#GT");
        txCosComp_03_05->Draw();

        cCosComp_03->cd(6);
        TH1D * hCosComp_03_6 = (TH1D *) hdummy_corr_00->Clone("hCosComp_03_6");
        hCosComp_03_6->GetXaxis()->SetTitleSize(0.07);
        hCosComp_03_6->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_03_6->Draw();
        cg_2p1_3p3_5p5->SetMarkerStyle(29);
        cg_2p1_3p3_5p5->Draw("same p");
        cg_EPJ_Fig12_cos2Phi1_3Phi3_5Phi5->Draw("same 3");
        TPaveText * txCosComp_03_06 = (TPaveText *) txCosComp_03_dum1->Clone();
        txCosComp_03_06->AddText("#LTcos(2#Phi_{1}* + 3#Phi_{3}* - 5#Phi_{5}*)#GT");
        txCosComp_03_06->Draw();

        cCosComp_03->cd(7);
        TH1D * hCosComp_03_7 = (TH1D *) hdummy_corr_00->Clone("hCosComp_03_7");
        hCosComp_03_7->GetXaxis()->SetTitleSize(0.07);
        hCosComp_03_7->GetXaxis()->SetTitleOffset(0.97);
        hCosComp_03_7->Draw();
        cg_4p1_9p3_5p5->SetMarkerStyle(34);
        cg_4p1_9p3_5p5->Draw("same p");
        cg_EPJ_Fig12_cos4Phi1_9Phi3_5Phi5->Draw("same 3");
        TPaveText * txCosComp_03_07 = (TPaveText *) txCosComp_03_dum1->Clone();
        txCosComp_03_07->AddText("#LTcos(4#Phi_{1}* - 9#Phi_{3}* + 5#Phi_{5}*)#GT");
        txCosComp_03_07->Draw();

        cCosComp_03->cd(8);
        TLegend * legCosComp_03 = new TLegend(0.12, 0.44, 0.57, 0.71);
        legCosComp_03->SetFillColor(0);
        legCosComp_03->SetBorderSize(0);
        legCosComp_03->SetTextFont(43);
        legCosComp_03->SetTextSize(18);
        legCosComp_03->SetHeader("Au+Au");
        legCosComp_03->AddEntry(cg_1p1_3p3_4p4," Toy Glauber","lp");
        legCosComp_03->AddEntry(cg_EPJ_Fig12_cos1Phi1_3Phi3_4Phi4," Eur.Phys.J. C73, 2510","f");
        legCosComp_03->Draw();

        if (print_plots && !minorAxes) cCosComp_03->Print("plots/majorAxes/cCosComp_03.png","png");
        if (close_plots) cCosComp_03->Close();

    }

}
