# include "TArrow.h"
# include "TCanvas.h"
# include "TDirectory.h"
# include "TEllipse.h"
# include "TF1.h"
# include "TFile.h"
# include "TGraphErrors.h"
# include "TH1.h"
# include "TH2.h"
# include "TH3.h"
# include "TLatex.h"
# include "TLegend.h"
# include "TLine.h"
# include "TMath.h"
# include "TPaveText.h"
# include "TRandom3.h"
# include "TStopwatch.h"
# include "TString.h"
# include "TStyle.h"
# include <fstream>
# include <iostream>
# include <stdio.h>

TRandom3 * ran;
TF1 * ranR;
TF1 * ranTheta;
TF1 * bconv = new TF1("bconv", "pol7", 0, 1);

TH1D * hws;
TH1D * hcent;
TH1D * hpsi[7][ncentbins];
TH2D * hA1;
TH2D * hA2;
TH2D * hOverlap;
TH3D * hA1_3D;
TH3D * hA2_3D;
TH3D * hCollision_3D;

// moment eccentricities
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

// naming convention: cg_k_pn_pm --> cos(k(Psi_{n} - Psi_{m}))
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
TGraphErrors * cg_1p1_4p2_5p5;
TGraphErrors * cg_3p1_2p2_5p5;
TGraphErrors * cg_1p1_6p3_5p5;
TGraphErrors * cg_2p1_3p3_5p5;
TGraphErrors * cg_4p1_9p3_5p5;

TFile * tfout;
int colors[7] = {kBlack, kBlue, kRed, kGreen+3, kMagenta, kCyan+2, kOrange+7};
TH2D * xy0[ncentbins];
TH2D * xy1[ncentbins];
TH2D * xy2[ncentbins];
TH2D * xy3[ncentbins];
TH2D * xy4[ncentbins];
TH2D * xy5[ncentbins];
TH2D * xy6[ncentbins];

TH2D * Psi2Psi1[ncentbins];
TH2D * Psi3Psi1[ncentbins];
TH2D * Psi3Psi2[ncentbins];
TH2D * Psi4Psi2[ncentbins];

TH1D * CosPsi2Psi1[ncentbins];
TH1D * CosPsi3Psi1[ncentbins];
TH1D * CosPsi3Psi2[ncentbins];
TH1D * CosPsi4Psi2[ncentbins];

TH2D * Psi2Psi1Diff[ncentbins];
TH2D * Psi3Psi1Diff[ncentbins];
TH2D * Psi5Psi1Diff[ncentbins];
TH2D * Psi3Psi2Diff[ncentbins];
TH2D * Psi4Psi2Diff[ncentbins];

TH1D * eccAng[7][ncentbins];
TH1D * eccDist[7][ncentbins];

double avNpart[ncentbins];
double avNpartCnt[ncentbins];
double avb[ncentbins];
double avbCnt[ncentbins];
double ecc[7][ncentbins];
double eccCnt[7][ncentbins];
double S[7][ncentbins];
double ecc_pow2[7][ncentbins];
double S_pow2[7][ncentbins];
double ecc_err[7][ncentbins];
double S_err[7][ncentbins];

double ecc_2_1[ncentbins];
double ecc_3_1[ncentbins];
double ecc_3_2[ncentbins];
double ecc_4_2[ncentbins];
double ecc_2_1_pow2[ncentbins];
double ecc_3_1_pow2[ncentbins];
double ecc_3_2_pow2[ncentbins];
double ecc_4_2_pow2[ncentbins];
double ecc_2_1_err[ncentbins];
double ecc_3_1_err[ncentbins];
double ecc_3_2_err[ncentbins];
double ecc_4_2_err[ncentbins];

double ecc_2_1_cnt[ncentbins];
double ecc_3_1_cnt[ncentbins];
double ecc_3_2_cnt[ncentbins];
double ecc_4_2_cnt[ncentbins];

double S_2_1[ncentbins];
double S_3_1[ncentbins];
double S_3_2[ncentbins];
double S_4_2[ncentbins];
double S_2_1_pow2[ncentbins];
double S_3_1_pow2[ncentbins];
double S_3_2_pow2[ncentbins];
double S_4_2_pow2[ncentbins];
double S_2_1_err[ncentbins];
double S_3_1_err[ncentbins];
double S_3_2_err[ncentbins];
double S_4_2_err[ncentbins];

double c_4_p4_p2[ncentbins];
double c_8_p4_p2[ncentbins];
double c_12_p4_p2[ncentbins];
double c_6_p3_p2[ncentbins];
double c_6_p2_p6[ncentbins];
double c_6_p3_p6[ncentbins];
double c_12_p3_p4[ncentbins];
double c_10_p2_p5[ncentbins];
double c_2p2_3p3_5p5[ncentbins];
double c_2p2_4p4_6p6[ncentbins];
double c_2p2_6p3_4p4[ncentbins];
double c_8p2_3p3_5p5[ncentbins];
double c_10p2_4p4_6p6[ncentbins];
double c_10p2_6p3_4p4[ncentbins];
double c_4_p4_p2_pow2[ncentbins];
double c_8_p4_p2_pow2[ncentbins];
double c_12_p4_p2_pow2[ncentbins];
double c_6_p3_p2_pow2[ncentbins];
double c_6_p2_p6_pow2[ncentbins];
double c_6_p3_p6_pow2[ncentbins];
double c_12_p3_p4_pow2[ncentbins];
double c_10_p2_p5_pow2[ncentbins];
double c_2p2_3p3_5p5_pow2[ncentbins];
double c_2p2_4p4_6p6_pow2[ncentbins];
double c_2p2_6p3_4p4_pow2[ncentbins];
double c_8p2_3p3_5p5_pow2[ncentbins];
double c_10p2_4p4_6p6_pow2[ncentbins];
double c_10p2_6p3_4p4_pow2[ncentbins];
double c_4_p4_p2_err[ncentbins];
double c_8_p4_p2_err[ncentbins];
double c_12_p4_p2_err[ncentbins];
double c_6_p3_p2_err[ncentbins];
double c_6_p2_p6_err[ncentbins];
double c_6_p3_p6_err[ncentbins];
double c_12_p3_p4_err[ncentbins];
double c_10_p2_p5_err[ncentbins];
double c_2p2_3p3_5p5_err[ncentbins];
double c_2p2_4p4_6p6_err[ncentbins];
double c_2p2_6p3_4p4_err[ncentbins];
double c_8p2_3p3_5p5_err[ncentbins];
double c_10p2_4p4_6p6_err[ncentbins];
double c_10p2_6p3_4p4_err[ncentbins];

double c_6p2_p3[ncentbins];
double c_4p2_p4[ncentbins];
double c_6p2_p6[ncentbins];
double c_6p3_p6[ncentbins];
double c_10p2_p5[ncentbins];
double c_15p3_p5[ncentbins];
double c_2p1_p2[ncentbins];
double c_3p1_p3[ncentbins];
double c_4p1_p4[ncentbins];
double c_5p1_p5[ncentbins];
double c_6p1_p6[ncentbins];
double c_1p1_2p2_3p3[ncentbins];
double c_4p1_2p2_6p3[ncentbins];
double c_2p1_4p2_6p3[ncentbins];
double c_5p1_2p2_3p3[ncentbins];
double c_2p1_2p2_4p4[ncentbins];
double c_2p1_6p2_4p4[ncentbins];
double c_2p1_6p2_8p4[ncentbins];
double c_1p1_3p3_4p4[ncentbins];
double c_2p1_6p3_4p4[ncentbins];
double c_2p1_4p2_5p5[ncentbins];
double c_3p1_2p2_5p5[ncentbins];
double c_1p1_6p3_5p5[ncentbins];
double c_2p1_3p3_5p5[ncentbins];
double c_4p1_9p3_5p5[ncentbins];
double c_6p2_p3_pow2[ncentbins];
double c_4p2_p4_pow2[ncentbins];
double c_6p2_p6_pow2[ncentbins];
double c_6p3_p6_pow2[ncentbins];
double c_10p2_p5_pow2[ncentbins];
double c_15p3_p5_pow2[ncentbins];
double c_2p1_p2_pow2[ncentbins];
double c_3p1_p3_pow2[ncentbins];
double c_4p1_p4_pow2[ncentbins];
double c_5p1_p5_pow2[ncentbins];
double c_6p1_p6_pow2[ncentbins];
double c_1p1_2p2_3p3_pow2[ncentbins];
double c_4p1_2p2_6p3_pow2[ncentbins];
double c_2p1_4p2_6p3_pow2[ncentbins];
double c_5p1_2p2_3p3_pow2[ncentbins];
double c_2p1_2p2_4p4_pow2[ncentbins];
double c_2p1_6p2_4p4_pow2[ncentbins];
double c_2p1_6p2_8p4_pow2[ncentbins];
double c_1p1_3p3_4p4_pow2[ncentbins];
double c_2p1_6p3_4p4_pow2[ncentbins];
double c_2p1_4p2_5p5_pow2[ncentbins];
double c_3p1_2p2_5p5_pow2[ncentbins];
double c_1p1_6p3_5p5_pow2[ncentbins];
double c_2p1_3p3_5p5_pow2[ncentbins];
double c_4p1_9p3_5p5_pow2[ncentbins];
double c_6p2_p3_err[ncentbins];
double c_4p2_p4_err[ncentbins];
double c_6p2_p6_err[ncentbins];
double c_6p3_p6_err[ncentbins];
double c_10p2_p5_err[ncentbins];
double c_15p3_p5_err[ncentbins];
double c_2p1_p2_err[ncentbins];
double c_3p1_p3_err[ncentbins];
double c_4p1_p4_err[ncentbins];
double c_5p1_p5_err[ncentbins];
double c_6p1_p6_err[ncentbins];
double c_1p1_2p2_3p3_err[ncentbins];
double c_4p1_2p2_6p3_err[ncentbins];
double c_2p1_4p2_6p3_err[ncentbins];
double c_5p1_2p2_3p3_err[ncentbins];
double c_2p1_2p2_4p4_err[ncentbins];
double c_2p1_6p2_4p4_err[ncentbins];
double c_2p1_6p2_8p4_err[ncentbins];
double c_1p1_3p3_4p4_err[ncentbins];
double c_2p1_6p3_4p4_err[ncentbins];
double c_2p1_4p2_5p5_err[ncentbins];
double c_3p1_2p2_5p5_err[ncentbins];
double c_1p1_6p3_5p5_err[ncentbins];
double c_2p1_3p3_5p5_err[ncentbins];
double c_4p1_9p3_5p5_err[ncentbins];

double c_cnt[ncentbins];

Double_t bounds(Double_t ord, Double_t ang) {
    while (ang >  TMath::Pi()/ord) ang-=TMath::TwoPi()/ord;
    while (ang < -TMath::Pi()/ord) ang+=TMath::TwoPi()/ord;
    return ang;
}

void wsFunc(int A, double * x, double * y, double * z){
    int i = 0;
    while (i<A) {
        double Rloc  = ranR->GetRandom();
        double th = ranTheta->GetRandom();
        double ph = ran->Uniform(0,2.*TMath::Pi());
        x[i] = Rloc*TMath::Sin(th)*TMath::Cos(ph);
        y[i] = Rloc*TMath::Sin(th)*TMath::Sin(ph);
        z[i] = Rloc*TMath::Cos(th);
        int sub = 0;
        for (int j = 0; j<i; j++) {
            if (sqrt( pow(x[j]-x[i],2) + pow(y[j]-y[i],2) + pow(z[j]-z[i],2) )<rHard) sub = 1;
        }
        if (sub>0) {
            // relocate 1
            th = ranTheta->GetRandom();
            ph = ran->Uniform(0,2.*TMath::Pi());
            x[i] = Rloc*TMath::Sin(th)*TMath::Cos(ph);
            y[i] = Rloc*TMath::Sin(th)*TMath::Sin(ph);
            z[i] = Rloc*TMath::Cos(th);
        }
        sub = 0;
        for (int j = 0; j<i; j++) {
            if (sqrt( pow(x[j]-x[i],2) + pow(y[j]-y[i],2) + pow(z[j]-z[i],2) )<rHard) sub = 1;
        }
        if (sub>0) {
            // relocate 2
            th = ranTheta->GetRandom();
            ph = ran->Uniform(0,2.*TMath::Pi());
            x[i] = Rloc*TMath::Sin(th)*TMath::Cos(ph);
            y[i] = Rloc*TMath::Sin(th)*TMath::Sin(ph);
            z[i] = Rloc*TMath::Cos(th);
        }
        sub = 0;
        for (int j = 0; j<i; j++) {
            if (sqrt( pow(x[j]-x[i],2) + pow(y[j]-y[i],2) + pow(z[j]-z[i],2) )<rHard) sub = 1;
        }
        if (sub>0) {
            // relocate 3
            th = ranTheta->GetRandom();
            ph = ran->Uniform(0,2.*TMath::Pi());
            x[i] = Rloc*TMath::Sin(th)*TMath::Cos(ph);
            y[i] = Rloc*TMath::Sin(th)*TMath::Sin(ph);
            z[i] = Rloc*TMath::Cos(th);
        }
        sub = 0;
        for (int j = 0; j<i; j++) {
            if (sqrt( pow(x[j]-x[i],2) + pow(y[j]-y[i],2) + pow(z[j]-z[i],2) )<rHard) sub = 1;
        }
        if (sub>0) {
            // relocate 4
            th = ranTheta->GetRandom();
            ph = ran->Uniform(0,2.*TMath::Pi());
            x[i] = Rloc*TMath::Sin(th)*TMath::Cos(ph);
            y[i] = Rloc*TMath::Sin(th)*TMath::Sin(ph);
            z[i] = Rloc*TMath::Cos(th);
        }
        sub = 0;
        for(int j = 0; j<i; j++) {
            if (sqrt( pow(x[j]-x[i],2) + pow(y[j]-y[i],2) + pow(z[j]-z[i],2) )<rHard) sub = 1;
        }
        if (sub>0) {
            // relocate 5
            th = ranTheta->GetRandom();
            ph = ran->Uniform(0,2.*TMath::Pi());
            x[i] = Rloc*TMath::Sin(th)*TMath::Cos(ph);
            y[i] = Rloc*TMath::Sin(th)*TMath::Sin(ph);
            z[i] = Rloc*TMath::Cos(th);
        }
        sub = 0;
        for (int j = 0; j<i; j++) {
            if (sqrt( pow(x[j]-x[i],2) + pow(y[j]-y[i],2) + pow(z[j]-z[i],2) )<rHard) sub = 1;
        }
        if (sub>0) {
            // relocate 6
            th = ranTheta->GetRandom();
            ph = ran->Uniform(0,2.*TMath::Pi());
            x[i] = Rloc*TMath::Sin(th)*TMath::Cos(ph);
            y[i] = Rloc*TMath::Sin(th)*TMath::Sin(ph);
            z[i] = Rloc*TMath::Cos(th);
        }
        sub = 0;
        for (int j = 0; j<i; j++) {
            if (sqrt( pow(x[j]-x[i],2) + pow(y[j]-y[i],2) + pow(z[j]-z[i],2) )<rHard) sub = 1;
        }

        i-=sub;
        i++;
    }
    return;
}

void collide(double * x1, double * y1, double * z1, double * x2, double * y2, double * z2, double * x, double * y, double * z, int &npart){
    double rmax = sqrt(sign/10./TMath::Pi());
    int hit1[500];
    int hit2[500];
    for (int i = 0; i<A1; i++) {hit1[i]=0;}
    for (int i = 0; i<A2; i++) {hit2[i]=0;}
    for (int i = 0; i<A1; i++) {
        for (int j = 0; j<A2; j++) {
            double d = sqrt( pow(x1[i]-x2[j],2) + pow(y1[i]-y2[j],2) );
            if (d<rmax) {
                hit1[i] = 1;
                hit2[j] = 1;
            }
        }
    }
    npart = 0;
    for (int i = 0; i<A1; i++) {
        if (hit1[i] == 1) {
            x[npart] = x1[i];
            y[npart] = y1[i];
            z[npart++] = z1[i];
        }
    }
    for (int i = 0; i<A2; i++) {
        if (hit2[i] == 1) {
            x[npart] = x2[i];
            y[npart] = y2[i];
            z[npart++] = z2[i];
        }
    }
}

double PsiN(double order, int np, double * x, double * y, double &eccang) {
    double sx = 0;
    double sy = 0;
    Double_t rs = order;
    if (order == 1) rs = 3;
    for (int i = 0; i<np; i++) {
        double r = sqrt( x[i]*x[i]+y[i]*y[i] );
        // double rn = pow(r,rs);
        double rn = pow(r,2); // for r^2 weighting
        double ang = TMath::ATan2(y[i],x[i]);
        sx += rn*TMath::Cos(order*ang);
        sy += rn*TMath::Sin(order*ang);
    }
    if (sx == 0 && sy == 0) return -10;
    double angEP;
    eccang = bounds(order,(1/order)*TMath::ATan2(sy,sx));
    angEP =  bounds(order,(1/order)*(TMath::ATan2(sy,sx)-TMath::Pi()));
    return angEP;
}

double eccentricity(double order, double psi, int npart, double * x, double * y, double &Sval) {
    double avx = 0;
    double avy = 0;
    double avx1 = 0;
    double avy1 = 0;
    double avx2 = 0;
    double avy2 = 0;
    double avr = 0;

    double rs = order;
    if (order == 1) rs = 3;
    for (int i = 0; i<npart; i++) {
        double r = sqrt( pow(x[i],2) + pow(y[i],2) );
        // double rn = pow(r, rs);
        double rn = pow(r, 2); // for r^2 weighting
        double ang = TMath::ATan2(y[i],x[i]);
        avx += rn*TMath::Cos(order*(ang-psi));
        avy += rn*TMath::Sin(order*(ang-psi));
        avx1 += r*TMath::Cos(order*(ang-psi));
        avy1 += r*TMath::Sin(order*(ang-psi));
        avx2 += pow(r*TMath::Cos(order*(ang-psi)),2);
        avy2 += pow(r*TMath::Sin(order*(ang-psi)),2);
        avr += rn;
    }

    avx/=(double)npart;
    avy/=(double)npart;
    avx1/=(double)npart;
    avy1/=(double)npart;
    avx2/=(double)npart;
    avy2/=(double)npart;
    avr/=(double)npart;

    double sigx = sqrt( avx2 - pow(avx1,2) );
    double sigy = sqrt( avy2 - pow(avy1,2) );
    Sval = TMath::Pi()*sigx*sigy;

    return avx/avr;
}

void recenter(int npart, double * x, double * y, double * z) {
    double xav = 0;
    double yav = 0;
    double zav = 0;
    for (int i = 0; i<npart; i++) {
        xav+=x[i];
        yav+=y[i];
        zav+=z[i];
    }
    xav/=npart;
    yav/=npart;
    zav/=npart;
    for (int i = 0; i<npart; i++) {
        x[i]-=xav;
        y[i]-=yav;
        z[i]-=zav;
    }
}

void PlotColl( double * xb1, double * y1, double * z1, double * xb2, double * y2, double * z2, double * x, double * y, double * z, int &npart, int cbin, double b )
{
    hOverlap = new TH2D(Form("Overlap_c%d",cbin), "Overlap", 400, -3*R, 3*R, 400, -3*R, 3*R);
    hA1 = new TH2D(Form("A1_c%d",cbin), "A1", 400, -3*R, 3*R, 400, -3*R, 3*R);
    hA2 = new TH2D(Form("A2_c%d",cbin), "A2", 400, -3*R, 3*R, 400, -3*R, 3*R);
    hA1_3D = new TH3D(Form("A1_3D_c%d",cbin), "A1_3D", 100, -3*R, 3*R, 100, -3*R, 3*R, 100, -3*R, 3*R);
    hA2_3D = new TH3D(Form("A2_3D_c%d",cbin), "A2_3D", 100, -3*R, 3*R, 100, -3*R, 3*R, 100, -3*R, 3*R);
    hCollision_3D = new TH3D(Form("Collision_3D_c%d",cbin), "Collision_3D", 100, -3*R, 3*R, 100, -3*R, 3*R, 100, -3*R, 3*R);
    hA1_3D->SetMarkerColor(kBlue);
    hA2_3D->SetMarkerColor(kRed);
    hCollision_3D->SetMarkerColor(kGreen);

    TEllipse * eA1[500];
    TEllipse * eA2[500];
    TEllipse * eOverlap[500];

    Double_t rnuc = sqrt( sign/10./TMath::Pi() )/2.;
    double scale = 50;
    double eccang1 = 0;
    double eccang2 = 0;
    double eccang3 = 0;
    double eccang4 = 0;
    double eccang5 = 0;
    double eccang6 = 0;
    double ph1 = PsiN(1., npart, x, y, eccang1);
    double ph2 = PsiN(2., npart, x, y, eccang2);
    double ph3 = PsiN(3., npart, x, y, eccang3);
    double ph4 = PsiN(4., npart, x, y, eccang4);
    double ph5 = PsiN(5., npart, x, y, eccang5);
    double ph6 = PsiN(6., npart, x, y, eccang6);

    TCanvas * cColl2D = new TCanvas(Form("cColl2D_c%d",cbin), "cColl2D", 800, 800);
    cColl2D->Range(0,0,1,1);
    cColl2D->cd();

    for (int i = 0; i<A1; i++) {
        hA1->Fill(xb1[i], y1[i]);
        hA1_3D->Fill(xb1[i], y1[i], z1[i]);
        eA1[i] = new TEllipse(xb1[i]/scale+0.5, y1[i]/scale+0.5, rnuc/scale, rnuc/scale);
        eA1[i]->SetLineColor(kBlue);
        eA1[i]->Draw();
    }

    for (int i = 0; i<A2; i++) {
        hA2->Fill(xb2[i], y2[i]);
        hA2_3D->Fill(xb2[i], y2[i], z2[i]);
        eA2[i] = new TEllipse(xb2[i]/scale+0.5, y2[i]/scale+0.5, rnuc/scale, rnuc/scale);
        eA2[i]->SetLineColor(kRed);
        eA2[i]->Draw();
    }

    for (int i = 0; i<npart; i++) {
        hOverlap->Fill(x[i],y[i]);
        hCollision_3D->Fill(x[i],y[i],z[i]);
        eOverlap[i] = new TEllipse(x[i]/scale+0.5, y[i]/scale+0.5, rnuc/scale, rnuc/scale);
        eOverlap[i]->SetFillColor(kGreen);
        eOverlap[i]->Draw();
    }

    double len = 0.3;

    if (ShowAngs) {
        TArrow * line0 = new TArrow(0.5, 0.5, 0.5+1.4*len, 0.5);
        line0->SetLineColor(kBlack);
        line0->SetLineWidth(2);
        line0->Draw();
        TLatex * lt0 = new TLatex(0.5+1.4*len, 0.5, "#Psi_{RP}");
        lt0->Draw();

        TArrow * line1 = new TArrow(0.5, 0.5, 0.5+len*cos(ph1), 0.5+len*sin(ph1));
        line1->SetLineColor(kBlue);
        line1->SetLineWidth(2);
        line1->Draw();
        TLatex * lt1 = new TLatex(0.5+1.06*len*cos(ph1), 0.5+1.06*len*sin(ph1), "#Psi_{1}");
        lt1->Draw();

        TArrow * line2 = new TArrow(0.5, 0.5, 0.5+len*cos(ph2), 0.5+len*sin(ph2));
        line2->SetLineColor(kRed);
        line2->SetLineWidth(2);
        line2->Draw();
        TLatex * lt2 = new TLatex(0.5+1.06*len*cos(ph2), 0.5+1.06*len*sin(ph2), "#Psi_{2}");
        lt2->Draw();

        TArrow * line3 = new TArrow(0.5, 0.5, 0.5+len*cos(ph3), 0.5+len*sin(ph3));
        line3->SetLineColor(kGreen+2);
        line3->SetLineWidth(2);
        line3->Draw();
        TLatex * lt3 = new TLatex(0.5+1.06*len*cos(ph3), 0.5+1.06*len*sin(ph3), "#Psi_{3}");
        lt3->Draw();

        TArrow * line4 = new TArrow(0.5, 0.5, 0.5+len*cos(ph4), 0.5+len*sin(ph4));
        line4->SetLineColor(kMagenta);
        line4->SetLineWidth(2);
        line4->Draw();
        TLatex * lt4 = new TLatex(0.5+1.06*len*cos(ph4), 0.5+1.06*len*sin(ph4), "#Psi_{4}");
        lt4->Draw();

        TArrow * line5 = new TArrow(0.5, 0.5, 0.5+len*cos(ph5), 0.5+len*sin(ph5));
        line5->SetLineColor(kCyan+2);
        line5->SetLineWidth(2);
        line5->Draw();
        TLatex * lt5 = new TLatex(0.5+1.06*len*cos(ph5), 0.5+1.06*len*sin(ph5), "#Psi_{5}");
        lt5->Draw();

        TArrow * line6 = new TArrow(0.5, 0.5, 0.5+len*cos(ph6), 0.5+len*sin(ph6));
        line6->SetLineColor(kOrange+7);
        line6->SetLineWidth(2);
        line6->Draw();
        TLatex * lt6 = new TLatex(0.5+1.06*len*cos(ph6), 0.5+1.06*len*sin(ph6), "#Psi_{6}");
        lt6->Draw();

    }

    hA1->SetMarkerStyle(24);
    hA2->SetMarkerStyle(24);
    hA1->SetMarkerColor(kBlue);
    hA2->SetMarkerColor(kRed);
    hOverlap->SetMarkerColor(kGreen);
    hA1->SetXTitle("x (fm)");
    hA1->SetYTitle("y (fm)");
    hA2->SetXTitle("x (fm)");
    hA2->SetYTitle("y (fm)");
    hOverlap->SetXTitle("x (fm)");
    hOverlap->SetYTitle("y (fm)");


    if (print_plots) {
        if (!fopen("plots/2DColl","r")) system("mkdir plots/2DColl");
        if (!fopen("plots/3DColl","r")) system("mkdir plots/3DColl");
    }

    TString ctag =Form("_%d-%d",(int)centbins[cbin],(int)centbins[cbin+1]);

    cColl2D->cd();
    TPaveText * txColl2D = new TPaveText(0.059, 0.84, 0.40, 0.95);
    txColl2D->SetFillColor(0);
    txColl2D->SetBorderSize(0);
    txColl2D->SetTextAlign(12);
    txColl2D->AddText(Form("%d - %d%% Centrality",(int)centbins[cbin],(int)centbins[cbin+1]));
    txColl2D->AddText(Form("b = %2.3f fm",b));
    txColl2D->Draw();
    if (print_plots) cColl2D->Print(Form("plots/2DColl/Collision2D%s.png",ctag.Data()),"png");
    if (close_plots) cColl2D->Close();


    TCanvas * cColl3D = new TCanvas(Form("cColl3D_c%d",cbin),"cColl3D",800,800);
    cColl3D->cd();
    hA1_3D->SetXTitle("x (fm)");
    hA1_3D->SetYTitle("y (fm)");
    hA1_3D->SetZTitle("z (fm)");
    hA2_3D->SetXTitle("x (fm)");
    hA2_3D->SetYTitle("y (fm)");
    hA2_3D->SetZTitle("z (fm)");
    hA1_3D->Draw();
    hA2_3D->Draw("same");
    hCollision_3D->Draw("same");
    if (print_plots) cColl3D->Print(Form("plots/3DColl/Collision3D%s.png",ctag.Data()),"png");
    if (close_plots) cColl3D->Close();

}
