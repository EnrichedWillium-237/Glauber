//-------------------------------------
// Configuration options
static const Bool_t PlotCent = kFALSE; // set to kFALSE to run over the full set of events
static const Bool_t ShowAngs = kTRUE;
static const Bool_t print_plots = kTRUE;
static const Bool_t close_plots = kTRUE;

static const Bool_t minorAxes = kFALSE; // chose minor or major axes for reaction plane angles

static const int nevents = 1e3;
static const int ncentbins = 20;
static const double centbins[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
// Pb
/*
static const double R = 6.62;
static const double a = 0.546;
static const int A1 = 208;
static const int A2 = 208;
static const double sign = 64;
static const double rHard = 0.5;
*/
// Au
static const double R = 6.38;
static const double a = 0.535;
static const int A1 = 197;
static const int A2 = 197;
static const double sign = 42;
static const double rHard = 0.5;
//-------------------------------------

# include "src/Glauber.h"
# include "src/Cumu.h"

using namespace std;

void wsFunc( int A, double * x, double * y, double * z );
void collide( double * x1, double * y1, double * z1,
              double * x2, double * y2, double * z2,
              double * x,  double * y,  double * z,  int &npart );
double PsiN( double order, int np, double * x, double * y, double &eccang );
double PsiNCumu( double order, int np, double * x, double * y, double &eccang );
double eccentricity( double order, double psi, int npart, double * x, double * y, double &S );
void SetupCumulant( int cbin, int lcmin, int lcmax );
double cumulantEcc( double order, int npart, double * x, double * y, double psi1, double psi2, double psi3, double psi4, double psi5, double psi6 );
void recenter( int npart, double * x, double * y, double * z );

void Glauber()
{
    TH1::SetDefaultSumw2();

    if (!fopen("results","r")) system("mkdir results");
    if (print_plots) {
        if (!fopen("plots","r")) system("mkdir plots");
    }

    TStopwatch * sw = new TStopwatch;
    sw->Start();

    bconv->SetParameters(0.469714, 84.373, -601.761, 2785.36, -7157.02,
                         10113.6, -7353.77, 2147.29);
    hcent = new TH1D("cent", "cent", ncentbins, centbins);
    if (!PlotCent) cout << "Running for " << nevents << " events" << endl;
    cout << " A1: " << A1 << "\tA2: " << A2 << endl;
    cout << " R:  " << R << " fm" << endl;
    cout << " a:  " << a << " fm  \n" << endl;
    TString wsDef = "4*3.1415926*x*x/(1 + exp((x - [0])/[1]))";
    ranR = new TF1("ranR", wsDef.Data(), 0, 2*R);
    ranR->FixParameter(0,R);
    ranR->FixParameter(1,a);
    ranTheta = new TF1("ranTheta", "sin(x)", 0, TMath::Pi());
    for (int cbin = 0; cbin<ncentbins; cbin++) {
        int lmin = (int)centbins[cbin];
        int lmax = (int)centbins[cbin+1];
        for (int iorder = 1; iorder<=6; iorder++) {
            hpsi[iorder][cbin] = new TH1D(Form("psi%d_%d-%d",iorder,lmin,lmax), "psi", 200, -TMath::Pi()/(double)iorder-0.2/(double)iorder, TMath::Pi()/(double)iorder+0.2/(double)iorder);
            hpsi[iorder][cbin]->SetMarkerColor(colors[iorder]);
            hpsi[iorder][cbin]->SetLineColor(colors[iorder]);
            hpsi[iorder][cbin]->SetXTitle(Form("#Psi_{%d}^{moment}",iorder));

            hpsiC[iorder][cbin] = new TH1D(Form("psiC%d_%d-%d",iorder,lmin,lmax), "psiC", 200, -TMath::Pi()/(double)iorder-0.2/(double)iorder, TMath::Pi()/iorder+0.2/(double)iorder);
            hpsiC[iorder][cbin]->SetMarkerColor(colors[iorder]);
            hpsiC[iorder][cbin]->SetLineColor(colors[iorder]);
            hpsiC[iorder][cbin]->SetXTitle(Form("#Psi_{%d}^{cumulant}",iorder));
        }
    }

    TH1D * hWS1 = new TH1D("WS1", "WS1", 400, 0, 3*R);
    TH1D * hWS2 = new TH1D("WS2", "WS2", 400, 0, 3*R);

    ran = new TRandom3(0);

    double x1[300], y1[300], z1[300];
    double x2[300], y2[300], z2[300];
    double  x[600],  y[600],  z[600];

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        avNpart[cbin] = 0;
        avNpartCnt[cbin] = 0;
        avb[cbin] = 0;
        avbCnt[cbin] = 0;
        Int_t lcmin = (int)(centbins[cbin]);
        Int_t lcmax = (int)(centbins[cbin+1]);
        xy0[cbin] = new TH2D(Form("xy0_%d-%d",lcmin,lcmax), "xy0", 100, -10, 10, 100, -10, 10);
        xy0[cbin]->Sumw2();
        xy0[cbin]->SetOption("colz");
        xy1[cbin] = (TH2D *) xy0[cbin]->Clone(Form("xy1_%d-%d",lcmin,lcmax));
        xy2[cbin] = (TH2D *) xy0[cbin]->Clone(Form("xy2_%d-%d",lcmin,lcmax));
        xy3[cbin] = (TH2D *) xy0[cbin]->Clone(Form("xy3_%d-%d",lcmin,lcmax));
        xy4[cbin] = (TH2D *) xy0[cbin]->Clone(Form("xy4_%d-%d",lcmin,lcmax));
        xy5[cbin] = (TH2D *) xy0[cbin]->Clone(Form("xy5_%d-%d",lcmin,lcmax));
        xy6[cbin] = (TH2D *) xy0[cbin]->Clone(Form("xy6_%d-%d",lcmin,lcmax));

        Psi2Psi1[cbin] = new TH2D(Form("Psi2Psi1_%d-%d",lcmin,lcmax), "Psi2Psi1", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi()/2., TMath::Pi()/2.);
        Psi3Psi1[cbin] = new TH2D(Form("Psi3Psi1_%d-%d",lcmin,lcmax), "Psi3Psi1", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi()/3., TMath::Pi()/3.);
        Psi3Psi2[cbin] = new TH2D(Form("Psi3Psi2_%d-%d",lcmin,lcmax), "Psi3Psi2", 100, -TMath::Pi()/2, TMath::Pi()/2, 100, -TMath::Pi()/3., TMath::Pi()/3.);
        Psi4Psi2[cbin] = new TH2D(Form("Psi4Psi2_%d-%d",lcmin,lcmax), "Psi4Psi2", 100, -TMath::Pi()/2, TMath::Pi()/2, 100, -TMath::Pi()/4., TMath::Pi()/4.);
        Psi2Psi1[cbin]->Sumw2();
        Psi3Psi1[cbin]->Sumw2();
        Psi3Psi2[cbin]->Sumw2();
        Psi4Psi2[cbin]->Sumw2();
        Psi2Psi1[cbin]->SetOption("colz");
        Psi3Psi1[cbin]->SetOption("colz");
        Psi3Psi2[cbin]->SetOption("colz");
        Psi4Psi2[cbin]->SetOption("colz");
        Psi2Psi1[cbin]->SetXTitle("#Psi_{1}");
        Psi2Psi1[cbin]->SetYTitle("#Psi_{2}");
        Psi3Psi1[cbin]->SetXTitle("#Psi_{1}");
        Psi3Psi1[cbin]->SetYTitle("#Psi_{3}");
        Psi3Psi2[cbin]->SetXTitle("#Psi_{2}");
        Psi3Psi2[cbin]->SetYTitle("#Psi_{3}");
        Psi4Psi2[cbin]->SetXTitle("#Psi_{2}");
        Psi4Psi2[cbin]->SetYTitle("#Psi_{4}");

        CosPsi2Psi1[cbin] = new TH1D(Form("CosPsi2Psi1_%d-%d",lcmin,lcmax), "CosPsi2Psi1", 100, -1, 1);
        CosPsi3Psi1[cbin] = new TH1D(Form("CosPsi3Psi1_%d-%d",lcmin,lcmax), "CosPsi3Psi1", 100, -1, 1);
        CosPsi3Psi2[cbin] = new TH1D(Form("CosPsi3Psi2_%d-%d",lcmin,lcmax), "CosPsi3Psi2", 100, -1, 1);
        CosPsi4Psi2[cbin] = new TH1D(Form("CosPsi4Psi2_%d-%d",lcmin,lcmax), "CosPsi4Psi2", 100, -1, 1);

        CosPsi2Psi1[cbin]->SetXTitle("cos(2#Psi_{2} - #Psi_{1})");
        CosPsi3Psi1[cbin]->SetXTitle("cos(3#Psi_{3} - #Psi_{1})");
        CosPsi3Psi2[cbin]->SetXTitle("cos(3#Psi_{3} - 2#Psi_{2})");
        CosPsi4Psi2[cbin]->SetXTitle("cos(4#Psi_{4} - 2#Psi_{2})");

        Psi2Psi1Diff[cbin] = new TH2D(Form("Psi2Psi1Diff_%d-%d",lcmin,lcmax), "Psi2Psi1Diff", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi()/2., TMath::Pi()/2.);
        Psi3Psi1Diff[cbin] = new TH2D(Form("Psi3Psi1Diff_%d-%d",lcmin,lcmax), "Psi3Psi1Diff", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi()/3., TMath::Pi()/3.);
        Psi5Psi1Diff[cbin] = new TH2D(Form("Psi5Psi1Diff_%d-%d",lcmin,lcmax), "Psi5Psi1Diff", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi()/5., TMath::Pi()/5.);
        Psi3Psi2Diff[cbin] = new TH2D(Form("Psi3Psi2Diff_%d-%d",lcmin,lcmax), "Psi3Psi2Diff", 100, -TMath::Pi()/2., TMath::Pi()/2., 100, -TMath::Pi()/3., TMath::Pi()/3.);
        Psi4Psi2Diff[cbin] = new TH2D(Form("Psi4Psi2Diff_%d-%d",lcmin,lcmax), "Psi4Psi2Diff", 100, -TMath::Pi()/2., TMath::Pi()/2., 100, -TMath::Pi()/4., TMath::Pi()/4.);

        Psi2Psi1Diff[cbin]->SetOption("colz");
        Psi3Psi1Diff[cbin]->SetOption("colz");
        Psi5Psi1Diff[cbin]->SetOption("colz");
        Psi3Psi2Diff[cbin]->SetOption("colz");
        Psi4Psi2Diff[cbin]->SetOption("colz");
        Psi2Psi1Diff[cbin]->SetXTitle("#Psi_{2} - #Psi_{1}");
        Psi2Psi1Diff[cbin]->SetYTitle("#Psi_{2}");
        Psi3Psi1Diff[cbin]->SetXTitle("#Psi_{3} - #Psi_{1}");
        Psi3Psi1Diff[cbin]->SetYTitle("#Psi_{3}}");
        Psi5Psi1Diff[cbin]->SetXTitle("#Psi_{5} - #Psi_{1}");
        Psi5Psi1Diff[cbin]->SetYTitle("#Psi_{5}");
        Psi3Psi2Diff[cbin]->SetXTitle("#Psi_{3} - #Psi_{2}");
        Psi3Psi2Diff[cbin]->SetYTitle("#Psi_{3}");
        Psi4Psi2Diff[cbin]->SetXTitle("#Psi_{4} - #Psi_{2}");
        Psi4Psi2Diff[cbin]->SetYTitle("#Psi_{4}");

        for (int iorder = 0; iorder<=6; iorder++) {
            eccAng[iorder][cbin] = new TH1D(Form("eccAng%d_%d-%d",iorder,lcmin,lcmax), "eccAng", 360, -TMath::Pi()/(double)iorder, TMath::Pi()/(double)iorder);
            eccDist[iorder][cbin] = new TH1D(Form("eccDist%d_%d-%d",iorder,lcmin,lcmax), "eccDist", 360, 0, 1.0);
        }

        SetupCumulant(cbin, lcmin, lcmax); // create cumulant histograms
    }

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        ecc_2_1[cbin] = 0;
        ecc_3_1[cbin] = 0;
        ecc_3_2[cbin] = 0;
        ecc_4_2[cbin] = 0;

        ecc_2_1_cnt[cbin] = 0;
        ecc_3_1_cnt[cbin] = 0;
        ecc_3_2_cnt[cbin] = 0;
        ecc_4_2_cnt[cbin] = 0;

        eccC_2_1[cbin] = 0;
        eccC_3_1[cbin] = 0;
        eccC_3_2[cbin] = 0;
        eccC_4_2[cbin] = 0;

        eccC_2_1_cnt[cbin] = 0;
        eccC_3_1_cnt[cbin] = 0;
        eccC_3_2_cnt[cbin] = 0;
        eccC_4_2_cnt[cbin] = 0;

        for (int j = 0; j<7; j++) {
            ecc[j][cbin] = 0;
            eccCnt[j][cbin] = 0;

            eccC[j][cbin] = 0;
            eccCCnt[j][cbin] = 0;
        }
        c_cnt[cbin] = 0;
    }

    double eccang = 0;
    double eccangC = 0;
    Int_t maxe = nevents;
    if (PlotCent) maxe = 1;

    //-- main event loop
    for (int ievent = 0; ievent<maxe; ievent++) {
        if (fmod(double(ievent+1), nevents/20.)==0) {
            sw->Continue();
            double elapse = sw->RealTime();
            cout<<(int)(100*(ievent/(double)nevents)+0.5)<<"  Elapsed: "<<elapse<<"\tTime per event: "<<elapse/(double)ievent<<endl;
        }
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            double cent = ran->Uniform(centbins[cbin],centbins[cbin+1]);
            double b = bconv->Eval(cent/100.);
            if (PlotCent && ievent == 0) {
                cout<<Form(" %0.2f%% Centrality",cent)<<"  bin: "<<centbins[cbin]<<"-"<<centbins[cbin+1]<<"%\t";
                cout<<" b = "<<b<<" fm \n"<<endl;
            }
            hcent->Fill(cent);

            wsFunc(A1,x1,y1,z1);
            wsFunc(A2,x2,y2,z2);
            for (int i = 0; i<A1; i++) {
                double r1 = sqrt( pow(x1[i],2) + pow(y1[i],2) + pow(z1[i],2) );
                hWS1->Fill(r1);
            }
            for (int i = 0; i<A2; i++) {
                double r2 = sqrt( pow(x2[i],2) + pow(y2[i],2) + pow(z2[i],2) );
                hWS2->Fill(r2);
            }

            double xb = b;
            double xb1[300];
            double xb2[300];
            for (int i = 0; i<A1; i++) {
                xb1[i] = x1[i] - xb/2.;
            }
            for (int i = 0; i<A2; i++) {
                xb2[i] = x2[i] + xb/2.;
            }

            int npart = 0;
            collide( xb1, y1, z1, xb2, y2, z2, x, y, z, npart );
            //-- draw single collision event
            if (PlotCent && ievent == 0) {
                PlotColl( xb1, y1, z1, xb2, y2, z2, x, y, z, npart, cbin, b );
            }

            avNpart[cbin]+=npart;
            ++avNpartCnt[cbin];
            avb[cbin]+=b;
            ++avbCnt[cbin];
            if (npart == 0) continue;
            recenter( npart, x, y, z );

            //-- moment calculation

            for (int iorder = 1; iorder<=6; iorder++) {
                double order = (double)iorder;
                double psi;
                if (minorAxes) psi = PsiN(order, npart, x, y, eccang);
                else psi = PsiN(order, npart, x, y, eccang) + TMath::Pi()/order;
                if (npart>5 && psi>-4) {
                    hpsi[iorder][cbin]->Fill(bounds(order,psi));
                    double Sval = 0;
                    double e = eccentricity(order, eccang, npart, x, y, Sval);

                    ecc[iorder][cbin] += e;
                    ecc_pow2[iorder][cbin] += pow(e,2);
                    S[iorder][cbin] += Sval;
                    S_pow2[iorder][cbin] += pow(Sval,2);
                    eccCnt[iorder][cbin] += 1;
                    eccDist[iorder][cbin]->Fill(e);

                    eccAng[iorder][cbin]->Fill(bounds(iorder,eccang/(double)iorder));
                }
            }

            double eccang1, eccang2, eccang3, eccang4, eccang5, eccang6;
            if (npart>5) {
                double psi1, psi2, psi3, psi4, psi5, psi6;
                // semiminor axes
                if (minorAxes) {
                    psi1 = PsiN(1., npart, x, y, eccang1);
                    psi2 = PsiN(2., npart, x, y, eccang2);
                    psi3 = PsiN(3., npart, x, y, eccang3);
                    psi4 = PsiN(4., npart, x, y, eccang4);
                    psi5 = PsiN(5., npart, x, y, eccang5);
                    psi6 = PsiN(6., npart, x, y, eccang6);
                } else { // semimajor axes
                    psi1 = PsiN(1., npart, x, y, eccang1) + TMath::Pi()/1.;
                    psi2 = PsiN(2., npart, x, y, eccang2) + TMath::Pi()/2.;
                    psi3 = PsiN(3., npart, x, y, eccang3) + TMath::Pi()/3.;
                    psi4 = PsiN(4., npart, x, y, eccang4) + TMath::Pi()/4.;
                    psi5 = PsiN(5., npart, x, y, eccang5) + TMath::Pi()/5.;
                    psi6 = PsiN(6., npart, x, y, eccang6) + TMath::Pi()/6.;
                }

                for (int k = 0; k<npart; k++) {
                    double r = sqrt( pow(x[k],2) + pow(y[k],2) );
                    double phi = TMath::ATan2(y[k],x[k]);
                    double xx = r*TMath::Cos(phi);
                    double yy = r*TMath::Sin(phi);
                    xy0[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi1);
                    yy = r*TMath::Sin(phi - psi1);
                    xy1[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi2);
                    yy = r*TMath::Sin(phi - psi2);
                    xy2[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi3);
                    yy = r*TMath::Sin(phi - psi3);
                    xy3[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi4);
                    yy = r*TMath::Sin(phi - psi4);
                    xy4[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi5);
                    yy = r*TMath::Sin(phi - psi5);
                    xy5[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi6);
                    yy = r*TMath::Sin(phi - psi6);
                    xy6[cbin]->Fill(xx, yy);
                }

                Psi2Psi1[cbin]->Fill(bounds(1,psi1), bounds(2,psi2));
                Psi3Psi1[cbin]->Fill(bounds(1,psi1), bounds(3,psi3));
                Psi3Psi2[cbin]->Fill(bounds(2,psi2), bounds(3,psi3));
                Psi4Psi2[cbin]->Fill(bounds(2,psi2), bounds(4,psi4));

                CosPsi2Psi1[cbin]->Fill(TMath::Cos(2*psi2 - psi1));
                CosPsi3Psi1[cbin]->Fill(TMath::Cos(3*psi3 - psi1));
                CosPsi3Psi2[cbin]->Fill(TMath::Cos(3*psi3 - 2*psi2));
                CosPsi4Psi2[cbin]->Fill(TMath::Cos(4*psi4 - 2*psi2));

                Psi2Psi1Diff[cbin]->Fill(psi1-psi2, psi2);
                Psi3Psi1Diff[cbin]->Fill(psi1-psi3, psi3);
                Psi5Psi1Diff[cbin]->Fill(psi1-psi5, psi5);
                Psi3Psi2Diff[cbin]->Fill(psi2-psi3, psi3);
                Psi4Psi2Diff[cbin]->Fill(psi2-psi4, psi4);

                double S21 = 0;
                double S31 = 0;
                double S32 = 0;
                double S42 = 0;

                double e21 = eccentricity(2., eccang1, npart, x, y, S21);
                double e31 = eccentricity(3., eccang1, npart, x, y, S31);
                double e32 = eccentricity(3., eccang2, npart, x, y, S32);
                double e42 = eccentricity(4., eccang2, npart, x, y, S42);

                ecc_2_1[cbin] += e21;
                ecc_3_1[cbin] += e31;
                ecc_3_2[cbin] += e32;
                ecc_4_2[cbin] += e42;

                S_2_1[cbin] += S21;
                S_3_1[cbin] += S31;
                S_3_2[cbin] += S32;
                S_4_2[cbin] += S42;

                ecc_2_1_pow2[cbin] += pow(e21,2);
                ecc_3_1_pow2[cbin] += pow(e31,2);
                ecc_3_2_pow2[cbin] += pow(e32,2);
                ecc_4_2_pow2[cbin] += pow(e42,2);

                S_2_1_pow2[cbin] += pow(S21,2);
                S_3_1_pow2[cbin] += pow(S31,2);
                S_3_2_pow2[cbin] += pow(S32,2);
                S_4_2_pow2[cbin] += pow(S42,2);

                ++ecc_2_1_cnt[cbin];
                ++ecc_3_1_cnt[cbin];
                ++ecc_3_2_cnt[cbin];
                ++ecc_4_2_cnt[cbin];

                c_4_p4_p2[cbin] += TMath::Cos(4*(psi4 - psi2));
                c_8_p4_p2[cbin] += TMath::Cos(8*(psi4 - psi2));
                c_12_p4_p2[cbin] += TMath::Cos(12*(psi4 - psi2));
                c_6_p3_p2[cbin] += TMath::Cos(6*(psi3 - psi2));
                c_6_p2_p6[cbin] += TMath::Cos(6*(psi2 - psi6));
                c_6_p3_p6[cbin] += TMath::Cos(6*(psi3 - psi6));
                c_12_p3_p4[cbin] += TMath::Cos(12*(psi3 - psi4));
                c_10_p2_p5[cbin] += TMath::Cos(10*(psi2 - psi5));

                c_2p2_3p3_5p5[cbin] += TMath::Cos(2*psi2 + 3*psi3 - 5*psi5);
                c_2p2_4p4_6p6[cbin] += TMath::Cos(2*psi2 + 4*psi4 - 6*psi6);
                c_2p2_6p3_4p4[cbin] += TMath::Cos(2*psi2 - 6*psi3 + 4*psi4);
                c_8p2_3p3_5p5[cbin] += TMath::Cos(-8*psi2 + 3*psi3 + 5*psi5);
                c_10p2_4p4_6p6[cbin] += TMath::Cos(-10*psi2 + 4*psi4 + 6*psi6);
                c_10p2_6p3_4p4[cbin] += TMath::Cos(-10*psi2 + 6*psi3 + 4*psi4);

                c_4_p4_p2_pow2[cbin] += pow(TMath::Cos(4*(psi4 - psi2)),2);
                c_8_p4_p2_pow2[cbin] += pow(TMath::Cos(8*(psi4 - psi2)),2);
                c_12_p4_p2_pow2[cbin] += pow(TMath::Cos(12*(psi4 - psi2)),2);
                c_6_p3_p2_pow2[cbin] += pow(TMath::Cos(6*(psi3 - psi2)),2);
                c_6_p2_p6_pow2[cbin] += pow(TMath::Cos(6*(psi2 - psi6)),2);
                c_6_p3_p6_pow2[cbin] += pow(TMath::Cos(6*(psi3 - psi6)),2);
                c_12_p3_p4_pow2[cbin] += pow(TMath::Cos(12*(psi3 - psi4)),2);
                c_10_p2_p5_pow2[cbin] += pow(TMath::Cos(10*(psi2 - psi5)),2);

                c_2p2_3p3_5p5_pow2[cbin] += pow(TMath::Cos(2*psi2 + 3*psi3 - 5*psi5),2);
                c_2p2_4p4_6p6_pow2[cbin] += pow(TMath::Cos(2*psi2 + 4*psi4 - 6*psi6),2);
                c_2p2_6p3_4p4_pow2[cbin] += pow(TMath::Cos(2*psi2 - 6*psi3 + 4*psi4),2);
                c_8p2_3p3_5p5_pow2[cbin] += pow(TMath::Cos(-8*psi2 + 3*psi3 + 5*psi5),2);
                c_10p2_4p4_6p6_pow2[cbin] += pow(TMath::Cos(-10*psi2 + 4*psi4 + 6*psi6),2);
                c_10p2_6p3_4p4_pow2[cbin] += pow(TMath::Cos(-10*psi2 + 6*psi3 + 4*psi4),2);

                c_6p2_p3[cbin] += TMath::Cos(6*(psi2 - psi3));
                c_4p2_p4[cbin] += TMath::Cos(4*(psi2 - psi4));
                c_6p2_p6[cbin] += TMath::Cos(6*(psi2 - psi6));
                c_6p3_p6[cbin] += TMath::Cos(6*(psi3 - psi6));
                c_10p2_p5[cbin] += TMath::Cos(10*(psi2 - psi5));
                c_15p3_p5[cbin] += TMath::Cos(15*(psi3 - psi5));

                c_2p1_p2[cbin] += TMath::Cos(2*(psi1 - psi2));
                c_3p1_p3[cbin] += TMath::Cos(3*(psi1 - psi3));
                c_4p1_p4[cbin] += TMath::Cos(4*(psi1 - psi4));
                c_5p1_p5[cbin] += TMath::Cos(5*(psi1 - psi5));
                c_6p1_p6[cbin] += TMath::Cos(6*(psi1 - psi6));

                c_1p1_2p2_3p3[cbin] += TMath::Cos(1*psi1 + 2*psi2 - 3*psi3);
                c_4p1_2p2_6p3[cbin] += TMath::Cos(4*psi1 + 2*psi2 - 6*psi3);
                c_2p1_4p2_6p3[cbin] += TMath::Cos(2*psi1 + 4*psi2 - 6*psi3);
                c_5p1_2p2_3p3[cbin] += TMath::Cos(5*psi1 - 2*psi2 - 3*psi3);
                c_2p1_2p2_4p4[cbin] += TMath::Cos(2*psi1 + 2*psi2 - 4*psi4);
                c_2p1_6p2_4p4[cbin] += TMath::Cos(2*psi1 - 6*psi2 + 4*psi4);
                c_2p1_6p2_8p4[cbin] += TMath::Cos(2*psi1 + 6*psi2 - 8*psi4);

                c_1p1_3p3_4p4[cbin] += TMath::Cos(1*psi1 + 3*psi3 - 4*psi4);
                c_2p1_6p3_4p4[cbin] += TMath::Cos(2*psi1 - 6*psi3 + 4*psi4);
                c_2p1_4p2_5p5[cbin] += TMath::Cos(1*psi1 + 4*psi2 - 5*psi5);
                c_3p1_2p2_5p5[cbin] += TMath::Cos(3*psi1 + 2*psi2 - 5*psi5);
                c_1p1_6p3_5p5[cbin] += TMath::Cos(1*psi1 - 6*psi3 + 5*psi5);
                c_2p1_3p3_5p5[cbin] += TMath::Cos(2*psi1 + 3*psi3 - 5*psi5);
                c_4p1_9p3_5p5[cbin] += TMath::Cos(4*psi1 - 9*psi3 + 5*psi5);


                c_6p2_p3_pow2[cbin] += pow(TMath::Cos(6*(psi2 - psi3)),2);
                c_4p2_p4_pow2[cbin] += pow(TMath::Cos(4*(psi2 - psi4)),2);
                c_6p2_p6_pow2[cbin] += pow(TMath::Cos(6*(psi2 - psi6)),2);
                c_6p3_p6_pow2[cbin] += pow(TMath::Cos(6*(psi3 - psi6)),2);
                c_10p2_p5_pow2[cbin] += pow(TMath::Cos(10*(psi2 - psi5)),2);
                c_15p3_p5_pow2[cbin] += pow(TMath::Cos(15*(psi3 - psi5)),2);

                c_2p1_p2_pow2[cbin] += pow(TMath::Cos(2*(psi1 - psi2)),2);
                c_3p1_p3_pow2[cbin] += pow(TMath::Cos(3*(psi1 - psi3)),2);
                c_4p1_p4_pow2[cbin] += pow(TMath::Cos(4*(psi1 - psi4)),2);
                c_5p1_p5_pow2[cbin] += pow(TMath::Cos(5*(psi1 - psi5)),2);
                c_6p1_p6_pow2[cbin] += pow(TMath::Cos(6*(psi1 - psi6)),2);

                c_1p1_2p2_3p3_pow2[cbin] += pow(TMath::Cos(1*psi1 + 2*psi2 - 3*psi3),2);
                c_4p1_2p2_6p3_pow2[cbin] += pow(TMath::Cos(4*psi1 + 2*psi2 - 6*psi3),2);
                c_2p1_4p2_6p3_pow2[cbin] += pow(TMath::Cos(2*psi1 + 4*psi2 - 6*psi3),2);
                c_5p1_2p2_3p3_pow2[cbin] += pow(TMath::Cos(5*psi1 - 2*psi2 - 3*psi3),2);
                c_2p1_2p2_4p4_pow2[cbin] += pow(TMath::Cos(2*psi1 + 2*psi2 - 4*psi4),2);
                c_2p1_6p2_4p4_pow2[cbin] += pow(TMath::Cos(2*psi1 - 6*psi2 + 4*psi4),2);
                c_2p1_6p2_8p4_pow2[cbin] += pow(TMath::Cos(2*psi1 + 6*psi2 - 8*psi4),2);

                c_1p1_3p3_4p4_pow2[cbin] += pow(TMath::Cos(1*psi1 + 3*psi3 - 4*psi4),2);
                c_2p1_6p3_4p4_pow2[cbin] += pow(TMath::Cos(2*psi1 - 6*psi3 + 4*psi4),2);
                c_2p1_4p2_5p5_pow2[cbin] += pow(TMath::Cos(1*psi1 + 4*psi2 - 5*psi5),2);
                c_3p1_2p2_5p5_pow2[cbin] += pow(TMath::Cos(3*psi1 + 2*psi2 - 5*psi5),2);
                c_1p1_6p3_5p5_pow2[cbin] += pow(TMath::Cos(1*psi1 - 6*psi3 + 5*psi5),2);
                c_2p1_3p3_5p5_pow2[cbin] += pow(TMath::Cos(2*psi1 + 3*psi3 - 5*psi5),2);
                c_4p1_9p3_5p5_pow2[cbin] += pow(TMath::Cos(4*psi1 - 9*psi3 + 5*psi5),2);

                c_cnt[cbin]++;
            }


            //-- cumulant calculation
            double psiC1, psiC2, psiC3, psiC4, psiC5, psiC6;
            if (minorAxes) {
                psiC1 = PsiNCumu(1., npart, x, y, eccangC);
                psiC2 = PsiNCumu(2., npart, x, y, eccangC);
                psiC3 = PsiNCumu(3., npart, x, y, eccangC);
                psiC4 = PsiNCumu(4., npart, x, y, eccangC);
                psiC5 = PsiNCumu(5., npart, x, y, eccangC);
                psiC6 = PsiNCumu(6., npart, x, y, eccangC);
            } else {
                psiC1 = PsiNCumu(1., npart, x, y, eccangC) + TMath::Pi()/1.;
                psiC2 = PsiNCumu(2., npart, x, y, eccangC) + TMath::Pi()/2.;
                psiC3 = PsiNCumu(3., npart, x, y, eccangC) + TMath::Pi()/3.;
                psiC4 = PsiNCumu(4., npart, x, y, eccangC) + TMath::Pi()/4.;
                psiC5 = PsiNCumu(5., npart, x, y, eccangC) + TMath::Pi()/5.;
                psiC6 = PsiNCumu(6., npart, x, y, eccangC) + TMath::Pi()/6.;
            }
            for (int iorder = 1; iorder<=6; iorder++) {
                double order = (double)iorder;
                double psi = PsiNCumu(iorder, npart, x, y, eccangC);
                if (npart>5 && psi>-4) {
                    hpsiC[iorder][cbin]->Fill(bounds(order,psi));
                    double e = cumulantEcc(order, npart, x, y, psiC1, psiC2, psiC3, psiC4, psiC5, psiC6);

                    eccC[iorder][cbin] += e;
                    eccC_pow2[iorder][cbin] += pow(e,2);
                    eccCCnt[iorder][cbin] += 1;
                    eccDistC[iorder][cbin]->Fill(e);

                    eccAngC[iorder][cbin]->Fill(bounds(iorder,eccangC/(double)iorder));
                }
            }

            double eccangC1, eccangC2, eccangC3, eccangC4, eccangC5, eccangC6;
            if (npart>5) {
                double psi1 = PsiNCumu(1., npart, x, y, eccangC1);
                double psi2 = PsiNCumu(2., npart, x, y, eccangC2);
                double psi3 = PsiNCumu(3., npart, x, y, eccangC3);
                double psi4 = PsiNCumu(4., npart, x, y, eccangC4);
                double psi5 = PsiNCumu(5., npart, x, y, eccangC5);
                double psi6 = PsiNCumu(6., npart, x, y, eccangC6);
                for (int k = 0; k<npart; k++) {
                    double r = sqrt( pow(x[k],2) + pow(y[k],2) );
                    double phi = TMath::ATan2(y[k],x[k]);
                    double xx = r*TMath::Cos(phi);
                    double yy = r*TMath::Sin(phi);
                    xyC0[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi1);
                    yy = r*TMath::Sin(phi - psi1);
                    xyC1[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi2);
                    yy = r*TMath::Sin(phi - psi2);
                    xyC2[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi3);
                    yy = r*TMath::Sin(phi - psi3);
                    xyC3[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi4);
                    yy = r*TMath::Sin(phi - psi4);
                    xyC4[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi5);
                    yy = r*TMath::Sin(phi - psi5);
                    xyC5[cbin]->Fill(xx, yy);

                    xx = r*TMath::Cos(phi - psi6);
                    yy = r*TMath::Sin(phi - psi6);
                    xyC6[cbin]->Fill(xx, yy);
                }

                Psi2Psi1C[cbin]->Fill(psi1, psi2);
                Psi3Psi1C[cbin]->Fill(psi1, psi3);
                Psi3Psi2C[cbin]->Fill(psi2, psi3);
                Psi4Psi2C[cbin]->Fill(psi2, psi4);

                CosPsi2Psi1C[cbin]->Fill(TMath::Cos(2*psi2 - psi1));
                CosPsi3Psi1C[cbin]->Fill(TMath::Cos(3*psi3 - psi1));
                CosPsi3Psi2C[cbin]->Fill(TMath::Cos(3*psi3 - 2*psi2));
                CosPsi4Psi2C[cbin]->Fill(TMath::Cos(4*psi4 - 2*psi2));

                Psi2Psi1DiffC[cbin]->Fill(psi1-psi2, psi2);
                Psi3Psi1DiffC[cbin]->Fill(psi1-psi3, psi3);
                Psi5Psi1DiffC[cbin]->Fill(psi1-psi5, psi5);
                Psi3Psi2DiffC[cbin]->Fill(psi2-psi3, psi3);
                Psi4Psi2DiffC[cbin]->Fill(psi2-psi4, psi4);

                double e21 = cumulantEcc(2., npart, x, y, psi1, psi2, psi3, psi4, psi5, psi6);
                double e31 = cumulantEcc(3., npart, x, y, psi1, psi2, psi3, psi4, psi5, psi6);
                double e32 = cumulantEcc(3., npart, x, y, psi1, psi2, psi3, psi4, psi5, psi6);
                double e42 = cumulantEcc(4., npart, x, y, psi1, psi2, psi3, psi4, psi5, psi6);

                eccC_2_1[cbin] += e21;
                eccC_3_1[cbin] += e31;
                eccC_3_2[cbin] += e32;
                eccC_4_2[cbin] += e42;

                eccC_2_1_pow2[cbin] += pow(e21,2);
                eccC_3_1_pow2[cbin] += pow(e31,2);
                eccC_3_2_pow2[cbin] += pow(e32,2);
                eccC_4_2_pow2[cbin] += pow(e42,2);

                ++eccC_2_1_cnt[cbin];
                ++eccC_3_1_cnt[cbin];
                ++eccC_3_2_cnt[cbin];
                ++eccC_4_2_cnt[cbin];

                cC_4_p4_p2[cbin] += TMath::Cos(4*(psi4 - psi2));
                cC_8_p4_p2[cbin] += TMath::Cos(8*(psi4 - psi2));
                cC_12_p4_p2[cbin] += TMath::Cos(12*(psi4 - psi2));
                cC_6_p3_p2[cbin] += TMath::Cos(6*(psi3 - psi2));
                cC_6_p2_p6[cbin] += TMath::Cos(6*(psi2 - psi6));
                cC_6_p3_p6[cbin] += TMath::Cos(6*(psi3 - psi6));
                cC_12_p3_p4[cbin] += TMath::Cos(12*(psi3 - psi4));
                cC_10_p2_p5[cbin] += TMath::Cos(10*(psi2 - psi5));

                cC_2p2_3p3_5p5[cbin] += TMath::Cos(2*psi2 + 3*psi3 - 5*psi5);
                cC_2p2_4p4_6p6[cbin] += TMath::Cos(2*psi2 + 4*psi4 - 6*psi6);
                cC_2p2_6p3_4p4[cbin] += TMath::Cos(2*psi2 - 6*psi3 + 4*psi4);
                cC_8p2_3p3_5p5[cbin] += TMath::Cos(-8*psi2 + 3*psi3 + 5*psi5);
                cC_10p2_4p4_6p6[cbin] += TMath::Cos(-10*psi2 + 4*psi4 + 6*psi6);
                cC_10p2_6p3_4p4[cbin] += TMath::Cos(-10*psi2 + 6*psi3 + 4*psi4);

                cC_4_p4_p2_pow2[cbin] += pow(TMath::Cos(4*(psi4 - psi2)),2);
                cC_8_p4_p2_pow2[cbin] += pow(TMath::Cos(8*(psi4 - psi2)),2);
                cC_12_p4_p2_pow2[cbin] += pow(TMath::Cos(12*(psi4 - psi2)),2);
                cC_6_p3_p2_pow2[cbin] += pow(TMath::Cos(6*(psi3 - psi2)),2);
                cC_6_p2_p6_pow2[cbin] += pow(TMath::Cos(6*(psi2 - psi6)),2);
                cC_6_p3_p6_pow2[cbin] += pow(TMath::Cos(6*(psi3 - psi6)),2);
                cC_12_p3_p4_pow2[cbin] += pow(TMath::Cos(12*(psi3 - psi4)),2);
                cC_10_p2_p5_pow2[cbin] += pow(TMath::Cos(10*(psi2 - psi5)),2);

                cC_2p2_3p3_5p5_pow2[cbin] += pow(TMath::Cos(2*psi2 + 3*psi3 - 5*psi5),2);
                cC_2p2_4p4_6p6_pow2[cbin] += pow(TMath::Cos(2*psi2 + 4*psi4 - 6*psi6),2);
                cC_2p2_6p3_4p4_pow2[cbin] += pow(TMath::Cos(2*psi2 - 6*psi3 + 4*psi4),2);
                cC_8p2_3p3_5p5_pow2[cbin] += pow(TMath::Cos(-8*psi2 + 3*psi3 + 5*psi5),2);
                cC_10p2_4p4_6p6_pow2[cbin] += pow(TMath::Cos(-10*psi2 + 4*psi4 + 6*psi6),2);
                cC_10p2_6p3_4p4_pow2[cbin] += pow(TMath::Cos(-10*psi2 + 6*psi3 + 4*psi4),2);

                cC_6p2_p3[cbin] += TMath::Cos(6*(psi2 - psi3));
                cC_4p2_p4[cbin] += TMath::Cos(4*(psi2 - psi4));
                cC_6p2_p6[cbin] += TMath::Cos(6*(psi2 - psi6));
                cC_6p3_p6[cbin] += TMath::Cos(6*(psi3 - psi6));
                cC_10p2_p5[cbin] += TMath::Cos(10*(psi2 - psi5));
                cC_15p3_p5[cbin] += TMath::Cos(15*(psi3 - psi5));

                cC_2p1_p2[cbin] += TMath::Cos(2*(psi1 - psi2));
                cC_3p1_p3[cbin] += TMath::Cos(3*(psi1 - psi3));
                cC_4p1_p4[cbin] += TMath::Cos(4*(psi1 - psi4));
                cC_5p1_p5[cbin] += TMath::Cos(5*(psi1 - psi5));
                cC_6p1_p6[cbin] += TMath::Cos(6*(psi1 - psi6));

                cC_1p1_2p2_3p3[cbin] += TMath::Cos(1*psi1 + 2*psi2 - 3*psi3);
                cC_4p1_2p2_6p3[cbin] += TMath::Cos(4*psi1 + 2*psi2 - 6*psi3);
                cC_2p1_4p2_6p3[cbin] += TMath::Cos(2*psi1 + 4*psi2 - 6*psi3);
                cC_5p1_2p2_3p3[cbin] += TMath::Cos(5*psi1 - 2*psi2 - 3*psi3);
                cC_2p1_2p2_4p4[cbin] += TMath::Cos(2*psi1 + 2*psi2 - 4*psi4);
                cC_2p1_6p2_4p4[cbin] += TMath::Cos(2*psi1 - 6*psi2 + 4*psi4);
                cC_2p1_6p2_8p4[cbin] += TMath::Cos(2*psi1 + 6*psi2 - 8*psi4);

                cC_1p1_3p3_4p4[cbin] += TMath::Cos(1*psi1 + 3*psi3 - 4*psi4);
                cC_2p1_6p3_4p4[cbin] += TMath::Cos(2*psi1 - 6*psi3 + 4*psi4);
                cC_2p1_4p2_5p5[cbin] += TMath::Cos(1*psi1 + 4*psi2 - 5*psi5);
                cC_3p1_2p2_5p5[cbin] += TMath::Cos(3*psi1 + 2*psi2 - 5*psi5);
                cC_1p1_6p3_5p5[cbin] += TMath::Cos(1*psi1 - 6*psi3 + 5*psi5);
                cC_2p1_3p3_5p5[cbin] += TMath::Cos(2*psi1 + 3*psi3 - 5*psi5);
                cC_4p1_9p3_5p5[cbin] += TMath::Cos(4*psi1 - 9*psi3 + 5*psi5);


                cC_6p2_p3_pow2[cbin] += pow(TMath::Cos(6*(psi2 - psi3)),2);
                cC_4p2_p4_pow2[cbin] += pow(TMath::Cos(4*(psi2 - psi4)),2);
                cC_6p2_p6_pow2[cbin] += pow(TMath::Cos(6*(psi2 - psi6)),2);
                cC_6p3_p6_pow2[cbin] += pow(TMath::Cos(6*(psi3 - psi6)),2);
                cC_10p2_p5_pow2[cbin] += pow(TMath::Cos(10*(psi2 - psi5)),2);
                cC_15p3_p5_pow2[cbin] += pow(TMath::Cos(15*(psi3 - psi5)),2);

                cC_2p1_p2_pow2[cbin] += pow(TMath::Cos(2*(psi1 - psi2)),2);
                cC_3p1_p3_pow2[cbin] += pow(TMath::Cos(3*(psi1 - psi3)),2);
                cC_4p1_p4_pow2[cbin] += pow(TMath::Cos(4*(psi1 - psi4)),2);
                cC_5p1_p5_pow2[cbin] += pow(TMath::Cos(5*(psi1 - psi5)),2);
                cC_6p1_p6_pow2[cbin] += pow(TMath::Cos(6*(psi1 - psi6)),2);

                cC_1p1_2p2_3p3_pow2[cbin] += pow(TMath::Cos(1*psi1 + 2*psi2 - 3*psi3),2);
                cC_4p1_2p2_6p3_pow2[cbin] += pow(TMath::Cos(4*psi1 + 2*psi2 - 6*psi3),2);
                cC_2p1_4p2_6p3_pow2[cbin] += pow(TMath::Cos(2*psi1 + 4*psi2 - 6*psi3),2);
                cC_5p1_2p2_3p3_pow2[cbin] += pow(TMath::Cos(5*psi1 - 2*psi2 - 3*psi3),2);
                cC_2p1_2p2_4p4_pow2[cbin] += pow(TMath::Cos(2*psi1 + 2*psi2 - 4*psi4),2);
                cC_2p1_6p2_4p4_pow2[cbin] += pow(TMath::Cos(2*psi1 - 6*psi2 + 4*psi4),2);
                cC_2p1_6p2_8p4_pow2[cbin] += pow(TMath::Cos(2*psi1 + 6*psi2 - 8*psi4),2);

                cC_1p1_3p3_4p4_pow2[cbin] += pow(TMath::Cos(1*psi1 + 3*psi3 - 4*psi4),2);
                cC_2p1_6p3_4p4_pow2[cbin] += pow(TMath::Cos(2*psi1 - 6*psi3 + 4*psi4),2);
                cC_2p1_4p2_5p5_pow2[cbin] += pow(TMath::Cos(1*psi1 + 4*psi2 - 5*psi5),2);
                cC_3p1_2p2_5p5_pow2[cbin] += pow(TMath::Cos(3*psi1 + 2*psi2 - 5*psi5),2);
                cC_1p1_6p3_5p5_pow2[cbin] += pow(TMath::Cos(1*psi1 - 6*psi3 + 5*psi5),2);
                cC_2p1_3p3_5p5_pow2[cbin] += pow(TMath::Cos(2*psi1 + 3*psi3 - 5*psi5),2);
                cC_4p1_9p3_5p5_pow2[cbin] += pow(TMath::Cos(4*psi1 - 9*psi3 + 5*psi5),2);

                cC_cnt[cbin]++;
            }
        } // end centrality loop
    } // end event loop

    cout << "Finished event loop... " << endl;

    for (int cbin = 0; cbin<ncentbins; cbin++) {
        avNpart[cbin]/=avNpartCnt[cbin];
        for (int iorder = 1; iorder<=6; iorder++) {
            ecc[iorder][cbin]/=eccCnt[iorder][cbin];
            ecc_pow2[iorder][cbin]/=eccCnt[iorder][cbin];
            ecc_err[iorder][cbin] = sqrt( ecc_pow2[iorder][cbin] - pow(ecc[iorder][cbin],2) )/sqrt( eccCnt[iorder][cbin] );
            S[iorder][cbin]/=eccCnt[iorder][cbin];
            S_pow2[iorder][cbin]/=eccCnt[iorder][cbin];
            S_err[iorder][cbin] = sqrt( S_pow2[iorder][cbin] - pow(S[iorder][cbin],2) )/sqrt( eccCnt[iorder][cbin] );

            eccC[iorder][cbin]/=eccCCnt[iorder][cbin];
            eccC_pow2[iorder][cbin]/=eccCCnt[iorder][cbin];
            eccC_err[iorder][cbin] = sqrt( eccC_pow2[iorder][cbin] - pow(eccC[iorder][cbin],2) )/sqrt( eccCCnt[iorder][cbin] );
        }
        //-- moment calculation
        ecc_2_1[cbin]/=ecc_2_1_cnt[cbin];
        ecc_2_1_pow2[cbin]/=ecc_2_1_cnt[cbin];
        ecc_3_1[cbin]/=ecc_3_1_cnt[cbin];
        ecc_3_1_pow2[cbin]/=ecc_3_1_cnt[cbin];
        ecc_3_2[cbin]/=ecc_3_2_cnt[cbin];
        ecc_3_2_pow2[cbin]/=ecc_3_2_cnt[cbin];
        ecc_4_2[cbin]/=ecc_4_2_cnt[cbin];
        ecc_4_2_pow2[cbin]/=ecc_4_2_cnt[cbin];

        S_2_1[cbin]/=ecc_2_1_cnt[cbin];
        S_2_1_pow2[cbin]/=ecc_2_1_cnt[cbin];
        S_3_1[cbin]/=ecc_3_1_cnt[cbin];
        S_3_1_pow2[cbin]/=ecc_3_1_cnt[cbin];
        S_3_2[cbin]/=ecc_3_2_cnt[cbin];
        S_3_2_pow2[cbin]/=ecc_3_2_cnt[cbin];
        S_4_2[cbin]/=ecc_4_2_cnt[cbin];
        S_4_2_pow2[cbin]/=ecc_4_2_cnt[cbin];

        ecc_2_1_err[cbin] = sqrt( ecc_2_1_pow2[cbin] - pow(ecc_2_1[cbin],2) )/sqrt( ecc_2_1_cnt[cbin] );
        ecc_3_1_err[cbin] = sqrt( ecc_3_1_pow2[cbin] - pow(ecc_3_1[cbin],2) )/sqrt( ecc_3_1_cnt[cbin] );
        ecc_3_2_err[cbin] = sqrt( ecc_3_2_pow2[cbin] - pow(ecc_3_2[cbin],2) )/sqrt( ecc_3_2_cnt[cbin] );
        ecc_4_2_err[cbin] = sqrt( ecc_4_2_pow2[cbin] - pow(ecc_4_2[cbin],2) )/sqrt( ecc_4_2_cnt[cbin] );

        S_2_1_err[cbin] = sqrt( S_2_1_pow2[cbin] - pow(S_2_1[cbin],2) )/sqrt( ecc_2_1_cnt[cbin] );
        S_3_1_err[cbin] = sqrt( S_3_1_pow2[cbin] - pow(S_3_1[cbin],2) )/sqrt( ecc_3_1_cnt[cbin] );
        S_3_2_err[cbin] = sqrt( S_3_2_pow2[cbin] - pow(S_3_2[cbin],2) )/sqrt( ecc_3_2_cnt[cbin] );
        S_4_2_err[cbin] = sqrt( S_4_2_pow2[cbin] - pow(S_4_2[cbin],2) )/sqrt( ecc_4_2_cnt[cbin] );

        c_4_p4_p2[cbin]/=c_cnt[cbin];
        c_8_p4_p2[cbin]/=c_cnt[cbin];
        c_12_p4_p2[cbin]/=c_cnt[cbin];
        c_6_p3_p2[cbin]/=c_cnt[cbin];
        c_6_p2_p6[cbin]/=c_cnt[cbin];
        c_6_p3_p6[cbin]/=c_cnt[cbin];
        c_12_p3_p4[cbin]/=c_cnt[cbin];
        c_10_p2_p5[cbin]/=c_cnt[cbin];
        c_2p2_3p3_5p5[cbin]/=c_cnt[cbin];
        c_2p2_4p4_6p6[cbin]/=c_cnt[cbin];
        c_2p2_6p3_4p4[cbin]/=c_cnt[cbin];
        c_8p2_3p3_5p5[cbin]/=c_cnt[cbin];
        c_10p2_4p4_6p6[cbin]/=c_cnt[cbin];
        c_10p2_6p3_4p4[cbin]/=c_cnt[cbin];

        c_4_p4_p2_pow2[cbin]/=c_cnt[cbin];
        c_8_p4_p2_pow2[cbin]/=c_cnt[cbin];
        c_12_p4_p2_pow2[cbin]/=c_cnt[cbin];
        c_6_p3_p2_pow2[cbin]/=c_cnt[cbin];
        c_6_p2_p6_pow2[cbin]/=c_cnt[cbin];
        c_6_p3_p6_pow2[cbin]/=c_cnt[cbin];
        c_12_p3_p4_pow2[cbin]/=c_cnt[cbin];
        c_10_p2_p5_pow2[cbin]/=c_cnt[cbin];
        c_2p2_3p3_5p5_pow2[cbin]/=c_cnt[cbin];
        c_2p2_4p4_6p6_pow2[cbin]/=c_cnt[cbin];
        c_2p2_6p3_4p4_pow2[cbin]/=c_cnt[cbin];
        c_8p2_3p3_5p5_pow2[cbin]/=c_cnt[cbin];
        c_10p2_4p4_6p6_pow2[cbin]/=c_cnt[cbin];
        c_10p2_6p3_4p4_pow2[cbin]/=c_cnt[cbin];

        c_4_p4_p2_err[cbin] = sqrt( c_4_p4_p2_pow2[cbin] - pow(c_4_p4_p2[cbin],2) )/sqrt( c_cnt[cbin] );
        c_8_p4_p2_err[cbin] = sqrt( c_8_p4_p2_pow2[cbin] - pow(c_8_p4_p2[cbin],2) )/sqrt( c_cnt[cbin] );
        c_12_p4_p2_err[cbin] = sqrt( c_12_p4_p2_pow2[cbin] - pow(c_12_p4_p2[cbin],2) )/sqrt( c_cnt[cbin] );
        c_6_p3_p2_err[cbin] = sqrt( c_6_p3_p2_pow2[cbin] - pow(c_6_p3_p2[cbin],2) )/sqrt( c_cnt[cbin] );
        c_6_p2_p6_err[cbin] = sqrt( c_6_p2_p6_pow2[cbin] - pow(c_6_p2_p6[cbin],2) )/sqrt( c_cnt[cbin] );
        c_6_p3_p6_err[cbin] = sqrt( c_6_p3_p6_pow2[cbin] - pow(c_6_p3_p6[cbin],2) )/sqrt( c_cnt[cbin] );
        c_12_p3_p4_err[cbin] = sqrt( c_12_p3_p4_pow2[cbin] - pow(c_12_p3_p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_10_p2_p5_err[cbin] = sqrt( c_10_p2_p5_pow2[cbin] - pow(c_10_p2_p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p2_3p3_5p5_err[cbin] = sqrt( c_2p2_3p3_5p5_pow2[cbin] - pow(c_2p2_3p3_5p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p2_4p4_6p6_err[cbin] = sqrt( c_2p2_4p4_6p6_pow2[cbin] - pow(c_2p2_4p4_6p6[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p2_6p3_4p4_err[cbin] = sqrt( c_2p2_6p3_4p4_pow2[cbin] - pow(c_2p2_6p3_4p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_8p2_3p3_5p5_err[cbin] = sqrt( c_8p2_3p3_5p5_pow2[cbin] - pow(c_8p2_3p3_5p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_10p2_4p4_6p6_err[cbin] = sqrt( c_10p2_4p4_6p6_pow2[cbin] - pow(c_10p2_4p4_6p6[cbin],2) )/sqrt( c_cnt[cbin] );
        c_10p2_6p3_4p4_err[cbin] = sqrt( c_10p2_6p3_4p4_pow2[cbin] - pow(c_10p2_6p3_4p4[cbin],2) )/sqrt( c_cnt[cbin] );

        c_6_p3_p2[cbin]*=10;
        c_6_p3_p2_err[cbin]*=10;
        c_10_p2_p5[cbin]*=10;
        c_10_p2_p5_err[cbin]*=10;

        c_6p2_p3[cbin]/=c_cnt[cbin];
        c_4p2_p4[cbin]/=c_cnt[cbin];
        c_6p2_p6[cbin]/=c_cnt[cbin];
        c_6p3_p6[cbin]/=c_cnt[cbin];
        c_10p2_p5[cbin]/=c_cnt[cbin];
        c_15p3_p5[cbin]/=c_cnt[cbin];
        c_2p1_p2[cbin]/=c_cnt[cbin];
        c_3p1_p3[cbin]/=c_cnt[cbin];
        c_4p1_p4[cbin]/=c_cnt[cbin];
        c_5p1_p5[cbin]/=c_cnt[cbin];
        c_6p1_p6[cbin]/=c_cnt[cbin];
        c_1p1_2p2_3p3[cbin]/=c_cnt[cbin];
        c_4p1_2p2_6p3[cbin]/=c_cnt[cbin];
        c_2p1_4p2_6p3[cbin]/=c_cnt[cbin];
        c_5p1_2p2_3p3[cbin]/=c_cnt[cbin];
        c_2p1_2p2_4p4[cbin]/=c_cnt[cbin];
        c_2p1_6p2_4p4[cbin]/=c_cnt[cbin];
        c_2p1_6p2_8p4[cbin]/=c_cnt[cbin];
        c_1p1_3p3_4p4[cbin]/=c_cnt[cbin];
        c_2p1_6p3_4p4[cbin]/=c_cnt[cbin];
        c_2p1_4p2_5p5[cbin]/=c_cnt[cbin];
        c_3p1_2p2_5p5[cbin]/=c_cnt[cbin];
        c_1p1_6p3_5p5[cbin]/=c_cnt[cbin];
        c_2p1_3p3_5p5[cbin]/=c_cnt[cbin];
        c_4p1_9p3_5p5[cbin]/=c_cnt[cbin];

        c_6p2_p3_pow2[cbin]/=c_cnt[cbin];
        c_4p2_p4_pow2[cbin]/=c_cnt[cbin];
        c_6p2_p6_pow2[cbin]/=c_cnt[cbin];
        c_6p3_p6_pow2[cbin]/=c_cnt[cbin];
        c_10p2_p5_pow2[cbin]/=c_cnt[cbin];
        c_15p3_p5_pow2[cbin]/=c_cnt[cbin];
        c_2p1_p2_pow2[cbin]/=c_cnt[cbin];
        c_3p1_p3_pow2[cbin]/=c_cnt[cbin];
        c_4p1_p4_pow2[cbin]/=c_cnt[cbin];
        c_5p1_p5_pow2[cbin]/=c_cnt[cbin];
        c_6p1_p6_pow2[cbin]/=c_cnt[cbin];
        c_1p1_2p2_3p3_pow2[cbin]/=c_cnt[cbin];
        c_4p1_2p2_6p3_pow2[cbin]/=c_cnt[cbin];
        c_2p1_4p2_6p3_pow2[cbin]/=c_cnt[cbin];
        c_5p1_2p2_3p3_pow2[cbin]/=c_cnt[cbin];
        c_2p1_2p2_4p4_pow2[cbin]/=c_cnt[cbin];
        c_2p1_6p2_4p4_pow2[cbin]/=c_cnt[cbin];
        c_2p1_6p2_8p4_pow2[cbin]/=c_cnt[cbin];
        c_1p1_3p3_4p4_pow2[cbin]/=c_cnt[cbin];
        c_2p1_6p3_4p4_pow2[cbin]/=c_cnt[cbin];
        c_2p1_4p2_5p5_pow2[cbin]/=c_cnt[cbin];
        c_3p1_2p2_5p5_pow2[cbin]/=c_cnt[cbin];
        c_1p1_6p3_5p5_pow2[cbin]/=c_cnt[cbin];
        c_2p1_3p3_5p5_pow2[cbin]/=c_cnt[cbin];
        c_4p1_9p3_5p5_pow2[cbin]/=c_cnt[cbin];

        c_6p2_p3_err[cbin] = sqrt( c_6p2_p3_pow2[cbin] - pow(c_6p2_p3[cbin],2) )/sqrt( c_cnt[cbin] );
        c_4p2_p4_err[cbin] = sqrt( c_4p2_p4_pow2[cbin] - pow(c_4p2_p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_6p2_p6_err[cbin] = sqrt( c_6p2_p6_pow2[cbin] - pow(c_6p2_p6[cbin],2) )/sqrt( c_cnt[cbin] );
        c_6p3_p6_err[cbin] = sqrt( c_6p3_p6_pow2[cbin] - pow(c_6p3_p6[cbin],2) )/sqrt( c_cnt[cbin] );
        c_10p2_p5_err[cbin] = sqrt( c_10p2_p5_pow2[cbin] - pow(c_10p2_p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_15p3_p5_err[cbin] = sqrt( c_15p3_p5_pow2[cbin] - pow(c_15p3_p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_p2_err[cbin] = sqrt( c_2p1_p2_pow2[cbin] - pow(c_2p1_p2[cbin],2) )/sqrt( c_cnt[cbin] );
        c_3p1_p3_err[cbin] = sqrt( c_3p1_p3_pow2[cbin] - pow(c_3p1_p3[cbin],2) )/sqrt( c_cnt[cbin] );
        c_4p1_p4_err[cbin] = sqrt( c_4p1_p4_pow2[cbin] - pow(c_4p1_p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_5p1_p5_err[cbin] = sqrt( c_5p1_p5_pow2[cbin] - pow(c_5p1_p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_6p1_p6_err[cbin] = sqrt( c_6p1_p6_pow2[cbin] - pow(c_6p1_p6[cbin],2) )/sqrt( c_cnt[cbin] );
        c_1p1_2p2_3p3_err[cbin] = sqrt( c_1p1_2p2_3p3_pow2[cbin] - pow(c_1p1_2p2_3p3[cbin],2) )/sqrt( c_cnt[cbin] );
        c_4p1_2p2_6p3_err[cbin] = sqrt( c_4p1_2p2_6p3_pow2[cbin] - pow(c_4p1_2p2_6p3[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_4p2_6p3_err[cbin] = sqrt( c_2p1_4p2_6p3_pow2[cbin] - pow(c_2p1_4p2_6p3[cbin],2) )/sqrt( c_cnt[cbin] );
        c_5p1_2p2_3p3_err[cbin] = sqrt( c_5p1_2p2_3p3_pow2[cbin] - pow(c_5p1_2p2_3p3[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_2p2_4p4_err[cbin] = sqrt( c_2p1_2p2_4p4_pow2[cbin] - pow(c_2p1_2p2_4p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_6p2_4p4_err[cbin] = sqrt( c_2p1_6p2_4p4_pow2[cbin] - pow(c_2p1_6p2_4p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_6p2_8p4_err[cbin] = sqrt( c_2p1_6p2_8p4_pow2[cbin] - pow(c_2p1_6p2_8p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_1p1_3p3_4p4_err[cbin] = sqrt( c_1p1_3p3_4p4_pow2[cbin] - pow(c_1p1_3p3_4p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_6p3_4p4_err[cbin] = sqrt( c_2p1_6p3_4p4_pow2[cbin] - pow(c_2p1_6p3_4p4[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_4p2_5p5_err[cbin] = sqrt( c_2p1_4p2_5p5_pow2[cbin] - pow(c_2p1_4p2_5p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_3p1_2p2_5p5_err[cbin] = sqrt( c_3p1_2p2_5p5_pow2[cbin] - pow(c_3p1_2p2_5p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_1p1_6p3_5p5_err[cbin] = sqrt( c_1p1_6p3_5p5_pow2[cbin] - pow(c_1p1_6p3_5p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_2p1_3p3_5p5_err[cbin] = sqrt( c_2p1_3p3_5p5_pow2[cbin] - pow(c_2p1_3p3_5p5[cbin],2) )/sqrt( c_cnt[cbin] );
        c_4p1_9p3_5p5_err[cbin] = sqrt( c_4p1_9p3_5p5_pow2[cbin] - pow(c_4p1_9p3_5p5[cbin],2) )/sqrt( c_cnt[cbin] );


        //-- cumulant calculation
        eccC_2_1[cbin]/=eccC_2_1_cnt[cbin];
        eccC_2_1_pow2[cbin]/=eccC_2_1_cnt[cbin];
        eccC_3_1[cbin]/=eccC_3_1_cnt[cbin];
        eccC_3_1_pow2[cbin]/=eccC_3_1_cnt[cbin];
        eccC_3_2[cbin]/=eccC_3_2_cnt[cbin];
        eccC_3_2_pow2[cbin]/=eccC_3_2_cnt[cbin];
        eccC_4_2[cbin]/=eccC_4_2_cnt[cbin];
        eccC_4_2_pow2[cbin]/=eccC_4_2_cnt[cbin];

        eccC_2_1_err[cbin] = sqrt( eccC_2_1_pow2[cbin] - pow(eccC_2_1[cbin],2) )/sqrt( eccC_2_1_cnt[cbin] );
        eccC_3_1_err[cbin] = sqrt( eccC_3_1_pow2[cbin] - pow(eccC_3_1[cbin],2) )/sqrt( eccC_3_1_cnt[cbin] );
        eccC_3_2_err[cbin] = sqrt( eccC_3_2_pow2[cbin] - pow(eccC_3_2[cbin],2) )/sqrt( eccC_3_2_cnt[cbin] );
        eccC_4_2_err[cbin] = sqrt( eccC_4_2_pow2[cbin] - pow(eccC_4_2[cbin],2) )/sqrt( eccC_4_2_cnt[cbin] );

        cC_4_p4_p2[cbin]/=cC_cnt[cbin];
        cC_8_p4_p2[cbin]/=cC_cnt[cbin];
        cC_12_p4_p2[cbin]/=cC_cnt[cbin];
        cC_6_p3_p2[cbin]/=cC_cnt[cbin];
        cC_6_p2_p6[cbin]/=cC_cnt[cbin];
        cC_6_p3_p6[cbin]/=cC_cnt[cbin];
        cC_12_p3_p4[cbin]/=cC_cnt[cbin];
        cC_10_p2_p5[cbin]/=cC_cnt[cbin];
        cC_2p2_3p3_5p5[cbin]/=cC_cnt[cbin];
        cC_2p2_4p4_6p6[cbin]/=cC_cnt[cbin];
        cC_2p2_6p3_4p4[cbin]/=cC_cnt[cbin];
        cC_8p2_3p3_5p5[cbin]/=cC_cnt[cbin];
        cC_10p2_4p4_6p6[cbin]/=cC_cnt[cbin];
        cC_10p2_6p3_4p4[cbin]/=cC_cnt[cbin];

        cC_4_p4_p2_pow2[cbin]/=cC_cnt[cbin];
        cC_8_p4_p2_pow2[cbin]/=cC_cnt[cbin];
        cC_12_p4_p2_pow2[cbin]/=cC_cnt[cbin];
        cC_6_p3_p2_pow2[cbin]/=cC_cnt[cbin];
        cC_6_p2_p6_pow2[cbin]/=cC_cnt[cbin];
        cC_6_p3_p6_pow2[cbin]/=cC_cnt[cbin];
        cC_12_p3_p4_pow2[cbin]/=cC_cnt[cbin];
        cC_10_p2_p5_pow2[cbin]/=cC_cnt[cbin];
        cC_2p2_3p3_5p5_pow2[cbin]/=cC_cnt[cbin];
        cC_2p2_4p4_6p6_pow2[cbin]/=cC_cnt[cbin];
        cC_2p2_6p3_4p4_pow2[cbin]/=cC_cnt[cbin];
        cC_8p2_3p3_5p5_pow2[cbin]/=cC_cnt[cbin];
        cC_10p2_4p4_6p6_pow2[cbin]/=cC_cnt[cbin];
        cC_10p2_6p3_4p4_pow2[cbin]/=cC_cnt[cbin];

        cC_4_p4_p2_err[cbin] = sqrt( cC_4_p4_p2_pow2[cbin] - pow(cC_4_p4_p2[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_8_p4_p2_err[cbin] = sqrt( cC_8_p4_p2_pow2[cbin] - pow(cC_8_p4_p2[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_12_p4_p2_err[cbin] = sqrt( cC_12_p4_p2_pow2[cbin] - pow(cC_12_p4_p2[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_6_p3_p2_err[cbin] = sqrt( cC_6_p3_p2_pow2[cbin] - pow(cC_6_p3_p2[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_6_p2_p6_err[cbin] = sqrt( cC_6_p2_p6_pow2[cbin] - pow(cC_6_p2_p6[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_6_p3_p6_err[cbin] = sqrt( cC_6_p3_p6_pow2[cbin] - pow(cC_6_p3_p6[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_12_p3_p4_err[cbin] = sqrt( cC_12_p3_p4_pow2[cbin] - pow(cC_12_p3_p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_10_p2_p5_err[cbin] = sqrt( cC_10_p2_p5_pow2[cbin] - pow(cC_10_p2_p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p2_3p3_5p5_err[cbin] = sqrt( cC_2p2_3p3_5p5_pow2[cbin] - pow(cC_2p2_3p3_5p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p2_4p4_6p6_err[cbin] = sqrt( cC_2p2_4p4_6p6_pow2[cbin] - pow(cC_2p2_4p4_6p6[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p2_6p3_4p4_err[cbin] = sqrt( cC_2p2_6p3_4p4_pow2[cbin] - pow(cC_2p2_6p3_4p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_8p2_3p3_5p5_err[cbin] = sqrt( cC_8p2_3p3_5p5_pow2[cbin] - pow(cC_8p2_3p3_5p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_10p2_4p4_6p6_err[cbin] = sqrt( cC_10p2_4p4_6p6_pow2[cbin] - pow(cC_10p2_4p4_6p6[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_10p2_6p3_4p4_err[cbin] = sqrt( cC_10p2_6p3_4p4_pow2[cbin] - pow(cC_10p2_6p3_4p4[cbin],2) )/sqrt( cC_cnt[cbin] );

        cC_6_p3_p2[cbin]*=10;
        cC_6_p3_p2_err[cbin]*=10;
        cC_10_p2_p5[cbin]*=10;
        cC_10_p2_p5_err[cbin]*=10;

        cC_6p2_p3[cbin]/=cC_cnt[cbin];
        cC_4p2_p4[cbin]/=cC_cnt[cbin];
        cC_6p2_p6[cbin]/=cC_cnt[cbin];
        cC_6p3_p6[cbin]/=cC_cnt[cbin];
        cC_10p2_p5[cbin]/=cC_cnt[cbin];
        cC_15p3_p5[cbin]/=cC_cnt[cbin];
        cC_2p1_p2[cbin]/=cC_cnt[cbin];
        cC_3p1_p3[cbin]/=cC_cnt[cbin];
        cC_4p1_p4[cbin]/=cC_cnt[cbin];
        cC_5p1_p5[cbin]/=cC_cnt[cbin];
        cC_6p1_p6[cbin]/=cC_cnt[cbin];
        cC_1p1_2p2_3p3[cbin]/=cC_cnt[cbin];
        cC_4p1_2p2_6p3[cbin]/=cC_cnt[cbin];
        cC_2p1_4p2_6p3[cbin]/=cC_cnt[cbin];
        cC_5p1_2p2_3p3[cbin]/=cC_cnt[cbin];
        cC_2p1_2p2_4p4[cbin]/=cC_cnt[cbin];
        cC_2p1_6p2_4p4[cbin]/=cC_cnt[cbin];
        cC_2p1_6p2_8p4[cbin]/=cC_cnt[cbin];
        cC_1p1_3p3_4p4[cbin]/=cC_cnt[cbin];
        cC_2p1_6p3_4p4[cbin]/=cC_cnt[cbin];
        cC_2p1_4p2_5p5[cbin]/=cC_cnt[cbin];
        cC_3p1_2p2_5p5[cbin]/=cC_cnt[cbin];
        cC_1p1_6p3_5p5[cbin]/=cC_cnt[cbin];
        cC_2p1_3p3_5p5[cbin]/=cC_cnt[cbin];
        cC_4p1_9p3_5p5[cbin]/=cC_cnt[cbin];

        cC_6p2_p3_pow2[cbin]/=cC_cnt[cbin];
        cC_4p2_p4_pow2[cbin]/=cC_cnt[cbin];
        cC_6p2_p6_pow2[cbin]/=cC_cnt[cbin];
        cC_6p3_p6_pow2[cbin]/=cC_cnt[cbin];
        cC_10p2_p5_pow2[cbin]/=cC_cnt[cbin];
        cC_15p3_p5_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_p2_pow2[cbin]/=cC_cnt[cbin];
        cC_3p1_p3_pow2[cbin]/=cC_cnt[cbin];
        cC_4p1_p4_pow2[cbin]/=cC_cnt[cbin];
        cC_5p1_p5_pow2[cbin]/=cC_cnt[cbin];
        cC_6p1_p6_pow2[cbin]/=cC_cnt[cbin];
        cC_1p1_2p2_3p3_pow2[cbin]/=cC_cnt[cbin];
        cC_4p1_2p2_6p3_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_4p2_6p3_pow2[cbin]/=cC_cnt[cbin];
        cC_5p1_2p2_3p3_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_2p2_4p4_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_6p2_4p4_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_6p2_8p4_pow2[cbin]/=cC_cnt[cbin];
        cC_1p1_3p3_4p4_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_6p3_4p4_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_4p2_5p5_pow2[cbin]/=cC_cnt[cbin];
        cC_3p1_2p2_5p5_pow2[cbin]/=cC_cnt[cbin];
        cC_1p1_6p3_5p5_pow2[cbin]/=cC_cnt[cbin];
        cC_2p1_3p3_5p5_pow2[cbin]/=cC_cnt[cbin];
        cC_4p1_9p3_5p5_pow2[cbin]/=cC_cnt[cbin];

        cC_6p2_p3_err[cbin] = sqrt( cC_6p2_p3_pow2[cbin] - pow(cC_6p2_p3[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_4p2_p4_err[cbin] = sqrt( cC_4p2_p4_pow2[cbin] - pow(cC_4p2_p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_6p2_p6_err[cbin] = sqrt( cC_6p2_p6_pow2[cbin] - pow(cC_6p2_p6[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_6p3_p6_err[cbin] = sqrt( cC_6p3_p6_pow2[cbin] - pow(cC_6p3_p6[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_10p2_p5_err[cbin] = sqrt( cC_10p2_p5_pow2[cbin] - pow(cC_10p2_p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_15p3_p5_err[cbin] = sqrt( cC_15p3_p5_pow2[cbin] - pow(cC_15p3_p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_p2_err[cbin] = sqrt( cC_2p1_p2_pow2[cbin] - pow(cC_2p1_p2[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_3p1_p3_err[cbin] = sqrt( cC_3p1_p3_pow2[cbin] - pow(cC_3p1_p3[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_4p1_p4_err[cbin] = sqrt( cC_4p1_p4_pow2[cbin] - pow(cC_4p1_p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_5p1_p5_err[cbin] = sqrt( cC_5p1_p5_pow2[cbin] - pow(cC_5p1_p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_6p1_p6_err[cbin] = sqrt( cC_6p1_p6_pow2[cbin] - pow(cC_6p1_p6[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_1p1_2p2_3p3_err[cbin] = sqrt( cC_1p1_2p2_3p3_pow2[cbin] - pow(cC_1p1_2p2_3p3[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_4p1_2p2_6p3_err[cbin] = sqrt( cC_4p1_2p2_6p3_pow2[cbin] - pow(cC_4p1_2p2_6p3[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_4p2_6p3_err[cbin] = sqrt( cC_2p1_4p2_6p3_pow2[cbin] - pow(cC_2p1_4p2_6p3[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_5p1_2p2_3p3_err[cbin] = sqrt( cC_5p1_2p2_3p3_pow2[cbin] - pow(cC_5p1_2p2_3p3[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_2p2_4p4_err[cbin] = sqrt( cC_2p1_2p2_4p4_pow2[cbin] - pow(cC_2p1_2p2_4p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_6p2_4p4_err[cbin] = sqrt( cC_2p1_6p2_4p4_pow2[cbin] - pow(cC_2p1_6p2_4p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_6p2_8p4_err[cbin] = sqrt( cC_2p1_6p2_8p4_pow2[cbin] - pow(cC_2p1_6p2_8p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_1p1_3p3_4p4_err[cbin] = sqrt( cC_1p1_3p3_4p4_pow2[cbin] - pow(cC_1p1_3p3_4p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_6p3_4p4_err[cbin] = sqrt( cC_2p1_6p3_4p4_pow2[cbin] - pow(cC_2p1_6p3_4p4[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_4p2_5p5_err[cbin] = sqrt( cC_2p1_4p2_5p5_pow2[cbin] - pow(cC_2p1_4p2_5p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_3p1_2p2_5p5_err[cbin] = sqrt( cC_3p1_2p2_5p5_pow2[cbin] - pow(cC_3p1_2p2_5p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_1p1_6p3_5p5_err[cbin] = sqrt( cC_1p1_6p3_5p5_pow2[cbin] - pow(cC_1p1_6p3_5p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_2p1_3p3_5p5_err[cbin] = sqrt( cC_2p1_3p3_5p5_pow2[cbin] - pow(cC_2p1_3p3_5p5[cbin],2) )/sqrt( cC_cnt[cbin] );
        cC_4p1_9p3_5p5_err[cbin] = sqrt( cC_4p1_9p3_5p5_pow2[cbin] - pow(cC_4p1_9p3_5p5[cbin],2) )/sqrt( cC_cnt[cbin] );

    }

    for (int iorder = 1; iorder<=6; iorder++) {
        g[iorder] = new TGraphErrors(ncentbins, avNpart, ecc[iorder], 0, ecc_err[iorder]);
        g[iorder]->SetMarkerStyle(21);
        g[iorder]->SetMarkerColor(colors[iorder]);
        g[iorder]->SetLineColor(colors[iorder]);

        sg[iorder] = new TGraphErrors(ncentbins, avNpart, S[iorder], 0, S_err[iorder]);
        sg[iorder]->SetMarkerStyle(21);
        sg[iorder]->SetMarkerColor(colors[iorder]);
        sg[iorder]->SetLineColor(colors[iorder]);

        gC[iorder] = new TGraphErrors(ncentbins, avNpart, eccC[iorder], 0, eccC_err[iorder]);
        gC[iorder]->SetMarkerStyle(21);
        gC[iorder]->SetMarkerColor(colors[iorder]);
        gC[iorder]->SetLineColor(colors[iorder]);
    }

    g_2_1 = new TGraphErrors(ncentbins, avNpart, ecc_2_1, 0, ecc_2_1_err);
    g_3_1 = new TGraphErrors(ncentbins, avNpart, ecc_3_1, 0, ecc_3_1_err);
    g_3_2 = new TGraphErrors(ncentbins, avNpart, ecc_3_2, 0, ecc_3_2_err);
    g_4_2 = new TGraphErrors(ncentbins, avNpart, ecc_4_2, 0, ecc_4_2_err);

    sg_2_1 = new TGraphErrors(ncentbins, avNpart, S_2_1, 0, S_2_1_err);
    sg_3_1 = new TGraphErrors(ncentbins, avNpart, S_3_1, 0, S_3_1_err);
    sg_3_2 = new TGraphErrors(ncentbins, avNpart, S_3_2, 0, S_3_2_err);
    sg_4_2 = new TGraphErrors(ncentbins, avNpart, S_4_2, 0, S_4_2_err);

    cg_4_p4_p2 = new TGraphErrors(ncentbins, avNpart, c_4_p4_p2, 0, c_4_p4_p2_err);
    cg_8_p4_p2 = new TGraphErrors(ncentbins, avNpart, c_8_p4_p2, 0, c_8_p4_p2_err);
    cg_12_p4_p2 = new TGraphErrors(ncentbins, avNpart, c_12_p4_p2, 0, c_12_p4_p2_err);
    cg_6_p3_p2 = new TGraphErrors(ncentbins, avNpart, c_6_p3_p2, 0, c_6_p3_p2_err);
    cg_6_p2_p6 = new TGraphErrors(ncentbins, avNpart, c_6_p2_p6, 0, c_6_p2_p6_err);
    cg_6_p3_p6 = new TGraphErrors(ncentbins, avNpart, c_6_p3_p6, 0, c_6_p3_p6_err);
    cg_12_p3_p4 = new TGraphErrors(ncentbins, avNpart, c_12_p3_p4, 0, c_12_p3_p4_err);
    cg_10_p2_p5 = new TGraphErrors(ncentbins, avNpart, c_10_p2_p5, 0, c_10_p2_p5_err);

    cg_2p2_3p3_5p5 = new TGraphErrors(ncentbins, avNpart, c_2p2_3p3_5p5, 0, c_2p2_3p3_5p5_err);
    cg_2p2_4p4_6p6 = new TGraphErrors(ncentbins, avNpart, c_2p2_4p4_6p6, 0, c_2p2_4p4_6p6_err);
    cg_2p2_6p3_4p4 = new TGraphErrors(ncentbins, avNpart, c_2p2_6p3_4p4, 0, c_2p2_6p3_4p4_err);
    cg_8p2_3p3_5p5 = new TGraphErrors(ncentbins, avNpart, c_8p2_3p3_5p5, 0, c_8p2_3p3_5p5_err);
    cg_10p2_4p4_6p6 = new TGraphErrors(ncentbins, avNpart, c_10p2_4p4_6p6, 0, c_10p2_4p4_6p6_err);
    cg_10p2_6p3_4p4 = new TGraphErrors(ncentbins, avNpart, c_10p2_6p3_4p4, 0, c_10p2_6p3_4p4_err);

    cg_6p2_p3 = new TGraphErrors(ncentbins, avNpart, c_6p2_p3, 0, c_6p2_p3_err);
    cg_4p2_p4 = new TGraphErrors(ncentbins, avNpart, c_4p2_p4, 0, c_4p2_p4_err);
    cg_6p2_p6 = new TGraphErrors(ncentbins, avNpart, c_6p2_p6, 0, c_6p2_p6_err);
    cg_6p3_p6 = new TGraphErrors(ncentbins, avNpart, c_6p3_p6, 0, c_6p3_p6_err);
    cg_10p2_p5 = new TGraphErrors(ncentbins, avNpart, c_10p2_p5, 0, c_10p2_p5_err);
    cg_15p3_p5 = new TGraphErrors(ncentbins, avNpart, c_15p3_p5, 0, c_15p3_p5_err);

    cg_2p1_p2 = new TGraphErrors(ncentbins, avNpart, c_2p1_p2, 0, c_2p1_p2_err);
    cg_3p1_p3 = new TGraphErrors(ncentbins, avNpart, c_3p1_p3, 0, c_3p1_p3_err);
    cg_4p1_p4 = new TGraphErrors(ncentbins, avNpart, c_4p1_p4, 0, c_4p1_p4_err);
    cg_5p1_p5 = new TGraphErrors(ncentbins, avNpart, c_5p1_p5, 0, c_5p1_p5_err);
    cg_6p1_p6 = new TGraphErrors(ncentbins, avNpart, c_6p1_p6, 0, c_6p1_p6_err);

    cg_1p1_2p2_3p3 = new TGraphErrors(ncentbins, avNpart, c_1p1_2p2_3p3, 0, c_1p1_2p2_3p3_err);
    cg_4p1_2p2_6p3 = new TGraphErrors(ncentbins, avNpart, c_4p1_2p2_6p3, 0, c_4p1_2p2_6p3_err);
    cg_2p1_4p2_6p3 = new TGraphErrors(ncentbins, avNpart, c_2p1_4p2_6p3, 0, c_2p1_4p2_6p3_err);
    cg_5p1_2p2_3p3 = new TGraphErrors(ncentbins, avNpart, c_5p1_2p2_3p3, 0, c_5p1_2p2_3p3_err);
    cg_2p1_2p2_4p4 = new TGraphErrors(ncentbins, avNpart, c_2p1_2p2_4p4, 0, c_2p1_2p2_4p4_err);
    cg_2p1_6p2_4p4 = new TGraphErrors(ncentbins, avNpart, c_2p1_6p2_4p4, 0, c_2p1_6p2_4p4_err);
    cg_2p1_6p2_8p4 = new TGraphErrors(ncentbins, avNpart, c_2p1_6p2_8p4, 0, c_2p1_6p2_8p4_err);

    cg_1p1_3p3_4p4 = new TGraphErrors(ncentbins, avNpart, c_1p1_3p3_4p4, 0, c_1p1_3p3_4p4_err);
    cg_2p1_6p3_4p4 = new TGraphErrors(ncentbins, avNpart, c_2p1_6p3_4p4, 0, c_2p1_6p3_4p4_err);
    cg_1p1_4p2_5p5 = new TGraphErrors(ncentbins, avNpart, c_2p1_4p2_5p5, 0, c_2p1_4p2_5p5_err);
    cg_3p1_2p2_5p5 = new TGraphErrors(ncentbins, avNpart, c_3p1_2p2_5p5, 0, c_3p1_2p2_5p5_err);
    cg_1p1_6p3_5p5 = new TGraphErrors(ncentbins, avNpart, c_1p1_6p3_5p5, 0, c_1p1_6p3_5p5_err);
    cg_2p1_3p3_5p5 = new TGraphErrors(ncentbins, avNpart, c_2p1_3p3_5p5, 0, c_2p1_3p3_5p5_err);
    cg_4p1_9p3_5p5 = new TGraphErrors(ncentbins, avNpart, c_4p1_9p3_5p5, 0, c_4p1_9p3_5p5_err);

    cgC_6p2_p3 = new TGraphErrors(ncentbins, avNpart, cC_6p2_p3, 0, cC_6p2_p3_err);
    cgC_4p2_p4 = new TGraphErrors(ncentbins, avNpart, cC_4p2_p4, 0, cC_4p2_p4_err);
    cgC_6p2_p6 = new TGraphErrors(ncentbins, avNpart, cC_6p2_p6, 0, cC_6p2_p6_err);
    cgC_6p3_p6 = new TGraphErrors(ncentbins, avNpart, cC_6p3_p6, 0, cC_6p3_p6_err);
    cgC_10p2_p5 = new TGraphErrors(ncentbins, avNpart, cC_10p2_p5, 0, cC_10p2_p5_err);
    cgC_15p3_p5 = new TGraphErrors(ncentbins, avNpart, cC_15p3_p5, 0, cC_15p3_p5_err);

    cgC_2p1_p2 = new TGraphErrors(ncentbins, avNpart, cC_2p1_p2, 0, cC_2p1_p2_err);
    cgC_3p1_p3 = new TGraphErrors(ncentbins, avNpart, cC_3p1_p3, 0, cC_3p1_p3_err);
    cgC_4p1_p4 = new TGraphErrors(ncentbins, avNpart, cC_4p1_p4, 0, cC_4p1_p4_err);
    cgC_5p1_p5 = new TGraphErrors(ncentbins, avNpart, cC_5p1_p5, 0, cC_5p1_p5_err);
    cgC_6p1_p6 = new TGraphErrors(ncentbins, avNpart, cC_6p1_p6, 0, cC_6p1_p6_err);

    cgC_1p1_2p2_3p3 = new TGraphErrors(ncentbins, avNpart, cC_1p1_2p2_3p3, 0, cC_1p1_2p2_3p3_err);
    cgC_4p1_2p2_6p3 = new TGraphErrors(ncentbins, avNpart, cC_4p1_2p2_6p3, 0, cC_4p1_2p2_6p3_err);
    cgC_2p1_4p2_6p3 = new TGraphErrors(ncentbins, avNpart, cC_2p1_4p2_6p3, 0, cC_2p1_4p2_6p3_err);
    cgC_5p1_2p2_3p3 = new TGraphErrors(ncentbins, avNpart, cC_5p1_2p2_3p3, 0, cC_5p1_2p2_3p3_err);
    cgC_2p1_2p2_4p4 = new TGraphErrors(ncentbins, avNpart, cC_2p1_2p2_4p4, 0, cC_2p1_2p2_4p4_err);
    cgC_2p1_6p2_4p4 = new TGraphErrors(ncentbins, avNpart, cC_2p1_6p2_4p4, 0, cC_2p1_6p2_4p4_err);
    cgC_2p1_6p2_8p4 = new TGraphErrors(ncentbins, avNpart, cC_2p1_6p2_8p4, 0, cC_2p1_6p2_8p4_err);

    cgC_1p1_3p3_4p4 = new TGraphErrors(ncentbins, avNpart, cC_1p1_3p3_4p4, 0, cC_1p1_3p3_4p4_err);
    cgC_2p1_6p3_4p4 = new TGraphErrors(ncentbins, avNpart, cC_2p1_6p3_4p4, 0, cC_2p1_6p3_4p4_err);
    cgC_1p1_4p2_5p5 = new TGraphErrors(ncentbins, avNpart, cC_2p1_4p2_5p5, 0, cC_2p1_4p2_5p5_err);
    cgC_3p1_2p2_5p5 = new TGraphErrors(ncentbins, avNpart, cC_3p1_2p2_5p5, 0, cC_3p1_2p2_5p5_err);
    cgC_1p1_6p3_5p5 = new TGraphErrors(ncentbins, avNpart, cC_1p1_6p3_5p5, 0, cC_1p1_6p3_5p5_err);
    cgC_2p1_3p3_5p5 = new TGraphErrors(ncentbins, avNpart, cC_2p1_3p3_5p5, 0, cC_2p1_3p3_5p5_err);
    cgC_4p1_9p3_5p5 = new TGraphErrors(ncentbins, avNpart, cC_4p1_9p3_5p5, 0, cC_4p1_9p3_5p5_err);


    //-- cumulants
    gC_2_1 = new TGraphErrors(ncentbins, avNpart, eccC_2_1, 0, eccC_2_1_err);
    gC_3_1 = new TGraphErrors(ncentbins, avNpart, eccC_3_1, 0, eccC_3_1_err);
    gC_3_2 = new TGraphErrors(ncentbins, avNpart, eccC_3_2, 0, eccC_3_2_err);
    gC_4_2 = new TGraphErrors(ncentbins, avNpart, eccC_4_2, 0, eccC_4_2_err);

    cgC_4_p4_p2 = new TGraphErrors(ncentbins, avNpart, cC_4_p4_p2, 0, cC_4_p4_p2_err);
    cgC_8_p4_p2 = new TGraphErrors(ncentbins, avNpart, cC_8_p4_p2, 0, cC_8_p4_p2_err);
    cgC_12_p4_p2 = new TGraphErrors(ncentbins, avNpart, cC_12_p4_p2, 0, cC_12_p4_p2_err);
    cgC_6_p3_p2 = new TGraphErrors(ncentbins, avNpart, cC_6_p3_p2, 0, cC_6_p3_p2_err);
    cgC_6_p2_p6 = new TGraphErrors(ncentbins, avNpart, cC_6_p2_p6, 0, cC_6_p2_p6_err);
    cgC_6_p3_p6 = new TGraphErrors(ncentbins, avNpart, cC_6_p3_p6, 0, cC_6_p3_p6_err);
    cgC_12_p3_p4 = new TGraphErrors(ncentbins, avNpart, cC_12_p3_p4, 0, cC_12_p3_p4_err);
    cgC_10_p2_p5 = new TGraphErrors(ncentbins, avNpart, cC_10_p2_p5, 0, cC_10_p2_p5_err);

    cgC_2p2_3p3_5p5 = new TGraphErrors(ncentbins, avNpart, cC_2p2_3p3_5p5, 0, cC_2p2_3p3_5p5_err);
    cgC_2p2_4p4_6p6 = new TGraphErrors(ncentbins, avNpart, cC_2p2_4p4_6p6, 0, cC_2p2_4p4_6p6_err);
    cgC_2p2_6p3_4p4 = new TGraphErrors(ncentbins, avNpart, cC_2p2_6p3_4p4, 0, cC_2p2_6p3_4p4_err);
    cgC_8p2_3p3_5p5 = new TGraphErrors(ncentbins, avNpart, cC_8p2_3p3_5p5, 0, cC_8p2_3p3_5p5_err);
    cgC_10p2_4p4_6p6 = new TGraphErrors(ncentbins, avNpart, cC_10p2_4p4_6p6, 0, cC_10p2_4p4_6p6_err);
    cgC_10p2_6p3_4p4 = new TGraphErrors(ncentbins, avNpart, cC_10p2_6p3_4p4, 0, cC_10p2_6p3_4p4_err);

    if (!PlotCent) {
        TString atag = "";
        if (minorAxes) atag+="_minorAxes";
        FILE * fout;
        FILE * foutS;
        FILE * foutc;
        FILE * foutNpart;
        fout = fopen(Form("results/glauberEcc%s.txt",atag.Data()),"w");
        foutS = fopen(Form("results/glauberS%s.txt",atag.Data()),"w");
        foutc = fopen(Form("results/glauberc%s.txt",atag.Data()),"w");
        foutNpart = fopen(Form("results/npart%s.txt",atag.Data()),"w");
        fprintf(fout,  "centmin\tcentmax\tnpart\te1\te2\te3\te4\te5\te6\te2_1\te3_1\te3_2\te4_2\n");
        fprintf(foutS, "centmin\tcentmax\tnpart\tS1\tS2\tS3\tS4\tS5\tS6\tS2_1\tS3_1\tS3_2\tS4_2\n");
        fprintf(foutc, "centmin\tcentmax\tnpart\t4_p4_p2\t8_p4_p2\t12_p4_p2\t6_p3_p2\t6_p2_p6\t6_p3_p6\t12_p3_p4\t10_p2_p5\t2p2_3p3_5p5\t2p2_4p4_6p6\t2p2_6p3_4p4\t8p2_3p3_5p5\t10p2_4p4_6p6\t10p2_6p3_4p4\n");
        fprintf(foutNpart, "centmin\tcentmax\tb\tnpart\n");
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            Int_t lmin = (int)centbins[cbin];
            Int_t lmax = (int)centbins[cbin+1];
            fprintf(fout, "%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n", lmin, lmax, avNpart[cbin], ecc[1][cbin], ecc[2][cbin], ecc[3][cbin], ecc[4][cbin], ecc[5][cbin], ecc[6][cbin], ecc_2_1[cbin], ecc_3_1[cbin], ecc_3_2[cbin], ecc_4_2[cbin]);
            fprintf(foutS, "%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n", lmin, lmax, avNpart[cbin], S[1][cbin], S[2][cbin], S[3][cbin], S[4][cbin], S[5][cbin], S[6][cbin], S_2_1[cbin], S_3_1[cbin], S_3_2[cbin], S_4_2[cbin]);
            fprintf(foutc, "%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n", lmin, lmax, avNpart[cbin], c_4_p4_p2[cbin], c_8_p4_p2[cbin], c_12_p4_p2[cbin], c_6_p3_p2[cbin], c_6_p2_p6[cbin], c_6_p3_p6[cbin],     c_12_p3_p4[cbin], c_10_p2_p5[cbin], c_2p2_3p3_5p5[cbin], c_2p2_4p4_6p6[cbin], c_2p2_6p3_4p4[cbin], c_8p2_3p3_5p5[cbin], c_10p2_4p4_6p6[cbin], c_10p2_6p3_4p4[cbin]);
            fprintf(foutNpart, "%d\t%d\t%f\t%f\n", lmin, lmax, avb[cbin], avNpart[cbin]);
        }
        fclose(fout);
        fclose(foutS);
        fclose(foutc);
        fclose(foutNpart);
    }


    // make plots WS plots
    TString tagWS = Form("A%d_A%d",A1,A2);
    TCanvas * cWS = new TCanvas("cWS","cWS",700,550);
    cWS->cd();
    hWS1->SetLineColor(kBlue);
    hWS2->SetLineColor(kRed);
    hWS1->Scale(1/(double)nevents);
    hWS2->Scale(1/(double)nevents);
    hWS1->SetXTitle("R (fm)");
    hWS2->SetXTitle("R (fm)");
    hWS1->Draw();
    hWS2->Draw("same");
    if (print_plots) cWS->Print(Form("plots/WS_%s.pdf",tagWS.Data()),"pdf");
    if (close_plots) cWS->Close();


    //-- write results
    if (!PlotCent) {
        TString atag = "";
        if (minorAxes) atag+="_minorAxes";
        tfout = new TFile(Form("results/glauber%s.root",atag.Data()),"recreate");
        TDirectory * tdmoment = (TDirectory *) tfout->mkdir("moments");
        TDirectory * tdcumulant = (TDirectory *) tfout->mkdir("cumulants");
        tdmoment->cd();
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * td = (TDirectory *) tdmoment->mkdir(Form("%d-%d",(int)centbins[cbin],(int)centbins[cbin+1]));
            td->cd();

            TDirectory * tdxy = (TDirectory *) td->mkdir("x_vs_y");
            tdxy->cd();
            xy0[cbin]->Write();
            xy1[cbin]->Write();
            xy2[cbin]->Write();
            xy3[cbin]->Write();
            xy4[cbin]->Write();
            xy5[cbin]->Write();
            xy6[cbin]->Write();

            TDirectory * tdpsi = (TDirectory *) td->mkdir("Psi");
            tdpsi->cd();
            for (int iorder = 1; iorder<=6; iorder++) {
                hpsi[iorder][cbin]->Write();
            }
            Psi2Psi1[cbin]->Write();
            Psi3Psi1[cbin]->Write();
            Psi3Psi2[cbin]->Write();
            Psi4Psi2[cbin]->Write();

            Psi2Psi1Diff[cbin]->Write();
            Psi3Psi1Diff[cbin]->Write();
            Psi5Psi1Diff[cbin]->Write();
            Psi3Psi2Diff[cbin]->Write();
            Psi4Psi2Diff[cbin]->Write();

            CosPsi2Psi1[cbin]->Write();
            CosPsi3Psi1[cbin]->Write();
            CosPsi3Psi2[cbin]->Write();
            CosPsi4Psi2[cbin]->Write();

            TDirectory * tdecc = (TDirectory *) td->mkdir("ecc");
            tdecc->cd();
            for (int iorder = 1; iorder<=6; iorder++) {eccAng[iorder][cbin]->Write();}
            for (int iorder = 1; iorder<=6; iorder++) {eccDist[iorder][cbin]->Write();}
        }

        for (int iorder = 1; iorder<=6; iorder++) {
            TDirectory * tdord = (TDirectory *) tdmoment->mkdir(Form("n_%d",iorder));
            tdord->cd();
            g[iorder]->SetName(Form("g%d",iorder));
            g[iorder]->SetTitle(Form("g%d",iorder));
            sg[iorder]->SetName(Form("sg%d",iorder));
            sg[iorder]->SetTitle(Form("sg%d",iorder));
            g[iorder]->Write();
            sg[iorder]->Write();
        }

        TDirectory * tdc = (TDirectory *) tdmoment->mkdir("correlations");
        tdc->cd();
        cg_4_p4_p2->SetName("<cos4(psi4-psi2)>");
        cg_4_p4_p2->SetTitle("<cos4(psi4-psi2)>");
        cg_8_p4_p2->SetName("<cos8(psi4-psi2)>");
        cg_8_p4_p2->SetTitle("<cos8(psi4-psi2)>");
        cg_12_p4_p2->SetName("<cos12(psi4-psi2)>");
        cg_12_p4_p2->SetTitle("<cos12(psi4-psi2)>");
        cg_6_p3_p2->SetName("<cos6(psi3-psi2)>");
        cg_6_p3_p2->SetTitle("<cos6(psi3-psi2)>");
        cg_6_p2_p6->SetName("<cos6(psi2-psi6)>");
        cg_6_p2_p6->SetTitle("<cos6(psi2-psi6)>");
        cg_6_p3_p6->SetName("<cos6(psi3-psi6)>");
        cg_6_p3_p6->SetTitle("<cos6(psi3-psi6)>");
        cg_12_p3_p4->SetName("<cos12(psi3-psi4)>");
        cg_12_p3_p4->SetTitle("<cos12(psi3-psi4)>");
        cg_10_p2_p5->SetName("<cos10(psi2-psi5)>");
        cg_10_p2_p5->SetTitle("<cos10(psi2-psi5)>");

        cg_2p2_3p3_5p5->SetName("<cos(2psi2+3psi3-5psi5)>");
        cg_2p2_3p3_5p5->SetTitle("<cos(2psi2+3psi3-5psi5)>");
        cg_2p2_4p4_6p6->SetName("<cos(2psi2+4psi4-6psi6)>");
        cg_2p2_4p4_6p6->SetTitle("<cos(2psi2+4psi4-6psi6)>");
        cg_2p2_6p3_4p4->SetName("<cos(2psi2-6psi3+4psi4)>");
        cg_2p2_6p3_4p4->SetTitle("<cos(2psi2-6psi3+4psi4)>");
        cg_8p2_3p3_5p5->SetName("<cos(-8psi2+3psi3+5psi5)>");
        cg_8p2_3p3_5p5->SetTitle("<cos(-8psi2+3psi3+5psi5)>");
        cg_10p2_4p4_6p6->SetName("<cos(-10psi2+4psi4+6psi6)>");
        cg_10p2_4p4_6p6->SetTitle("<cos(-10psi2+4psi4+6psi6)>");
        cg_10p2_6p3_4p4->SetName("<cos(-10psi2+6psi3+4psi4)>");
        cg_10p2_6p3_4p4->SetTitle("<cos(-10psi2+6psi3+4psi4)>");

        cg_6p2_p3->SetName("<cos6(psi2-psi3)>");
        cg_6p2_p3->SetTitle("<cos6(psi2-psi3)>");
        cg_4p2_p4->SetName("<cos4(psi2-psi4)>");
        cg_4p2_p4->SetTitle("<cos4(psi2-psi4)>");
        cg_6p2_p6->SetName("<cos6(psi2-psi6)>");
        cg_6p2_p6->SetTitle("<cos6(psi2-psi6)>");
        cg_6p3_p6->SetName("<cos6(psi3-psi6)>");
        cg_6p3_p6->SetTitle("<cos6(psi3-psi6)>");
        cg_10p2_p5->SetName("<cos10(psi2-psi5)>");
        cg_10p2_p5->SetTitle("<cos10(psi2-psi5)>");
        cg_15p3_p5->SetName("<cos15(psi3-psi5)>");
        cg_15p3_p5->SetTitle("<cos15(psi3-psi5)>");

        cg_2p1_p2->SetName("<cos2(psi1-psi2)>");
        cg_2p1_p2->SetTitle("<cos2(psi1-psi2)>");
        cg_3p1_p3->SetName("<cos3(psi1-psi3)>");
        cg_3p1_p3->SetTitle("<cos3(psi1-psi3)>");
        cg_4p1_p4->SetName("<cos4(psi1-psi4)>");
        cg_4p1_p4->SetTitle("<cos4(psi1-psi4)>");
        cg_5p1_p5->SetName("<cos5(psi1-psi5)>");
        cg_5p1_p5->SetTitle("<cos5(psi1-psi5)>");
        cg_6p1_p6->SetName("<cos6(psi1-psi6)>");
        cg_6p1_p6->SetTitle("<cos6(psi1-psi6)>");

        cg_1p1_2p2_3p3->SetName("<cos(1psi1+2psi2-3psi3)>");
        cg_1p1_2p2_3p3->SetTitle("<cos(1psi1+2psi2-3psi3)>");
        cg_4p1_2p2_6p3->SetName("<cos(4psi1+2psi2-6psi3)>");
        cg_4p1_2p2_6p3->SetTitle("<cos(4psi1+2psi2-6psi3)>");
        cg_2p1_4p2_6p3->SetName("<cos(2psi1+4psi2-6psi3)>");
        cg_2p1_4p2_6p3->SetTitle("<cos(2psi1+4psi2-6psi3)>");
        cg_5p1_2p2_3p3->SetName("<cos(5psi1-2psi2-3psi3)>");
        cg_5p1_2p2_3p3->SetTitle("<cos(5psi1-2psi2-3psi3)>");
        cg_2p1_2p2_4p4->SetName("<cos(2psi1+2psi2-4psi4)>");
        cg_2p1_2p2_4p4->SetTitle("<cos(2psi1+2psi2-4psi4)>");
        cg_2p1_6p2_4p4->SetName("<cos(2psi1-6psi2+4psi4)>");
        cg_2p1_6p2_4p4->SetTitle("<cos(2psi1-6psi2+4psi4)>");
        cg_2p1_6p2_8p4->SetName("<cos(2psi1+6psi2-8psi4)>");
        cg_2p1_6p2_8p4->SetTitle("<cos(2psi1+6psi2-8psi4)>");

        cg_1p1_3p3_4p4->SetName("<cos(1psi1+3psi3-4psi4)>");
        cg_1p1_3p3_4p4->SetTitle("<cos(1psi1+3psi3-4psi4)>");
        cg_2p1_6p3_4p4->SetName("<cos(2psi1-6psi3+4psi4)>");
        cg_2p1_6p3_4p4->SetTitle("<cos(2psi1-6psi3+4psi4)>");
        cg_1p1_4p2_5p5->SetName("<cos(1psi1+4psi2-5psi5)>");
        cg_1p1_4p2_5p5->SetTitle("<cos(1psi1+4psi2-5psi5)>");
        cg_3p1_2p2_5p5->SetName("<cos(3psi1+2psi2-5psi5)>");
        cg_3p1_2p2_5p5->SetTitle("<cos(3psi1+2psi2-5psi5)>");
        cg_1p1_6p3_5p5->SetName("<cos(1psi1-6psi3+5psi5)>");
        cg_1p1_6p3_5p5->SetTitle("<cos(1psi1-6psi3+5psi5)>");
        cg_2p1_3p3_5p5->SetName("<cos(2psi1+3psi3-5psi5)>");
        cg_2p1_3p3_5p5->SetTitle("<cos(2psi1+3psi3-5psi5)>");
        cg_4p1_9p3_5p5->SetName("<cos(4psi1-9psi3+5psi5)>");
        cg_4p1_9p3_5p5->SetTitle("<cos(4psi1-9psi3+5psi5)>");

        cg_4_p4_p2->Write();
        cg_8_p4_p2->Write();
        cg_12_p4_p2->Write();
        cg_6_p3_p2->Write();
        cg_6_p2_p6->Write();
        cg_6_p3_p6->Write();
        cg_12_p3_p4->Write();
        cg_10_p2_p5->Write();
        cg_2p2_3p3_5p5->Write();
        cg_2p2_4p4_6p6->Write();
        cg_2p2_6p3_4p4->Write();
        cg_8p2_3p3_5p5->Write();
        cg_10p2_4p4_6p6->Write();
        cg_10p2_6p3_4p4->Write();

        cg_6p2_p3->Write();
        cg_4p2_p4->Write();
        cg_6p2_p6->Write();
        cg_6p3_p6->Write();
        cg_10p2_p5->Write();
        cg_15p3_p5->Write();
        cg_2p1_p2->Write();
        cg_3p1_p3->Write();
        cg_4p1_p4->Write();
        cg_5p1_p5->Write();
        cg_6p1_p6->Write();
        cg_1p1_2p2_3p3->Write();
        cg_4p1_2p2_6p3->Write();
        cg_2p1_4p2_6p3->Write();
        cg_5p1_2p2_3p3->Write();
        cg_2p1_2p2_4p4->Write();
        cg_2p1_6p2_4p4->Write();
        cg_2p1_6p2_8p4->Write();
        cg_1p1_3p3_4p4->Write();
        cg_2p1_6p3_4p4->Write();
        cg_1p1_4p2_5p5->Write();
        cg_3p1_2p2_5p5->Write();
        cg_1p1_6p3_5p5->Write();
        cg_2p1_3p3_5p5->Write();
        cg_4p1_9p3_5p5->Write();

        tdmoment->cd();
        g_2_1->SetName("g_2_1");
        g_2_1->SetTitle("g_2_1");
        g_3_1->SetName("g_3_1");
        g_3_1->SetTitle("g_3_1");
        g_3_2->SetName("g_3_2");
        g_3_2->SetTitle("g_3_2");
        g_4_2->SetName("g_4_2");
        g_4_2->SetTitle("g_4_2");

        g_2_1->Write();
        g_3_1->Write();
        g_3_2->Write();
        g_4_2->Write();

        sg_2_1->SetName("sg_2_1");
        sg_2_1->SetTitle("sg_2_1");
        sg_3_1->SetName("sg_3_1");
        sg_3_1->SetTitle("sg_3_1");
        sg_3_2->SetName("sg_3_2");
        sg_3_2->SetTitle("sg_3_2");
        sg_4_2->SetName("sg_4_2");
        sg_4_2->SetTitle("sg_4_2");

        sg_2_1->Write();
        sg_3_1->Write();
        sg_3_2->Write();
        sg_4_2->Write();

        tdcumulant->cd();
        for (int cbin = 0; cbin<ncentbins; cbin++) {
            TDirectory * td = (TDirectory *) tdcumulant->mkdir(Form("%d-%d",(int)centbins[cbin],(int)centbins[cbin+1]));
            td->cd();

            TDirectory * tdxy = (TDirectory *) td->mkdir("x_vs_y");
            tdxy->cd();
            xyC0[cbin]->Write();
            xyC1[cbin]->Write();
            xyC2[cbin]->Write();
            xyC3[cbin]->Write();
            xyC4[cbin]->Write();
            xyC5[cbin]->Write();
            xyC6[cbin]->Write();

            TDirectory * tdpsi = (TDirectory *) td->mkdir("Psi");
            tdpsi->cd();
            for (int iorder = 1; iorder<=6; iorder++) {
                hpsiC[iorder][cbin]->Write();
            }
            Psi2Psi1C[cbin]->Write();
            Psi3Psi1C[cbin]->Write();
            Psi3Psi2C[cbin]->Write();
            Psi4Psi2C[cbin]->Write();

            Psi2Psi1DiffC[cbin]->Write();
            Psi3Psi1DiffC[cbin]->Write();
            Psi5Psi1DiffC[cbin]->Write();
            Psi3Psi2DiffC[cbin]->Write();
            Psi4Psi2DiffC[cbin]->Write();

            CosPsi2Psi1C[cbin]->Write();
            CosPsi3Psi1C[cbin]->Write();
            CosPsi3Psi2C[cbin]->Write();
            CosPsi4Psi2C[cbin]->Write();

            TDirectory * tdecc = (TDirectory *) td->mkdir("ecc");
            tdecc->cd();
            for (int iorder = 1; iorder<=6; iorder++) {eccAngC[iorder][cbin]->Write();}
            for (int iorder = 1; iorder<=6; iorder++) {eccDistC[iorder][cbin]->Write();}
        }

        for (int iorder = 1; iorder<=6; iorder++) {
            TDirectory * tdord = (TDirectory *) tdcumulant->mkdir(Form("n_%d",iorder));
            tdord->cd();
            gC[iorder]->SetName(Form("g%d",iorder));
            gC[iorder]->SetTitle(Form("g%d",iorder));
            gC[iorder]->Write();
        }

        TDirectory * tdcC = (TDirectory *) tdcumulant->mkdir("correlations");
        tdcC->cd();
        cgC_4_p4_p2->SetName("<cos4(psi4-psi2)>");
        cgC_4_p4_p2->SetTitle("<cos4(psi4-psi2)>");
        cgC_8_p4_p2->SetName("<cos8(psi4-psi2)>");
        cgC_8_p4_p2->SetTitle("<cos8(psi4-psi2)>");
        cgC_12_p4_p2->SetName("<cos12(psi4-psi2)>");
        cgC_12_p4_p2->SetTitle("<cos12(psi4-psi2)>");
        cgC_6_p3_p2->SetName("<cos6(psi3-psi2)>");
        cgC_6_p3_p2->SetTitle("<cos6(psi3-psi2)>");
        cgC_6_p2_p6->SetName("<cos6(psi2-psi6)>");
        cgC_6_p2_p6->SetTitle("<cos6(psi2-psi6)>");
        cgC_6_p3_p6->SetName("<cos6(psi3-psi6)>");
        cgC_6_p3_p6->SetTitle("<cos6(psi3-psi6)>");
        cgC_12_p3_p4->SetName("<cos12(psi3-psi4)>");
        cgC_12_p3_p4->SetTitle("<cos12(psi3-psi4)>");
        cgC_10_p2_p5->SetName("<cos10(psi2-psi5)>");
        cgC_10_p2_p5->SetTitle("<cos10(psi2-psi5)>");

        cgC_2p2_3p3_5p5->SetName("<cos(2psi2+3psi3-5psi5)>");
        cgC_2p2_3p3_5p5->SetTitle("<cos(2psi2+3psi3-5psi5)>");
        cgC_2p2_4p4_6p6->SetName("<cos(2psi2+4psi4-6psi6)>");
        cgC_2p2_4p4_6p6->SetTitle("<cos(2psi2+4psi4-6psi6)>");
        cgC_2p2_6p3_4p4->SetName("<cos(2psi2-6psi3+4psi4)>");
        cgC_2p2_6p3_4p4->SetTitle("<cos(2psi2-6psi3+4psi4)>");
        cgC_8p2_3p3_5p5->SetName("<cos(-8psi2+3psi3+5psi5)>");
        cgC_8p2_3p3_5p5->SetTitle("<cos(-8psi2+3psi3+5psi5)>");
        cgC_10p2_4p4_6p6->SetName("<cos(-10psi2+4psi4+6psi6)>");
        cgC_10p2_4p4_6p6->SetTitle("<cos(-10psi2+4psi4+6psi6)>");
        cgC_10p2_6p3_4p4->SetName("<cos(-10psi2+6psi3+4psi4)>");
        cgC_10p2_6p3_4p4->SetTitle("<cos(-10psi2+6psi3+4psi4)>");

        cgC_6p2_p3->SetName("<cos6(psi2-psi3)>");
        cgC_6p2_p3->SetTitle("<cos6(psi2-psi3)>");
        cgC_4p2_p4->SetName("<cos4(psi2-psi4)>");
        cgC_4p2_p4->SetTitle("<cos4(psi2-psi4)>");
        cgC_6p2_p6->SetName("<cos6(psi2-psi6)>");
        cgC_6p2_p6->SetTitle("<cos6(psi2-psi6)>");
        cgC_6p3_p6->SetName("<cos6(psi3-psi6)>");
        cgC_6p3_p6->SetTitle("<cos6(psi3-psi6)>");
        cgC_10p2_p5->SetName("<cos10(psi2-psi5)>");
        cgC_10p2_p5->SetTitle("<cos10(psi2-psi5)>");
        cgC_15p3_p5->SetName("<cos15(psi3-psi5)>");
        cgC_15p3_p5->SetTitle("<cos15(psi3-psi5)>");

        cgC_2p1_p2->SetName("<cos2(psi1-psi2)>");
        cgC_2p1_p2->SetTitle("<cos2(psi1-psi2)>");
        cgC_3p1_p3->SetName("<cos3(psi1-psi3)>");
        cgC_3p1_p3->SetTitle("<cos3(psi1-psi3)>");
        cgC_4p1_p4->SetName("<cos4(psi1-psi4)>");
        cgC_4p1_p4->SetTitle("<cos4(psi1-psi4)>");
        cgC_5p1_p5->SetName("<cos5(psi1-psi5)>");
        cgC_5p1_p5->SetTitle("<cos5(psi1-psi5)>");
        cgC_6p1_p6->SetName("<cos6(psi1-psi6)>");
        cgC_6p1_p6->SetTitle("<cos6(psi1-psi6)>");

        cgC_1p1_2p2_3p3->SetName("<cos(1psi1+2psi2-3psi3)>");
        cgC_1p1_2p2_3p3->SetTitle("<cos(1psi1+2psi2-3psi3)>");
        cgC_4p1_2p2_6p3->SetName("<cos(4psi1+2psi2-6psi3)>");
        cgC_4p1_2p2_6p3->SetTitle("<cos(4psi1+2psi2-6psi3)>");
        cgC_2p1_4p2_6p3->SetName("<cos(2psi1+4psi2-6psi3)>");
        cgC_2p1_4p2_6p3->SetTitle("<cos(2psi1+4psi2-6psi3)>");
        cgC_5p1_2p2_3p3->SetName("<cos(5psi1-2psi2-3psi3)>");
        cgC_5p1_2p2_3p3->SetTitle("<cos(5psi1-2psi2-3psi3)>");
        cgC_2p1_2p2_4p4->SetName("<cos(2psi2+2psi2-4psi4)>");
        cgC_2p1_2p2_4p4->SetTitle("<cos(2psi2+2psi2-4psi4)>");
        cgC_2p1_6p2_4p4->SetName("<cos(2psi1-6psi2+4psi4)>");
        cgC_2p1_6p2_4p4->SetTitle("<cos(2psi1-6psi2+4psi4)>");
        cgC_2p1_6p2_8p4->SetName("<cos(2psi1+6psi2-8psi4)>");
        cgC_2p1_6p2_8p4->SetTitle("<cos(2psi1+6psi2-8psi4)>");

        cgC_1p1_3p3_4p4->SetName("<cos(1psi1+3psi3-4psi4)>");
        cgC_1p1_3p3_4p4->SetTitle("<cos(1psi1+3psi3-4psi4)>");
        cgC_2p1_6p3_4p4->SetName("<cos(2psi1-6psi3+4psi4)>");
        cgC_2p1_6p3_4p4->SetTitle("<cos(2psi1-6psi3+4psi4)>");
        cgC_1p1_4p2_5p5->SetName("<cos(1psi1+4psi2-5psi5)>");
        cgC_1p1_4p2_5p5->SetTitle("<cos(1psi1+4psi2-5psi5)>");
        cgC_3p1_2p2_5p5->SetName("<cos(3psi1+2psi2-5psi5)>");
        cgC_3p1_2p2_5p5->SetTitle("<cos(3psi1+2psi2-5psi5)>");
        cgC_1p1_6p3_5p5->SetName("<cos(1psi1-6psi3+4psi4)>");
        cgC_1p1_6p3_5p5->SetTitle("<cos(1psi1-6psi3+4psi4)>");
        cgC_2p1_3p3_5p5->SetName("<cos(2psi1+3psi3-5psi5)>");
        cgC_2p1_3p3_5p5->SetTitle("<cos(2psi1+3psi3-5psi5)>");
        cgC_4p1_9p3_5p5->SetName("<cos(4psi1-9psi3+5psi5)>");
        cgC_4p1_9p3_5p5->SetTitle("<cos(4psi1-9psi3+5psi5)>");

        cgC_4_p4_p2->Write();
        cgC_8_p4_p2->Write();
        cgC_12_p4_p2->Write();
        cgC_6_p3_p2->Write();
        cgC_6_p2_p6->Write();
        cgC_6_p3_p6->Write();
        cgC_12_p3_p4->Write();
        cgC_10_p2_p5->Write();
        cgC_2p2_3p3_5p5->Write();
        cgC_2p2_4p4_6p6->Write();
        cgC_2p2_6p3_4p4->Write();
        cgC_8p2_3p3_5p5->Write();
        cgC_10p2_4p4_6p6->Write();
        cgC_10p2_6p3_4p4->Write();

        cgC_6p2_p3->Write();
        cgC_4p2_p4->Write();
        cgC_6p2_p6->Write();
        cgC_6p3_p6->Write();
        cgC_10p2_p5->Write();
        cgC_15p3_p5->Write();
        cgC_2p1_p2->Write();
        cgC_3p1_p3->Write();
        cgC_4p1_p4->Write();
        cgC_5p1_p5->Write();
        cgC_6p1_p6->Write();
        cgC_1p1_2p2_3p3->Write();
        cgC_4p1_2p2_6p3->Write();
        cgC_2p1_4p2_6p3->Write();
        cgC_5p1_2p2_3p3->Write();
        cgC_2p1_2p2_4p4->Write();
        cgC_2p1_6p2_4p4->Write();
        cgC_2p1_6p2_8p4->Write();
        cgC_1p1_3p3_4p4->Write();
        cgC_2p1_6p3_4p4->Write();
        cgC_1p1_4p2_5p5->Write();
        cgC_3p1_2p2_5p5->Write();
        cgC_1p1_6p3_5p5->Write();
        cgC_2p1_3p3_5p5->Write();
        cgC_4p1_9p3_5p5->Write();

        tdcumulant->cd();
        gC_2_1->SetName("gC_2_1");
        gC_2_1->SetTitle("gC_2_1");
        gC_3_1->SetName("gC_3_1");
        gC_3_1->SetTitle("gC_3_1");
        gC_3_2->SetName("gC_3_2");
        gC_3_2->SetTitle("gC_3_2");
        gC_4_2->SetName("gC_4_2");
        gC_4_2->SetTitle("gC_4_2");

        gC_2_1->Write();
        gC_3_1->Write();
        gC_3_2->Write();
        gC_4_2->Write();

        tfout->Close();
    }
}
