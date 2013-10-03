// Cumulant vn code
// Hui Wang
// April, 2012

// Includes
#ifndef __CINT__
#include "TApplication.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TStyle.h>
#include <TComplex.h>

using namespace std;
#endif

#define PARTICLE_WEIGHT

//  Canvas and pads
TCanvas *c;
//  Histograms
TH1D *ZDCHist;
TH1D *refMultHist;

TProfile *C2Hist;
TProfile *C3Hist;
TProfile *T16Hist;
TProfile *C7Hist;
TProfile *C8Hist;
TProfile *C9Hist;
TProfile *C10Hist;
TProfile *T18Hist;

enum CorrelationTypes {
    C2 = 0, C3, T16, C7, C8, C9, C10, T18,
    C14, C15, T28, C17a, C17b, C18a, C18b, C19a, C19b, T32,
    nRefCorr = 8, nDiffCorr = 10
};

// Globals
int    EventNum;
int    RunNum;
int    refMult;
int    numberOfPrimaryVertices;
int    zdcCoincidenceRate;
int    trigger;
double zdcUnAttenuatedEast;
double zdcUnAttenuatedWest;
float  VertexX;
float  VertexY;
float  VertexZ;
int    numFiles;
float  vpdVz;
int    TOFmult;
int    Multiplicity;


//  Centrality
int  myCentral;
int  zdcHold[10];
char zdcName[10][200];
int  iRun;
char name[200];

// Prototypes
void balance_main(void);
//void initCanvas(void);
void initHists(void);
//void showItAll();
void cal_v2(void);

void balance_main(void) {
    FILE *fNames;
    FILE *fd;
    char dataFile[200];
    int iFile;

    int totalNumberOfEvents;
    int myDisplay;
    int myLocal;

#ifdef PARTICLE_WEIGHT
    float weight;
    const int nMax = 5;
    const int kMax = 5;
    const int pMax = 5;
    
    TComplex Q[nMax][kMax];     // Q_{n, k} = \sum w_i^k exp(i n \phi_i)
    TComplex QStar[nMax][kMax];     // complex conjugate of Q
    double   S[pMax][kMax];     // S_{p, k} = (\sum w_i^k)^p
#else    
    double sin2phiSum;
    double cos2phiSum;
    double sin4phiSum;
    double cos4phiSum;
#endif // m    
    //  Numbers for each event
    int num_tracks;

    myDisplay=1000000;
    myLocal=0;
    totalNumberOfEvents=0;

    TFile saveFile("v2_2.root","RECREATE","Check of Regular Events");

    //	setMyCenCuts();
    //	initCanvas();
    initHists();

    fNames=fopen("../good_runs_raw.txt","r");
    fscanf(fNames,"%d",&numFiles);
    cout << "We will read " << numFiles << " files" << endl;

    for(iFile=0; iFile<numFiles; iFile++) {
        fscanf(fNames,"%s",&dataFile);
        cout << "Begin file " << iFile+1 << ", " << dataFile << endl;
        fd=fopen(dataFile,"rb");
        cout << "Files open" << endl;

        fread(&EventNum,sizeof(int),1,fd);
        while(!feof(fd)) { //read each event
            fread(&RunNum,sizeof(int),1,fh);
            fread(&zdcUnAttenuatedEast,sizeof(double),1,fh);
            fread(&refMult,sizeof(int),1,fh);
            fread(&trigger,sizeof(int),1,fh);
            fread(&NtracksGlobal,sizeof(int),1,fh);
            fread(&numberOfPrimaryVertices,sizeof(int),1,fh);
            fread(&zdcCoincidenceRate,sizeof(int),1,fh);
            fread(&zdcUnAttenuatedWest,sizeof(double),1,fh);
            fread(&VertexX,sizeof(float),1,fh);
            fread(&VertexY,sizeof(float),1,fh);
            fread(&VertexZ,sizeof(float),1,fh);
            fread(&TOFmult,sizeof(int),1,fh);
            fread(&TOFmatchedTr,sizeof(int),1,fh);
            fread(&vpdVz,sizeof(float),1,fh);
            fread(&goodMult,sizeof(int),1,fh);

            refMultCorr=(int)(0.5+getRefMultCorr(refMult,VertexZ,RunNum,zdcCoincidenceRate));
            myCentral=cenBin[refMultCorr];

            fread(&EventNumd,sizeof(int),1,fd);
            fread(&RunNumd,sizeof(int),1,fd);
            //  		      cout << "EventNum ==  " <<EventNum << " EventNumd == " << EventNumd << "RunNum ==  " <<RunNum << " RunNumd == " << RunNumd << endl;
            if(EventNum != EventNumd) {
                cout<<"****************** You are in big trouble ***********" << endl;
                cout << "Header and track do not match,exiting............." << endl;
                cout << " EventNum = " << EventNum << "EventNumd = " << EventNumd << endl;
                exit(0);
            }
            // cout << "RunNum = " << RunNum << ", EventNum = " << EventNum << ", refMult = " << refMult << ", num_tracks = " << num_tracks <<", cos2phiSum = " << cos2phiSum << endl;
            //				myCentral=cenBin[refMult];

            myLocal++;
            totalNumberOfEvents++;
            if(myLocal == myDisplay) {
                cout << "Event " << totalNumberOfEvents << endl;
                // cout << "RunNum = " << RunNum << ", EventNum = " << EventNum << ", refMult = " << refMult << ", num_pos = " << num_pos << ", num_neg = " << num_neg << endl;
                myLocal=0;
                //				showItAll();
            }

            Multiplicity=0;
            num_tracks=0;

#ifdef PARTICLE_WEIGHT
            for (int n = 0; n < nMax; ++n) {
                for (int k = 0; k < kMax; ++k) {
                    Q[n][k] = TComplex(0, 0);
                }
            }

            for (int p = 0; p < pMax; ++n) {
                for (int k = 0; k < kMax; ++k) {
                    S[p][k] = 0.0;
                }
            }
#else
            sin2phiSum=0.0;
            cos2phiSum=0.0;
            sin4phiSum=0.0;
            cos4phiSum=0.0;
#endif // PARTICLE_WEIGHT 
            for(iEvent=0; iEvent<goodMult; iEvent++) { // loop of tracks in each event

                fread(&PID,sizeof(short),1,fd);
                fread(&pt,sizeof(float),1,fd);
                fread(&eta,sizeof(float),1,fd);
                fread(&phi,sizeof(float),1,fd);


                //          cout << "pt = " << pt<< " ,eta = " << eta<< " ,phi = " << phi << endl;

                if(fabs(eta)<1.0) Multiplicity++;

                if((pt > 0.2) && (pt < 2.0) && (fabs(eta) < 1.0)) {

#ifdef PARTICLE_WEIGHT
                    weight = getWeight(pt,eta,phi,mySagita,myCentral); // read each track's weight

                    for (int n = 0; n < nMax; ++n) {
                        for (int k = 0; k < kMax; ++k) {
                            Q[n][k] += TMath::Power(weight, k) * TComplex(1, n*phi); // TODO, optimize speed
                        }
                    }

                    for (int p = 0; p < pMax; ++n) { // S matrix is not ready
                        for (int k = 0; k < kMax; ++k) {
                            S[p][k] += TMath::Power(weight, k); // TODO, optimize speed
                        }
                    }
                    
#else
                    sin2phiSum+=sin(2.0*phi);
                    cos2phiSum+=cos(2.0*phi);
                    sin4phiSum+=sin(4.0*phi);
                    cos4phiSum+=cos(4.0*phi);
#endif // PARTICLE_WEIGHT

                    num_tracks++;
                }

            } // track loop end

#ifdef PARTICLE_WEIGHT
            for (int n = 0; n < nMax; ++n) {
                for (int k = 0; k < kMax; ++k) {
                    QStar[n][k] = TComplex::Conjugate( Q[n][k] );
                }
            }

            for (int p = 0; p < pMax; ++n) { // S matrix is ready
                for (int k = 0; k < kMax; ++k) {
                    S[p][k] = TMath::Power(S[p][k], p);
                }
            }
#endif // PARTICLE_WEIGHT
            
            if((fabs(VertexZ)<30)&&(num_tracks>3)) {

#ifdef PARTICLE_WEIGHT
                double M1    = S[1][1];
                double M11   = S[2][1] - S[1][2];
                double M111  = S[3][1] - 3*S[1][2]*S[1][1] + 2*S[1]*[3];
                double M1111 = S[4][1] - 6*S[1][2]*S[2][1] + 8*S[1][3]*S[1][1] + 3*S[2]*[2] - 6*S[1][4];

                double EventWeight1    = 1.0;
                double EventWeight11   = 1.0;
                double EventWeight111  = 1.0;
                double EventWeight1111 = 1.0;

                const int n = 2; // flow harmonics
                
                double B8 = ( (Q[n][1]).Rho2() - S[1][2] ) / M11; 
                double& coor22 = B8; // this is a wrong name, should use eq. number

                double B9 = ( TMath::Power( Q[n][1].Rho(), 4 ) + Q[2*n][2].Rho2() - 2 * (Q[2*n][2]*QStar[n][1]*QStar[n][1]).Re()
                              + 8 * (Q[n][3]*QStar[n][1]).Re() - 4 * S[1][2] * Q[n][1].Rho()
                              - 6 * S[1][4] - 2 * S[2][2] ) / M1111;
                double& coor24 = B9;

                // number of conbination have been divided
                // apply event weights, just in case non-unit weigts are needed
                T16Hist->Fill(refMult, coor22, EventWeight11); 
                T18Hist->Fill(refMult, coor24, EventWeight1111);
                
#else
                // ref cumulants w/o particle weight
                Double_t M       = num_tracks;
                Double_t nComb1  = M;
                Double_t nComb2  = M * (M-1);
                Double_t nComb3  = M * (M-1) * (M-2);
                Double_t nComb4  = M * (M-1) * (M-2) * (M-3);

                TComplex Q2(cos2phiSum, sin2phiSum);
                TComplex Q4(cos4phiSum, sin4phiSum);
                TComplex Q2Star   = TComplex::Conjugate(Q2);
                TComplex Q4Star   = TComplex::Conjugate(Q4);

                Double_t Q2Square = Q2.Rho2();
                Double_t Q4Square = Q4.Rho2();
                Double_t ReQQQ    = (Q4 * Q2Star * Q2Star).Re();

                Double_t coor22   = (Q2Square - M);
                Double_t coor24   = (Q2Square*Q2Square + Q4Square - 2*ReQQQ
                                     - 4*(M-2)*Q2Square + 2*M*(M-3));
                TComplex C7_8     = Q2 * Q2 - Q4;
                TComplex C9_10    = Q2*Q2Star*Q2Star - Q2*Q4Star - 2*(M-1)*Q2Star;

                C2Hist->Fill(refMult,Q2.Re()/nComb1);
                C3Hist->Fill(refMult,Q2.Im()/nComb1);
                T16Hist->Fill(refMult,coor22/nComb2);
                C7Hist->Fill(refMult,C7_8.Re()/nComb2);
                C8Hist->Fill(refMult,C7_8.Im()/nComb2);
                C9Hist->Fill(refMult,C9_10.Re()/nComb3);
                C10Hist->Fill(refMult,C9_10.Im()/nComb3);
                T18Hist->Fill(refMult,coor24/nComb4);
#endif  // PARTICLE_WEIGHT

                refMultHist->Fill(refMult);
                ZDCHist->Fill(zdcUnAttenuatedEast+zdcUnAttenuatedWest);

            }

            fread(&EventNum,sizeof(int),1,fd);
        }
        fclose(fd);
    }
    cout << "Total number of events read = " << (float)totalNumberOfEvents/1000000 <<"M"<< endl;

    //	showItAll();
    cal_v2();

    //	c->Print("v2_2.eps");
    //	c->Print("V2_2.png");
    saveFile.Write();
    saveFile.Close();


}

void cal_v2() {

    FILE *fout;
    int i,iRef;
    float c2_2[80],c2_2_err[80],c2_4[80],c2_4_err[80];
    float v1,v2,v3,v4,v5,v6;



    Double_t refMean[nRefCorr]        = {0};
    Double_t refMeanError2[nRefCorr]  = {0};


    for(iRef=0; iRef<80; iRef++) {
        //		v1=refMultHist->GetBinCenter(i+1);
        //		v2=refMultHist->GetBinContent(i+1);

        refMean[C2]=C2Hist->GetBinContent(iRef+1);
        refMean[C3]=C3Hist->GetBinContent(iRef+1);
        refMean[T16]=T16Hist->GetBinContent(iRef+1);
        refMean[C7]=C7Hist->GetBinContent(iRef+1);
        refMean[C8]=C8Hist->GetBinContent(iRef+1);
        refMean[C9]=C9Hist->GetBinContent(iRef+1);
        refMean[C10]=C10Hist->GetBinContent(iRef+1);
        refMean[T18]=T18Hist->GetBinContent(iRef+1);

        refMeanError2[C2]=TMath::Power(C2Hist->GetBinError(iRef+1), 2);
        refMeanError2[C3]=TMath::Power(C3Hist->GetBinError(iRef+1), 2);
        refMeanError2[T16]=TMath::Power(T16Hist->GetBinError(iRef+1), 2);
        refMeanError2[C7]=TMath::Power(C7Hist->GetBinError(iRef+1), 2);
        refMeanError2[C8]=TMath::Power(C8Hist->GetBinError(iRef+1), 2);
        refMeanError2[C9]=TMath::Power(C9Hist->GetBinError(iRef+1), 2);
        refMeanError2[C10]=TMath::Power(C10Hist->GetBinError(iRef+1), 2);
        refMeanError2[T18]=TMath::Power(T18Hist->GetBinError(iRef+1), 2);



        Double_t c22     = refMean[T16] - refMean[C2]*refMean[C2] - refMean[C3]*refMean[C3];
        Double_t c22Err2 = 0;

        Double_t dC22dx2[3];
        dC22dx2[C2]  = 4*refMean[C2]*refMean[C2];
        dC22dx2[C3]  = 4*refMean[C3]*refMean[C3];
        dC22dx2[T16] = 1.0;
        for ( int i = 0; i < 3; ++i ) {
            c22Err2 += dC22dx2[i] * refMeanError2[i];
        }

        Double_t c24 = refMean[T18] - 2*refMean[T16]*refMean[T16]
            - 4*refMean[C2]*refMean[C9]
            + 4*refMean[C3]*refMean[C10]
            - refMean[C7]*refMean[C7] - refMean[C8]*refMean[C8]
            + 4*refMean[C7]*(refMean[C2]*refMean[C2] - refMean[C3]*refMean[C3])
            + 8*refMean[C8]*refMean[C3]*refMean[C2]
            + 8*refMean[T16]*(refMean[C2]*refMean[C2] + refMean[C3]*refMean[C3])
            - 6*pow(refMean[C2]*refMean[C2] + refMean[C3]*refMean[C3], 2);

        // c24 = -TMath::Abs(c24);

        Double_t dC24dx2[8];
        dC24dx2[C2]  = 16 * TMath::Power(6*pow(refMean[C2], 3) - 2*refMean[C3]*refMean[C8] + refMean[C9]
                                         - 2*refMean[C2]*(-3*refMean[C3]*refMean[C3] + refMean[C7] + 2*refMean[T16]), 2);
        dC24dx2[C3]  = 16 * TMath::Power(refMean[C10] - 6*refMean[C3]*(refMean[C2]*refMean[C2]+refMean[C3]*refMean[C3])
                                         - 2*refMean[C3]*refMean[C7] + 2*refMean[C2]*refMean[C8] + 4*refMean[C3]*refMean[T16], 2);
        dC24dx2[T16] = TMath::Power(8*(refMean[C2]*refMean[C2] + refMean[C3]*refMean[C3]) - 4*refMean[T16], 2);
        dC24dx2[C7]  = TMath::Power(4*refMean[C2]*refMean[C2] - 2*(2*refMean[C3]*refMean[C3]+refMean[C7]), 2);
        dC24dx2[C8]  = TMath::Power(8*refMean[C2]*refMean[C3] -2*refMean[C8], 2);
        dC24dx2[C9]  = 16 * refMean[C2] * refMean[C2];
        dC24dx2[C10] = 16 * refMean[C3] * refMean[C3];
        dC24dx2[T18] = 1;

        Double_t c24Err2 = 0;
        for ( int i = 0; i < nRefCorr; ++i ) {
            c24Err2 += dC24dx2[i] * refMeanError2[i];
        }

        c2_2[iRef]=c22;
        c2_2_err[iRef]=sqrt(fabs(c22Err2));
        c2_4[iRef]=c24;
        c2_4_err[iRef]=sqrt(fabs(c24Err2));


    }

    fout=fopen("v2_2_parameters.txt","w");
    fprintf(fout,"bin#	bin_center #_of_events v2_2	err v2_4 err\n");

    for(i=0; i<80; i++) {
        v1=refMultHist->GetBinCenter(i+1);
        v2=refMultHist->GetBinContent(i+1);

        v3=sqrt(fabs(c2_2[i]));
        if(c2_2[i]<0) v3=-v3;
        if(c2_2[i]!=0) {
            v4=0.5*v3*(c2_2_err[i]/c2_2[i]);
        } else {
            v4=0;
        }


        v5=sqrt(sqrt(fabs(c2_4[i])));
        if(c2_4[i]>0) v5=-v5;
        if(c2_4[i]!=0) {
            v6=-0.25*v5*(c2_4_err[i]/c2_4[i]);
        } else {
            v6=0;
        }


        fprintf(fout,"%d %e %e %e %e %e %e\n",i,v1,v2,v3,v4,v5,v6);
    }

    fclose(fout);

    cout <<"v2 calculation finished" << endl;

}

/*
void initCanvas(){
	c = new TCanvas("c","STAR Data",10,30,1200,800);
	c->SetFillColor(0); // White canvas
//	c->GetFrame()->SetFillColor(0); // White frame
//	c->GetFrame()->SetBorderSize(12);
	c->Divide(3,2);

	cout << "Canvas initialized" << endl;
}
*/

void initHists() {

    ZDCHist = new TH1D("ZDC","ZDC",1000,0,4000);
    ZDCHist->GetXaxis()->SetTitle("ZDC");
    ZDCHist->GetYaxis()->SetTitle("Counts");
    ZDCHist->SetFillColor(2);

    refMultHist = new TH1D("refMult","refMult",80,0,800);
    refMultHist->GetXaxis()->SetTitle("refMult");
    refMultHist->GetYaxis()->SetTitle("Counts");
    refMultHist->SetFillColor(3);

    C2Hist= new TProfile("C2","C2",80,0,800);
    C2Hist->GetXaxis()->SetTitle("Mult");
    C2Hist->GetYaxis()->SetTitle("C2");

    C3Hist= new TProfile("C3","C3",80,0,800);
    C3Hist->GetXaxis()->SetTitle("Mult");
    C3Hist->GetYaxis()->SetTitle("C3");

    T16Hist= new TProfile("T16","T16",80,0,800);
    T16Hist->GetXaxis()->SetTitle("Mult");
    T16Hist->GetYaxis()->SetTitle("T16");

    C7Hist= new TProfile("C7","C7",80,0,800);
    C7Hist->GetXaxis()->SetTitle("Mult");
    C7Hist->GetYaxis()->SetTitle("C7");

    C8Hist= new TProfile("C8","C8",80,0,800);
    C8Hist->GetXaxis()->SetTitle("Mult");
    C8Hist->GetYaxis()->SetTitle("C8");

    C9Hist= new TProfile("C9","C9",80,0,800);
    C9Hist->GetXaxis()->SetTitle("Mult");
    C9Hist->GetYaxis()->SetTitle("C9");

    C10Hist= new TProfile("C10","C10",80,0,800);
    C10Hist->GetXaxis()->SetTitle("Mult");
    C10Hist->GetYaxis()->SetTitle("C10");

    T18Hist= new TProfile("T18","T18",80,0,800);
    T18Hist->GetXaxis()->SetTitle("Mult");
    T18Hist->GetYaxis()->SetTitle("T18");


    //	gStyle->SetOptStat(1);
    //	gStyle->SetPalette(1,0);

    cout << "Histograms initialized" << endl;
}



#ifndef __CINT__
int main(int argc,char **argv) {
    TApplication theApp("App",&argc,argv);
    balance_main();
    theApp.Run();
    exit(0);
}
#endif
