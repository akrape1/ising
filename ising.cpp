#include <cmath>
#include "TCanvas.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"

static const int NX = 64;
static const int NY = 64;

static int spin[NX+2][NY+2];
static int ntherm = 1000;


// different random guy instead of that drand48 
inline double myrand() {
  return gRandom->Rndm();
}

void update_spin(int *s, double env){
  int oldspin = *s;
  int newspin = ( myrand() < 0.5 ? 1 : -1 );
  double DeltaBetaE = -(newspin - oldspin) * env;

  if (DeltaBetaE <= 0 || myrand() < exp(-DeltaBetaE))
    *s = newspin;
}


void sweep(double beta, double h){
  for (int nx = 1; nx <= NX; nx++) {
    for (int ny = 1; ny <= NY; ny++) {
      double environment =
        beta * (spin[nx][ny-1] + spin[nx][ny+1] +
                spin[nx-1][ny] + spin[nx+1][ny])
        + h;

      update_spin(&spin[nx][ny], environment);
    }
  }
}


void InitializeHot(){
  for (int nx = 0; nx <= NX+1; nx++) {
    for (int ny = 0; ny <= NY+1; ny++) {
      if (nx == 0 || nx == NX+1 || ny == 0 || ny == NY+1)
        spin[nx][ny] = 0;
      else
        spin[nx][ny] = ( myrand() < 0.5 ? 1 : -1 );
    }
  }
}

double Magnetization(){
  int nmag = 0;
  for (int nx = 1; nx <= NX; nx++)
    for (int ny = 1; ny <= NY; ny++)
      nmag += spin[nx][ny];

  return (double)nmag / (NX * NY);  
}

// Energy per spin 
double Energy(double beta, double h){
  const int Nsites = NX * NY;

  double H = 0.0;
  if (beta > 0.0)
    H = h / beta;

  double E = 0.0;

  // interaction term
  for (int nx = 1; nx <= NX; nx++) {
    for (int ny = 1; ny <= NY; ny++) {
      int s = spin[nx][ny];

      if (nx < NX) E += - s * spin[nx+1][ny];
      if (ny < NY) E += - s * spin[nx][ny+1];
    }
  }

  // magnetic field term
  int M = 0;
  for (int nx = 1; nx <= NX; nx++)
    for (int ny = 1; ny <= NY; ny++)
      M += spin[nx][ny];

  E += -H * M;

  return E / Nsites;   
}


void ising(int nsweep = 3000, double h = 0.0, double Tmax = 4.0, int ntemp = 30){
  gStyle->SetOptStat(0);
  gRandom->SetSeed(0);

  double *T = new double[ntemp];
  double *Mval = new double[ntemp];
  double *Eval = new double[ntemp];
  double *Cval = new double[ntemp];

  InitializeHot();

  for (int itemp = ntemp; itemp > 0; --itemp) {

    int idx = itemp - 1;
    T[idx] = Tmax * itemp / ntemp;
    double beta = 1.0 / T[idx];

    //thermalize
    for (int n = 0; n < ntherm; n++)
      sweep(beta, h);

    double sumM = 0.0;
    double sumE = 0.0;
    double sumE2 = 0.0;

    // measurements
    for (int n = 0; n < nsweep; n++) {
      sweep(beta, h);

      double m = Magnetization();
      double e = Energy(beta, h);

      sumM  += m;
      sumE  += e;
      sumE2 += e * e;
    }

    // averages
    Mval[idx] = sumM / nsweep;
    Eval[idx] = sumE / nsweep;

    double Eavg  = Eval[idx];
    double E2avg = sumE2 / nsweep;

    //C = 1/N^2 (<E^2> - <E>^2)/T^2 but E=E/N so the N's should cancel out
    Cval[idx] = (E2avg - Eavg * Eavg) / (T[idx] * T[idx]);
  }


  TCanvas *c = new TCanvas("c", "2D Ising: E(T), M(T), C(T)", 900, 900);
  c->Divide(1,3);

  c->cd(1);
  TGraph *gE = new TGraph(ntemp, T, Eval);
  gE->SetTitle("Energy;T;E(T)");
  gE->SetLineWidth(2);
  gE->Draw("AL");

  c->cd(2);
  TGraph *gM = new TGraph(ntemp, T, Mval);
  gM->SetTitle("Magnetization;T;M(T)");
  gM->SetLineWidth(2);
  gM->Draw("AL");

  c->cd(3);
  TGraph *gC = new TGraph(ntemp, T, Cval);
  gC->SetTitle("Specific Heat;T;C(T)");
  gC->SetLineWidth(2);
  gC->Draw("AL");

  c->SaveAs("ising.pdf");

  delete [] T;
  delete [] Mval;
  delete [] Eval;
  delete [] Cval;
}
