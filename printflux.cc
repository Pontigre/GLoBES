#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include "math.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"

char * fInputName = "nuprism_spectra.root";

int main(int argc, char *argv[])
{
  TFile *finput;
  TH2D *numu_oa_enu;
  TH2D *numub_oa_enu;
  TH2D *nue_oa_enu;
  TH2D *nueb_oa_enu;
  TH2D *anumu_oa_enu;
  TH2D *anumub_oa_enu;
  TH2D *anue_oa_enu;
  TH2D *anueb_oa_enu;

  if(fInputName!=NULL){
    finput = new TFile(fInputName);
    numu_oa_enu = (TH2D*)finput->Get("oa_enu_numode_numu");
    numub_oa_enu = (TH2D*)finput->Get("oa_enu_numode_numub");
    nue_oa_enu = (TH2D*)finput->Get("oa_enu_numode_nue");
    nueb_oa_enu = (TH2D*)finput->Get("oa_enu_numode_nueb");
    anumu_oa_enu = (TH2D*)finput->Get("oa_enu_anumode_numu");
    anumub_oa_enu = (TH2D*)finput->Get("oa_enu_anumode_numub");
    anue_oa_enu = (TH2D*)finput->Get("oa_enu_anumode_nue");
    anueb_oa_enu = (TH2D*)finput->Get("oa_enu_anumode_nueb");
  }
  
  TH1D *numu_enu = (TH1D*)numu_oa_enu->ProjectionX()->Clone("");
  TH1D *numub_enu = (TH1D*)numub_oa_enu->ProjectionX()->Clone("");
  TH1D *nue_enu = (TH1D*)nue_oa_enu->ProjectionX()->Clone("");
  TH1D *nueb_enu = (TH1D*)nueb_oa_enu->ProjectionX()->Clone("");
  TH1D *anumu_enu = (TH1D*)anumu_oa_enu->ProjectionX()->Clone("");
  TH1D *anumub_enu = (TH1D*)anumub_oa_enu->ProjectionX()->Clone("");
  TH1D *anue_enu = (TH1D*)anue_oa_enu->ProjectionX()->Clone("");
  TH1D *anueb_enu = (TH1D*)anueb_oa_enu->ProjectionX()->Clone("");

  numu_enu->Rebin();
  std::cout<<numu_enu->GetXaxis()->GetNbins();
  FILE *fp1 = fopen("JHF_nuPRISM_plus.dat","w");
  
  if (fp1 != NULL){
    for(int i = 1; i<=numu_enu->GetXaxis()->GetNbins(); i++){
      double E = numu_enu->GetBinCenter(i);
      double numu = numu_enu->GetBinContent(i);
      double numub = numub_enu->GetBinContent(i);
      double nue = nue_enu->GetBinContent(i);
      double nueb = nueb_enu->GetBinContent(i);
      fprintf(fp1,"%f %f %f 0 %f %f 0\n",E, numu, numub, nue, nueb);
    }
  }

  FILE *fp2 = fopen("JHF_nuPRISM_minus.dat","w");

  if (fp2 != NULL){
    for(int i = 1; i<=numu_enu->GetXaxis()->GetNbins(); i++){
      double E = anumu_enu->GetBinCenter(i);
      double numu = anumu_enu->GetBinContent(i);
      double numub = anumub_enu->GetBinContent(i);
      double nue = anue_enu->GetBinContent(i);
      double nueb = anueb_enu->GetBinContent(i);
      fprintf(fp2,"%f %f %f 0 %f %f 0\n",E, numu, numub, nue, nueb);
    }
  }
}
