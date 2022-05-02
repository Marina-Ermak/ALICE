#include "TCanvas.h"
#include "TFrame.h"
#include "TBenchmark.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TError.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TMath.h"
#include <cstring>

void eta_LHC_E100_E300() {
   TCanvas *c1 = new TCanvas("c1_fit2","The Fit Canvas",1200,800);
   TCanvas *c2 = new TCanvas("c2_fit2", "The Fit Canvas", 1200, 800);
   gStyle->SetOptFit(0);
   c1->Divide(2,2);
   c2->Divide(2,2);
   c1->SetGridx();
   c1->SetGridy();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderMode(-1);
   c1->GetFrame()->SetBorderSize(5);

     c2->SetGridx();
     c2->SetGridy();
     c2->GetFrame()->SetFillColor(21);
     c2->GetFrame()->SetBorderMode(-1);
     c2->GetFrame()->SetBorderSize(5);

   TFile *_file0 = TFile::Open("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/LHC18.root");//watch the name!
   THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEtaPHI7;1");
   THashList *_l1 =(THashList*)_file0->FindObjectAny("PHOSEtaMB;1");

   TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");
   TH2D *_h_all_cut=(TH2D*)_l->FindObject("InvMass_All_Rejection_E300");
   TH2D *_h_both=(TH2D*)_l->FindObject("InvMass_Both_E300");
   TH2D *_h_both_cut=(TH2D*)_l->FindObject("InvMass_Both_Rejection_E300");

  TH2D *_h_all_cut_mb=(TH2D*)_l1->FindObject("InvMass_All_Rejection_E300");
  TH2D *_h_all_mb=(TH2D*)_l1->FindObject("InvMass_All_E300");
  TH2D *_h_both_mb=(TH2D*)_l1->FindObject("InvMass_Both_E300");
  TH2D *_h_both_cut_mb=(TH2D*)_l1->FindObject("InvMass_Both_Rejection_E300");

   double pTbins[20]={1.5,3.5,5.5,7.5,9.5};// !!! to be corrected
/*   for (int i=0; i<16; i++)
          {
              pTbins[i]=1.5+i*1.5;
              cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
          }
*/

   TH1D* etaspectrum_both_cut = new TH1D("etaspectrum_both_cut","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_all_cut = new TH1D("etaspectrum_all_cut","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_both = new TH1D("etaspectrum_both","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_all = new TH1D("etaspectrum_all","eta spectrum", 4, pTbins);

   TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
   Double_t events=(Double_t)_ev->GetBinContent(2);

    cout<< "!!! Number of events - " << events<< "!!!"<<endl;

   Int_t bin1=0;
   Int_t bin2=0;

   TH1D* eta_eff_both = new TH1D("eta_eff_both","eta efficiency both", 4, pTbins);
   TH1D* eta_eff_both_cut = new TH1D("eta_eff_both_cut","eta efficiency both cut", 4, pTbins);
   TH1D* eta_eff_all_cut = new TH1D("eta_eff_all_cut","eta efficiency all cut", 4, pTbins);
   TH1D* eta_eff_all_cut_mb = new TH1D("eta_eff_all_cut_mb","eta efficiency all cut", 4, pTbins);
   TH1D* efficiency = new TH1D("efficiency","both cut/both", 4, pTbins);
   TH1D* efficiency_mb = new TH1D("efficiency_mb","both cut/both MB", 4, pTbins);
  Double_t x[4];
  Double_t y[4];
  Double_t exl[4];
  Double_t eyl[4];
  Double_t exh[4];
  Double_t eyh[4];

  TH1D* signal_bg_all = new TH1D("signal/background all","signal/bg all", 4, pTbins);
  TH1D* signal_bg_all_cut = new TH1D("signal/background all cut","signal/bg all cut", 4, pTbins);
  TH1D* signal_bg_both = new TH1D("signal/background both","signal/bg both", 4, pTbins);
  TH1D* signal_bg_both_cut = new TH1D("signal/background both cut","signal/bg both cut", 4, pTbins);

  TH1D* signal_bg_all_mb = new TH1D("signal/background all mb","signal/bg all", 4, pTbins);
  TH1D* signal_bg_all_cut_mb = new TH1D("signal/background all cut mb","signal/bg all cut", 4, pTbins);
  TH1D* signal_bg_both_mb = new TH1D("signal/background both mb","signal/bg both", 4, pTbins);
  TH1D* signal_bg_both_cut_mb = new TH1D("signal/background both cut mb","signal/bg both cut", 4, pTbins);

   int count=1;
   TString hist_name;
  for(int i=0; i<4; i++)
    {

        bin1 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i]);
        bin2 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i+1]);

            hist_name.Form("Projection %.2f < p_{T} < %.2f GeV/c",pTbins[i],pTbins[i+1]);
            cout<<"Hist_name - "<< hist_name<<endl;
            TH1D *_hf_all=_h_all->ProjectionX(hist_name,bin1,bin2);
            TH1D *_hf_all_cut=_h_all_cut->ProjectionX(hist_name+"cut",bin1,bin2);
            TH1D *_hf_both=_h_both->ProjectionX(hist_name+"both",bin1,bin2);
            TH1D *_hf_both_cut=_h_both_cut->ProjectionX(hist_name+"bothcut",bin1,bin2);
            TH1D *_hf_all_cut_mb=_h_all_cut_mb->ProjectionX(hist_name+"cut mb",bin1,bin2);
            TH1D *_hf_all_mb=_h_all_mb->ProjectionX(hist_name+"mb",bin1,bin2);
            TH1D *_hf_both_mb=_h_both_mb->ProjectionX(hist_name+"both mb",bin1,bin2);
            TH1D *_hf_both_cut_mb=_h_both_cut_mb->ProjectionX(hist_name+"bothcut mb",bin1,bin2);

            _hf_all->Rebin(10);
            _hf_all_cut->Rebin(10);
            _hf_both->Rebin(10);//6
            _hf_both_cut->Rebin(10);//6
            _hf_all_cut_mb->Rebin(10);
            _hf_all_mb->Rebin(10);
            _hf_both_mb->Rebin(10);//6
            _hf_both_cut_mb->Rebin(10);//6


            _hf_all_cut_mb->Sumw2();
            _hf_all_mb->Sumw2();
            _hf_all->Sumw2();
            _hf_all_cut->Sumw2();
            _hf_both->Sumw2();
            _hf_both_cut->Sumw2();
            _hf_both_mb->Sumw2();//6
            _hf_both_cut_mb->Sumw2();

           c1->cd(count);

    TF1  *f1 = new TF1("f1","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f2 = new TF1("f2","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f3 = new TF1("f3","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f4 = new TF1("f4","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f5g = new TF1("f5g","[0]*TMath::Exp(-(x-[1])*(x-[1])/2./[2]/[2])", 0.4, 0.7);
    TF1  *f_mb = new TF1("f_mb","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);//!
    TF1  *f_mb_cut = new TF1("f_mb_cut","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);//!
    TF1  *f3_mb = new TF1("f3_mb","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f4_mb = new TF1("f4_mb","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);

              f1->SetParameters(3.5*_hf_both->GetBinContent(66),0.55,0.02,-2e5/(i+1),3e4/(i+1),2.5*_hf_both->GetBinContent(51)/(i+1)) ;
              f2->SetParameters(2.5*_hf_both->GetBinContent(66),0.55,0.017,-2e5/(i+1),3e4/(i+1),2.5*_hf_both->GetBinContent(51)/(i+1)) ;
//              f2->SetParameters(1.*_hf_both->GetBinContent(66),0.55,0.017,1.,3.,2.*_hf_both->GetBinContent(51)) ;
              f3->SetParameters(1.*_hf_both->GetBinContent(66),0.55,0.017,1.,3.,2.*_hf_both->GetBinContent(51)) ;
              f4->SetParameters(1.*_hf_both->GetBinContent(66),0.55,0.017,1.,3.,2.*_hf_both->GetBinContent(51)) ;
              f3_mb->SetParameters(1.*_hf_both_mb->GetBinContent(66),0.55,0.02,1.,3.,2.*_hf_both_mb->GetBinContent(51)) ;//
              f4_mb->SetParameters(1.*_hf_both_mb->GetBinContent(66),0.55,0.02,1.,3.,2.*_hf_both_mb->GetBinContent(51)) ;//
              f_mb->SetParameters(1.*_hf_all_mb->GetBinContent(66),0.55,0.02,-1e4/(i+1)/(i+1),2e3/(i+1)/(i+1),1.7e4/(i+1)/(i+1)) ;
              f_mb_cut->SetParameters(1.*_hf_all_cut_mb->GetBinContent(66),0.55,0.02,-1e4/(i+1)/(i+1),1e3/(i+1)/(i+1),1.7e4/(i+1)/(i+1)) ;


              f1->SetParLimits(1, 0.5, 0.6);
              f1->SetParLimits(2, 0.01, 0.02);
              f_mb->SetParLimits(1, 0.5, 0.6);
              f_mb->SetParLimits(2, 0.01, 0.02);
              f2->SetParLimits(1, 0.5, 0.6);
              f2->SetParLimits(2, 0.01, 0.02);
              f_mb_cut->SetParLimits(1, 0.5, 0.6);
              f_mb_cut->SetParLimits(2, 0.01, 0.02);
              f3->SetParLimits(1, 0.5, 0.6);
              f3->SetParLimits(2, 0.01, 0.02);
              f3_mb->SetParLimits(1, 0.5, 0.6);
              f3_mb->SetParLimits(2, 0.01, 0.02);
              f4->SetParLimits(1, 0.5, 0.6);
              f4->SetParLimits(2, 0.01, 0.02);
              f4_mb->SetParLimits(1, 0.5, 0.6);
              f4_mb->SetParLimits(2, 0.01, 0.02);
              f1->SetParLimits(0, 0., 1e6);
              f2->SetParLimits(0, 0., 1e6);
              f_mb->SetParLimits(0, 0., 1e6);
              f_mb_cut->SetParLimits(0, 0., 1e6);

TFitResultPtr r_all = _hf_all->Fit("f1", "RMS+");
TFitResultPtr r_all_cut = _hf_all_cut->Fit("f2", "RMS+");
TFitResultPtr r_both = _hf_both->Fit("f3", "RMS+");
TFitResultPtr r_both_cut = _hf_both_cut->Fit("f4", "RMS+");
TFitResultPtr r_all_cut_mb = _hf_all_cut_mb->Fit("f_mb_cut", "RMS+");
TFitResultPtr r_all_mb = _hf_all_mb->Fit("f_mb", "RMS+");
TFitResultPtr r_both_mb = _hf_both_mb->Fit("f3_mb", "RMS+");//
TFitResultPtr r_both_cut_mb = _hf_both_cut_mb->Fit("f4_mb", "RMS+");//



             _hf_all->SetName(Form("all%d",i));
              _hf_all->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_all->SetTitle(hist_name);
             _hf_all->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV/c^{2}");
             _hf_all->SetLineColor(kBlue);

             _hf_all_mb->SetName(Form("all mb %d",i));
              _hf_all_mb->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_all_mb->SetTitle(hist_name);
             _hf_all_mb->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV/c^{2}");
             _hf_all_mb->SetLineColor(kBlue);

             _hf_all_cut->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_all_cut->SetName(Form("all_cut%d",i));
             _hf_all_cut->SetTitle(hist_name);
             _hf_all_cut->SetLineColor(kRed);

             _hf_all_cut_mb->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_all_cut_mb->SetName(Form("all_cut mb %d",i));
             _hf_all_cut_mb->SetTitle(hist_name);
             _hf_all_cut_mb->SetLineColor(kRed);

             _hf_both->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_both->SetName(Form("both%i",i));
             _hf_both->SetTitle(hist_name);
             _hf_both->SetLineColor(kBlue);

             _hf_both_cut->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_both_cut->SetName(Form("both_cut%i",i));
             _hf_both_cut->SetTitle(hist_name);
             _hf_both_cut->SetLineColor(kRed);

             _hf_both_mb->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_both_mb->SetName(Form("both%i mb",i));
             _hf_both_mb->SetTitle(hist_name);
             _hf_both_mb->SetLineColor(kBlue);

             _hf_both_cut_mb->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_both_cut_mb->SetName(Form("both_cut%i mb",i));
             _hf_both_cut_mb->SetTitle(hist_name);
             _hf_both_cut_mb->SetLineColor(kRed);

              _hf_all->SetMarkerColor(kBlue);
              _hf_all_cut->SetMarkerColor(kRed);
              _hf_all_mb->SetMarkerColor(kBlue);
              _hf_all_cut_mb->SetMarkerColor(kRed);
              _hf_both->SetMarkerColor(kBlue);
              _hf_both_cut->SetMarkerColor(kRed);
              _hf_both_mb->SetMarkerColor(kBlue);
              _hf_both_cut_mb->SetMarkerColor(kRed);


             Double_t p1[3],ep1[3], p2[3], ep2[3], p3[3], ep3[3], p4[3], ep4[3], p5[3], ep5[3], p6[3], ep6[3],  p1_mb[3],ep1_mb[3],  p2_mb[3], ep2_mb[3], p3_mb[3], ep3_mb[3], p4_mb[3], ep4_mb[3];
             TF1 *Number_eta_both = new TF1("number_eta_both_gaus","gaus", 0.4, 0.7);
             TF1 *Number_eta_both_cut = new TF1("number_eta_both_cut","gaus", 0.4, 0.7);
             TF1 *Number_eta_all_cut = new TF1("number_eta_all_cut","gaus", 0.4, 0.7);
             TF1 *Number_eta_all = new TF1("number_eta_all","gaus", 0.4, 0.7);
              TF1 *Number_eta_all_mb = new TF1("number_eta_all mb","gaus", 0.4, 0.7);
              TF1 *Number_eta_all_cut_mb = new TF1("number_eta_all_cut mb","gaus", 0.4, 0.7);
              TF1 *Number_eta_both_mb = new TF1("number_eta_both_gaus mb ","gaus", 0.4, 0.7);
              TF1 *Number_eta_both_cut_mb = new TF1("number_eta_both_cut mb","gaus", 0.4, 0.7);

              TMatrixDSym covrfit2 = r_all_cut->GetCovarianceMatrix();
              TMatrixTSym<double> dataBG2;
              covrfit2.GetSub(0,4,dataBG2);

              for(int j=0; j<3; j++)
              {
                  p2[j]=f2->GetParameter(j);
                  ep2[j]=f2->GetParError(j);
                  Number_eta_all_cut->SetParameter(j,p2[j]);

              }
              double num_all_cut = Number_eta_all_cut->Integral(0.4, 0.7)/_hf_all_cut->GetXaxis()->GetBinWidth(0);//num->f2
              double err_all_cut = Number_eta_all_cut->IntegralError(0.4, 0.7,r_all_cut->GetParams(), r_all_cut->GetCovarianceMatrix().GetMatrixArray())/_hf_all_cut->GetXaxis()->GetBinWidth(0);
              cout<<"Integral all cut - "<< num_all_cut  <<endl;
              cout<<"Integral all error - "<<err_all_cut <<endl;

              double signal_all_cut = f2->Integral(0.4, 0.7)/_hf_all_cut->GetXaxis()->GetBinWidth(0);//num
              double signal_all_cut_error = f2->IntegralError(0.4, 0.7,r_all_cut->GetParams(), dataBG2.GetMatrixArray())/_hf_all_cut->GetXaxis()->GetBinWidth(0);
              cout<<"signal_all_cut_error: "<<signal_all_cut_error<<endl;

              TMatrixDSym covrfit1 = r_all->GetCovarianceMatrix();
              TMatrixTSym<double> dataBG1;
              covrfit1.GetSub(0,4,dataBG1);

            for(int j=0; j<3; j++)
            {
                p1[j]=f1->GetParameter(j);
                ep1[j]=f1->GetParError(j);
                Number_eta_all->SetParameter(j,p1[j]);

            }
            cout<<"par: "<<f1->GetParameter(0)<<endl;
            double num_all = Number_eta_all->Integral(0.4, 0.7)/_hf_all->GetXaxis()->GetBinWidth(0);
            double err_all = Number_eta_all->IntegralError(0.4, 0.7,r_all->GetParams(), r_all->GetCovarianceMatrix().GetMatrixArray())/_hf_all->GetXaxis()->GetBinWidth(0);
            cout<<"Integral all - "<< num_all  <<endl;
            cout<<"Integral all error - "<<err_all <<endl;

            double signal_all = f1->Integral(0.4, 0.7)/_hf_all->GetXaxis()->GetBinWidth(0);
            double signal_all_error = f1->IntegralError(0.4, 0.7,r_all->GetParams(), dataBG1.GetMatrixArray())/_hf_all->GetXaxis()->GetBinWidth(0);



            TMatrixDSym covrfit1_mb = r_all_mb->GetCovarianceMatrix();
            TMatrixTSym<double> dataBG1_mb;
            covrfit1_mb.GetSub(0,4,dataBG1_mb);

            for(int j=0; j<3; j++)
            {
                p1_mb[j]=f_mb->GetParameter(j);
                ep1_mb[j]=f_mb->GetParError(j);
                Number_eta_all_mb->SetParameter(j,p1_mb[j]);
            }

            cout<<"par: "<<f_mb->GetParameter(0)<<endl;
            double num_all_mb = Number_eta_all_mb->Integral(0.4, 0.7)/_hf_all_mb->GetXaxis()->GetBinWidth(0);//not working
            double err_all_mb = Number_eta_all_mb->IntegralError(0.4, 0.7,r_all_mb->GetParams(), r_all_mb->GetCovarianceMatrix().GetMatrixArray())/_hf_all_mb->GetXaxis()->GetBinWidth(0);
            cout<<"Integral all mb - "<< num_all_mb  <<endl;
            cout<<"Integral all mb error - "<<err_all_mb <<endl;

            double signal_all_mb = f_mb->Integral(0.4, 0.7)/_hf_all_mb->GetXaxis()->GetBinWidth(0);
            double signal_all_error_mb = f_mb->IntegralError(0.4, 0.7,r_all_mb->GetParams(), dataBG1_mb.GetMatrixArray())/_hf_all_mb->GetXaxis()->GetBinWidth(0);


            TMatrixDSym covrfit2_mb = r_all_cut_mb->GetCovarianceMatrix();
            TMatrixTSym<double> dataBG2_mb;
            covrfit2_mb.GetSub(0,4,dataBG2_mb);

            for(int j=0; j<3; j++)
             {
                 p2_mb[j]=f_mb_cut->GetParameter(j);
                 ep2_mb[j]=f_mb_cut->GetParError(j);
                 Number_eta_all_cut_mb->SetParameter(j,p2_mb[j]);

            }
            Double_t err_true_cut_h;
            _hf_all_cut_mb->IntegralAndError(0, 19, err_true_cut_h);

            double num_all_cut_mb = Number_eta_all_cut_mb->Integral(0.4, 0.7)/_hf_all_cut_mb->GetXaxis()->GetBinWidth(0);///&&&&&
            double err_all_cut_mb = Number_eta_all_cut_mb->IntegralError(0.4, 0.7, r_all_cut_mb->GetParams(), r_all_cut_mb->GetCovarianceMatrix().GetMatrixArray())/_hf_all_cut_mb->GetXaxis()->GetBinWidth(0);
            cout<<"Integral all cut mb - "<< num_all_cut_mb  <<endl;
            cout<<"Integral all mb error - "<< err_true_cut_h <<endl;

            double signal_all_cut_mb = f_mb_cut->Integral(0.4, 0.7)/_hf_all_cut_mb->GetXaxis()->GetBinWidth(0);
            double signal_all_cut_error_mb = f_mb_cut->IntegralError(0.4, 0.7,r_all_cut_mb->GetParams(), dataBG2_mb.GetMatrixArray())/_hf_all_cut_mb->GetXaxis()->GetBinWidth(0);

                  TMatrixDSym covrfit3 = r_both->GetCovarianceMatrix();
                   TMatrixTSym<double> dataBG3;
                   covrfit3.GetSub(0,4,dataBG3);

            for(int j=0; j<3; j++)
            {
                  p3[j]=f3->GetParameter(j);
                  ep3[j]=f3->GetParError(j);
                  Number_eta_both->SetParameter(j,p3[j]);

            }

              double num_both = Number_eta_both->Integral(0.4, 0.7)/_hf_both->GetXaxis()->GetBinWidth(0);
              double err_both = Number_eta_both->IntegralError(0.4, 0.7,r_both->GetParams(), r_both->GetCovarianceMatrix().GetMatrixArray())/_hf_both->GetXaxis()->GetBinWidth(0);
              cout<<"Integral both - "<< num_both  <<endl;
              cout<<"Integral both error - "<< err_both <<endl;

              double signal_both = f3->Integral(0.4, 0.7)/_hf_both->GetXaxis()->GetBinWidth(0);
              double signal_both_error = f3->IntegralError(0.4, 0.7,r_both->GetParams(), dataBG3.GetMatrixArray())/_hf_both->GetXaxis()->GetBinWidth(0);

              TMatrixDSym covrfit4 = r_both_cut->GetCovarianceMatrix();
                 TMatrixTSym<double> dataBG4;
                 covrfit4.GetSub(0,4,dataBG4);

            for(int j=0; j<3; j++)
            {
                   p4[j]=f4->GetParameter(j);
                   ep4[j]=f4->GetParError(j);
                   Number_eta_both_cut->SetParameter(j,p4[j]);

            }
               double num_both_cut = Number_eta_both_cut->Integral(0.4, 0.7)/_hf_both_cut->GetXaxis()->GetBinWidth(0);//!
               double err_both_cut = Number_eta_both_cut->IntegralError(0.4, 0.7,r_both_cut->GetParams(), r_both_cut->GetCovarianceMatrix().GetMatrixArray())/_hf_both_cut->GetXaxis()->GetBinWidth(0);//!
               cout<<"Integral both cut - "<< num_both_cut  <<endl;
               cout<<"Integral both cut error - "<< err_both_cut <<endl;

               double signal_both_cut = f4->Integral(0.4, 0.7)/_hf_both_cut->GetXaxis()->GetBinWidth(0);
               double signal_both_cut_error = f4->IntegralError(0.4, 0.7,r_both_cut->GetParams(), dataBG4.GetMatrixArray())/_hf_both_cut->GetXaxis()->GetBinWidth(0);

               TMatrixDSym covrfit3_mb = r_both_mb->GetCovarianceMatrix();
                TMatrixTSym<double> dataBG3_mb;
                covrfit3_mb.GetSub(0,4,dataBG3_mb);

               for(int j=0; j<3; j++)
               {
                     p3_mb[j]=f3_mb->GetParameter(j);
                     ep3_mb[j]=f3_mb->GetParError(j);
                     Number_eta_both_mb->SetParameter(j,p3_mb[j]);

               }
       //        cout << "par 0 both - " << p3[1] <<endl;
                 double num_both_mb = Number_eta_both_mb->Integral(0.4, 0.7)/_hf_both_mb->GetXaxis()->GetBinWidth(0);//smth wrong
                 double err_both_mb = Number_eta_both_mb->IntegralError(0.4, 0.7,r_both_mb->GetParams(), r_both_mb->GetCovarianceMatrix().GetMatrixArray())/_hf_both_mb->GetXaxis()->GetBinWidth(0);
                 cout<<"Integral both mb - "<< num_both_mb  <<endl;
                 cout<<"Integral both mb error - "<< err_both_mb <<endl;

                 double signal_both_mb = f3_mb->Integral(0.4, 0.7)/_hf_both_mb->GetXaxis()->GetBinWidth(0);
                 double signal_both_error_mb = f3_mb->IntegralError(0.4, 0.7,r_both_mb->GetParams(), dataBG3_mb.GetMatrixArray())/_hf_both_mb->GetXaxis()->GetBinWidth(0);


                 TMatrixDSym covrfit4_mb = r_both_cut_mb->GetCovarianceMatrix();
                    TMatrixTSym<double> dataBG4_mb;
                    covrfit4_mb.GetSub(0,4,dataBG4_mb);
               for(int j=0; j<3; j++)
               {
                 p4_mb[j]=f4_mb->GetParameter(j);
                      ep4_mb[j]=f4_mb->GetParError(j);
                      Number_eta_both_cut_mb->SetParameter(j,p4_mb[j]);
               }
                  double num_both_cut_mb = Number_eta_both_cut_mb->Integral(0.4, 0.7)/_hf_both_cut_mb->GetXaxis()->GetBinWidth(0);//!
                  double err_both_cut_mb = Number_eta_both_cut_mb->IntegralError(0.4, 0.7,r_both_cut_mb->GetParams(), r_both_cut_mb->GetCovarianceMatrix().GetMatrixArray())/_hf_both_cut_mb->GetXaxis()->GetBinWidth(0);//!
                  cout<<"Integral both cut mb - "<< num_both_cut_mb  <<endl;
                  cout<<"Integral both cut mb error - "<< err_both_cut_mb <<endl;


                  double signal_both_cut_mb = f4_mb->Integral(0.4, 0.7)/_hf_both_cut_mb->GetXaxis()->GetBinWidth(0);
                  double signal_both_cut_error_mb = f4_mb->IntegralError(0.4, 0.7,r_both_cut_mb->GetParams(), dataBG4_mb.GetMatrixArray())/_hf_both_cut_mb->GetXaxis()->GetBinWidth(0);

                double dpt = pTbins[i+1]-pTbins[i];
                cout<<"DPT"<<dpt<<endl;
                etaspectrum_both_cut->SetBinContent(i+1, num_both_cut / dpt);
                etaspectrum_both_cut->SetBinError(i, err_both_cut / dpt);
                cout<<"!!BOTH CUT!!"<< num_both_cut/ dpt<<endl;
                  etaspectrum_both->SetBinContent(i+1, num_both / dpt);
                 etaspectrum_both->SetBinError(i, err_both / dpt);
                 cout<<"!!BOTH !!"<< num_both/ dpt<<endl;
                 etaspectrum_all_cut->SetBinContent(i+1, num_all_cut / dpt);
                 etaspectrum_all_cut->SetBinError(i, err_all_cut / dpt);
                 cout<<"!!ALL CUT !!"<< num_all_cut/ dpt<<endl;
                 etaspectrum_all->SetBinContent(i+1, num_all / dpt);
                 etaspectrum_all->SetBinError(i, err_all / dpt);
                 cout<<"!!ALL = !!"<< num_all/ dpt<<endl;


                 if(i==0){
                     signal_bg_all->SetBinContent(i+1, -1000);
                     signal_bg_all->SetBinError(i+1, 0.1);
                     signal_bg_all_cut->SetBinContent(i+1, 1000);
                     signal_bg_all_cut->SetBinError(i+1, 0.1);

                 }
                 else{
                     signal_bg_all->SetBinContent(i+1, num_all/(signal_all+num_all)/ dpt);
                     double errS1_all = err_all/(num_all+signal_all);
                     double errS2_all = (signal_all_error+err_all)*num_all/(signal_all+num_all)/(signal_all+num_all);
                     signal_bg_all->SetBinError(i+1, TMath::Sqrt(errS1_all*errS1_all+errS2_all*errS2_all)/ dpt);

                     signal_bg_all_cut->SetBinContent(i+1, num_all_cut/(signal_all_cut+num_all_cut)/ dpt);
                     double errS1_all_cut = err_all_cut/(signal_all_cut+num_all_cut);
                     double errS2_all_cut = (signal_all_cut_error+err_all_cut)*num_all_cut/(signal_all_cut+num_all_cut)/(signal_all_cut+num_all_cut);
                     signal_bg_all_cut->SetBinError(i+1, TMath::Sqrt(errS1_all_cut*errS1_all_cut+errS2_all_cut*errS2_all_cut) / dpt);
                 }


                 signal_bg_both->SetBinContent(i+1, num_both/(signal_both+num_both)/ dpt);

                  double errS1_both = err_both/(signal_both+num_both);
                  double errS2_both = (signal_both_error+err_both)*num_both/(signal_both+num_both)/(signal_both+num_both);

                 signal_bg_both->SetBinError(i+1, TMath::Sqrt(errS1_both*errS1_both+errS2_both*errS2_both)/ dpt);


                 signal_bg_both_cut->SetBinContent(i+1, num_both_cut/(signal_both_cut+num_both_cut)/ dpt);
                 double errS1_both_cut = err_both_cut/(signal_both_cut+num_both_cut);
                 double errS2_both_cut = (signal_both_cut_error+err_both_cut)*num_both_cut/(signal_both_cut+num_both_cut)/(signal_both_cut+num_both_cut);
                 signal_bg_both_cut->SetBinError(i+1, TMath::Sqrt(errS1_both*errS1_both+errS2_both*errS2_both) / dpt);
///////////////////////////////////////////////////////////////////

                 signal_bg_all_mb->SetBinContent(i+1, num_all_mb/(signal_all_mb+num_all_mb)/ dpt);

                 double errS1_all_mb = err_all_mb/(num_all_mb+signal_all_mb);
                 double errS2_all_mb = (signal_all_error_mb+err_all_mb)*num_all_mb/(signal_all_mb+num_all_mb)/(signal_all_mb+num_all_mb);


                 signal_bg_all_mb->SetBinError(i+1, TMath::Sqrt(errS1_all_mb*errS1_all_mb+errS2_all_mb*errS2_all_mb)/ dpt);

                 signal_bg_all_cut_mb->SetBinContent(i+1, num_all_cut_mb/(signal_all_cut_mb+num_all_cut_mb)/ dpt);

                 double errS1_all_cut_mb = err_all_mb/(signal_all_cut_mb+num_all_cut_mb);
                 double errS2_all_cut_mb = (signal_all_cut_error_mb+err_all_cut_mb)*num_all_cut_mb/(signal_all_cut_mb+num_all_cut_mb)/(signal_all_cut_mb+num_all_cut_mb);

                 signal_bg_all_cut_mb->SetBinError(i+1, TMath::Sqrt(errS1_all_cut_mb*errS1_all_cut_mb+errS2_all_cut_mb*errS2_all_cut_mb) / dpt);

                 signal_bg_both_mb->SetBinContent(i+1, num_both_mb/(signal_both_mb+num_both_mb)/ dpt);

                  double errS1_both_mb = err_both_mb/(signal_both_mb+num_both_mb);
                  double errS2_both_mb = (signal_both_error_mb+err_both_mb)*num_both_mb/(signal_both_mb+num_both_mb)/(signal_both_mb+num_both_mb);

                 signal_bg_both_mb->SetBinError(i+1, TMath::Sqrt(errS1_both_mb*errS1_both_mb+errS2_both_mb*errS2_both_mb)/ dpt);


                 signal_bg_both_cut_mb->SetBinContent(i+1, num_both_cut_mb/(signal_both_cut_mb+num_both_cut_mb)/ dpt);
                 double errS1_both_cut_mb = err_both_cut_mb/(signal_both_cut_mb+num_both_cut_mb);
                 double errS2_both_cut_mb = (signal_both_cut_error_mb+err_both_cut_mb)*num_both_cut_mb/(signal_both_cut_mb+num_both_cut_mb)/(signal_both_cut_mb+num_both_cut_mb);
                 signal_bg_both_cut_mb->SetBinError(i+1, TMath::Sqrt(errS1_both_mb*errS1_both_mb+errS2_both_mb*errS2_both_mb) / dpt);

            TLatex* t=new TLatex();
                auto legend_1 = new TLegend(0.14,0.13,0.88,0.25);
             //  legend_1->SetHeader("Efficiency E300","C"); // optio
              _hf_all_mb->SetMinimum(0.);
             _hf_all_mb->Draw();
            _hf_all_cut_mb->Draw("same");
          //   _hf_both_mb->SetMinimum(0.);
            // _hf_both_mb->Draw();
          //   _hf_both_cut_mb->Draw("same");
           legend_1->AddEntry(_hf_all_mb,Form("Integral All  MB: %2.1f #pm %2.1f",num_all_mb,err_all_mb),"lep");
           legend_1->AddEntry(_hf_all_cut_mb,Form("Integral All Rejection  MB: %2.1f #pm %2.1f",num_all_cut_mb,err_all_cut_mb),"lep");
          //  legend_1->AddEntry(_hf_all,Form("Integral Both MB: %2.1f #pm %2.1f",num_both_mb,err_both_mb),"lep");
          //  legend_1->AddEntry(_hf_all_cut,Form("Integral Both Rejection MB: %2.1f #pm %2.1f",num_both_cut_mb,err_both_cut_mb),"lep");
             legend_1->Draw("same");
             auto legend_2 = new TLegend(0.14,0.13,0.88,0.25);

              c2->cd(count);
              count++;
              _hf_both->SetMinimum(0.);
             _hf_both->Draw();
             _hf_both_cut->Draw("same");
         _hf_all->SetMinimum(0.);
          _hf_all->Draw();
         _hf_all_cut->Draw("same");
          legend_2->AddEntry(_hf_all,Form("Integral All PHI7: %2.1f #pm %2.1f",num_all,err_all),"lep");
          legend_2->AddEntry(_hf_all_cut,Form("Integral All Rejection PHI7: %2.1f #pm %2.1f",num_all_cut,err_all_cut),"lep");
            //  legend_2->AddEntry(_hf_both,Form("Integral Both PHI7: %2.1f #pm %2.1f",num_both,err_both),"lep");
          //    legend_2->AddEntry(_hf_both_cut,Form("Integral Both Rejection PHI7: %2.1f #pm %2.1f", num_both_cut, err_both_cut),"lep");
              legend_2->Draw("same");
          //    auto legend_2 = new TLegend(0.098,0.1,0.289,0.216);

               if(num_all>0)
               {
                eta_eff_all_cut->SetBinContent(i+1, num_all_cut/num_all);
                double errA1 = err_all_cut/num_all;
                double errA2 = err_all*num_all_cut/num_all/num_all;
                cout << " eta_eff_all_cut " << num_all_cut/num_all<<endl;

                eta_eff_both->SetBinContent(i+1, num_both/num_all);
                double errB1 = err_both/num_all;
                double errB2 = err_all*num_both/num_all/num_all;
                cout << " eta_eff_both " << num_both/num_all<<endl;


                eta_eff_both_cut->SetBinContent(i+1, num_both_cut/num_all);
                double errBC1 = err_both_cut/num_all;
                double errBC2 = err_all*num_both_cut/num_all/num_all;

                cout << " eta_eff_both_cut " << num_both_cut/num_all<<endl;

                eta_eff_both->SetBinError(i+1, TMath::Sqrt(errB1*errB1+errB2*errB2));
                eta_eff_both_cut->SetBinError(i+1, TMath::Sqrt(errBC1*errBC1+errBC2*errBC2));
                eta_eff_all_cut->SetBinError(i+1, TMath::Sqrt(errA1*errA1+errA2*errA2));


                double err_bc=TMath::Sqrt(errBC1*errBC1+errBC2*errBC2);
                double err_b = TMath::Sqrt(errB1*errB1+errB2*errB2);

                efficiency->SetBinContent(i+1, (num_both_cut/num_both));

                double erref1 = err_both_cut/num_both;
                double erref2 = err_both*num_both_cut/num_both/num_both;
                efficiency->SetBinError(i+1, TMath::Sqrt(erref1*erref1+erref2*erref2));
                cout<<"cycle number i = "<<i<<endl;
                x[i]=pTbins[i+1]-1.5/2;
                y[i]=num_both_cut/num_both;
                exl[i]=exh[i]=(pTbins[i+1]-pTbins[i])/2;
                eyl[i]=TMath::Sqrt(erref1*erref1+erref2*erref2);
                eyh[i]=(1-num_both_cut/num_both);
                cout<<x[i]<<endl;
                cout<<y[i]<<endl;
                cout<<exl[i]<<endl;
                cout<<eyl[i]<<endl;
                cout<<exh[i]<<endl;
                cout<<eyh[i]<<endl;

             }
             if(num_both_mb>0){
             efficiency_mb->SetBinContent(i+1, (num_both_cut_mb/num_both_mb));
             double err_mb1 = err_both_cut_mb/num_both_mb;
             double err_mb2 = err_both_mb*num_both_cut_mb/num_both_mb/num_both_mb;
             efficiency_mb->SetBinError(i+1, TMath::Sqrt(err_mb1*err_mb1+err_mb2*err_mb2));

           }
             if(num_all_cut_mb>0)
             {
               eta_eff_all_cut_mb->SetBinContent(i+1, num_all_cut_mb/num_all_mb);
               cout << " eta_eff_all_cut " << num_all_cut_mb/num_all_mb<<endl;
               double errA1_mb = err_all_cut_mb/num_all_mb;
               double errA2_mb = err_all_mb*num_all_cut_mb/num_all_mb/num_all_mb;
               cout << " eta_eff_all_cut " << num_all_cut_mb/num_all_mb<<endl;
               eta_eff_all_cut_mb->SetBinError(i+1, TMath::Sqrt(errA1_mb*errA1_mb+errA2_mb*errA2_mb));
             }
    }

 TCanvas *c3 = new TCanvas("c3_fit2", "The Fit Canvas", 1200, 800);
 TLatex* t2=new TLatex();

 gStyle->SetOptStat(0);
 gStyle->SetEndErrorSize(2);
 eta_eff_all_cut_mb->SetMinimum(0);
 eta_eff_all_cut_mb->SetMaximum(1);
 eta_eff_all_cut->SetMinimum(0);
 eta_eff_all_cut->SetMaximum(1);
 efficiency->SetMinimum(0);
 efficiency->SetMaximum(1.1);
 efficiency_mb->SetMinimum(0);
 efficiency_mb->SetMaximum(1.1);
 eta_eff_both->GetXaxis()->SetRangeUser(0., 10.);
 eta_eff_both->GetXaxis()->SetTitle("p_{T}, GeV/c");
 eta_eff_both->GetYaxis()->SetTitle("#epsilon");
 efficiency->GetXaxis()->SetTitle("p_{T}, GeV/c");
 efficiency->GetYaxis()->SetTitle("#Epsilon");
 eta_eff_both->SetMarkerStyle(21);
 eta_eff_both_cut->SetMarkerStyle(25);
 eta_eff_all_cut->SetMarkerStyle(41);
 efficiency->SetMarkerStyle(21);
 efficiency_mb->SetMarkerStyle(20);
 eta_eff_all_cut_mb->SetMarkerStyle(31);


  eta_eff_both->SetMarkerColor(kGreen);
  eta_eff_both_cut->SetMarkerColor(kGreen);
  eta_eff_all_cut->SetMarkerColor(kBlue);
  efficiency->SetMarkerColor(kBlack);
  eta_eff_all_cut_mb->SetMarkerColor(kRed);
  efficiency_mb->SetMarkerColor(kRed);

  eta_eff_both->SetMarkerSize(2);
  eta_eff_both_cut->SetMarkerSize(2);
  eta_eff_all_cut->SetMarkerSize(2);
 efficiency->SetMarkerSize(2);
   eta_eff_all_cut_mb->SetMarkerSize(2);
   efficiency_mb->SetMarkerSize(2);

 eta_eff_both->SetLineColor(kGreen);
 eta_eff_both_cut->SetLineColor(kGreen);
 eta_eff_all_cut->SetLineColor(kBlue);
 efficiency->SetLineColor(kBlack);
 eta_eff_all_cut_mb->SetLineColor(kRed);
 efficiency_mb->SetLineColor(kRed);
// eta_eff_both->Draw("p e1");
 //efficiency->Draw("p e1 same");

 auto legend = new TLegend(0.5,0.16,0.88,0.35);
  legend->SetHeader("Efficiency (E > 100 MeV)","C"); // option "C" allows to center the header
// legend->AddEntry(eta_eff_both,"Efficiency both", "lep");
//  legend->AddEntry(eta_eff_both_cut,"Efficiency both rejection","lep");
//  legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection PHI7","lep");
  //legend->AddEntry(eta_eff_all_cut_mb,"Efficiency all rejection MB","lep");
  legend->AddEntry(efficiency,"(both rejection)/both PHI7","lep");
  legend->AddEntry(efficiency_mb,"(both rejection)/both MB","lep");
 ///eta_eff_both_cut->Draw("p e1 same");
// eta_eff_all_cut->Draw("p e1");
// eta_eff_all_cut_mb->Draw("p e1 same");
 //efficiency->Draw("p e1");
// efficiency_mb->Draw("p e1 same");
// legend->Draw("same");
//////////////////////////////////////////////////////////////////////////////////////////////////////
   TCanvas *c_eff_mb_phi7 = new TCanvas("c4_fit2", "The Fit Canvas", 1200, 800);
   c_eff_mb_phi7->cd(1);
   cout<<"war!"<<endl;
   //auto gr_eff_PHI7 = new TGraphAsymmErrors(efficiency);
   auto gr_eff_PHI7 = new TGraphAsymmErrors(4,x,y,exl,exh,eyl,eyh);
   auto gr_eff_MB = new TGraphAsymmErrors(efficiency_mb);
   gr_eff_PHI7->SetMinimum(0.);
   gr_eff_PHI7->SetMaximum(1.);
   gr_eff_PHI7->SetMarkerColor(4);
   gr_eff_PHI7->SetMarkerStyle(21);
   gr_eff_PHI7->Draw();
//   gr_eff_MB->Draw("same");
//   c_eff_mb_phi7->SaveAs("Correct_errors_efficiency.pdf");
////////////////////////////////////////////////////////////////////////////////////////////////////
 TCanvas *c4 = new TCanvas("c4_spectrum", "The Fit Canvas", 1200, 800);
 gPad->SetLogy();
 etaspectrum_both_cut->Scale(1./events);
 etaspectrum_both_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
 etaspectrum_both_cut->GetYaxis()->SetTitle("dN/dp_{T}");
 etaspectrum_both_cut->SetMarkerStyle(21);

  etaspectrum_both->Scale(1./events);
  etaspectrum_both->GetXaxis()->SetTitle("p_{T}, GeV/c");
  etaspectrum_both->GetYaxis()->SetTitle("dN/dp_{T}");
  etaspectrum_both->SetMarkerStyle(31);
  etaspectrum_both->SetLineColor(kGreen);

  etaspectrum_all_cut->Scale(1./events);
  etaspectrum_all_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
  etaspectrum_all_cut->GetYaxis()->SetTitle("dN/dp_{T}");
  etaspectrum_all_cut->SetMarkerStyle(41);
  etaspectrum_all_cut->SetLineColor(kOrange);

  etaspectrum_all->Scale(1./events);
  etaspectrum_all->GetXaxis()->SetTitle("p_{T}, GeV/c");
  etaspectrum_all->GetYaxis()->SetTitle("dN/dp_{T}");
  etaspectrum_all->SetMarkerStyle(51);
  etaspectrum_all->SetLineColor(kRed);


//  auto legend_sp = new TLegend(0.1,0.1,0.38,0.3);
 // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
//    legend_sp->AddEntry(etaspectrum_both_cut,"spectrum both cut","lep");
//    legend_sp->AddEntry(etaspectrum_both,"spectrum both","lep");
  //  legend_sp->AddEntry(etaspectrum_all_cut,"spectrum all cut","lep");
  //  legend_sp->AddEntry(etaspectrum_all,"spectrum all","lep");
// TF1 *f6 = new TF1("f6","[0]*TMath::Exp(-[1]*x)", 6, 16);
// etaspectrum->Fit("f6", "RMS");
 //etaspectrum_both->Draw("p e");
// etaspectrum_both_cut->Draw("p e same");
// etaspectrum_all_cut->Draw("p e same");
 //etaspectrum_all->Draw("p e same");
// legend_sp->Draw("p e same");
// TF1 *f6 = new TF1("f6","[0]*TMath::Exp(-[1]*x)", 6, 16);
// etaspectrum_both_cut->Fit("f6", "RMS");

TCanvas *c_s_bg = new TCanvas("signal/background trig", "The Fit Canvas", 1200, 800);
c_s_bg->cd(1);
  signal_bg_all->SetStats(0);
  signal_bg_all_cut->SetStats(0);
  signal_bg_both->SetStats(0);
//  signal_bg_all->SetMaximum(0.062);
//  signal_bg_all->SetMinimum(0);
  signal_bg_all->GetXaxis()->SetRangeUser(0., 10.);
  signal_bg_all->GetXaxis()->SetTitle("p_{T}, GeV/c");
  signal_bg_all->GetYaxis()->SetTitle("s/(s+bg)");


  signal_bg_all->SetMarkerStyle(21);
  signal_bg_all_cut->SetMarkerStyle(21);
  signal_bg_both->SetMarkerStyle(20);
  signal_bg_both_cut->SetMarkerStyle(20);
  signal_bg_all->SetMarkerSize(2);
  signal_bg_all_cut->SetMarkerSize(2);
  signal_bg_both->SetMarkerSize(2);
  signal_bg_both_cut->SetMarkerSize(2);

  signal_bg_all->SetLineColor(kGreen+2);
  signal_bg_all_cut->SetLineColor(kBlack+2);
  signal_bg_both->SetLineColor(kBlue+2);
  signal_bg_both_cut->SetLineColor(kRed+2);
  signal_bg_all->SetMarkerColor(kGreen+2);
  signal_bg_all_cut->SetMarkerColor(kBlack+2);
  signal_bg_both->SetMarkerColor(kBlue+2);
  signal_bg_both_cut->SetMarkerColor(kRed+2);

  signal_bg_all->SetMinimum(-0.02);
  signal_bg_all->SetMaximum(0.15);

  signal_bg_all->SetTitle("");

  signal_bg_all->Draw("p e1");
  signal_bg_all_cut->Draw("p e1 same");
  signal_bg_both->Draw("p e1 same");
  signal_bg_both_cut->Draw("p e1 same");

  auto legend2 = new TLegend(0.12, 0.6,0.5,0.89);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend2->AddEntry(signal_bg_all,"Trigger: PHI7, PID: All","p");
  legend2->AddEntry(signal_bg_all_cut,"Trigger: PHI7, PID: All + rejection","p");
  legend2->AddEntry(signal_bg_both,"Trigger: PHI7, PID: Both","p");
  legend2->AddEntry(signal_bg_both_cut,"Trigger: PHI7, PID: Both + rejection","p");
  legend2->Draw("same");


TCanvas *c_s_bg_mb = new TCanvas("signal/background mb", "The Fit Canvas", 1200, 800);
c_s_bg_mb->cd(1);
signal_bg_all_mb->SetStats(0);
signal_bg_all_cut_mb->SetStats(0);
signal_bg_both_mb->SetStats(0);
signal_bg_all_mb->SetMaximum(0.062);
signal_bg_all_mb->SetMinimum(0);
signal_bg_all_mb->GetXaxis()->SetRangeUser(0., 10.);
signal_bg_all_mb->GetXaxis()->SetTitle("p_{T}, GeV/c");
signal_bg_all_mb->GetYaxis()->SetTitle("s/(s+bg)");
signal_bg_all_mb->SetMarkerStyle(21);
signal_bg_all_cut_mb->SetMarkerStyle(21);
signal_bg_both_mb->SetMarkerStyle(20);
signal_bg_both_cut_mb->SetMarkerStyle(20);
signal_bg_all_mb->SetMarkerSize(2);
signal_bg_all_cut_mb->SetMarkerSize(2);
signal_bg_both_mb->SetMarkerSize(2);
signal_bg_both_cut_mb->SetMarkerSize(2);

signal_bg_all_mb->SetLineColor(kGreen+2);
signal_bg_all_cut_mb->SetLineColor(kBlack+2);
signal_bg_both_mb->SetLineColor(kBlue+2);
signal_bg_both_cut_mb->SetLineColor(kRed+2);
signal_bg_all_mb->SetMarkerColor(kGreen+2);
signal_bg_all_cut_mb->SetMarkerColor(kBlack+2);
signal_bg_both_mb->SetMarkerColor(kBlue+2);
signal_bg_both_cut_mb->SetMarkerColor(kRed+2);

signal_bg_all_mb->SetMinimum(-0.02);
signal_bg_all_mb->SetMaximum(0.15);

signal_bg_all_mb->SetTitle("");

signal_bg_all_mb->Draw("p e1");
signal_bg_all_cut_mb->Draw("p e1 same");
signal_bg_both_mb->Draw("p e1 same");
signal_bg_both_cut_mb->Draw("p e1 same");

auto legend3 = new TLegend(0.12, 0.6,0.5,0.89);
// legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend3->AddEntry(signal_bg_all_mb,"Trigger: MB, PID: All","p");
  legend3->AddEntry(signal_bg_all_cut_mb,"Trigger: MB, PID: All + rejection","p");
  legend3->AddEntry(signal_bg_both_mb,"Trigger: MB, PID: Both","p");
  legend3->AddEntry(signal_bg_both_cut_mb,"Trigger: MB, PID: Both + rejection","p");
  legend3->Draw("same");


  TCanvas *c_s_bg_mb_ratio = new TCanvas("signal/background double ratio", "The Fit Canvas", 1200, 800);


  TH1D* ratio_all_mb = (TH1D*)signal_bg_all_cut_mb->Clone("ratio_all_mb");
  ratio_all_mb->Divide(signal_bg_all_mb);
  TH1D* ratio_both_mb = (TH1D*)signal_bg_both_cut_mb->Clone("ratio_both_mb");
  ratio_both_mb->Divide(signal_bg_both_mb);
  TH1D* ratio_all = (TH1D*)signal_bg_all_cut->Clone("ratio_all");
  ratio_all->Divide(signal_bg_all);
  TH1D* ratio_both = (TH1D*)signal_bg_both_cut->Clone("ratio_both");
  ratio_both->Divide(signal_bg_both);
  ratio_both->SetMarkerStyle(24);
  ratio_all->SetMarkerStyle(25);

  ratio_all->SetMarkerColor(kBlue+2);
  ratio_both->SetMarkerColor(kMagenta+2);
  ratio_all->SetLineColor(kBlue+2);
  ratio_both->SetLineColor(kMagenta+2);

  ratio_all_mb->SetMinimum(0.5);
  ratio_all_mb->SetMaximum(5.);

  ratio_all_mb->GetXaxis()->SetRangeUser(0., 10.);
  ratio_all_mb->GetXaxis()->SetTitle("p_{T}, GeV/c");
  ratio_all_mb->GetYaxis()->SetTitle("double ratio: rejection / no rejection");

  ratio_all_mb->SetTitle("");

  ratio_all_mb->Draw("p e1 ");
  ratio_both_mb->Draw("p e1 same");
  ratio_both->Draw("p e1 same");
  ratio_all->Draw("p e1 same");

  auto legend4 = new TLegend(0.12, 0.6,0.5,0.89);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend4->AddEntry(ratio_all_mb,"Trigger: MB, PID: All","p");
  legend4->AddEntry(ratio_both_mb,"Trigger: MB, PID: Both","p");
  legend4->AddEntry(ratio_all,"Trigger: PHI7, PID: All","p");
  legend4->AddEntry(ratio_both,"Trigger: PHI7, PID: Both","p");
  legend4->Draw("same");
c_s_bg_mb_ratio->Update();
c1->SaveAs("eta_hist_all_LHC_E300_MB.pdf");
c2->SaveAs("eta_hist_both_LHC_E300_MB.pdf");
c4->SaveAs("hist_eta_spectrum_model_E300_MB.pdf");
c3->SaveAs("eta_efficiency_model__everything_E300_MB.pdf");

//c_eff_mb_phi7->SaveAs("eta_efficiency_model_assym.pdf");
}
/*

           c2->cd(count);
           auto legend_2 = new TLegend(0.14,0.13,0.88,0.25);
           _hf_both_mb->SetMinimum(0.);
           _hf_both_mb->Draw();
           _hf_both_cut_mb->Draw("same");
           legend_2->AddEntry(_hf_both_mb,Form("Integral Both  MB: %2.1f #pm %2.1f",num_both_mb,err_both_mb),"lep");
           legend_2->AddEntry(_hf_both_cut_mb,Form("Integral Both Rejection  MB: %2.1f #pm %2.1f",num_both_cut_mb,err_both_cut_mb),"lep");
          legend_2->Draw("same");

           c3->cd(count);
           auto legend_3 = new TLegend(0.14,0.13,0.88,0.25);
           _hf_all_phi7->SetMinimum(0.);
           _hf_all_phi7->Draw();
           _hf_all_cut_phi7->Draw("same");
           legend_3->AddEntry(_hf_all_phi7,Form("Integral All phi7: %2.1f #pm %2.1f",num_all_phi7,err_all_phi7),"lep");
           legend_3->AddEntry(_hf_all_cut_phi7,Form("Integral All Rejection phi7: %2.1f #pm %2.1f",num_all_cut_phi7,err_all_cut_phi7),"lep");
           legend_3->Draw("same");

           c4->cd(count);
           auto legend_4 = new TLegend(0.14,0.13,0.88,0.25);
           _hf_both_phi7->SetMinimum(0.);
           _hf_both_phi7->Draw();
           _hf_both_cut_phi7->Draw("same");
           legend_4->AddEntry(_hf_both_phi7,Form("Integral Both phi7: %2.1f #pm %2.1f",num_both_phi7,err_both_phi7),"lep");
           legend_4->AddEntry(_hf_both_cut_phi7,Form("Integral Both Rejection  phi7: %2.1f #pm %2.1f",num_both_cut_mb,err_both_cut_phi7),"lep");
           legend_4->Draw("same");

           count++;
           */
