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
using namespace std;

Double_t CB1(Double_t * x, Double_t * par){
    Double_t kMean=0.136;
    Double_t m=par[1];
    Double_t s=par[2];
    Double_t n=par[3];
    Double_t a=par[4];
    Double_t dx=(x[0]-m)/s;
    if(dx>-a)
        return par[0]*exp(-dx*dx/2.)+par[5]+par[6]*(x[0]-kMean)+par[7]*(x[0]-kMean)*(x[0]-kMean);
        else{
            Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2);
            Double_t B=n/TMath::Abs(a)-TMath::Abs(a);
            return par[0]*A*TMath::Power((B-dx),-n)+par[5]+par[6]*(x[0]-kMean)+par[7]*(x[0]-kMean)*(x[0]-kMean);
            }
}

Double_t CB1S(Double_t * x, Double_t * par){
    Double_t kMean=0.136;
    Double_t m=par[1];
    Double_t s=par[2];
    Double_t n=par[3];
    Double_t a=par[4];
    Double_t dx=(x[0]-m)/s;
    if(dx>-a)
        return par[0]*exp(-dx*dx/2.);
        else{
            Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2);
            Double_t B=n/TMath::Abs(a)-TMath::Abs(a);
            return par[0]*A*TMath::Power((B-dx),-n);
            }
}

void pi0_model() {
   TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",1200,800);
   gStyle->SetOptFit(1);
   c1->Divide(3,6);//! watch
   c1->SetGridx();
   c1->SetGridy();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderMode(-1);
   c1->GetFrame()->SetBorderSize(5);

   TFile *_file0 = TFile::Open("histos_MC_21092020.root");//watch the name!
   THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEta;1");
   TH2D *_h_both_CB1=(TH2D*)_l->FindObject("InvMass_Both_E300");
//   TH2D *_h_both_gaus=(TH2D*)_l->FindObject("Both_Inv_mass");
   TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");

   TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
   Double_t events=(Double_t)_ev->GetBinContent(2);
   cout<< "!!! Number of events - " << events<< "!!!"<<endl;

   Int_t bin1=0;
   Int_t bin2=0;
   double pTbins[32];
   pTbins[0]=0.5;// {0.5, 0.7, 0.9, 1.2, 1.5, 2, 2.5, 3.2, 4., 5., 7., 10.}
    for (int i=1; i<=32; i++)
        {
          if(i<=10)
            pTbins[i]=pTbins[i-1]+0.5;
          else
            pTbins[i]=pTbins[i-1]+1;
            cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
        }

   TString hist_name;

   TF1  *f_gaus = new TF1("f_gaus","gaus+[3]*x+[4]+[5]*x*x", 0.08, 0.25);
   TF1  *f4_gaus = new TF1("f4_gaus","gaus+[3]*x+[4]+[5]*x*x", 0.08, 0.25);
   TF1  *f_CB1 = new TF1("f_CB1",CB1, 0.08, 0.2,8);
   TF1  *f_CB_all = new TF1("f_CB_all",CB1, 0.08, 0.2,8);

   TF1 *Number_pi0_both_CB1 = new TF1("number_pi0_both_CB1", CB1S, 0.1, 0.25, 5);//why gaus instead of CB1?
   TF1 *Number_pi0_both_gaus = new TF1("number_pi0_both_gaus","gaus", 0.1, 0.25);
   TF1 *Number_pi0_all = new TF1("number_pi0_2","gaus", 0.1, 0.25);
   TF1 *Number_pi0_all_CB = new TF1("number_pi0_2","gaus", 0.1, 0.25);


   TH1D* pi0_eff_both = new TH1D("pi0_eff_both","pi0 efficiency both", 32, pTbins);
   TH1D* pi0_eff_both_cut = new TH1D("pi0_eff_both_cut","pi0 efficiency both cut", 32, pTbins);
   TH1D* pi0_eff_all_cut = new TH1D("pi0_eff_all_cut","pi0 efficiency all cut", 32, pTbins);
   TH1D* pi0_eff_CB = new TH1D("pi0_eff_CB","pi0 efficiency CB", 32, pTbins);
   TH1D* pi0_eff_CB_m = new TH1D("pi0_eff_CB model","pi0 efficiency CB model", 32, pTbins);

   TH1D* pi0spectrum = new TH1D("pi0spectrum","pi0 spectrum", 32, pTbins);//11
   TH1D* pi0position = new TH1D("pi0position","pi0 position ", 32, pTbins);//11
   TH1D* pi0width = new TH1D("pi0width", "pi0 width", 32, pTbins);//watch the number of bins!

   TH1D* pi0spectrumCB = new TH1D("pi0spectrumCB","pi0 spectrum CB", 32, pTbins);//11
   TH1D* pi0positionCB = new TH1D("pi0positionCB","pi0 position CB", 32, pTbins);//11
   TH1D* pi0widthCB = new TH1D("pi0widthCB", "pi0 width CB", 32, pTbins);//watch the number of bins!

  Double_t p[8],ep[8], p2[3], ep2[3], p3[3], ep3[3], p4[3], ep4[3], p5[3], ep5[3], p6[3], ep6[3];

//  f_CB1->SetParLimits(2, 0.0, 0.025);
  for(int i=0; i<18; i++)
    {
            bin1 = (Int_t) _h_both_CB1->GetYaxis()->FindBin(pTbins[i]);
            bin2 = (Int_t) _h_both_CB1->GetYaxis()->FindBin(pTbins[i+1]);

            hist_name.Form("Projection %.2f < pT < %.2f",pTbins[i],pTbins[i+1]);
            cout<<"Hist_name - "<< hist_name<<endl;
            TH1D *_hf_both_CB1=_h_both_CB1->ProjectionX(hist_name+"CB1 both",bin1,bin2);
            TH1D *_hf_both_gaus=_h_both_CB1->ProjectionX(hist_name+"both",bin1,bin2);
            TH1D *_hf_all=_h_all->ProjectionX(hist_name+"all",bin1,bin2);
            TH1D *_hf_all_CB=_h_all->ProjectionX(hist_name+"all CB",bin1,bin2);


            _hf_both_CB1->GetXaxis()->SetRangeUser(0., 0.3);
            _hf_both_CB1->SetName(Form("projection both CB1 %i",i));
            _hf_both_CB1->SetTitle(hist_name);


             _hf_both_gaus->GetXaxis()->SetRangeUser(0., 0.3);
             _hf_both_gaus->SetName(Form("projection both gaus %i",i));
             _hf_both_gaus->SetTitle(hist_name);

            _hf_all->GetXaxis()->SetRangeUser(0., 0.3);
            _hf_all->SetName(Form("projection all %i",i));
            _hf_all->SetTitle(hist_name);
            _hf_all->SetLineColor(kBlue);

            _hf_all_CB->GetXaxis()->SetRangeUser(0., 0.3);
            _hf_all_CB->SetName(Form("projection all CB!!!!!!!! %i",i));
            _hf_all_CB->SetTitle(hist_name);
            _hf_all_CB->SetLineColor(kBlue);

               /* Comparison between CB1 and gaus*/

                f_CB1->SetParameters(1.*_hf_both_CB1->GetBinContent(66),0.135,0.0055,1.,3.,2.*_hf_both_CB1->GetBinContent(51),-0.001,0.0001) ;
                f_CB_all->SetParameters(1.*_hf_all_CB->GetBinContent(66),0.135,0.0055,1.,3.,2.*_hf_all_CB->GetBinContent(51),-0.001,0.0001) ;
                f_gaus->SetParameters(1.*_hf_both_CB1->GetBinContent(66),0.135,0.01,1.,3.,2.*_hf_both_CB1->GetBinContent(51)) ;
                f4_gaus->SetParameters(1.*_hf_both_CB1->GetBinContent(66),0.135,0.01,1.,3.,2.*_hf_both_CB1->GetBinContent(51)) ;

                f_gaus->SetLineStyle(2);
                f4_gaus->SetLineStyle(2);

                f_CB1->SetLineColor(kGreen);
                f_CB_all->SetLineColor(kGreen);
               c1->cd(i+1);

                if(i==17){
                 //  f_CB1->SetParameters(5.248, 0.138, -0.003591);
                  f_gaus->SetParameter(1, 0.138);
                }

              TFitResultPtr rCB = _hf_both_CB1->Fit("f_CB1", "RMS+");//hf

              TMatrixDSym covrfit1 = rCB->GetCovarianceMatrix();
 //      covrfit1.Print() ;
       TMatrixTSym<double> dataBG1;
       covrfit1.GetSub(0,4,dataBG1);

               for(int j=0; j<8; j++)
               {
                p[j]=f_CB1->GetParameter(j);
                ep[j]=f_CB1->GetParError(j);
                Number_pi0_both_CB1->SetParameter(j,p[j]);
                Number_pi0_both_CB1->SetParError(j,ep[j]);
               }

               double numpi0 = Number_pi0_both_CB1->Integral(0.1, 0.25)/_hf_both_CB1->GetXaxis()->GetBinWidth(0);
               double errpi0 = Number_pi0_both_CB1->IntegralError(0.1, 0.25,rCB->GetParams(), dataBG1.GetMatrixArray())/_hf_both_CB1->GetXaxis()->GetBinWidth(0);

               double dpt = pTbins[i+1]-pTbins[i];
               cout<<"Integral - "<< numpi0  <<endl;
               cout<<"Integral error - "<<errpi0 <<endl;

             cout<<"work2"<<endl;

               TFitResultPtr r_gaus = _hf_both_gaus->Fit("f_gaus", "RMS+");

                TMatrixDSym covrfit_g1 = r_gaus->GetCovarianceMatrix();
 //      covrfit1.Print() ;

       TMatrixTSym<double> dataBG_g1;
       covrfit_g1.GetSub(0,2,dataBG_g1);


                for(int j=0; j<3; j++)
                {
                 p2[j]=f_gaus->GetParameter(j);
                 ep2[j]=f_gaus->GetParError(j);
                 Number_pi0_both_gaus->SetParameter(j,p2[j]);
                 Number_pi0_both_gaus->SetParError(j,ep2[j]);
                }

                double num2pi0 = Number_pi0_both_gaus->Integral(0.1, 0.25)/_hf_both_gaus->GetXaxis()->GetBinWidth(0);
                double err2pi0 = Number_pi0_both_gaus->IntegralError(0.1, 0.25,r_gaus->GetParams(), dataBG_g1.GetMatrixArray())/_hf_both_gaus->GetXaxis()->GetBinWidth(0);

                 cout<<"Integral both - "<< num2pi0  <<endl;
                 cout<<"Integral both  error - "<<err2pi0 <<endl;

             /*   if(i==0)
               {
                   f_CB_all->SetParameters(2.535e+04, 0.1379, -0.008292, -2.562e+07, 725.6, 6.418e+04, -5.396e+05, 1.304e+06);
                   f4_gaus->SetParameters(3.112e+04, 0.1375, 0.008617, -2.012e+06, 3.631e+05, 3.032e+06);
               }
                if (i==13)
                {
                    f_CB_all->SetParameters(25.03, 0.1329, -0.0055, -5.84e+06, 140.9, 2.104, 5.078, 24);
                }
               */
                 TFitResultPtr r_all = _hf_all->Fit("f4_gaus", "RMS+");



                 TMatrixDSym covrfit_g2 = r_all->GetCovarianceMatrix();
 //      covrfit1.Print() ;

       TMatrixTSym<double> dataBG_g2;
       covrfit_g2.GetSub(0,2,dataBG_g2);

                   for(int j=0; j<3; j++)
                   {
                    p5[j]=f4_gaus->GetParameter(j);
                    ep5[j]=f4_gaus->GetParError(j);
                    Number_pi0_all->SetParameter(j,p5[j]);
                    Number_pi0_all->SetParError(j,ep5[j]);
                  }

                   double num5pi0 = Number_pi0_all->Integral(0.1, 0.25)/_hf_all->GetXaxis()->GetBinWidth(0);
                   double err5pi0 = Number_pi0_all->IntegralError(0.1, 0.25,r_all->GetParams(), dataBG_g2.GetMatrixArray())/_hf_all->GetXaxis()->GetBinWidth(0);

                    cout<<"Integral all  - "<< num5pi0  <<endl;
                    cout<<"Integral all   error - "<<err5pi0 <<endl;

                TFitResultPtr r_all_CB = _hf_all_CB->Fit("f_CB_all", "RMS+");


                TMatrixDSym covrfit2 = r_all_CB->GetCovarianceMatrix();
                 //      covrfit1.Print() ;

       TMatrixTSym<double> dataBG2;
       covrfit2.GetSub(0,4,dataBG2);

                    for(int j=0; j<3; j++)
                    {
                     p6[j]=f_CB_all->GetParameter(j);
                     ep6[j]=f_CB_all->GetParError(j);
                     Number_pi0_all_CB->SetParameter(j,p6[j]);
                     Number_pi0_all_CB->SetParError(j,ep6[j]);
                   }

                    double num_all_CB = Number_pi0_all_CB->Integral(0.1, 0.25)/_hf_all_CB->GetXaxis()->GetBinWidth(0);
                    double err_all_CB = Number_pi0_all_CB->IntegralError(0.1, 0.25,r_all_CB->GetParams(), dataBG2.GetMatrixArray())/_hf_all_CB->GetXaxis()->GetBinWidth(0);

                     cout<<"Integral all CB - "<< num_all_CB  <<endl;
                     cout<<"Integral all CB  error - "<<err_all_CB <<endl;


               pi0spectrum->SetBinContent(i+1, num2pi0 / (dpt));
               pi0spectrum->SetBinError(i+1,err2pi0 / dpt);

               pi0position->SetBinContent(i+1, p2[1]);
               pi0position->SetBinError(i+1, ep2[1]);

               pi0width->SetBinContent(i+1, abs(p2[2]));
               pi0width->SetBinError(i+1, abs(ep2[2]));

               pi0spectrumCB->SetBinContent(i+1, numpi0 / (dpt));
               pi0spectrumCB->SetBinError(i+1,errpi0 / dpt);

               pi0positionCB->SetBinContent(i+1, p[1]);
               pi0positionCB->SetBinError(i+1, ep[1]);

               pi0widthCB->SetBinContent(i+1, abs(p[2]));
               pi0widthCB->SetBinError(i+1, abs(ep[2]));


               _hf_both_CB1->Sumw2();
               _hf_both_gaus->Sumw2();
               _hf_all->Sumw2();
               _hf_all_CB->Sumw2();
             //  _hf_all_CB->Draw();
             //  _hf_all->Draw("same");
             //  _hf_all_CB->Draw("same");

               _hf_both_CB1->Draw();
               _hf_both_gaus->Draw("same");
              if(num5pi0>0){
              pi0_eff_both->SetBinContent(i+1, num2pi0/num5pi0);
              double errB1 = err2pi0/num5pi0;
              double errB2 = err5pi0*num2pi0/num5pi0/num5pi0;
              pi0_eff_both->SetBinError(i+1, TMath::Sqrt(errB1*errB1+errB2*errB2));
              }
              if(num_all_CB>0){
              pi0_eff_CB->SetBinContent(i+1, numpi0/num_all_CB);
              double errA1 = errpi0/num_all_CB;
              double errA2 = errpi0*numpi0/num_all_CB/num_all_CB;

             // pi0_eff_both->SetBinError(i+1, TMath::Sqrt(errB1*errB1+errB2*errB2));
              pi0_eff_CB->SetBinError(i+1, TMath::Sqrt(errA1*errA1+errA2*errA2));
        }
    }

TCanvas *c_eff = new TCanvas("c_eff","The Fit Canvas",1200,800);
pi0_eff_both->SetLineColor(kGreen);
//pi0_eff_both_cut->SetLineColor(kRed);
//pi0_eff_all_cut->SetLineColor(kBlue);
pi0_eff_both->GetXaxis()->SetTitle("p_{T}, GeV/c");
//pi0_eff_both->SetMaximum(1);
pi0_eff_both->SetMinimum(0);

pi0_eff_CB->Draw();
//pi0_eff_both->Draw("p e");


//pi0_eff_both_cut->Draw("p e same");
//pi0_eff_all_cut->Draw("p e same");

TCanvas *c2 = new TCanvas("c2_results","The Fit Canvas",1200,800);

c2->Divide(1,3);
c2->cd(1);
gPad->SetLogy();
//gPad->SetLogx();
pi0spectrum->Scale(1./events);
pi0spectrumCB->Scale(1./events);

pi0spectrum->GetXaxis()->SetTitle("p_{T}, GeV/c");
pi0spectrum->GetYaxis()->SetTitle("dN/dp_{T}");
pi0spectrum->SetMarkerStyle(21);
//TF1 *f3 = new TF1("f3","[0]*TMath::Exp(-[1]*x)", 0.5, 10);
//TF1 *f3 = new TF1("f3","[0]*TMath::Exp(-[1]*x)+TMath::Power([2]*(1+[3]*x*x), 13)", 0.5, 20);
//pi0spectrum->Fit("f3", "RMS");
pi0spectrum->Draw("p e ");
pi0spectrumCB->SetLineColor(kGreen);
pi0spectrumCB->SetMarkerColor(kGreen);
pi0spectrumCB->SetMarkerStyle(22);
pi0spectrumCB->Draw("p e same");

c2->cd(2);
pi0position->GetXaxis()->SetTitle("p_{T}, GeV/c");
pi0position->GetYaxis()->SetTitle("Peak position,9,pTbins, Gev");
pi0position->SetMarkerStyle(21);
pi0position->SetMinimum(0.13);
pi0position->SetMaximum(0.14);

TF1  *f4  = new TF1("f4","[0]", 0.5, 10);
pi0position->Fit("f4", "RMS");

pi0position->Draw("p e");
pi0positionCB->SetLineColor(kGreen);
pi0positionCB->SetMarkerStyle(22);
pi0positionCB->SetMarkerColor(kGreen);

pi0positionCB->Draw("p e same");

c2->cd(3);
TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*TMath::Power(x, [5])", 0.5, 10);
//TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*TMath::Power(x, [5])", 0.5, 10);
//TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*x+[5]*x*x", 0.5, 10);
//TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x)", 0.5, 20);
pi0width->GetXaxis()->SetTitle("p_{T}, GeV/c");
pi0width->GetYaxis()->SetTitle("Width, Gev");
pi0width->SetMinimum(0);
pi0width->SetMaximum(0.02);

pi0width->Fit("f5", "RMS");
pi0width->SetMarkerStyle(21);
pi0widthCB->SetLineColor(kGreen);
pi0widthCB->SetMarkerColor(kGreen);
pi0widthCB->SetMarkerStyle(22);

pi0width->Draw("p e");
pi0widthCB->Draw("p e same");

c1->SaveAs("fits_gaus+CB1_model_E300.pdf");
c2->SaveAs("pi0spectrum_CB1_model_E300.pdf");
c_eff->SaveAs("efficiency_pi0_model_E300.pdf");


}
