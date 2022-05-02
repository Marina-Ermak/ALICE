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
const Double_t kMean=0.136 ;
Double_t CB1 (Double_t * x, Double_t * par){
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

Double_t CB1S (Double_t * x, Double_t * par){
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

void pi0_LHC() {
   TCanvas *c1 = new TCanvas("c1_fit1","The Fit Canvas",1200,800);
   gStyle->SetOptFit(1);
   c1->Divide(8,4);//! watch
   c1->SetGridx();
   c1->SetGridy();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderMode(-1);
   c1->GetFrame()->SetBorderSize(5);

   TFile *_file0 = TFile::Open("LHC16.root");//watch the name! sum.root or AnalysisResults.root
   THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEta;1");
   TH2D *_h_both_CB1=(TH2D*)_l->FindObject("Both_Inv_mass");
   TH2D *_h_all = (TH2D*)_l->FindObject("A_Inv_mass");
   
    TFile *_file1 = TFile::Open("histos_sum.root");//watch the name!
    THashList *_l1 =(THashList*)_file1->FindObjectAny("PHOSEta;1");
    TH2D *_h_both_CB1_m=(TH2D*)_l1->FindObject("Both_Inv_mass");    
    TH2D *_h_all_m=(TH2D*)_l1->FindObject("A_Inv_mass");
    TH2D *_h_all_CB_m=(TH2D*)_l1->FindObject("A_Inv_mass");
    
   TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
   Double_t events=(Double_t)_ev->GetBinContent(2);
   cout<< "!!! Number of events - " << events<< "!!!"<<endl;
   
   Int_t bin1=0;
   Int_t bin2=0;
   double pTbins[32];// {0.5, 0.7, 0.9, 1.2, 1.5, 2, 2.5, 3.2, 4., 5., 7., 10.}
    for (int i=0; i<=32; i++)
        {
            pTbins[i]=0.5+i*0.5;
            cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
        }

   TString hist_name;

    TF1  *f_gaus = new TF1("f_gaus","gaus+[3]*x+[4]+[5]*x*x", 0.08, 0.25);
    TF1  *f_gaus_m = new TF1("f_gaus_m","gaus+[3]*x+[4]+[5]*x*x", 0.08, 0.25);
  
    TF1  *f4_gaus = new TF1("f4_gaus","gaus+[3]*x+[4]+[5]*x*x", 0.08, 0.25);
    TF1  *f4_gaus_m = new TF1("f4_gaus_m","gaus+[3]*x+[4]+[5]*x*x", 0.08, 0.25);
    TF1  *f_CB1 = new TF1("f_CB1",CB1, 0.08, 0.2,8);
    TF1  *f_CB_all = new TF1("f_CB_all",CB1, 0.08, 0.2,8);
    TF1  *f_CB_all_m = new TF1("f_CB_all_m",CB1, 0.08, 0.2,8);
    TF1  *f_CB_both_m = new TF1("f_CB_both_m",CB1, 0.08, 0.2,8);


   TF1 *Number_pi0_both_CB1 = new TF1("number_pi0_both_CB1", CB1S, 0.1, 0.25, 5);//why gaus instead of CB1?
   TF1 *Number_pi0_both_gaus = new TF1("number_pi0_both_gaus","gaus", 0.1, 0.25);
  
   TF1 *Number_pi0_all = new TF1("number_pi0_2","gaus", 0.1, 0.25);
   TF1 *Number_pi0_both_CB1_m = new TF1("number_pi0_both_CB1 model", CB1S, 0.1, 0.25, 5);
   TF1 *Number_pi0_both_gaus_m = new TF1("number_pi0_both_gaus model","gaus", 0.1, 0.25);
   TF1 *Number_pi0_all_m = new TF1("number_pi0_2_m","gaus", 0.1, 0.25);
   TF1 *Number_pi0_all_CB_m = new TF1("number_pi0_all_CB_m",CB1S, 0., 1.,5);
   TF1 *Number_pi0_all_CB = new TF1("number_pi0_all_CB",CB1S, 0.1, 0.25,5);
   TF1 *Number_pi0_both_CB_m = new TF1("number_pi0_both_CB_m",CB1S, 0.1, 0.25,5);
   
   TH1D* pi0_eff_both = new TH1D("pi0_eff_both","pi0 efficiency both", 32, pTbins);
   TH1D* pi0_eff_both_m = new TH1D("pi0_eff_both_m","pi0 efficiency both", 32, pTbins);
   
   TH1D* pi0_eff_both_CB = new TH1D("pi0_eff_both_CB","pi0 efficiency both CB", 32, pTbins);
   TH1D* pi0_eff_both_CB_m = new TH1D("pi0_eff_both_CB_m","pi0 efficiency both CB model", 32, pTbins);

   TH1D* pi0spectrum = new TH1D("pi0spectrum","pi0 spectrum", 32, pTbins);//11
   TH1D* pi0position = new TH1D("pi0position","pi0 position ", 32, pTbins);//11
   TH1D* pi0width = new TH1D("pi0width", "pi0 width", 32, pTbins);//watch the number of bins!

   TH1D* pi0spectrumCB = new TH1D("pi0spectrumCB","pi0 spectrum CB", 32, pTbins);//11
   TH1D* pi0positionCB = new TH1D("pi0positionCB","pi0 position CB", 32, pTbins);//11
   TH1D* pi0widthCB = new TH1D("pi0widthCB", "pi0 width CB", 32, pTbins);//watch the number of bins!

  Double_t p[8],ep[8], p2[3], ep2[3], p3[3], ep3[3], p4[3], ep4[3], p5[3], ep5[3], p6[8], ep6[8],p7[8], ep7[8], p8[3], ep8[3], p_cb_m[3], ep_cb_m[3], p_all_m[3], ep_all_m[3];

//  f_CB1->SetParLimits(2, 0.0, 0.025);
  for(int i=0; i<31; i++)
    {
            bin1 = (Int_t) _h_both_CB1->GetYaxis()->FindBin(pTbins[i]);
            bin2 = (Int_t) _h_both_CB1->GetYaxis()->FindBin(pTbins[i+1]);

            hist_name.Form("Projection %.2f < pT < %.2f",pTbins[i],pTbins[i+1]); 
            cout<<"Hist_name - "<< hist_name<<endl;
            TH1D *_hf_both_CB1=_h_both_CB1->ProjectionX(hist_name+"CB1 both",bin1,bin2);
            TH1D *_hf_both_CB1_m=_h_both_CB1_m->ProjectionX(hist_name+"CB1 both model",bin1,bin2);
            TH1D *_hf_both_gaus=_h_both_CB1->ProjectionX(hist_name+"both",bin1,bin2);
            TH1D *_hf_both_gaus_m=_h_both_CB1_m->ProjectionX(hist_name+"both model",bin1,bin2);
           
            TH1D *_hf_all=_h_all->ProjectionX(hist_name+"all",bin1,bin2);
            TH1D *_hf_all_m=_h_all_m->ProjectionX(hist_name+"all model",bin1,bin2);
            TH1D *_hf_all_CB=_h_all->ProjectionX(hist_name+"all CB",bin1,bin2);
            TH1D *_hf_all_CB_m=_h_all_CB_m->ProjectionX(hist_name+"all CB model",bin1,bin2);

            _hf_both_CB1->GetXaxis()->SetRangeUser(0., 0.3);
            _hf_both_CB1->SetName(Form("projection both CB1 %i",i));
            _hf_both_CB1->SetTitle(hist_name);

             _hf_both_CB1_m->GetXaxis()->SetRangeUser(0., 0.3);
             _hf_both_CB1_m->SetName(Form("projection both CB1 model %i",i));
             _hf_both_CB1_m->SetTitle(hist_name);

             _hf_both_gaus->GetXaxis()->SetRangeUser(0., 0.3);
             _hf_both_gaus->SetName(Form("projection both gaus %i",i));
             _hf_both_gaus->SetTitle(hist_name);

              _hf_both_gaus_m->GetXaxis()->SetRangeUser(0., 0.3);
              _hf_both_gaus_m->SetName(Form("projection both gaus %i",i));
              _hf_both_gaus_m->SetTitle(hist_name);
            
            _hf_all->GetXaxis()->SetRangeUser(0., 0.3);
            _hf_all->SetName(Form("projection all %i",i));
            _hf_all->SetTitle(hist_name);
            _hf_all->SetLineColor(kRed);
            
            _hf_all_CB->GetXaxis()->SetRangeUser(0., 0.3);
            _hf_all_CB->SetName(Form("projection all CB %i",i));
            _hf_all_CB->SetTitle(hist_name);
            _hf_all_CB->SetLineColor(kGreen);

            _hf_all_m->GetXaxis()->SetRangeUser(0., 0.3);
            _hf_all_m->SetName(Form("projection all model gaus %i",i));
            _hf_all_m->SetTitle(hist_name);
            _hf_all->SetLineColor(38);
	              
             _hf_all_CB_m->GetXaxis()->SetRangeUser(0., 0.3);
             _hf_all_CB_m->SetName(Form("projection all CB model %i",i));
             _hf_all_CB_m->SetTitle(hist_name+"model");
             _hf_all_CB_m->SetLineColor(kBlue);

               /* Comparison between CB1 and gaus*/

                f_CB1->SetParameters(1.*_hf_both_CB1->GetBinContent(66),0.136,0.005,1.,3.,2.*_hf_both_CB1->GetBinContent(51),-0.001,0.0001) ;
                f_CB_all->SetParameters(1.*_hf_all_CB->GetBinContent(66),0.136,0.005,1.,3.,2.*_hf_all_CB->GetBinContent(51),-0.001,0.0001) ;      
                
                f_CB_all_m->SetParameters(1.*_hf_all_CB_m->GetBinContent(66),0.136,0.005,1.,3.,2.*_hf_all_CB_m->GetBinContent(51),-0.001,0.0001) ;
                
                f_CB_both_m->SetParameters(1.*_hf_all_CB_m->GetBinContent(66),0.136,0.005,1.,3.,2.*_hf_all_CB_m->GetBinContent(51),-0.001,0.0001) ;

                f_gaus->SetParameters(1.*_hf_both_CB1->GetBinContent(66),0.136,0.01,1.,3.,2.*_hf_both_CB1->GetBinContent(51)) ;
                f_gaus_m->SetParameters(1.*_hf_both_CB1->GetBinContent(66),0.136,0.01,1.,3.,2.*_hf_both_CB1->GetBinContent(51)) ;
                f4_gaus->SetParameters(1.*_hf_both_CB1->GetBinContent(66),0.136,0.01,1.,3.,2.*_hf_both_CB1->GetBinContent(51)) ;
            //    f_CB_both_m->SetParameters(1.*_hf_both_CB1_m->GetBinContent(66),0.136,0.01,1.,3.,2.*_hf_both_CB1_m->GetBinContent(51)) ;

                f_gaus->SetLineStyle(2);
                f_gaus_m->SetLineStyle(4);
                f4_gaus->SetLineStyle(4);
                
                f_CB1->SetLineColor(kGreen);
                f_CB_all->SetLineColor(kGreen);
                f_CB_all_m->SetLineColor(kGreen);
                f_CB_both_m->SetLineColor(kGreen);


               c1->cd(i+1);

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
               cout<<"Integral both CB - "<< numpi0  <<endl;
               cout<<"Integral both CB error - "<<errpi0 <<endl;           
	       
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
                double err2pi0 = Number_pi0_both_gaus->IntegralError(0.1, 0.25, r_gaus->GetParams(), dataBG_g1.GetMatrixArray())/_hf_both_gaus->GetXaxis()->GetBinWidth(0);
                
                 cout<<"Integral both - "<< num2pi0  <<endl;
                 cout<<"Integral both  error - "<<err2pi0 <<endl;
                 
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
                    
                    if(i==1)
                    {
                        f_CB_all->SetParameter(1, 0.137);
                    } else if(i==23) {
                        f_CB_all->SetParameter(0, 5040);
                        f_CB_all->SetParameter(1, 0.137);
                        f_CB_all->SetParameter(2, 0.004886);
                        f_CB_all->SetParameter(3, -4.304e+05);
                        f_CB_all->SetParameter(4, 5.984e+04);
                        f_CB_all->SetParameter(5, 807.1);
                        f_CB_all->SetParameter(6, -366);
                        f_CB_all->SetParameter(7, -8.861e+04);
                    }
                        
                    TFitResultPtr r_all_CB = _hf_all_CB->Fit("f_CB_all", "RMS+");
 
                          TMatrixDSym covrfit2 = r_all_CB->GetCovarianceMatrix(); 
//      covrfit1.Print() ;
      
      TMatrixTSym<double> dataBG2;
      covrfit2.GetSub(0,4,dataBG2); 
                    
                    for(int j=0; j<8; j++)
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


		       if(i<=19)
                {

                    f_gaus_m->SetParameters(1.*_hf_both_CB1_m->GetBinContent(66),0.136,0.01,1.,3.,2.*_hf_both_CB1_m->GetBinContent(51)) ;
                    f4_gaus_m->SetParameters(1.*_hf_both_CB1_m->GetBinContent(66),0.136,0.01,1.,3.,2.*_hf_both_CB1_m->GetBinContent(51)) ;
                    
                   f_gaus_m->SetParameter(1, 0.137);
                   f4_gaus_m->SetParameter(1, 0.137); 
                    

               TFitResultPtr r_gaus_m = _hf_both_gaus_m->Fit("f_gaus_m", "RMS+");
               
                 TMatrixDSym covrfit_g3 = r_gaus_m->GetCovarianceMatrix(); 
//      covrfit1.Print() ;
      
      TMatrixTSym<double> dataBG_g3;
      covrfit_g3.GetSub(0,2,dataBG_g3);               
               
                for(int j=0; j<8; j++)
                {
                 p_cb_m[j]=f_gaus_m->GetParameter(j);
                 ep_cb_m[j]=f_gaus_m->GetParError(j);
                 Number_pi0_both_gaus_m->SetParameter(j,p_cb_m[j]);
                 Number_pi0_both_gaus_m->SetParError(j,ep_cb_m[j]);
                }

                double numpi0_m = Number_pi0_both_gaus_m->Integral(0.1, 0.25)/_hf_both_gaus_m->GetXaxis()->GetBinWidth(0);
                double errpi0_m = Number_pi0_both_gaus_m->IntegralError(0.1, 0.25,r_gaus_m->GetParams(), dataBG_g3.GetMatrixArray())/_hf_both_gaus_m->GetXaxis()->GetBinWidth(0);

                cout<<"Integral gaus both model- "<< numpi0_m  <<endl;
                cout<<"Integral gaus both model error - "<<errpi0_m <<endl;
                
                TFitResultPtr r_all_m = _hf_all_m->Fit("f4_gaus_m", "RMS+");
                 TMatrixDSym covrfit_g4 = r_all_m->GetCovarianceMatrix(); 
//      covrfit1.Print() ;
      
      TMatrixTSym<double> dataBG_g4;
      covrfit_g4.GetSub(0,2,dataBG_g4);  
      
                     for(int j=0; j<3; j++)
                     {
                      p_all_m[j]=f4_gaus_m->GetParameter(j);
                      ep_all_m[j]=f4_gaus_m->GetParError(j);
                      Number_pi0_all_m->SetParameter(j,p_all_m[j]);
                      Number_pi0_all_m->SetParError(j,ep_all_m[j]);
                    }
 
                     double num5pi0_m = Number_pi0_all_m->Integral(0.1, 0.25)/_hf_all_m->GetXaxis()->GetBinWidth(0);
                     double err5pi0_m = Number_pi0_all_m->IntegralError(0.1, 0.25,r_all_m->GetParams(), dataBG_g4.GetMatrixArray())/_hf_all_m->GetXaxis()->GetBinWidth(0);
                      cout<<"Integral all  model - "<< num5pi0_m <<endl;
                      cout<<"Integral all  model - "<< err5pi0_m <<endl;
                    
                pi0_eff_both_m->SetBinContent(i+1, numpi0_m/num5pi0_m);
               double errB1_m = errpi0_m/num5pi0_m;
               double errB2_m = err5pi0_m*numpi0_m/num5pi0_m/num5pi0_m;
	       
               pi0_eff_both_m->SetBinError(i+1, TMath::Sqrt(errB1_m*errB1_m+errB2_m*errB2_m));
                    
                    	       
                 f_CB_all_m->SetParameters(1*_hf_all_CB_m->GetBinContent(55),0.136,0.015,1.,3.,1.*_hf_all_CB_m->GetBinContent(71),-0.001,0.0001) ;
            /*     if(i==1)
                 {
                     f_CB_all_m->SetParameters(1.658e+05, 0.1372, -0.006325, -9.794e+06, 333.3, 2.223e+05, -1.112e+06, -3.051e+07);
                 }
            */
                  TFitResultPtr r_all_CB_m = _hf_all_CB_m->Fit("f_CB_all_m", "RMS+");

                                            TMatrixDSym covrfit3 = r_all_CB_m->GetCovarianceMatrix(); 
//      covrfit1.Print() ;
      
      TMatrixTSym<double> dataBG3;
      covrfit3.GetSub(0,4,dataBG3);
      
                  for(int j=0; j<8; j++)
                     {
                      p7[j]=f_CB_all_m->GetParameter(j);
                      ep7[j]=f_CB_all_m->GetParError(j);
                      Number_pi0_all_CB_m->SetParameter(j,p7[j]);
                      Number_pi0_all_CB_m->SetParError(j,ep7[j]);
                    }
                     cout<<p7[0]<<endl;
                 //   cout<<Number_pi0_all_CB_m->Integral(0.1, 0.25)<<endl;
                     double num_all_CB_m = Number_pi0_all_CB_m->Integral(0.1, 0.25)/_hf_all_CB_m->GetXaxis()->GetBinWidth(0);
                     cout<<"work f3 -  "<<i<<endl;
                     double err_all_CB_m = Number_pi0_all_CB_m->IntegralError(0.1, 0.25,r_all_CB_m->GetParams(), dataBG3.GetMatrixArray())/_hf_all_CB_m->GetXaxis()->GetBinWidth(0);
                     cout<<"work f3 -  "<<i<<endl;
                      cout<<"Integral all  CB model - "<<  num_all_CB_m  <<endl;
                      cout<<"Integral all  CB model  error - "<< err_all_CB_m <<endl;

                    
                   TFitResultPtr r_both_CB_m = _hf_both_CB1_m->Fit("f_CB_both_m", "RMS+");

                                                               TMatrixDSym covrfit4 = r_both_CB_m->GetCovarianceMatrix(); 
//      covrfit1.Print() ;
      
      TMatrixTSym<double> dataBG4;
      covrfit4.GetSub(0,4,dataBG4);
      
                   for(int j=0; j<8; j++)
                      {
                       p8[j]=f_CB_both_m->GetParameter(j);
                       ep8[j]=f_CB_both_m->GetParError(j);
                       Number_pi0_both_CB_m->SetParameter(j,p8[j]);
                       Number_pi0_both_CB_m->SetParError(j,ep8[j]);
                     }
                   //    Int_t binPi0 = _hf_both_CB1_m->FindBin(kMean);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   //     _hf_both_CB1_m->Fit("f_CB_both_m", "Q");
                      double num_both_CB_m = Number_pi0_both_CB_m->Integral(0.1, 0.25)/_hf_both_CB1_m->GetXaxis()->GetBinWidth(0);
                      double err_both_CB_m = Number_pi0_both_CB_m->IntegralError(0.1, 0.25, r_both_CB_m->GetParams(), dataBG4.GetMatrixArray())/_hf_both_CB1_m->GetXaxis()->GetBinWidth(0);
 
                       cout<<"Integral all  CB model - "<< num_both_CB_m  <<endl;
                       cout<<"Integral all  CB model  error - "<< err_both_CB_m  <<endl;
                        
                        pi0_eff_both_CB_m->SetBinContent(i+1, num_both_CB_m/num_all_CB_m);
                double errBC1_m = err_both_CB_m/num_all_CB_m;
                double errBC2_m = err_all_CB_m*num_both_CB_m/num_all_CB_m/num_all_CB_m;

             pi0_eff_both_CB_m->SetBinError(i+1, TMath::Sqrt(errBC1_m*errBC1_m+errBC2_m*errBC2_m));
        }
             pi0_eff_both_CB->SetBinContent(i+1, numpi0/num_all_CB);
               double errBC1 = errpi0/num_all_CB;
               double errBC2 = err_all_CB*numpi0/num_all_CB/num_all_CB;
 
             pi0_eff_both_CB->SetBinError(i+1, TMath::Sqrt(errBC1*errBC1+errBC2*errBC2));
                

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
               _hf_both_gaus_m->Sumw2();
               _hf_all_m->Sumw2();
                _hf_all_CB->Sumw2();
                _hf_all_CB->Draw();
               _hf_all->Sumw2();
               _hf_all->Draw("same");
              // _hf_both_CB1->Draw();
              // _hf_both_gaus->Draw("same");
           //    _hf_both_gaus_m->Draw("same");//! same
           //    _hf_all_m->Draw("same");
              pi0_eff_both->SetBinContent(i+1, num2pi0/num5pi0);
              double errB1 = err2pi0/num5pi0;
              double errB2 = err5pi0*num2pi0/num5pi0/num5pi0;
              
              pi0_eff_both->SetBinError(i+1, TMath::Sqrt(errB1*errB1+errB2*errB2));

        }

TCanvas *c_eff = new TCanvas("c_eff","The Fit Canvas",1200,800);
pi0_eff_both->SetLineColor(kGreen);
pi0_eff_both->GetXaxis()->SetTitle("p_{T}, GeV/c");
pi0_eff_both->SetStats(0);
pi0_eff_both->SetMaximum(1);
pi0_eff_both->SetMinimum(0);
pi0_eff_both_m->SetLineColor(kRed);
pi0_eff_both->Draw("p e");
pi0_eff_both_m->Draw("p e same");
 auto legend1 = new TLegend(0.1,0.1,0.48,0.2);
   // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
      legend1->AddEntry(pi0_eff_both,"LHC","lep");
      legend1->AddEntry(pi0_eff_both_m,"Monte-Carlo","lep");
 //pi0_eff_both_cut->Draw("p e same");
 //pi0_eff_all_cut->Draw("p e same");
 legend1->Draw("same");


TCanvas *c_eff_m = new TCanvas("c_eff_m","The Fit Canvas",1200,800);
 pi0_eff_both_CB->SetLineColor(kGreen);
 pi0_eff_both_CB_m->SetLineColor(kRed);
 
 pi0_eff_both_CB->SetStats(0);
  pi0_eff_both_CB_m->SetStats(0);
 pi0_eff_both_CB->GetXaxis()->SetTitle("p_{T}, GeV/c");
 pi0_eff_both_CB_m->SetMaximum(1);
 pi0_eff_both_CB_m->SetMinimum(0);


 pi0_eff_both_CB->SetMaximum(1);
 pi0_eff_both_CB->SetMinimum(0);
 pi0_eff_both_CB->Draw("p e");
 pi0_eff_both_CB_m->Draw("p e same");


 auto legend = new TLegend(0.1,0.1,0.48,0.2);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
     legend->AddEntry(pi0_eff_both_CB,"LHC","lep");
     legend->AddEntry(pi0_eff_both_CB_m,"Monte-Carlo","lep");
//pi0_eff_both_cut->Draw("p e same");
//pi0_eff_all_cut->Draw("p e same");
legend->Draw("same");
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
TF1 *f3 = new TF1("f3","[0]*TMath::Exp(-[1]*x)+TMath::Power([2]*(1+[3]*x*x), 13)", 0.5, 20);
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

TF1  *f4  = new TF1("f4","[0]", 0.5, 20);
pi0position->Fit("f4", "RMS");
pi0position->Draw("p e");
pi0positionCB->SetLineColor(kGreen);
pi0positionCB->SetMarkerStyle(22);
pi0positionCB->SetMarkerColor(kGreen);

pi0positionCB->Draw("p e same");

c2->cd(3);
TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*TMath::Power(x, [5])", 0.5, 16);
//TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*x+[5]*x*x", 0.5, 16);
//TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x)+[1]", 0.5, 20);
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

c1->SaveAs("fits_gaus+CB1_LHC.pdf"); 
c2->SaveAs("pi0spectrum_CB1_LHC.pdf");
c_eff->SaveAs("efficiency_pi0_LHC_&_MC.pdf");
c_eff_m->SaveAs("efficiency_pi0_LHC_&_MC crystall Ball.pdf");
    
}
