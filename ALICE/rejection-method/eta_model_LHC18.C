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


void eta_model() {
   TCanvas *c1 = new TCanvas("c1_fit2","The Fit Canvas",1200,800);
   TCanvas *c2 = new TCanvas("c2_fit2", "The Fit Canvas", 1200, 800);
   gStyle->SetOptFit(1);
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

   TFile *_file0 = TFile::Open("LHC18.root");//watch the name!
   THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEta;1");

   TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E100");
   TH2D *_h_all_cut=(TH2D*)_l->FindObject("InvMass_All_Rejection_E100");
   TH2D *_h_true=(TH2D*)_l->FindObject("True_eta_Inv_mass");
   TH2D *_h_both=(TH2D*)_l->FindObject("InvMass_Both_E100");
   TH2D *_h_both_cut=(TH2D*)_l->FindObject("InvMass_Both_Rejection_E100");
   TH2D *_h_true_cut=(TH2D*)_l->FindObject("True_eta_Inv_mass_Rejection");


   double pTbins[20];// !!! to be corrected
   for (int i=0; i<16; i++)
          {
              pTbins[i]=1.5+i*1.5;
              cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
          }

   TH1D* etaspectrum_both_cut = new TH1D("etaspectrum_both_cut","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_all_cut = new TH1D("etaspectrum","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_both = new TH1D("etaspectrum","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_all = new TH1D("etaspectrum","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_true = new TH1D("etaspectrum","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_true_cut = new TH1D("etaspectrum","eta spectrum", 4, pTbins);

   TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
   Double_t events=(Double_t)_ev->GetBinContent(2);

    cout<< "!!! Number of events - " << events<< "!!!"<<endl;

   Int_t bin1=0;
   Int_t bin2=0;

   TH1D* eta_eff_both = new TH1D("eta_eff_both","eta efficiency both", 4, pTbins);
   TH1D* eta_eff_both_cut = new TH1D("eta_eff_both_cut","eta efficiency both cut", 4, pTbins);
   TH1D* eta_eff_all_cut = new TH1D("eta_eff_all_cut","eta efficiency all cut", 4, pTbins);
   TH1D* efficiency = new TH1D("efficiency","both cut/both", 4, pTbins);
   TH1D* eta_eff_true = new TH1D("eta_eff_true","eta efficiency true", 4, pTbins);

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
             TH1D *_hf_true=(TH1D *)_h_true->ProjectionX(hist_name+"true",bin1,bin2);
             TH1D *_hf_true_cut=(TH1D *)_h_true_cut->ProjectionX(hist_name+"true cut",bin1,bin2);

            _hf_all->Rebin(4);
            _hf_all_cut->Rebin(4);
            _hf_both->Rebin(4);//6
            _hf_both_cut->Rebin(4);//6
            _hf_true->Rebin(4);//6
            _hf_true_cut->Rebin(4);

            _hf_all->Sumw2();
            _hf_all_cut->Sumw2();
            _hf_both->Sumw2();
            _hf_both_cut->Sumw2();
            _hf_true->Sumw2();
            _hf_true_cut->Sumw2(4);
           c1->cd(i+1);

    TF1  *f1 = new TF1("f1","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f2 = new TF1("f2","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f3 = new TF1("f3","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f4 = new TF1("f4","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    TF1  *f_true = new TF1("f_true","[0]*TMath::Exp(-(x-[1])*(x-[1])/2./[2]/[2])", 0.4, 0.7);
    TF1  *f5g = new TF1("f5g","[0]*TMath::Exp(-(x-[1])*(x-[1])/2./[2]/[2])", 0.4, 0.7);

            f5g->SetParameters(30., 0.5456, 0.015);

            if(i==0)
                f5g->SetParameters(84., 0.55, 0.023);

            if(i==3){
                f5g->SetParameter(0, 4.45627e+00);
                f5g->SetParameter(1, 5.57616e-01);
              //  f5g->SetParLimits(2, 0, 0.05);
                f5g->SetParameter(2, 2.03748e-02);
              }
               TFitResultPtr r_true = _hf_true->Fit("f5g", "RMS+");

                 f1->SetParameters(2.5*_hf_both->GetBinContent(66),0.55,0.017,-2e5/(i+1),3e4/(i+1),2.5*_hf_both->GetBinContent(51)/(i+1)) ;
                 f2->SetParameters(1.*_hf_both->GetBinContent(66),0.5456,0.015,1.,3.,2.*_hf_both->GetBinContent(51)) ;
                 f3->SetParameters(1.*_hf_both->GetBinContent(66),0.5456,0.015,1.,3.,2.*_hf_both->GetBinContent(51)) ;
                 f4->SetParameters(30.,0.5456,0.015,-2e3/(i+1),3e2/(i+1),2.5*_hf_both->GetBinContent(51)/(i+1)) ;//both_cut or just both?
                 f_true->SetParameters(1.*_hf_true_cut->GetBinContent(66),0.5456,0.015,1.,3.,2.*_hf_true_cut->GetBinContent(51)) ;

                 double mean = f5g->GetParameter(1);
                 double sigma = f5g->GetParameter(2);

                 f1->FixParameter(1,mean);
                 f1->FixParameter(2,sigma);
                 f2->FixParameter(1,mean);
                 f2->FixParameter(2,sigma);
                 f3->FixParameter(1,mean);
                 f3->FixParameter(2,sigma);
                 f4->FixParameter(1,mean);
                 f4->FixParameter(2,sigma);
                 f_true->FixParameter(1,mean);
                 f_true->FixParameter(2,sigma);

TFitResultPtr r_all = _hf_all->Fit("f1", "RMS+");
TFitResultPtr r_all_cut = _hf_all_cut->Fit("f2", "RMS");
TFitResultPtr r_both = _hf_both->Fit("f3", "RMS+");
TFitResultPtr r_both_cut = _hf_both_cut->Fit("f4", "RMS+");
TFitResultPtr r_true_cut = _hf_true_cut->Fit("f_true", "RMS+");

              _hf_true->SetName(Form("true %d",i));
              _hf_true->GetXaxis()->SetRangeUser(0.4, 0.7);
              _hf_true->SetTitle(hist_name+"true");
              _hf_true->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV/c^{2}");

              _hf_true_cut->SetName(Form("true cut %d",i));
              _hf_true_cut->GetXaxis()->SetRangeUser(0.4, 0.7);
              _hf_true_cut->SetTitle(hist_name+"true cut");
              _hf_true_cut->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV/c^{2}");


             _hf_all->SetName(Form("all%d",i));
              _hf_all->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_all->SetTitle(hist_name);
             _hf_all->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV/c^{2}");

             _hf_all_cut->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_all_cut->SetName(Form("all_cut%d",i));
             _hf_all_cut->SetTitle(hist_name);
             _hf_all_cut->SetLineColor(kRed);
           /*
              _hf_true->GetXaxis()->SetRangeUser(0.4, 0.7);
              _hf_true->SetName(Form("true",i));
              _hf_true->SetTitle(hist_name);
              _hf_true->SetLineColor(kGreen);
            */
             _hf_both->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_both->SetName(Form("both%i",i));
             _hf_both->SetTitle(hist_name);
             _hf_both->SetLineColor(kOrange);

             _hf_both_cut->GetXaxis()->SetRangeUser(0.4, 0.7);
             _hf_both_cut->SetName(Form("both_cut%i",i));
             _hf_both_cut->SetTitle(hist_name);
             _hf_both_cut->SetLineColor(kBlack);

              _hf_true_cut->SetLineColor(kGreen);
              _hf_true->SetLineColor(kBlack);



             Double_t p1[3],ep1[3], p2[3], ep2[3], p3[3], ep3[3], p4[3], ep4[3], p5[3], ep5[3], p6[3], ep6[3];
   TF1 *Number_eta_both = new TF1("number_eta_both_gaus","gaus", 0.4, 0.7);
   TF1 *Number_eta_both_cut = new TF1("number_eta_both_cut","gaus", 0.4, 0.7);
   TF1 *Number_eta_all_cut = new TF1("number_eta_all_cut","gaus", 0.4, 0.7);
   TF1 *Number_eta_all = new TF1("number_eta_all","gaus", 0.4, 0.7);
   TF1 *Number_eta_true = new TF1("number_eta_true","gaus", 0.4, 0.7);
   TF1 *Number_eta_true_cut = new TF1("number_eta_true cut","gaus", 0.4, 0.7);

        for(int j=0; j<3; j++)
             {
                 p5[j]=f5g->GetParameter(j);
                 ep5[j]=f5g->GetParError(j);
                 Number_eta_true->SetParameter(j,p5[j]);
                 Number_eta_true->SetParError(j, ep5[j]);
             }
             cout<<"par: "<<f1->GetParameter(0)<<endl;
             double num_true = Number_eta_true->Integral(0.4, 0.7)/_hf_true->GetXaxis()->GetBinWidth(0);
             double err_true = Number_eta_true->IntegralError(0.4, 0.7,r_true->GetParams(), r_true->GetCovarianceMatrix().GetMatrixArray())/_hf_true->GetXaxis()->GetBinWidth(0);
             cout<<"Integral true - "<< num_true  <<endl;
             cout<<"Integral true error - "<<err_true <<endl;

             for(int j=0; j<3; j++)
                  {
                      p6[j]=f_true->GetParameter(j);
                      ep6[j]=f_true->GetParError(j);
                      Number_eta_true_cut->SetParameter(j,p6[j]);
                      Number_eta_true_cut->SetParError(j, ep6[j]);
                  }
                  cout<<"par: "<<f1->GetParameter(0)<<endl;
                  double num_true_cut = Number_eta_true_cut->Integral(0.4, 0.7)/_hf_true_cut->GetXaxis()->GetBinWidth(0);
                  double err_true_cut = Number_eta_true_cut->IntegralError(0.4, 0.7,r_true_cut->GetParams(), r_true_cut->GetCovarianceMatrix().GetMatrixArray())/_hf_true_cut->GetXaxis()->GetBinWidth(0);
                  cout<<"Integral true cut - "<< num_true_cut  <<endl;
                  cout<<"Integral true cut error - "<<err_true_cut <<endl;

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



           for(int j=0; j<3; j++)
             {
                 p2[j]=f2->GetParameter(j);
                 ep2[j]=f2->GetParError(j);
                 Number_eta_all_cut->SetParameter(j,p2[j]);

            }
            double num_all_cut = Number_eta_all_cut->Integral(0.4, 0.7)/_hf_all_cut->GetXaxis()->GetBinWidth(0);
            double err_all_cut = Number_eta_all_cut->IntegralError(0.4, 0.7,r_all_cut->GetParams(), r_all_cut->GetCovarianceMatrix().GetMatrixArray())/_hf_all_cut->GetXaxis()->GetBinWidth(0);
            cout<<"Integral all cut - "<< num_all_cut  <<endl;
            cout<<"Integral all error - "<<err_all_cut <<endl;


            for(int j=0; j<3; j++)
            {
                  p3[j]=f3->GetParameter(j);
                  ep3[j]=f3->GetParError(j);
                  Number_eta_both->SetParameter(j,p3[j]);

            }
    //        cout << "par 0 both - " << p3[1] <<endl;
              double num_both = Number_eta_both->Integral(0.4, 0.7)/_hf_both->GetXaxis()->GetBinWidth(0);
              double err_both = Number_eta_both->IntegralError(0.4, 0.7,r_both->GetParams(), r_both->GetCovarianceMatrix().GetMatrixArray())/_hf_both->GetXaxis()->GetBinWidth(0);
              cout<<"Integral both - "<< num_both  <<endl;
              cout<<"Integral both error - "<< err_both <<endl;

            for(int j=0; j<3; j++)
            {
                   p4[j]=f4->GetParameter(j);
                   ep4[j]=f4->GetParError(j);
                   Number_eta_both_cut->SetParameter(j,p4[j]);

            }
               cout << "par 0 both cut - " << p4[0] <<endl;
               cout << "par 1 both cut - " << p4[1] <<endl;
               double num_both_cut = Number_eta_both_cut->Integral(0.4, 0.7)/_hf_both_cut->GetXaxis()->GetBinWidth(0);//!
               double err_both_cut = Number_eta_both_cut->IntegralError(0.4, 0.7,r_both_cut->GetParams(), r_both_cut->GetCovarianceMatrix().GetMatrixArray())/_hf_both_cut->GetXaxis()->GetBinWidth(0);//!
               cout<<"Integral both cut - "<< num_both_cut  <<endl;
               cout<<"Integral both cut error - "<< err_both_cut <<endl;

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
                 etaspectrum_true->SetBinContent(i+1, num_true / dpt);
                 etaspectrum_true->SetBinError(i, err_true / dpt);
                 cout<<"!!TRUE = !!"<< num_true/ dpt<<endl;
                 etaspectrum_true_cut->SetBinContent(i+1, num_true_cut/ dpt);
                 etaspectrum_true_cut->SetBinError(i, err_true_cut / dpt);
                 cout<<"!!TRUE CUT = !!"<< num_true_cut/ dpt<<endl;

            TLatex* t=new TLatex();

              _hf_all->SetMinimum(0.);
              _hf_all->Draw();
              t->DrawLatexNDC(0.12,0.82,Form("Integral (All): %2.1f #pm %2.1f",num_all,err_all));
              _hf_all_cut->Draw("same");
              t->DrawLatexNDC(0.12,0.42,Form("Integral with rejection (All): %2.1f #pm %2.1f",num_all_cut,err_all_cut));
             _hf_true->Draw("same");
              t->DrawLatexNDC(0.12,0.22,Form("Integral true: %2.1f #pm %2.1f",num_true,err_true));
              t->Draw("same");
              _hf_true_cut->Draw("same");
             t->DrawLatexNDC(0.12,0.12,Form("Integral true cut: %2.1f #pm %2.1f",num_true_cut,err_true_cut));
             t->Draw("same");

              c2->cd(i+1);
              _hf_both->SetMinimum(0.);
             _hf_both->Draw();
              t->DrawLatexNDC(0.19,0.72,Form("Integral both: %2.1f #pm %2.1f",num_both,err_both));
              _hf_both_cut->Draw("same");
             t->DrawLatexNDC(0.12,0.42,Form("Integral with rejection (Both): %2.1f #pm %2.1f", num_both_cut, err_both_cut));
              _hf_true->Draw("same");
             t->DrawLatexNDC(0.12,0.22,Form("Integral true: %2.1f #pm %2.1f",num_true,err_true));
             t->Draw("same");
             _hf_true_cut->Draw("same");
            t->DrawLatexNDC(0.05,0.12,Form("Integral true cut: %2.1f #pm %2.1f",num_true_cut,err_true_cut));
            t->Draw("same");

               if(num_all>0)
               {
                cout<<"true eff"<<num_true_cut/num_true<<endl;

                eta_eff_true->SetBinContent(i+1, num_true_cut/num_true);
                double errT1 = err_true_cut/num_true;
                double errT2 = err_true*num_true_cut/num_true/num_true;
                cout << " eta_eff_all_cut " << num_true_cut/num_true<<endl;
                eta_eff_true->SetBinError(i+1, TMath::Sqrt(errT1*errT1+errT2*errT2));



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

                efficiency->SetBinContent(i+1, (num_both_cut/num_all)/(num_both/num_all));

                double erref1 = err_bc/(num_both/num_all);
                double erref2 = err_b*(num_both_cut/num_all)/(num_both/num_all)/(num_both/num_all);
                efficiency->SetBinError(i+1, TMath::Sqrt(erref1*erref1+erref2*erref2));
             }
    }

 TCanvas *c3 = new TCanvas("c3_fit2", "The Fit Canvas", 1200, 800);
 TLatex* t2=new TLatex();

 gStyle->SetOptStat(0);
 eta_eff_both->GetXaxis()->SetRangeUser(0., 10.);
 eta_eff_both->GetXaxis()->SetTitle("p_{T}, GeV/c");
 eta_eff_both->GetYaxis()->SetTitle("Efficiency for E100");
 eta_eff_both->SetMarkerStyle(21);
 eta_eff_both_cut->SetMarkerStyle(31);
 eta_eff_all_cut->SetMarkerStyle(41);
 efficiency->SetMarkerStyle(5);
 eta_eff_true->SetMarkerStyle(34);

 eta_eff_true->SetMaximum(1.);

  eta_eff_both->SetMarkerColor(kGreen);
  eta_eff_both_cut->SetMarkerColor(kRed);
  eta_eff_all_cut->SetMarkerColor(kBlue);
  efficiency->SetMarkerColor(6);
  eta_eff_true->SetMarkerColor(kBlack);

  eta_eff_both->SetMarkerSize(2);
  eta_eff_both_cut->SetMarkerSize(2);
  eta_eff_all_cut->SetMarkerSize(2);
 efficiency->SetMarkerSize(2);
 eta_eff_true->SetMarkerSize(2);

 eta_eff_both->SetLineColor(kRed);
 eta_eff_both_cut->SetLineColor(kGreen);
 eta_eff_all_cut->SetLineColor(6);
 efficiency->SetLineColor(kBlue);
 eta_eff_true->SetLineColor(kBlack);
// eta_eff_both->Draw("p e");
//efficiency->Draw("p e same");
 eta_eff_true->Draw("p e");

 auto legend = new TLegend(0.099,0.1,0.369,0.3);
  legend->SetHeader("Efficiency E100","C"); // option "C" allows to center the header
/*  legend->AddEntry(eta_eff_both,"Efficiency both", "lep");
  legend->AddEntry(eta_eff_both_cut,"Efficiency both rejection","lep");
  legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection","lep");
  legend->AddEntry(efficiency,"(both rejection)/both","lep");
*/  legend->AddEntry(eta_eff_true,"(true rejection)/true","lep");
 // t2->DrawLatexNDC(0.12,0.72,"All cut");
 //eta_eff_both_cut->Draw("p e same");
 // t2->DrawLatexNDC(0.12,0.62, "Both");
 //eta_eff_all_cut->Draw("p e same");
 legend->Draw("same");
 //t2->DrawLatexNDC(0.12,0.52, "Both cut");
 TCanvas *c4 = new TCanvas("c4_fit2", "The Fit Canvas", 1200, 800);

 gPad->SetLogy();
 etaspectrum_both_cut->Scale(1./events);
 etaspectrum_both_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
 etaspectrum_both_cut->GetYaxis()->SetTitle("dN/dp_{T}");
 etaspectrum_both_cut->SetMarkerStyle(21);

 etaspectrum_true_cut->Scale(1./events);
 etaspectrum_true_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
 etaspectrum_true_cut->GetYaxis()->SetTitle("dN/dp_{T}");
 etaspectrum_true_cut->SetMarkerStyle(71);

 etaspectrum_true->Scale(1./events);
 etaspectrum_true->GetXaxis()->SetTitle("p_{T}, GeV/c");
 etaspectrum_true->GetYaxis()->SetTitle("dN/dp_{T}");
 etaspectrum_true->SetMarkerStyle(91);

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


  auto legend_sp = new TLegend(0.1,0.1,0.38,0.3);
 // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend_sp->AddEntry(etaspectrum_both_cut,"spectrum both cut","lep");
    legend_sp->AddEntry(etaspectrum_both,"spectrum both","lep");
    legend_sp->AddEntry(etaspectrum_all_cut,"spectrum all cut","lep");
    legend_sp->AddEntry(etaspectrum_all,"spectrum all","lep");
    legend_sp->AddEntry(etaspectrum_true,"spectrum true","lep");
    legend_sp->AddEntry(etaspectrum_true_cut,"spectrum true cut","lep");
// TF1 *f6 = new TF1("f6","[0]*TMath::Exp(-[1]*x)", 6, 16);
// etaspectrum->Fit("f6", "RMS");
 etaspectrum_both->Draw("p e");
// etaspectrum_both_cut->Draw("p e same");
 etaspectrum_all_cut->Draw("p e same");
 etaspectrum_all->Draw("p e same");
 etaspectrum_true->Draw("p e same");
 etaspectrum_true_cut->Draw("p e same");
 legend_sp->Draw("p e same");
// TF1 *f6 = new TF1("f6","[0]*TMath::Exp(-[1]*x)", 6, 16);
// etaspectrum_both_cut->Fit("f6", "RMS");


c1->SaveAs("eta_hist_all_model_E100.pdf");
c2->SaveAs("eta_hist_both_model_E100.pdf");
c4->SaveAs("hist_eta_spectrum_model_E100.pdf");
c3->SaveAs("eta_efficiency_model_E100.pdf");
}
