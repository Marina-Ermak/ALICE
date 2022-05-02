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
#include <ctime>


void eta_model_mult_LHC() {

   TCanvas *c1 = new TCanvas("c1_fit2","The Fit Canvas",1200,800);
   TCanvas *c2 = new TCanvas("c2_fit2", "The Fit Canvas", 1200, 800);
   TCanvas *c3 = new TCanvas("c3_fit2", "The Fit Canvas", 1200, 800);




   gStyle->SetOptFit(1);
   c1->Divide(2,3);
   c2->Divide(2,5);
   c3->Divide(3,2);

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

     c3->SetGridx();
     c3->SetGridy();
     c3->GetFrame()->SetFillColor(21);
     c3->GetFrame()->SetBorderMode(-1);
     c3->GetFrame()->SetBorderSize(5);

   TFile *_file0 = TFile::Open("LHC16.root");//watch the name!
   THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEtaMB;1");

   TH2D *_h_all;
   TH2D *_h_all_cut;
   TH2D *_h_both;
   TH2D *_h_both_cut;

   string All="InvMass_All_E100";
   string Both="InvMass_Both_E100";
   string All_R="InvMass_All_Rejection_E100";
   string Both_R="InvMass_Both_Rejection_E100";

   TH1D *_hf_all;
   TH1D *_hf_all_cut;
   TH1D *_hf_both;
   TH1D *_hf_both_cut;

   double pTbins[5];// !!! to be corrected
   for (int i=0; i<5; i++)
          {
              pTbins[i]=1.5+i*2;
              cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
          }

   TH1D* etaspectrum_both_cut = new TH1D("etaspectrum_both_cut","eta spectrum", 4, pTbins);//WATCH NUMBER OF BINS!!!
   TH1D* etaspectrum_all_cut = new TH1D("etaspectrum","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_both = new TH1D("etaspectrum","eta spectrum", 4, pTbins);
   TH1D* etaspectrum_all = new TH1D("etaspectrum","eta spectrum", 4, pTbins);

   TH1D* eta_eff_both;
   TH1D* eta_eff_both_cut;
   TH1D* eta_eff_all_cut;
   TH1D* efficiency;

   TF1  *f1 = new TF1("f1","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
   TF1  *f2 = new TF1("f2","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
   TF1  *f3 = new TF1("f3","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
   TF1  *f4 = new TF1("f4","gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);

   TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
   Double_t events=(Double_t)_ev->GetBinContent(2);

    cout<< "!!! Number of events - " << events<< "!!!"<<endl;

    Double_t p1[3],ep1[3], p2[3], ep2[3], p3[3], ep3[3], p4[3], ep4[3], p5[3], ep5[3], p6[3], ep6[3];


   Int_t bin1=0;
   Int_t bin2=0;
   TString hist_name;
   TString s;


   double mean;
   double sigma;

   TFitResultPtr r_all;
   TFitResultPtr r_all_cut;
   TFitResultPtr r_both;
   TFitResultPtr r_both_cut;

   TF1 *Number_eta_both = new TF1("number_eta_both_gaus","gaus", 0.4, 0.7);
   TF1 *Number_eta_both_cut = new TF1("number_eta_both_cut","gaus", 0.4, 0.7);
   TF1 *Number_eta_all_cut = new TF1("number_eta_all_cut","gaus", 0.4, 0.7);
   TF1 *Number_eta_all = new TF1("number_eta_all","gaus", 0.4, 0.7);
   auto legend = new TLegend(0.66,0.1,0.9,0.9);
   int count=1;
   for(int j=0; j<5; j++)///& j=0 or j=1?
   {
     s=j | 0x30;//int to char
     eta_eff_both = new TH1D("eta_eff_both","#eta efficiency M"+s, 4, pTbins);//WATCH NUMBER OF BINS!!!
     eta_eff_both_cut = new TH1D("eta_eff_both_cut","#eta efficiency both rejection M"+s, 4, pTbins);
     eta_eff_all_cut = new TH1D("eta_eff_all_cut","#eta efficiency all rejection M"+s, 4, pTbins);
     efficiency = new TH1D("efficiency","both rejection/both M"+s, 4, pTbins);

     if(j==0){
      eta_eff_both = new TH1D("eta_eff_both","0 < M < 5 PHOS clusters", 4, pTbins);
      eta_eff_both_cut = new TH1D("eta_eff_both_rejection","0 < M < 5 PHOS clusters", 4, pTbins);
      eta_eff_all_cut = new TH1D("eta_eff_all_rejection","0 < M < 5 PHOS clusters", 4, pTbins);
      efficiency = new TH1D("eta_eff","0 < M < 5 PHOS clusters", 4, pTbins);
    }else if(j==1){
      eta_eff_both = new TH1D("eta_eff_both","Efficiency for both", 4, pTbins);
      eta_eff_both_cut = new TH1D("eta_eff_both_rejection","Efficiency for both rejection", 4, pTbins);
      eta_eff_all_cut = new TH1D("eta_eff_all_rejection","Efficiency for all rejection", 4, pTbins);
      efficiency = new TH1D("eta_eff","both rejection/all", 4, pTbins);
    }else if(j==2){
      eta_eff_both = new TH1D("eta_eff_both","10 < M < 20 PHOS clusters", 4, pTbins);
      eta_eff_both_cut = new TH1D("eta_eff_both_rejection","10 < M < 20 PHOS clusters", 4, pTbins);
      eta_eff_all_cut = new TH1D("eta_eff_all_rejection","10 < M < 20 PHOS clusters", 4, pTbins);
      efficiency = new TH1D("eta_eff","10 < M < 20 PHOS clusters", 4, pTbins);
    }else if(j==3){
       eta_eff_both = new TH1D("eta_eff_both","20 < M < 30 PHOS clusters", 4, pTbins);
       eta_eff_both_cut = new TH1D("eta_eff_both_rejection","20 < M < 30 PHOS clusters", 4, pTbins);
       eta_eff_all_cut = new TH1D("eta_eff_all_rejection","20 < M < 30 PHOS clusters", 4, pTbins);
       efficiency = new TH1D("eta_eff","20 < M < 30 PHOS clusters", 4, pTbins);
    }else if(j==4){
        eta_eff_both = new TH1D("eta_eff_both","M > 30 PHOS clusters", 4, pTbins);
        eta_eff_both_cut = new TH1D("eta_eff_both_rejection","M > 30 PHOS clusters", 4, pTbins);
        eta_eff_all_cut = new TH1D("eta_eff_all_rejection","M > 30 PHOS clusters", 4, pTbins);
        efficiency = new TH1D("eta_eff","M > 30 PHOS clusters", 4, pTbins);
    }
    _h_all=(TH2D*)_l->FindObject(All+"_M"+s);
    _h_all_cut=(TH2D*)_l->FindObject(All_R+"_M"+s);
    _h_both=(TH2D*)_l->FindObject(Both+"_M"+s);
    _h_both_cut=(TH2D*)_l->FindObject(Both_R+"_M"+s);
     for(int i=0; i<1; i++)
       {
           bin1 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i]);
           bin2 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i+1]);

               hist_name.Form("Projection %.2f < p_{T} < %.2f GeV/c M"+s+" ",pTbins[i],pTbins[i+1]);
               cout<<"Hist_name - "<< hist_name<<endl;
               _hf_all=_h_all->ProjectionX(hist_name,bin1,bin2);
               _hf_all_cut=_h_all_cut->ProjectionX(hist_name+"cut",bin1,bin2);
               _hf_both=_h_both->ProjectionX(hist_name+"both",bin1,bin2);
               _hf_both_cut=_h_both_cut->ProjectionX(hist_name+"bothcut",bin1,bin2);

               _hf_all->Rebin(4);
               _hf_all_cut->Rebin(4);
               _hf_both->Rebin(4);
               _hf_both_cut->Rebin(4);

               _hf_all->Sumw2();
               _hf_all_cut->Sumw2();
               _hf_both->Sumw2();
               _hf_both_cut->Sumw2();

              _hf_all->SetName(Form("all%d",i));
              _hf_all->GetXaxis()->SetRangeUser(0.4, 0.7);
              _hf_all->SetTitle(hist_name);
              _hf_all->GetXaxis()->SetTitle("m_{#gamma#gamma}, GeV/c^{2}");
              _hf_all->SetLineColor(kBlue);

              _hf_all_cut->GetXaxis()->SetRangeUser(0.4, 0.7);
              _hf_all_cut->SetName(Form("all_cut%d",i));
              _hf_all_cut->SetTitle(hist_name);
              _hf_all_cut->SetLineColor(kRed);

              _hf_both->GetXaxis()->SetRangeUser(0.4, 0.7);
              _hf_both->SetName(Form("both%i",i));
              _hf_both->SetTitle(hist_name);
              _hf_both->SetLineColor(kOrange);

              _hf_both_cut->GetXaxis()->SetRangeUser(0.4, 0.7);
              _hf_both_cut->SetName(Form("both_cut%i",i));
              _hf_both_cut->SetTitle(hist_name);
              _hf_both_cut->SetLineColor(kBlack);

                 f1->SetParameters(2.5*_hf_both->GetBinContent(66),0.55,0.017,-2e5/(i+1),3e4/(i+1),2.5*_hf_both->GetBinContent(51)/(i+1)) ;
                 f2->SetParameters(1.*_hf_both->GetBinContent(66),0.5456,0.015,1.,3.,2.*_hf_both->GetBinContent(51)) ;
                 f3->SetParameters(1.*_hf_both->GetBinContent(66),0.5456,0.015,1.,3.,2.*_hf_both->GetBinContent(51)) ;
                 f4->SetParameters(30.,0.5456,0.015,-2e3/(i+1),3e2/(i+1),2.5*_hf_both_cut->GetBinContent(51)/(i+1)) ;//both_cut or just both?

                 /*
                 f1->FixParameter(1,mean);
                 f1->FixParameter(2,sigma);
                 f2->FixParameter(1,mean);
                 f2->FixParameter(2,sigma);
                 f3->FixParameter(1,mean);
                 f3->FixParameter(2,sigma);
                 f4->FixParameter(1,mean);
                 f4->FixParameter(2,sigma);
                 */
   r_all = _hf_all->Fit("f1", "RMS+");
   r_all_cut = _hf_all_cut->Fit("f2", "RMS");
   r_both = _hf_both->Fit("f3", "RMS+");
   r_both_cut = _hf_both_cut->Fit("f4", "RMS+");

               for(int k=0; k<3; k++)
               {
                   p1[k]=f1->GetParameter(k);
                   ep1[k]=f1->GetParError(k);
                   Number_eta_all->SetParameter(k,p1[k]);

               }
               cout<<"par: "<<f1->GetParameter(0)<<endl;
               double num_all = Number_eta_all->Integral(0.4, 0.7)/_hf_all->GetXaxis()->GetBinWidth(0);
               double err_all = Number_eta_all->IntegralError(0.4, 0.7,r_all->GetParams(), r_all->GetCovarianceMatrix().GetMatrixArray())/_hf_all->GetXaxis()->GetBinWidth(0);
               cout<<"Integral all - "<< num_all  <<endl;
               cout<<"Integral all error - "<<err_all <<endl;

               for(int k=0; k<3; k++)
                {
                    p2[k]=f2->GetParameter(k);
                    ep2[k]=f2->GetParError(k);
                    Number_eta_all_cut->SetParameter(k,p2[k]);

               }

               double num_all_cut = Number_eta_all_cut->Integral(0.4, 0.7)/_hf_all_cut->GetXaxis()->GetBinWidth(0);
               double err_all_cut = Number_eta_all_cut->IntegralError(0.4, 0.7,r_all_cut->GetParams(), r_all_cut->GetCovarianceMatrix().GetMatrixArray())/_hf_all_cut->GetXaxis()->GetBinWidth(0);
               cout<<"Integral all cut - "<< num_all_cut  <<endl;
               cout<<"Integral all error - "<<err_all_cut <<endl;

               for(int k=0; k<3; k++)
               {
                     p3[k]=f3->GetParameter(k);
                     ep3[k]=f3->GetParError(k);
                     Number_eta_both->SetParameter(k,p3[k]);

               }
       //        cout << "par 0 both - " << p3[1] <<endl;
                 double num_both = Number_eta_both->Integral(0.4, 0.7)/_hf_both->GetXaxis()->GetBinWidth(0);
                 double err_both = Number_eta_both->IntegralError(0.4, 0.7,r_both->GetParams(), r_both->GetCovarianceMatrix().GetMatrixArray())/_hf_both->GetXaxis()->GetBinWidth(0);
                 cout<<"Integral both - "<< num_both  <<endl;
                 cout<<"Integral both error - "<< err_both <<endl;

               for(int k=0; k<3; k++)
               {
                      p4[k]=f4->GetParameter(k);
                      ep4[k]=f4->GetParError(k);
                      Number_eta_both_cut->SetParameter(k,p4[k]);

               }
                  cout << "par 0 both cut - " << p4[0] <<endl;
                  cout << "par 1 both cut - " << p4[1] <<endl;
                  double num_both_cut = Number_eta_both_cut->Integral(0.4, 0.7)/_hf_both_cut->GetXaxis()->GetBinWidth(0);//!
                  double err_both_cut = Number_eta_both_cut->IntegralError(0.4, 0.7,r_both_cut->GetParams(), r_both_cut->GetCovarianceMatrix().GetMatrixArray())/_hf_both_cut->GetXaxis()->GetBinWidth(0);//!
                  cout<<"Integral both cut - "<< num_both_cut  <<endl;
                  cout<<"Integral both cut error - "<< err_both_cut <<endl;

                   double dpt = pTbins[i+1]-pTbins[i];

                   etaspectrum_both_cut->SetBinContent(i+1, num_both_cut / dpt);
                   etaspectrum_both_cut->SetBinError(i, err_both_cut / dpt);

                    etaspectrum_both->SetBinContent(i+1, num_both / dpt);
                    etaspectrum_both->SetBinError(i, err_both / dpt);

                    etaspectrum_all_cut->SetBinContent(i+1, num_all_cut / dpt);
                    etaspectrum_all_cut->SetBinError(i, err_all_cut / dpt);

                    etaspectrum_all->SetBinContent(i+1, num_all / dpt);
                    etaspectrum_all->SetBinError(i, err_all / dpt);

                    TLatex* t=new TLatex();
                   c1->cd(count);
                  _hf_all->Draw();
                  _hf_both->Draw("same");
                  _hf_all_cut->Draw("same");
                  _hf_both_cut->Draw("same");
              /*  t->DrawLatexNDC(0.12,0.22,Form("Integral true: %2.1f #pm %2.1f",num_true,err_true));
                t->Draw("same");
                t->DrawLatexNDC(0.07,0.12,Form("Integral true cut: %2.1f #pm %2.1f",num_true_cut,err_true_cut));
                t->Draw("same");
*/
                 count++;
                 cout<<"count - "<<count<<endl;

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

                 efficiency->SetBinContent(i+1, (num_both_cut/num_all)/(num_both/num_all));

                 double erref1 = err_bc/(num_both/num_all);
                 double erref2 = err_b*(num_both_cut/num_all)/(num_both/num_all)/(num_both/num_all);
                 efficiency->SetBinError(i+1, TMath::Sqrt(erref1*erref1+erref2*erref2));

                }
       }
       c3->cd(0);
       TLatex* t2=new TLatex();
       gStyle->SetOptStat(0);

       eta_eff_all_cut->SetMinimum(0.);
       eta_eff_all_cut->SetMaximum(1);
       eta_eff_both_cut->SetMinimum(0.);
       eta_eff_both_cut->SetMaximum(1);
       eta_eff_both->SetMinimum(0.);
       eta_eff_both->SetMaximum(1);
       efficiency->SetMaximum(1);

       eta_eff_both->GetXaxis()->SetRangeUser(0., 10.);
       eta_eff_both->GetXaxis()->SetTitle("p_{T}, GeV/c");
       eta_eff_both_cut->GetXaxis()->SetRangeUser(0., 10.);
       eta_eff_both_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
       eta_eff_all_cut->GetXaxis()->SetRangeUser(0., 10.);
       eta_eff_all_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");

       eta_eff_both->GetYaxis()->SetTitle("Efficiency for M"+s);

        eta_eff_both->SetMarkerSize(2);
        eta_eff_both_cut->SetMarkerSize(2);
        eta_eff_all_cut->SetMarkerSize(2);
       efficiency->SetMarkerSize(2);


       if(j==1){
      //  eta_eff_true->Draw("p e");
        // eta_eff_both_cut->Draw("p e");
        //eta_eff_all_cut->Draw("p e");
        //  eta_eff_both->Draw("p e");
          efficiency->Draw("p e");
     }else{
         //eta_eff_true->Draw("p e same");true_cut
          //eta_eff_both_cut->Draw("p e same");
          //eta_eff_all_cut->Draw("p e same");
      //  eta_eff_both->Draw("p e same");
            efficiency->Draw("p e same");
  }
      //  legend->SetHeader("Efficiency M"+s,"C"); // option "C" allows to center the header         //legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection M"+s,"lep");

       if(j==0){


        eta_eff_both->SetMarkerStyle(21);
        eta_eff_both_cut->SetMarkerStyle(31);
        eta_eff_all_cut->SetMarkerStyle(41);
        efficiency->SetMarkerStyle(51);

         eta_eff_both->SetMarkerColor(kRed);
         eta_eff_both_cut->SetMarkerColor(kBlue);
         eta_eff_all_cut->SetMarkerColor(kBlack);
         efficiency->SetMarkerColor(kGreen);

         eta_eff_both->SetLineColor(kBlue);
         eta_eff_both_cut->SetLineColor(kGreen);
         eta_eff_all_cut->SetLineColor(kRed);
         efficiency->SetLineColor(kBlack);

      //  legend->AddEntry(eta_eff_both_cut,"Efficiency both rejection 0 < M < 5 PHOS clusters","lep");
      //  legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection 0 < M < 5 PHOS clusters","lep");
      // legend->AddEntry(eta_eff_both,"Efficiency both 0 < M < 5 PHOS clusters", "lep");
         legend->AddEntry(efficiency,"(both rejection)/all 0 < M < 5 PHOS clusters" ,"lep");

      }else if(j==1){
        eta_eff_both->SetMarkerStyle(20);
        eta_eff_both_cut->SetMarkerStyle(21);
        eta_eff_all_cut->SetMarkerStyle(22);
        efficiency->SetMarkerStyle(23);

         eta_eff_both->SetMarkerColor(kBlue);
         eta_eff_both_cut->SetMarkerColor(kBlack);
         eta_eff_all_cut->SetMarkerColor(kGreen);
         efficiency->SetMarkerColor(kOrange);

         eta_eff_both->SetLineColor(kBlue);
         eta_eff_both_cut->SetLineColor(kBlack);
         eta_eff_all_cut->SetLineColor(kGreen);
         efficiency->SetLineColor(kOrange);

          // legend->AddEntry(eta_eff_both_cut,"Efficiency both rejection 5 < M < 10 PHOS clusters","lep");
      //   legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection 5 < M < 10 PHOS clusters","lep");
        //legend->AddEntry(eta_eff_both,"Efficiency both 5 < M < 10 PHOS clusters", "lep");
           legend->AddEntry(efficiency,"(both rejection)/all 5 < M < 10 PHOS clusters" ,"lep");
      }else if(j==2){
        eta_eff_both->SetMarkerStyle(21);
        eta_eff_both_cut->SetMarkerStyle(22);
        eta_eff_all_cut->SetMarkerStyle(23);
        efficiency->SetMarkerStyle(34);

         eta_eff_both->SetMarkerColor(kBlack);
         eta_eff_both_cut->SetMarkerColor(kGreen);
         eta_eff_all_cut->SetMarkerColor(kOrange);
         efficiency->SetMarkerColor(kRed);

         eta_eff_both->SetLineColor(kBlack);
         eta_eff_both_cut->SetLineColor(kGreen);
         eta_eff_all_cut->SetLineColor(kOrange);
         efficiency->SetLineColor(kRed);

          //   legend->AddEntry(eta_eff_both_cut,"Efficiency both rejection 10 < M < 20 PHOS clusters","lep");
          // legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection 10 < M < 20 PHOS clusters","lep");
          // legend->AddEntry(eta_eff_both,"Efficiency both 10 < M < 20 PHOS clusters", "lep");
          legend->AddEntry(efficiency,"(both rejection)/all 10 < M < 20 PHOS clusters" ,"lep");

      }else if(j==3){
        eta_eff_both->SetMarkerStyle(22);
        eta_eff_both_cut->SetMarkerStyle(23);
        eta_eff_all_cut->SetMarkerStyle(34);
        efficiency->SetMarkerStyle(20);

         eta_eff_both->SetMarkerColor(kGreen);
         eta_eff_both_cut->SetMarkerColor(kOrange);
         eta_eff_all_cut->SetMarkerColor(kRed);
         efficiency->SetMarkerColor(kBlue);

         eta_eff_both->SetLineColor(kGreen);
         eta_eff_both_cut->SetLineColor(kOrange);
         eta_eff_all_cut->SetLineColor(kRed);
         efficiency->SetLineColor(kBlue);

          // legend->AddEntry(eta_eff_both_cut,"Efficiency both rejection 20 < M < 30 PHOS clusters","lep");
        //  legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection 20 < M < 30 PHOS clusters","lep");
        //legend->AddEntry(eta_eff_both,"Efficiency both 20 < M < 30 PHOS clusters", "lep");
           legend->AddEntry(efficiency,"(both rejection)/all 20 < M < 30 PHOS clusters" ,"lep");

      }else if(j==4){
          eta_eff_both->SetMarkerStyle(23);
          eta_eff_both_cut->SetMarkerStyle(34);
          eta_eff_all_cut->SetMarkerStyle(20);
          efficiency->SetMarkerStyle(21);

           eta_eff_both->SetMarkerColor(kOrange);
           eta_eff_both_cut->SetMarkerColor(kRed);
           eta_eff_all_cut->SetMarkerColor(kBlue);
           efficiency->SetMarkerColor(kBlack);

           eta_eff_both->SetLineColor(kOrange);
           eta_eff_both_cut->SetLineColor(kRed);
           eta_eff_all_cut->SetLineColor(kBlue);
           efficiency->SetLineColor(kBlack);

          //  legend->AddEntry(eta_eff_both_cut,"Efficiency both rejection M > 30 PHOS clusters","lep");
            // legend->AddEntry(eta_eff_all_cut,"Efficiency all rejection M > 30 PHOS clusters","lep");
          //  legend->AddEntry(eta_eff_both,"Efficiency both M > 30 PHOS clusters", "lep");
            legend->AddEntry(efficiency,"(both rejection)/all M > 30 PHOS clusters" ,"lep");
      }
        legend->Draw("same");


       TCanvas *c4 = new TCanvas("c4_fit2", "The Fit Canvas", 1200, 800);
       c4->cd(0);
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

        etaspectrum_all_cut->Scale(1./events);c1->cd(count);
        etaspectrum_all_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
        etaspectrum_all_cut->GetYaxis()->SetTitle("dN/dp_{T}");
        etaspectrum_all_cut->SetMarkerStyle(41);
        etaspectrum_all_cut->SetLineColor(kOrange);

        etaspectrum_all->Scale(1./events);
        etaspectrum_all->GetXaxis()->SetTitle("p_{T}, GeV/c");
        etaspectrum_all->GetYaxis()->SetTitle("dN/dp_{T}");
        etaspectrum_all->SetMarkerStyle(51);
        etaspectrum_all->SetLineColor(kRed);


        auto legend_sp = new TLegend(0.1,0.1,0.48,0.2);
       // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
          legend_sp->AddEntry(etaspectrum_both_cut,"spectrum both cut","lep");
          legend_sp->AddEntry(etaspectrum_both,"spectrum both","lep");
          legend_sp->AddEntry(etaspectrum_all_cut,"spectrum all cut","lep");
          legend_sp->AddEntry(etaspectrum_all,"spectrum all","lep");
      // TF1 *f6 = new TF1("f6","[0]*TMath::Exp(-[1]*x)", 6, 16);
      // etaspectrum->Fit("f6", "RMS");
       etaspectrum_both->Draw("p e");
      // etaspectrum_both_cut->Draw("p e same");
       etaspectrum_all_cut->Draw("p e same");
       etaspectrum_all->Draw("p e same");
       legend_sp->Draw("p e same");
      // TF1 *f6 = new TF1("f6","[0]*TMath::Exp(-[1]*x)", 6, 16);
      // etaspectrum_both_cut->Fit("f6", "RMS");
      c2->SaveAs("eta_hist_both_model_E100_M"+s+".pdf");
      c4->SaveAs("hist_eta_spectrum_model_E100_M"+s+".pdf");
     }

     c3->SaveAs("eta_efficiency_E100.pdf");
     c1->SaveAs("eta_hist_all_model_E100_M"+s+"1.5_3.5pT.pdf");
 }
