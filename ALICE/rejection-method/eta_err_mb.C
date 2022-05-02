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

TH1D* fill_hist(TH2D *_h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i);
//std::tuple<double, double>
void IntFit(TCanvas *c,TH1D* _hf, char* name, TH1D* signal, TH1D* signal_bg, double dpt, int i, int count);
void beauty_draw(TH1D* ( &h )[2][2], int i, TString name, TString hist_name);
void draw_err(TH1D* err, TH1D* signal, TH1D* cut_err, TH1D* cut_signal, TH1D* err2, TH1D* signal2, TH1D* cut_err2, TH1D* cut_signal2, TString hist_name, TString hist_pid, TString hist_name2, TString hist_pid2);

void eta_err_mb() {
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("all ","all ",1200,800);

   c1->Divide(2);
   TCanvas *c2 = new TCanvas("Disp&CPV ", "Disp&CPV ", 1200, 800);
   c2->Divide(2);

   TFile *_file0 = TFile::Open("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/LHC16_17.root");//watch the name!
   THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEtaMB;1");

   TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");
   TH2D *_h_all_cut=(TH2D*)_l->FindObject("InvMass_All_Rejection_E300");
   TH2D *_h_both_=(TH2D*)_l->FindObject("InvMass_Both_E300");
   TH2D *_h_both_cut=(TH2D*)_l->FindObject("InvMass_Both_Rejection_E300");

   double pTbins[20];// !!! to be corrected 1.5,3.5,5.5,7.5,9.5

   for (int i=0; i<17; i++)
          {
              pTbins[i]=1.5+i*1.5;
              cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
          }

   TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
   Double_t events=(Double_t)_ev->GetBinContent(2);
    cout<< "!!! Number of events  - " << events<< "!!!"<<endl;

   Int_t bin1=0;
   Int_t bin2=0;

  TH1D* s_all = new TH1D("s+bg all "," all", 16, pTbins);
  TH1D* s_all_cut = new TH1D("s+bg  all cut "," all cut", 16, pTbins);
  TH1D* s_both = new TH1D("s+bg  Disp&CPV "," Disp&CPV", 16, pTbins);
  TH1D* s_both_cut = new TH1D("s+bg  Disp&CPV cut "," Disp&CPV cut", 16, pTbins);

  TH1D* s_all_err = new TH1D("signal all error"," all", 16, pTbins);
  TH1D* s_all_cut_err = new TH1D("signal all cut error"," all cut", 16, pTbins);
  TH1D* s_both_err = new TH1D("signal Disp&CPV error"," Disp&CPV", 16, pTbins);
  TH1D* s_both_cut_err = new TH1D("signal Disp&CPV cut error"," Disp&CPV cut", 16, pTbins);

  TH1D *_hf_all;
  TH1D *_hf_all_cut;
  TH1D *_hf_both_;
  TH1D *_hf_both_cut;


double num_all, err_all;
double num_all_cut, err_all_cut;
double num_both_, err_both_;
double num_both_cut, err_both_cut;

   int count=0;
   TString hist_name;

TH1D * h[2][2];

  for(int i=0; i<6; i++)
    {

        bin1 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i]);
        bin2 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i+1]);

            hist_name.Form("%.2f < p_{T} < %.2f GeV/c",pTbins[i],pTbins[i+1]);
            cout<<"Hist_name - "<< hist_name<<endl;

            _hf_all=fill_hist(_h_all, bin1, bin2, "all", hist_name, i);
            _hf_all_cut=fill_hist(_h_all_cut, bin1, bin2, "all_cut", hist_name, i);
            _hf_both_=fill_hist(_h_both_, bin1, bin2, "Disp&CPV_", hist_name, i);
            _hf_both_cut=fill_hist(_h_both_cut, bin1, bin2, "Disp&CPV_cut", hist_name, i);

            _hf_all->SetLineColor(kBlue);
            _hf_all_cut->SetLineColor(kRed);
            _hf_both_->SetLineColor(kBlue);
            _hf_both_cut->SetLineColor(kRed);


            _hf_all->SetMarkerColor(kBlue);
            _hf_all_cut->SetMarkerColor(kRed);
            _hf_both_->SetMarkerColor(kBlue);
            _hf_both_cut->SetMarkerColor(kRed);

            double dpt = pTbins[i+1]-pTbins[i];

            c1->cd(i+1);
            count++;// 2 or 5
            IntFit(c1,_hf_all, "all", s_all, s_all_err, dpt, i, count);// считает интегралы, ошибки
            IntFit(c1,_hf_all_cut, "all_cut", s_all_cut, s_all_cut_err, dpt, i, count); // считает интегралы, ошибки

            c2->cd(i+1);
            count++;// 3 or 6
            IntFit(c2,_hf_both_, "Disp&CPV_", s_both, s_both_err, dpt, i, count); // считает интегралы, ошибки
            IntFit(c2,_hf_both_cut, "Disp&CPV_cut", s_both_cut, s_both_cut_err, dpt, i, count); // считает интегралы, ошибки

            h[0][0]=_hf_all;
            h[0][1]=_hf_all_cut;
            beauty_draw(h, 0, "All", hist_name);//!
            h[1][0]=_hf_both_;
            h[1][1]=_hf_both_cut;
            beauty_draw(h, 1, "Disp&CPV", hist_name);

    }
    draw_err(s_all_err, s_all, s_all_cut_err, s_all_cut, s_both_err, s_both, s_both_cut_err, s_both_cut, "All signal error", "PID: All signal", "Disp&CPV signal error", "PID: Disp&CPV signal");
}

TH1D* fill_hist(TH2D* _h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i)
{
  TH1D*_hf=(TH1D*)_h->ProjectionX(name+Form("_%d",i),bin1,bin2);
  _hf->Rebin(10);
  _hf->Sumw2();
//  _hf->SetName(name+Form("%d",i));//all SetName-s are different
  _hf->GetXaxis()->SetRangeUser(0.4, 0.7);
  _hf->SetTitle(hist_name);
  _hf->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
  return _hf;
}

void IntFit(TCanvas *c, TH1D* _hf, char* name, TH1D* signal, TH1D* signal_err, double dpt, int i, int count)
{
  //  c->cd(i+1);
    char * fname = Form("f1_%s_%d",name,i);
    char * fname2 = Form("N_%s_%d",name,i);
    TF1  *f1 = new TF1(fname,"gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    f1->SetParameters(1.*_hf->GetBinContent(56),0.552,0.013,0,0,0);//-1e4/(i+1)/(i+1),1e3/(i+1)/(i+1),1.7e4/(i+1)/(i+1));//-2e5/(i+1),3e4/(i+1),2.5*_hf->GetBinContent(51)/(i+1)) ;
    f1->SetParLimits(1, 0.52, 0.58);
    if (count==1 or count==4){
 //     f1->SetParameter(1, 0.55);
//      f1->SetParLimits(2, 0.016, 0.025);
 //     f1->SetParameter(2, 0.2);
    }
    f1->SetParLimits(0, 0,2*_hf->GetBinContent(56));
    _hf->Fit(fname, "0");
  TFitResultPtr r = _hf->Fit(fname, "RMS+");//check if f1 works

  Double_t p[3],ep[3];
  TMatrixDSym covrfit = r->GetCovarianceMatrix();
  TMatrixTSym<double> dataBG;
  covrfit.GetSub(0,2,dataBG);///
  TF1 *Number= new TF1(fname2,"gaus", 0.4, 0.7);
  for(int j=0; j<3; j++)
  {
      p[j]=f1->GetParameter(j);
      ep[j]=f1->GetParError(j);
      Number->SetParameter(j,p[j]);//to be corrected later
  }
  double mean = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);

  double num = f1->Integral(mean-2.*sigma, mean+2.*sigma)/_hf->GetXaxis()->GetBinWidth(0);
  double num_err = f1->IntegralError(mean-2.*sigma, mean+2.*sigma,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/_hf->GetXaxis()->GetBinWidth(0);

  cout<<"Integral: "<< num  <<endl;
  cout<<"Integral error: "<< num_err <<endl;

//  cout<<"DPT"<<dpt<<endl;
//  etaspectrum->SetBinContent(i+1, num / dpt);
//  etaspectrum->SetBinError(i+1, num_err / dpt);

  double signal_v = Number->Integral(mean-2.*sigma, mean+2.*sigma)/_hf->GetXaxis()->GetBinWidth(0);
  double signal_v_error = Number->IntegralError(mean-2.*sigma, mean+2.*sigma,r->GetParams(), dataBG.GetMatrixArray())/_hf->GetXaxis()->GetBinWidth(0);
  cout<<"Signal: "<<signal_v<<endl;
  cout<<"Signal_error: "<<signal_v_error<<endl;

  signal->SetBinContent(i+1, signal_v/ dpt);//только сигнал
  signal_err->SetBinContent(i+1, signal_v_error/ dpt);

//  signal_bg->SetBinContent(i+1, num/ dpt - signal_v/ dpt);///
//  signal_bg_err->SetBinContent(i+1, sqrt((num_err/ dpt)*(num_err/ dpt) + (signal_v_error/ dpt)*(signal_v_error/ dpt)));///!!!
  //  return  std::make_tuple(num, num_err);
}


void beauty_draw(TH1D* ( &h )[2][2], int i, TString name, TString hist_name)
{
  TCanvas * c1DP = new TCanvas("Sginal1","Signal1",600,800) ;
  gPad->SetBorderMode(0);
  gPad->SetBorderSize(2);
  gPad->SetFrameBorderMode(0);
  gPad->SetLeftMargin(0.18);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(0.03);
  cout<<"Beaty work!"<<endl;
  TH1D * box = new TH1D("box","",100,0.401,0.699) ;
  box->SetMinimum(0.);

  box->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})") ;
  box->SetYTitle("Counts") ;
  box->GetXaxis()->SetTitleOffset(0.9) ;
  box->SetStats(0) ;
  box->GetXaxis()->SetTitleSize(0.055) ;
  box->GetXaxis()->SetLabelSize(0.045) ;
  box->GetYaxis()->SetTitleSize(0.055) ;
  box->GetYaxis()->SetLabelSize(0.045) ;
  box->Draw() ;

  h[i][1]->SetMarkerStyle(20) ;
  h[i][1]->SetMarkerColor(1) ;
  h[i][1]->SetLineColor(1) ;
  h[i][1]->SetLineWidth(2) ;

  h[i][0]->SetMarkerStyle(24) ;
  h[i][0]->SetMarkerColor(kPink) ;
  h[i][0]->SetLineColor(kPink) ;
  h[i][0]->SetLineWidth(2) ;

  h[i][1]->Draw("same") ;//&&
  h[i][0]->Draw("same") ;

  TLegend * leg = new TLegend(0.67,0.87,0.95,0.95) ;
  leg->AddEntry(h[i][0],name+" pairs","p") ;
  leg->AddEntry(h[i][1],"rejected #pi^{0}","p") ;
  leg->SetLineColor(0) ;
  leg->Draw() ;
  TLatex t1 = TLatex();
  t1.SetNDC();
  t1.DrawLatex(0.23, 0.9, "#splitline{ALICE LHC}{work in progress}");//12.01.2021
  t1.SetTextSize(0.018) ;
  t1.Draw() ;
  TLatex t2 = TLatex();
  t2.SetNDC();
  t2.DrawLatex(0.227, 0.825, "12.01.2021");
  t2.SetTextSize(0.018) ;
  t2.Draw() ;
  TLatex t3 = TLatex();
  t3.SetNDC();
  t3.DrawLatex(0.229, 0.775, hist_name);
  t3.SetTextSize(0.018) ;
  t3.Draw() ;
  c1DP->SaveAs("LHC_16_17/MB/"+name+"/"+Form("_%d",i)+".pdf");
}

void draw_err(TH1D* err, TH1D* signal, TH1D* cut_err, TH1D* cut_signal, TH1D* err2, TH1D* signal2, TH1D* cut_err2, TH1D* cut_signal2, TString hist_name, TString hist_pid, TString hist_name2, TString hist_pid2)
{
  TCanvas *c = new TCanvas(hist_name, hist_name,1200,800);
  err->Divide(err, signal,1,1,"B");
  cut_err->Divide(cut_err, cut_signal,1,1,"B");

  err->SetLineColor(kBlue);
  cut_err->SetLineColor(kRed);

  err->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  cut_err->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  err->GetXaxis()->SetRangeUser(0, 10.5);
  cut_err->GetXaxis()->SetRangeUser(0, 10.5);

  err->GetYaxis()->SetTitle("Rel. error");
  cut_err->GetYaxis()->SetTitle("Rel. error");

  err2->Divide(err2, signal2,1,1,"B");
  cut_err2->Divide(cut_err2, cut_signal2,1,1,"B");

  err2->SetLineColor(kBlack);
  cut_err2->SetLineColor(kGreen);

  err2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  cut_err2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  err2->GetXaxis()->SetRangeUser(0, 10.5);
  cut_err2->GetXaxis()->SetRangeUser(0, 10.5);

  err2->GetYaxis()->SetTitle("Rel. error");
  cut_err2->GetYaxis()->SetTitle("Rel. error");

  err->Draw("p x0 e1");
  cut_err->Draw("same x0 p e1");

  err2->Draw("same p x0 e1");
  cut_err2->Draw("same x0 p e1");

  err->Draw("HIST");
  cut_err->Draw("same HIST");

  err2->Draw("same HIST");
  cut_err2->Draw("same HIST");


  err->SetMarkerColor(kBlue);
  cut_err->SetMarkerColor(kRed);

  auto legend1 = new TLegend(0.4, 0.17,0.66,0.32);

  legend1->AddEntry(err, hist_pid,"l");
  legend1->AddEntry(cut_err, hist_pid+" Rejection","l");
  legend1->AddEntry(err2, hist_pid2,"l");
  legend1->AddEntry(cut_err2, hist_pid2+" Rejection","l");
  legend1->Draw("same");

  c->SaveAs("LHC_16_17/MB/"+hist_name+" "+hist_name2+".pdf");

TCanvas *c2 = new TCanvas(hist_name, hist_name,1200,800);

  cut_err->Divide(cut_err, err,1,1,"B");
  cut_err2->Divide(cut_err2, err2,1,1,"B");

  cut_err->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  cut_err2->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  cut_err2->GetYaxis()->SetTitle("err_{rej.#pi^{0}} / err_{all}");
  cut_err->GetYaxis()->SetTitle("err_{rej.#pi^{0}} / err_{all}");

  cut_err2->SetLineColor(kBlue);
  cut_err->SetLineColor(kRed);

  cut_err->Draw("p x0 e1");
  cut_err2->Draw("same p x0 e1");

  cut_err->Draw("HIST");
  cut_err2->Draw("same HIST");

  cut_err->GetYaxis()->SetRangeUser(0.7, 1.);

  auto legend2 = new TLegend(0.4, 0.17,0.66,0.32);

  legend2->AddEntry(cut_err, hist_pid, "l");
  legend2->AddEntry(cut_err2, hist_pid2, "l");

  legend2->Draw("same");

  c2->SaveAs("LHC_16_17/MB/double_err.pdf");

}
