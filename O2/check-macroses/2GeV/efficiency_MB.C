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
#include <tuple>

TH1D* fill_hist(TH2D *_h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i);
//std::tuple<double, double>
void IntFit(TCanvas *c,TH1D* _hf, char* name, TH1D* signal, TH1D* signal_bg, double dpt, int i, int count, TH1D * etaspectrum);
void drawing(TCanvas *c, int i, TH1D *_hf1, TH1D *_hf2,  TH1D *_hf3, TH1D *_hf4,  TH1D *s1, TH1D *s2, TH1D *s3, TH1D *s4, TString name1, TString name2, TString name3, TString name4, int count);
void draw_sg(double up_lim, double down_lim, TH1D*s_all, TH1D*s_all_cut, TH1D* s_both, TH1D* s_both_cut);//bad variant
void draw_eff(double up_lim, double down_lim, TH1D*s_all, TH1D*s_both, TH1D* s_true, Double_t events);
void draw_spec(TH1D * true_h, TH1D * true_cut, TH1D * all, TH1D * all_cut, TH1D * both, TH1D * both_cut, Double_t events, TH1D* _h_gen, Double_t events2);
TH1D* fill_eff(TH1D* eff, double num1, double err1, double num2, double err2, int i);
TH1D* draw_sg_bg(TH1D* s, TH1D* s_bg);
void beauty_draw(TH1D* ( &h )[2][4], int i, TString name);

void efficiency_MC(){

gStyle->SetOptFit(0);
gStyle->SetOptStat(0);

TCanvas *c1 = new TCanvas("all ","all ",1200,800);
c1->Divide(2);
TCanvas *c2 = new TCanvas("Disp&CPV ", "Disp&CPV ", 1200, 800);
c2->Divide(2);
TCanvas *c3 = new TCanvas("true","true",1200,800);
c3->Divide(2);

TFile *_file0 = TFile::Open("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/LHC16_17.root");//watch the name!
THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEtaMB;1");
TFile *_file1 = TFile::Open("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/histos.root");//watch the name!
THashList *_l1 =(THashList*)_file1->FindObjectAny("PHOSEta;1");

TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");
TH2D *_h_all_cut=(TH2D*)_l->FindObject("InvMass_All_Rejection_E300");
TH2D *_h_both_=(TH2D*)_l->FindObject("InvMass_Both_E300");
TH2D *_h_both_cut=(TH2D*)_l->FindObject("InvMass_Both_Rejection_E300");

TH2D *_h_true=(TH2D*)_l->FindObject("True_eta_Inv_mass");
TH2D *_h_true_cut=(TH2D*)_l->FindObject("True_eta_Inv_mass_Rejection");


TH1D* _h_gen = (TH1D*)_l1->FindObject("EtaGeneratedSpectra");
//cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Bins - "<<_h_gen->GetNbinsY()<<endl;
//_h_gen->Rebin(36);
//_h_gen->Rebin(2);
//cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Bins - "<<_h_gen->GetNbinsX()<<endl;

_h_gen->Sumw2();
_h_gen->GetXaxis()->SetRangeUser(0.4, 0.7);
_h_gen->SetTitle("gen");
_h_gen->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");

double pTbins[20]={2.5, 4.5, 9};// !!! to be corrected 1.5,3.5,5.5,7.5,9.5

TH1D* etaspectrum_both_cut = new TH1D("etaspectrum_both_cut ","eta spectrum both_cut ", 2, pTbins);
TH1D* etaspectrum_all_cut = new TH1D("etaspectrum_all_cut ","eta spectrum all_cut ", 2, pTbins);
TH1D* etaspectrum_both = new TH1D("etaspectrum_both ","eta spectrum both ", 2, pTbins);
TH1D* etaspectrum_all = new TH1D("etaspectrum_all ","eta spectrum all ", 2, pTbins);

TH1D* etaspectrum_true = new TH1D("etaspectrum_true","eta spectrum true", 2, pTbins);
TH1D* etaspectrum_true_cut = new TH1D("etaspectrum_true cut","eta spectrum true cut", 2, pTbins);

TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
Double_t events=(Double_t)_ev->GetBinContent(2);
 cout<< "!!! Number of events  - " << events<< "!!!"<<endl;

 TH1D *_ev1=(TH1D*)_l1->FindObject("hSelEvents");
 Double_t events1=(Double_t)_ev1->GetBinContent(2);
  cout<< "!!! Number of events  - " << events1<< "!!!"<<endl;

Int_t bin1=0;
Int_t bin2=0;

Double_t x[4];
Double_t y[4];
Double_t exl[4];
Double_t eyl[4];
Double_t exh[4];
Double_t eyh[4];

TH1D* signal_bg_all = new TH1D("signal/background bg all "," all", 2, pTbins);
TH1D* signal_bg_all_cut = new TH1D("signal/background bg all cut "," all cut", 2, pTbins);
TH1D* signal_bg_both_ = new TH1D("signal/background bg Disp&CPV "," Disp&CPV", 2, pTbins);
TH1D* signal_bg_both_cut = new TH1D("signal/background bg Disp&CPV cut "," Disp&CPV cut", 2, pTbins);
TH1D* signal_bg_true = new TH1D("signal/background bg true"," true", 2, pTbins);
TH1D* signal_bg_true_cut = new TH1D("signal/background bg true cut"," true cut", 2, pTbins);

TH1D* signal_all = new TH1D("signal/background all "," all", 2, pTbins);
TH1D* signal_all_cut = new TH1D("signal/background all cut "," all cut", 2, pTbins);
TH1D* signal_both_ = new TH1D("signal/background Disp&CPV "," Disp&CPV",2, pTbins);
TH1D* signal_both_cut = new TH1D("signal/background Disp&CPV cut "," Disp&CPV cut", 2, pTbins);
TH1D* signal_true = new TH1D("signal/background true","true", 2, pTbins);
TH1D* signal_true_cut = new TH1D("signal/background true cut"," true cut", 2, pTbins);

TH1D *_hf_all;
TH1D *_hf_all_cut;
TH1D *_hf_both_;
TH1D *_hf_both_cut;

TH1D *_hf_true;
TH1D *_hf_true_cut;

double num_all, err_all;
double num_all_cut, err_all_cut;
double num_both_, err_both_;
double num_both_cut, err_both_cut;

double num_true, err_true;
double num_true_cut, err_true_cut;

int count=0;
TString hist_name;

TH1D * h[2][4];

for(int i=0; i<2; i++)
 {

     bin1 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i]);
     bin2 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i+1]);






         hist_name.Form("%.2f < p_{T} < %.2f GeV/c",pTbins[i],pTbins[i+1]);
         cout<<"Hist_name - "<< hist_name<<endl;

         _hf_all=fill_hist(_h_all, bin1, bin2, "all", hist_name, i);
         _hf_all_cut=fill_hist(_h_all_cut, bin1, bin2, "all_cut", hist_name, i);
         _hf_both_=fill_hist(_h_both_, bin1, bin2, "Disp&CPV_", hist_name, i);
         _hf_both_cut=fill_hist(_h_both_cut, bin1, bin2, "Disp&CPV_cut", hist_name, i);
         _hf_true=fill_hist(_h_true, bin1, bin2, "true", hist_name, i);
         _hf_true_cut=fill_hist(_h_true_cut, bin1, bin2, "true_cut", hist_name, i);

         _hf_all->SetLineColor(kBlue);
         _hf_all_cut->SetLineColor(kRed);
         _hf_both_->SetLineColor(kBlue);
         _hf_both_cut->SetLineColor(kRed);

         _hf_true->SetLineColor(kCyan-3);
         _hf_true_cut->SetLineColor(kBlue+4);

         _hf_all->SetMarkerColor(kBlue);
         _hf_all_cut->SetMarkerColor(kRed);
         _hf_both_->SetMarkerColor(kBlue);
         _hf_both_cut->SetMarkerColor(kRed);
         _hf_true->SetMarkerColor(kCyan-3);
         _hf_true_cut->SetMarkerColor(kBlue+4);

         double dpt = pTbins[i+1]-pTbins[i];

         c3->cd(i+1);
         count++;//1 or 4
         cout<<"cout = "<<count<<endl;
         IntFit(c3,_hf_true, "true", signal_true, signal_bg_true, dpt, i, count, etaspectrum_true); // считает интегралы, ошибки
         IntFit(c3,_hf_true_cut, "true_cut", signal_true_cut, signal_bg_true_cut, dpt, i, count, etaspectrum_true_cut); // считает интегралы, ошибки

         drawing(c3, i+1, _hf_true, _hf_true_cut, _hf_true, _hf_true_cut, signal_true, signal_true_cut, signal_true, signal_true_cut, "true","True rejection","true","True rejection", count);//рисовка пар гистограмм

         c1->cd(i+1);
         count++;// 2 or 5
         IntFit(c1,_hf_all, "all", signal_all, signal_bg_all, dpt, i, count, etaspectrum_all);// считает интегралы, ошибки
         IntFit(c1,_hf_all_cut, "all_cut", signal_all_cut, signal_bg_all_cut, dpt, i, count, etaspectrum_all_cut); // считает интегралы, ошибки

         drawing(c1, i+1, _hf_all, _hf_all_cut, _hf_true, _hf_true_cut, signal_all, signal_all_cut, signal_true, signal_true_cut,"ALL ", "ALL REJECTION", "TRUE", "TRUE REJECTION", count);

         c2->cd(i+1);
         count++;// 3 or 6
         IntFit(c2,_hf_both_, "Disp&CPV_", signal_both_, signal_bg_both_, dpt, i, count, etaspectrum_both); // считает интегралы, ошибки
         IntFit(c2,_hf_both_cut, "Disp&CPV_cut", signal_both_cut, signal_bg_both_cut, dpt, i, count, etaspectrum_both_cut); // считает интегралы, ошибки
         drawing(c2, i+1, _hf_both_, _hf_both_cut, _hf_true, _hf_true_cut, signal_both_, signal_both_cut, signal_true, signal_true_cut, "Disp&CPV ", "Disp&CPV REJECTION", "TRUE", "TRUE REJECTION", count);

         h[i][0]=_hf_true;
         h[i][1]=_hf_all;
         h[i][2]=_hf_true_cut;
         h[i][3]=_hf_all_cut;
         beauty_draw(h, i, "All");
         h[i][1]=_hf_both_;
         h[i][3]=_hf_both_cut;
         beauty_draw(h, i, "Disp&CPV");
      }
  draw_spec(etaspectrum_true,  etaspectrum_true_cut, etaspectrum_all, etaspectrum_all_cut, etaspectrum_both, etaspectrum_both_cut, events, _h_gen);

 }

void draw_spec(TH1D * true_h, TH1D * true_cut, TH1D * all, TH1D * all_cut, TH1D * both, TH1D * both_cut, Double_t events, TH1D* _h_gen, Double_t events2){
  TCanvas *c2 = new TCanvas(true_h->GetName(),true_h->GetName(),1200,800);


  //c2->cd(1)->SetLogy();

  true_h->Scale(1./events);
  true_cut->Scale(1./events);
  all->Scale(1./events);
  all_cut->Scale(1./events);
  both->Scale(1./events);
  both_cut->Scale(1./events);
  _h_gen->Scale(1./events2);

  _h_gen->GetXaxis()->SetLimits(2., 9.);
  _h_gen->Rebin(36);
  cout<<"Number of bins X "<<_h_gen->GetNbinsX()<<endl;
  cout<<"Number of bins Y "<<_h_gen->GetNbinsY()<<endl;
  //_h_gen->GetYaxis()->Rebin(2);
  true_h->GetXaxis()->SetLimits(2., 9);
  true_cut->GetXaxis()->SetLimits(2., 9);
  all->GetXaxis()->SetLimits(2., 9);
  all_cut->GetXaxis()->SetLimits(2., 9);
  both->GetXaxis()->SetLimits(2., 9);
  both_cut->GetXaxis()->SetLimits(2., 9);
/*
  _h_gen->GetYaxis()->SetLimits(1e-12, 0.001);
  true_h->GetYaxis()->SetLimits(1e-12, 0.001);
  true_cut->GetYaxis()->SetLimits(1e-12, 0.001);
  all->GetYaxis()->SetLimits(1e-12, 0.001);
  all_cut->GetYaxis()->SetLimits(1e-12, 0.001);
  both->GetYaxis()->SetLimits(1e-12, 0.001);
  both_cut->GetYaxis()->SetLimits(1e-12, 0.001);
*/
//  cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Number of bins X - "<<true_h->GetNbinsX()<<endl;
//  cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Number of bins Y - "<<true_h->GetNbinsY()<<endl;

  TH1D* eff_true = (TH1D*)true_h->Clone("eff_true");
  TH1D* eff_true_cut = (TH1D*) true_cut->Clone("eff_true_cut ");

  eff_true->Divide(_h_gen);
  eff_true_cut->Divide(_h_gen);

  TH1D* eff_all = (TH1D*)all->Clone("eff_all");
  TH1D* eff_all_cut = (TH1D*) all_cut->Clone("eff_all_cut ");

  eff_all->Divide(_h_gen);
  eff_all_cut->Divide(_h_gen);

  TH1D* eff_both = (TH1D*)all->Clone("eff_both");
  TH1D* eff_both_cut = (TH1D*) all_cut->Clone("eff_both_cut ");

  eff_both->Divide(_h_gen);
  eff_both_cut->Divide(_h_gen);


  eff_true->GetXaxis()->SetTitle("p_{T}, GeV/c");
  eff_true_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");

  eff_true->GetXaxis()->SetRangeUser(2., 9.);
  eff_true->GetYaxis()->SetRangeUser(0., 0.12);

  eff_all->GetXaxis()->SetTitle("p_{T}, GeV/c");
  eff_all_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");

  eff_both->GetXaxis()->SetTitle("p_{T}, GeV/c");
  eff_both_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");

  eff_true->SetMarkerStyle(21);
  eff_true_cut->SetMarkerStyle(25);

  eff_all->SetLineColor(kRed);
  eff_all->SetMarkerColor(kRed);
  eff_all->SetMarkerStyle(23);

  eff_all_cut->SetLineColor(kRed);
  eff_all_cut->SetMarkerColor(kRed);
  eff_all_cut->SetMarkerStyle(32);

  eff_both->SetLineColor(kGreen);
  eff_both->SetMarkerColor(kGreen);
  eff_both->SetMarkerStyle(20);

  eff_both_cut->SetLineColor(kGreen);
  eff_both_cut->SetMarkerColor(kGreen);
  eff_both_cut->SetMarkerStyle(24);

  auto legend1 = new TLegend(0.11, 0.75,0.34,0.9);
  legend1->AddEntry(eff_true, "true", "lep");
  legend1->AddEntry(eff_true_cut, "true cut", "lep");
  legend1->AddEntry(eff_all, "all", "lep");
  legend1->AddEntry(eff_all_cut, "all cut", "lep");
  legend1->AddEntry(eff_both, "both", "lep");
  legend1->AddEntry(eff_both_cut, "both cut", "lep");

  eff_true->Draw("p e ");
  legend1->Draw();
  eff_true_cut->Draw("p e same");
  eff_both->Draw("p e same");
  eff_both_cut->Draw("p e same");
  eff_all->Draw("p e same");
  eff_all_cut->Draw("p e same");

  c2->SaveAs("eta_eff_MC.root");

  TCanvas *c3 = new TCanvas(true_cut->GetName(),true_cut->GetName(),1200,800);


  c3->cd(1)->SetLogy();

  true_h->GetXaxis()->SetTitle("p_{T}, GeV/c");
  true_h->GetXaxis()->SetRangeUser(2., 10);
  true_h->SetMaximum(0.001);
  true_h->GetYaxis()->SetTitle("dN/dp_{T}");

  true_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
  true_cut->GetXaxis()->SetRangeUser(2., 10);
  true_cut->GetYaxis()->SetTitle("dN/dp_{T}");


  all->GetXaxis()->SetTitle("p_{T}, GeV/c");
  all->GetXaxis()->SetRangeUser(2., 10);
  all->GetYaxis()->SetTitle("dN/dp_{T}");


  all_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
  all_cut->GetXaxis()->SetRangeUser(2., 10);
  all_cut->GetYaxis()->SetTitle("dN/dp_{T}");


  both->GetXaxis()->SetTitle("p_{T}, GeV/c");
  both->GetXaxis()->SetRangeUser(2., 10);
  both->GetYaxis()->SetTitle("dN/dp_{T}");

  both_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
  both_cut->GetXaxis()->SetRangeUser(2., 10);
  both_cut->GetYaxis()->SetTitle("dN/dp_{T}");

  true_h->SetMarkerStyle(21);

  //TF1 *f3 = new TF1("f3","[0]*TMath::Exp(-[1]*x)+TMath::Power([2]*(1+[3]*x*x), 13)", 0.5, 10);
  //pi0spectrum->Fit("f3", "RMS");

  true_cut->SetMarkerStyle(25);

  all->SetLineColor(kRed);
  all->SetMarkerColor(kRed);
  all->SetMarkerStyle(23);

  all_cut->SetLineColor(kRed);
  all_cut->SetMarkerColor(kRed);
  all_cut->SetMarkerStyle(32);

  both->SetLineColor(kBlack);
  both->SetMarkerColor(kBlack);
  both->SetMarkerStyle(20);

  both_cut->SetLineColor(kBlack);
  both_cut->SetMarkerColor(kBlack);
  both_cut->SetMarkerStyle(24);

  auto legend = new TLegend(0.25, 0.77,0.9,0.9);
  legend->AddEntry(true_h, "true", "lep");
  legend->AddEntry(true_cut, "true cut", "lep");
  legend->AddEntry(all, "all", "lep");
  legend->AddEntry(all_cut, "all cut", "lep");
  legend->AddEntry(both, "both", "lep");
  legend->AddEntry(both_cut, "both cut", "lep");

  true_h->Draw("p e ");
  legend->Draw();
  true_cut->Draw("p e same");
  both->Draw("p e same");
  both_cut->Draw("p e same");
  all->Draw("p e same");
  all_cut->Draw("p e same");


  c3->SaveAs("eta_spec_MC.root");
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

void IntFit(TCanvas *c, TH1D* _hf, char* name, TH1D* signal, TH1D* signal_bg, double dpt, int i, int count, TH1D * etaspectrum)
{
  //  c->cd(i+1);
    char * fname = Form("f1_%s_%d",name,i);
    char * fname2 = Form("N_%s_%d",name,i);
    TF1  *f1 = new TF1(fname,"gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    f1->SetParameters(1.*_hf->GetBinContent(56),0.55,0.013,0,0,0);//-1e4/(i+1)/(i+1),1e3/(i+1)/(i+1),1.7e4/(i+1)/(i+1));//-2e5/(i+1),3e4/(i+1),2.5*_hf->GetBinContent(51)/(i+1)) ;
    f1->SetParLimits(1, 0.52, 0.58);
    if (count==1 or count==4){
      f1->SetParameter(1, 0.55);
      f1->SetParLimits(2, 0.016, 0.025);
      f1->SetParameter(2, 0.2);
    }
    f1->SetParLimits(0, 0,2*_hf->GetBinContent(56));
    _hf->Fit(fname, "0");
//    _hf->Fit(fname, "0");
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
  cout<<fname<<endl;
  cout<<"Integral: "<< num  <<endl;
  cout<<"Integral error: "<< num_err <<endl;

//  cout<<"DPT"<<dpt<<endl;
  etaspectrum->SetBinContent(i+1, num / dpt);
  etaspectrum->SetBinError(i+1, num_err / dpt);

  double signal_v = Number->Integral(mean-2.*sigma, mean+2.*sigma)/_hf->GetXaxis()->GetBinWidth(0);
  double signal_v_error = Number->IntegralError(mean-2.*sigma, mean+2.*sigma,r->GetParams(), dataBG.GetMatrixArray())/_hf->GetXaxis()->GetBinWidth(0);
  cout<<"Signal: "<<signal_v<<endl;
  cout<<"Signal_error: "<<signal_v_error<<endl;

  signal->SetBinContent(i+1, signal_v/ dpt);//только сигнал
  signal->SetBinError(i+1, signal_v_error/ dpt);

  signal_bg->SetBinContent(i+1, num/ dpt - signal_v/ dpt);///
  signal_bg->SetBinError(i+1, sqrt((num_err/ dpt)*(num_err/ dpt) + (signal_v_error/ dpt)*(signal_v_error/ dpt)));///!!!
  //  return  std::make_tuple(num, num_err);
}

void drawing(TCanvas *c, int i, TH1D *_hf1, TH1D *_hf2,  TH1D *_hf3, TH1D *_hf4,  TH1D *s1, TH1D *s2, TH1D *s3, TH1D *s4, TString name1, TString name2, TString name3, TString name4, int count)
{
    gStyle->SetOptFit(0);

  auto legend = new TLegend(0.25, 0.77,0.9,0.9);//TO BE CORRECTED

  _hf1->SetMinimum(0.);
  _hf2->SetMinimum(0.);
  _hf3->SetMinimum(0.);
  _hf4->SetMinimum(0.);
  if(count==2){
    _hf1->SetMaximum(1350);
    _hf2->SetMaximum(1350);
    _hf3->SetMaximum(1350);
    _hf4->SetMaximum(1350);
  }
  if(count==5){
    _hf1->SetMaximum(350);
    _hf2->SetMaximum(350);
    _hf3->SetMaximum(350);
    _hf4->SetMaximum(350);
  }
  if(count==3){
    _hf1->SetMaximum(655);
    _hf2->SetMaximum(655);
    _hf3->SetMaximum(655);
    _hf4->SetMaximum(655);
  }
  if(count==6){
   _hf1->SetMaximum(165);
   _hf2->SetMaximum(165);
   _hf3->SetMaximum(165);
   _hf4->SetMaximum(165);
 }
  _hf1->Draw();
  _hf2->Draw("same");
  _hf3->Draw("same");
  _hf4->Draw("same");

  double num1 = s1->GetBinContent(i);
  double num2 = s2->GetBinContent(i);
  double num3 = s3->GetBinContent(i);
  double num4 = s4->GetBinContent(i);

  double err1 = s1->GetBinError(i);
  double err2 = s2->GetBinError(i);
  double err3 = s3->GetBinError(i);
  double err4 = s4->GetBinError(i);

  legend->AddEntry(_hf1,Form("Integral: %2.1f #pm %2.1f ",num1,err1)+name1,"lep");
  legend->AddEntry(_hf2,Form("Integral: %2.1f #pm %2.1f ",num2,err2)+name2,"lep");
  legend->AddEntry(_hf3,Form("Integral: %2.1f #pm %2.1f ",num3,err3)+name3,"lep");
  legend->AddEntry(_hf4,Form("Integral: %2.1f #pm %2.1f ",num4,err4)+name4,"lep");

  legend->Draw("same");
}

void draw_eff(double up_lim, double down_lim, TH1D*s_all, TH1D*s_both, TH1D* s_true)
{
  TCanvas *c_s_bg = new TCanvas("eff all", "Efficiency All", 1200, 800);

   s_all->SetStats(0);
   s_both->SetStats(0);
   s_all->SetMaximum(up_lim);
   s_all->SetMinimum(down_lim);
   s_all->GetXaxis()->SetRangeUser(0., 10.);
    cout<<"0"<<endl;
   s_all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   s_all->GetYaxis()->SetTitle("#epsilon");
   s_both->GetXaxis()->SetRangeUser(0., 10.);
    cout<<"01"<<endl;
   s_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   s_both->GetYaxis()->SetTitle("N_{#eta}^{Rejection}/N_{#eta}");

   s_all->SetMarkerStyle(21);
   s_both->SetMarkerStyle(20);
    s_true->SetMarkerStyle(20);
   s_all->SetMarkerSize(2);
   s_both->SetMarkerSize(2);
    s_true->SetMarkerSize(2);

   s_all->SetLineColor(kGreen+2);
   s_both->SetLineColor(kBlue+2);
    s_true->SetLineColor(kRed+2);
   s_all->SetMarkerColor(kGreen+2);
   s_both->SetMarkerColor(kBlue+2);
    s_true->SetMarkerColor(kRed+2);

   s_all->SetMinimum(0.);
   s_all->SetMaximum(1.2);
   s_all->SetTitle("");
   s_both->SetTitle("");
   s_all->Draw("p e1");
   s_both->Draw("same p e1");
    cout<<"02"<<endl;
    auto legend = new TLegend(0.25, 0.77,0.9,0.9);
    cout<<"03"<<endl;
    legend->AddEntry(s_all,"efficiency  PID: All Rejection/All","p");
    legend->AddEntry(s_both,"efficiency PID: Disp&CPV Rejection/Disp&CPV","p");
    legend->Draw("same");
    cout<<"04"<<endl;
    c_s_bg->SaveAs("MC/both_efficiency_mc.root");
    cout<<"06"<<endl;
}


void draw_sg(double up_lim, double down_lim, TH1D*s_all, TH1D*s_all_cut, TH1D* s_both, TH1D* s_both_cut)
{
    TCanvas *c_s_bg = new TCanvas("s/bg", "s/bg", 1200, 800);
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

    cout<<"1"<<endl;
   s_all->SetStats(0);
   s_all_cut->SetStats(0);
    s_both->SetStats(0);
   s_all->SetMaximum(0.8);
   s_all->SetMinimum(down_lim);
   s_all->GetXaxis()->SetRangeUser(0., 10.);
    cout<<"2"<<endl;
  s_all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  s_all->GetYaxis()->SetTitle("(S/Bg)");

   s_all->SetMarkerStyle(21);
   s_all_cut->SetMarkerStyle(25);
    s_both->SetMarkerStyle(20);
    s_both_cut->SetMarkerStyle(24);
   s_all->SetMarkerSize(2);
   s_all_cut->SetMarkerSize(2);
    s_both->SetMarkerSize(2);
    s_both_cut->SetMarkerSize(2);

   s_all->SetLineColor(kGreen+2);
   s_all_cut->SetLineColor(kBlack+2);
    s_both->SetLineColor(kBlue+2);
    s_both_cut->SetLineColor(kRed+2);
   s_all->SetMarkerColor(kGreen+2);
   s_all_cut->SetMarkerColor(kBlack+2);
    s_both->SetMarkerColor(kBlue+2);
    s_both_cut->SetMarkerColor(kRed+2);

   s_all->SetTitle("");
   s_all->Draw("p e1");
   s_all_cut->Draw("p e1 same");
    s_both->Draw("p e1 same");
    s_both_cut->Draw("p e1 same");

    auto legend = new TLegend(0.2, 0.6,0.45,0.78);
    legend->AddEntry(s_all,"PID: All","p");
    legend->AddEntry(s_all_cut,"PID: All Rejection","p");
    legend->AddEntry(s_both,"PID: Disp&CPV","p");
    legend->AddEntry(s_both_cut,"PID: Disp&CPV Rejection","p");
    legend->SetLineColor(0);
    legend->Draw("same");

    TLatex * t1 =new TLatex(-0.13, 0.73,"ALICE simulations") ;
    t1->SetTextSize(0.04) ;
    t1->Draw() ;
    TLatex * t2 =new TLatex(-0.14, 0.69,"work in progress") ;
    t2->SetTextSize(0.04) ;
    t2->Draw() ;
    TLatex * t3 =new TLatex(0.35, 0.65,"10.11.2020") ;
    t3->SetTextSize(0.04) ;
    t3->Draw() ;

    TLatex * t4 =new TLatex(9.8, 0.73,"pp, #sqrt{s}=13 TeV") ;
    t4->SetTextSize(0.04) ;
    t4->Draw() ;

    c_s_bg->SaveAs("MC/sg_bg.pdf");

    TCanvas *c_dr = new TCanvas("doubleratio", "doubleratio", 1200, 800);

    c_dr->cd(1);
    TH1D* dr_all = (TH1D*)s_all_cut->Clone("dr_all");
    TH1D* dr_both = (TH1D*) s_both_cut->Clone("dr_both");
    dr_all->Divide(s_all);
    dr_both->Divide(s_both);
    dr_all->SetMinimum(0);
    dr_all->SetMaximum(2.7);
    dr_all->SetTitle("");

    dr_all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    dr_all->GetYaxis()->SetTitle("(S/Bg)_{all} / (S/Bg)_{rej.#pi^{0}}");
    dr_all->GetXaxis()->SetRangeUser(0., 10.);

    dr_all->Draw("p e1");
    dr_both->Draw("p e1 same");

    TLatex * t_1 =new TLatex(0.045, 2.5,"ALICE simulations") ;
    t_1->SetTextSize(0.04) ;
    t_1->Draw() ;
    TLatex * t_2 =new TLatex(0.31, 2.27,"work in progress") ;
    t_2->SetTextSize(0.04) ;
    t_2->Draw() ;
    TLatex * t_3 =new TLatex(0.78, 2.,"10.11.2020") ;
    t_3->SetTextSize(0.04) ;
    t_3->Draw() ;

    TLatex * t_4 =new TLatex(5.3, 2.5,"pp, #sqrt{s}=13 TeV") ;
    t_4->SetTextSize(0.04) ;
    t_4->Draw() ;

    auto legend2 = new TLegend(0.71, 0.72,0.89,0.89);
    legend2->AddEntry(dr_all,"PID: All","p");
    legend2->AddEntry(dr_both,"PID: Disp&CPV","p");
    legend2->SetLineColor(0);
    legend2->Draw("same");

    c_dr->SaveAs("MC/doubleratio.pdf");
}

void beauty_draw(TH1D* ( &h )[2][4], int i, TString name)
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

  if(i==0 and name=="All"){
  box->SetMaximum(1999.) ;
  }
  if(i==1 and name=="All"){
    box->SetMaximum(500.) ;
  }
    if(i==0 and name=="Disp&CPV"){
        box->SetMaximum(900.) ;
    }
    if(i==1 and name=="Disp&CPV")
    {
        box->SetMaximum(300.) ;
    }


  h[i][3]->SetMarkerStyle(21) ;
  h[i][3]->SetMarkerColor(kAzure+9) ;
  h[i][3]->SetLineColor(kAzure+9) ;
  h[i][3]->SetLineWidth(2) ;

  h[i][0]->SetMarkerStyle(24) ;
  h[i][0]->SetMarkerColor(kPink) ;
  h[i][0]->SetLineColor(kPink) ;
  h[i][0]->SetLineWidth(2) ;

  h[i][2]->SetMarkerStyle(25) ;
  h[i][2]->SetMarkerColor(kGreen+3) ;
  h[i][2]->SetLineColor(kGreen+3) ;
  h[i][2]->SetLineWidth(2) ;

  h[i][1]->Draw("same") ;
  h[i][0]->Draw("same") ;
  h[i][3]->Draw("same") ;
  h[i][2]->Draw("same") ;

  TLegend * leg = new TLegend(0.56,0.60,0.96,0.79) ;
  leg->AddEntry(h[i][1],name+" pairs","p") ;
  leg->AddEntry(h[i][3]," rejected #pi^{0} ","p") ;
  leg->AddEntry(h[i][0], name+" pairs, PID: True","p") ;
  leg->AddEntry(h[i][2],"rejected #pi^{0}, PID: True","p") ;
  leg->SetLineColor(0) ;
  leg->Draw() ;
  if(i==0 and name=="All"){
  TLatex * t1 =new TLatex(0.415,1900.,"ALICE simulations") ;
  t1->SetTextSize(0.04) ;
  t1->Draw() ;
  TLatex * t2 =new TLatex(0.415,1790.,"work in progress") ;
  t2->SetTextSize(0.04) ;
  t2->Draw() ;
  TLatex * t3 =new TLatex(0.415,1680.,"10.11.2020") ;
  t3->SetTextSize(0.04) ;
  t3->Draw() ;

  TLatex * t4 =new TLatex(0.596,1900.,"pp, #sqrt{s}=13 TeV") ;
  t4->SetTextSize(0.04) ;
  t4->Draw() ;

  TLatex * t5 =new TLatex(0.58,1790.,"2.5<p_{T}<4.5 GeV/c") ;
  t5->SetTextSize(0.04) ;
  t5->Draw() ;
  }

  if(i==1 and name=="All"){
  TLatex * t1 =new TLatex(0.415,475.,"ALICE simulations") ;
  t1->SetTextSize(0.04) ;
  t1->Draw() ;
  TLatex * t2 =new TLatex(0.415,448.,"work in progress") ;
  t2->SetTextSize(0.04) ;
  t2->Draw() ;
  TLatex * t3 =new TLatex(0.415, 421.,"10.11.2020") ;
  t3->SetTextSize(0.04) ;
  t3->Draw() ;

  TLatex * t4 =new TLatex(0.596,475.,"pp, #sqrt{s}=13 TeV") ;
  t4->SetTextSize(0.04) ;
  t4->Draw() ;

    TLatex * t5 =new TLatex(0.58, 448.,"4.5<p_{T}<9 GeV/c") ;
    t5->SetTextSize(0.04) ;
    t5->Draw() ;
  }

  if(i==0 and name=="Disp&CPV"){
  TLatex * t1 =new TLatex(0.415,855.,"ALICE simulations") ;
  t1->SetTextSize(0.04) ;
  t1->Draw() ;
  TLatex * t2 =new TLatex(0.415,805.,"work in progress") ;
  t2->SetTextSize(0.04) ;
  t2->Draw() ;
  TLatex * t3 =new TLatex(0.415, 758.,"10.11.2020") ;
  t3->SetTextSize(0.04) ;
  t3->Draw() ;

  TLatex * t4 =new TLatex(0.596, 855.,"pp, #sqrt{s}=13 TeV") ;
  t4->SetTextSize(0.04) ;
  t4->Draw() ;

  TLatex * t5 =new TLatex(0.58,805.,"2.5<p_{T}<4.5 GeV/c") ;
  t5->SetTextSize(0.04) ;
  t5->Draw() ;
  }

  if(i==1 and name=="Disp&CPV"){
  TLatex * t1 =new TLatex(0.415, 280.,"ALICE simulations") ;
  t1->SetTextSize(0.04) ;
  t1->Draw() ;
  TLatex * t2 =new TLatex(0.415, 265.,"work in progress") ;
  t2->SetTextSize(0.04) ;
  t2->Draw() ;
  TLatex * t3 =new TLatex(0.415, 252.,"10.11.2020") ;
  t3->SetTextSize(0.04) ;
  t3->Draw() ;

  TLatex * t4 =new TLatex(0.596, 280.,"pp, #sqrt{s}=13 TeV") ;
  t4->SetTextSize(0.04) ;
  t4->Draw() ;

  TLatex * t5 =new TLatex(0.58,265.,"4.5<p_{T}<9 GeV/c") ;
  t5->SetTextSize(0.04) ;
  t5->Draw() ;
  }


  c1DP->SaveAs("MC/"+name+"/"+Form("_%d",i)+".pdf");
}
