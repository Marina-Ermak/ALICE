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
void draw_eff(double up_lim, double down_lim, TH1D*s_all, TH1D*s_both);
TH1D* fill_eff(TH1D* eff, double num1, double err1, double num2, double err2, int i);
TH1D* draw_sg_bg(TH1D* s, TH1D* s_bg);
void beauty_draw(TH1D* ( &h )[16][2], int i, TString name, TString hist_name);
void draw_spec(TH1D* all, TH1D* all_cut, TH1D* both, TH1D* both_cut, Double_t events);

void mod_eta_LHC_phi7_2016_2017() {
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("all ","all ",1200,800);

   c1->Divide(2);
   TCanvas *c2 = new TCanvas("Disp&CPV ", "Disp&CPV ", 1200, 800);
   c2->Divide(2);

   TFile *_file0 = TFile::Open("histos_MC_09102020.root");//watch the name!
   THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEta;1");

   TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");
   TH2D *_h_all_cut=(TH2D*)_l->FindObject("InvMass_All_Rejection_E300");
   TH2D *_h_both_=(TH2D*)_l->FindObject("InvMass_Both_E300");
   TH2D *_h_both_cut=(TH2D*)_l->FindObject("InvMass_Both_Rejection_E300");

   Int_t bins=11;
   double pTbins[11]={};
   pTbins[0]=0.5;// {0.5, 0.7, 0.9, 1.2, 1.5, 2, 2.5, 3.2, 4., 5., 7., 10.}
    for (int i=1; i<11; i++)
        {
             pTbins[i]=pTbins[i-1]+4;
            cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
        }

   TH1D* etaspectrum_both_cut = new TH1D("etaspectrum_both_cut ","eta spectrum both_cut ", bins, pTbins);
   TH1D* etaspectrum_all_cut = new TH1D("etaspectrum_all_cut ","eta spectrum all_cut ", bins, pTbins);
   TH1D* etaspectrum_both = new TH1D("etaspectrum_both ","eta spectrum both ", bins, pTbins);
   TH1D* etaspectrum_all = new TH1D("etaspectrum_all ","eta spectrum all ", bins, pTbins);

   TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
   Double_t events=(Double_t)_ev->GetBinContent(2);
    cout<< "!!! Number of events  - " << events<< "!!!"<<endl;

   Int_t bin1=0;
   Int_t bin2=0;

  TH1D* signal_bg_all = new TH1D("signal/background bg all "," all", bins, pTbins);
  TH1D* signal_bg_all_cut = new TH1D("signal/background bg all cut "," all cut", bins, pTbins);
  TH1D* signal_bg_both_ = new TH1D("signal/background bg Disp&CPV "," Disp&CPV", bins, pTbins);
  TH1D* signal_bg_both_cut = new TH1D("signal/background bg Disp&CPV cut "," Disp&CPV cut", bins, pTbins);

  TH1D* signal_all = new TH1D("signal/background all "," all", bins, pTbins);
  TH1D* signal_all_cut = new TH1D("signal/background all cut "," all cut", bins, pTbins);
  TH1D* signal_both_ = new TH1D("signal/background Disp&CPV "," Disp&CPV", bins, pTbins);
  TH1D* signal_both_cut = new TH1D("signal/background Disp&CPV cut "," Disp&CPV cut", bins, pTbins);


  TH1D *_hf_all;
  TH1D *_hf_all_cut;
  TH1D *_hf_both_;
  TH1D *_hf_both_cut;

   int count=0;
   TString hist_name;

TH1D * h[16][2];

  for(int i=0; i<3; i++)
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
            IntFit(c1,_hf_all, "all", signal_all, signal_bg_all, dpt, i, count, etaspectrum_all);// считает интегралы, ошибки
            IntFit(c1,_hf_all_cut, "all_cut", signal_all_cut, signal_bg_all_cut, dpt, i, count, etaspectrum_all_cut); // считает интегралы, ошибки

            c2->cd(i+1);
            count++;// 3 or 6
            IntFit(c2,_hf_both_, "Disp&CPV_", signal_both_, signal_bg_both_, dpt, i, count, etaspectrum_both); // считает интегралы, ошибки
            IntFit(c2,_hf_both_cut, "Disp&CPV_cut", signal_both_cut, signal_bg_both_cut, dpt, i, count,  etaspectrum_both_cut); // считает интегралы, ошибки

            h[i][0]=_hf_all;
            h[i][1]=_hf_all_cut;
            beauty_draw(h, i, "All", hist_name);//!
            h[i][0]=_hf_both_;
            h[i][1]=_hf_both_cut;
            beauty_draw(h, i, "Disp&CPV", hist_name);

    }

    TH1D* efficiency_both=(TH1D*)signal_both_cut->Clone("efficiency");

    efficiency_both->Divide(efficiency_both, signal_both_,1,1,"B");

    TH1D* efficiency_all=(TH1D*)signal_all_cut->Clone("efficiency_all");

    efficiency_all->Divide(efficiency_all,signal_all,1,1,"B");

    draw_eff(1., 0., efficiency_all, efficiency_both);

    TH1D * sbg_all = (TH1D*)signal_all->Clone("sbg_all");
    TH1D * sbg_all_cut = (TH1D*)signal_all_cut->Clone("sbg_all_cut");
    TH1D * sbg_both = (TH1D*)signal_both_->Clone("sbg_both");
    TH1D * sbg_both_cut = (TH1D*)signal_both_cut->Clone("sbg_both_cut");

    sbg_all->Divide(sbg_all,signal_bg_all,1,1,"B");
    sbg_all_cut->Divide(sbg_all_cut,signal_bg_all_cut,1,1,"B");
    sbg_both->Divide(sbg_both,signal_bg_both_,1,1,"B");
    sbg_both_cut->Divide(sbg_both_cut,signal_bg_both_cut,1,1,"B");

    draw_sg(0.25, 0., sbg_all, sbg_all_cut, sbg_both, sbg_both_cut);//!!рисует отношения сигнала к фону и сохраняет в пдф
    draw_spec(etaspectrum_all, etaspectrum_all_cut, etaspectrum_both, etaspectrum_both_cut, events);
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
  signal->SetBinError(i+1, signal_v_error/ dpt);

  signal_bg->SetBinContent(i+1, num/ dpt - signal_v/ dpt);///
  signal_bg->SetBinError(i+1, sqrt((num_err/ dpt)*(num_err/ dpt) + (signal_v_error/ dpt)*(signal_v_error/ dpt)));///!!!

  etaspectrum->SetBinContent(i+1, signal_v/ dpt);
  etaspectrum->SetBinError(i+1, signal_v_error/ dpt);
  //  return  std::make_tuple(num, num_err);
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

    TH1D * box = new TH1D("box","",150,0.400,0.690) ;

    box->SetMinimum(0.);
    box->SetMaximum(2.5);
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
   s_all->SetMaximum(2.5);
   s_all->SetMinimum(down_lim);
   s_all->GetXaxis()->SetRangeUser(0., 20.);
    cout<<"2"<<endl;
  s_all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  s_all->GetYaxis()->SetTitle("(S/Bg)");

   s_all->SetMarkerStyle(21);
   s_all_cut->SetMarkerStyle(21);
    s_both->SetMarkerStyle(20);
    s_both_cut->SetMarkerStyle(20);
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

    auto legend = new TLegend(0.19, 0.77,0.44,0.95);
    legend->AddEntry(s_all,"PID: All","p");
    legend->AddEntry(s_all_cut,"PID: All Rejection","p");
    legend->AddEntry(s_both,"PID: Disp&CPV","p");
    legend->AddEntry(s_both_cut,"PID: Disp&CPV Rejection","p");
    legend->SetLineColor(0);
    legend->Draw("same");

    TLatex * t__1 =new TLatex(0.93, 1.5,"ALICE LHC PHI7") ;
    t__1->SetTextSize(0.04) ;
    t__1->Draw() ;
    TLatex * t__2 =new TLatex(0.93, 1.35,"work in progress") ;
    t__2->SetTextSize(0.04) ;
    t__2->Draw() ;
    TLatex * t__3 =new TLatex(0.93, 1.2,"10.11.2020") ;
    t__3->SetTextSize(0.04) ;
    t__3->Draw() ;

    TLatex * t__4 =new TLatex(0.93, 1.66,"pp, #sqrt{s}=13 TeV") ;
    t__4->SetTextSize(0.04) ;
    t__4->Draw("same") ;

    c_s_bg->SaveAs("LHC_16_17/PHI7/sg_bg_PHI7.pdf");

    TCanvas *c_dr = new TCanvas("doubleratio", "doubleratio", 1200, 800);

    c_dr->cd(1);
    TH1D* dr_all = (TH1D*)s_all_cut->Clone("dr_all");
    TH1D* dr_both = (TH1D*) s_both_cut->Clone("dr_both");
    dr_all->Divide(s_all);
    dr_both->Divide(s_both);
    dr_all->SetMinimum(0);
    dr_all->SetMaximum(4.);
    dr_all->SetTitle("");

    dr_all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    dr_all->GetYaxis()->SetTitle("(S/Bg)_{rej.#pi^{0}} / (S/Bg)_{all}");
    dr_all->GetXaxis()->SetRangeUser(0., 20.);

    dr_all->Draw("p e1");
    dr_both->Draw("p e1 same");
    TLatex * t_1 =new TLatex(1.56, 2.4,"ALICE LHC PHI7 Trigger") ;
    t_1->SetTextSize(0.04) ;
    t_1->Draw() ;
    TLatex * t_2 =new TLatex(1.56, 2.2,"work in progress") ;
    t_2->SetTextSize(0.04) ;
    t_2->Draw() ;
    TLatex * t_3 =new TLatex(1.56, 2.0,"10.11.2020") ;
    t_3->SetTextSize(0.04) ;
    t_3->Draw() ;

    TLatex * t_4 =new TLatex(1.56, 2.6,"pp, #sqrt{s}=13 TeV") ;
    t_4->SetTextSize(0.04) ;
    t_4->Draw() ;

    auto legend2 = new TLegend(0.13, 0.71,0.32,0.88);
    legend2->AddEntry(dr_all,"PID: All","p");
    legend2->AddEntry(dr_both,"PID: Disp&CPV","p");
    legend2->SetLineColor(0);
    legend2->Draw("same");
    c_dr->SaveAs("LHC_16_17/PHI7/doubleratio.pdf");
}

void beauty_draw(TH1D* ( &h )[16][2], int i, TString name, TString hist_name)
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
  box->SetMaximum(145000);
  if(name=="All")
{
  if(i==0 or i==1 or i==2 or i==3){
    box->SetMaximum(145000);
  }
  else if (i==4 and name=="All")
  {
    box->SetMaximum(100000.) ;
  }

  else if (i==5 or i==6 and name=="All")
  {
    box->SetMaximum(40000.) ;
  }

  else if (i==7 or i==8 and name=="All")
  {
    box->SetMaximum(20000.) ;
  }

  else if (i==9 or i==10 and name=="All")
  {
    box->SetMaximum(8000.) ;
  }

  else if (i==11 or i==12 or i==13 and name=="All")
  {
    box->SetMaximum(5000.) ;
  }

  else if(i==14 and name=="All")
  {
    box->SetMaximum(2000.) ;
  }

  else if(name=="All")
  {
    box->SetMaximum(600.) ;
  }
}
else{
  if(i==0 or i==1 or i==2 or i==3){//!! надо смотреть на картинки
    box->SetMaximum(70000);
  }
  else if (i==4)
  {
    box->SetMaximum(40000.) ;
  }

  else if (i==5 or i==6)
  {
    box->SetMaximum(30000.) ;
  }

  else if (i==7 or i==8)
  {
    box->SetMaximum(10000.) ;
  }

  else if (i==9 or i==10)
  {
    box->SetMaximum(5000.) ;
  }

  else if (i==11 or i==12 or i==13)
  {
    box->SetMaximum(2000.) ;
  }

  else if(i==14)
  {
    box->SetMaximum(600.) ;
  }

  else if(name=="Disp&CPV")
  {
    box->SetMaximum(600.) ;
  }
}
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
  c1DP->SaveAs("LHC_16_17/PHI7/"+name+"/"+Form("_%d",i)+".pdf");
}

void draw_eff(double up_lim, double down_lim, TH1D*s_all, TH1D*s_both)
{
  TCanvas *c_s_bg = new TCanvas("eff all", "Efficiency All", 1200, 800);

   s_all->SetStats(0);
   s_both->SetStats(0);
   s_all->SetMaximum(up_lim);
   s_all->SetMinimum(down_lim);
   s_all->GetXaxis()->SetRangeUser(0., 20.);
    cout<<"0"<<endl;
   s_all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   s_all->GetYaxis()->SetTitle("#epsilon");
   s_both->GetXaxis()->SetRangeUser(0., 20.);
    cout<<"01"<<endl;
   s_both->GetXaxis()->SetTitle("p_{T} (GeV/c)");
   s_both->GetYaxis()->SetTitle("N_{#eta}^{Rejection}/N_{#eta}");

   s_all->SetMarkerStyle(21);
   s_both->SetMarkerStyle(20);
   s_all->SetMarkerSize(2);
   s_both->SetMarkerSize(2);

   s_all->SetLineColor(kGreen+2);
   s_both->SetLineColor(kBlue+2);
   s_all->SetMarkerColor(kGreen+2);
   s_both->SetMarkerColor(kBlue+2);

   s_all->SetMinimum(0.);
   s_all->SetMaximum(1.2);
   s_all->SetTitle("");
   s_both->SetTitle("");
   s_all->Draw("p e1");
   s_both->Draw("same p e1");
    auto legend = new TLegend(0.23, 0.18,0.89,0.32);
    legend->AddEntry(s_all,"efficiency  PID: All Rejection/All","p");
    legend->AddEntry(s_both,"efficiency PID: Disp&CPV Rejection/Disp&CPV","p");
    legend->Draw("same");
    c_s_bg->SaveAs("LHC_16_17/PHI7/eficiency_PHI7.pdf");

}

void draw_spec(TH1D* all, TH1D* all_cut, TH1D* both, TH1D* both_cut, Double_t events)
{

  TCanvas *c4 = new TCanvas("c4_fit2", "The Fit Canvas", 1200, 800);

  gPad->SetLogy();
  both_cut->Scale(1./events);
  both_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
  both_cut->GetYaxis()->SetTitle("dN/dp_{T}");
  both_cut->SetMarkerStyle(21);

   both->Scale(1./events);
   both->GetXaxis()->SetTitle("p_{T}, GeV/c");
   both->GetYaxis()->SetTitle("dN/dp_{T}");
   both->SetMarkerStyle(31);
   both->SetLineColor(kGreen);

   all_cut->Scale(1./events);
   all_cut->GetXaxis()->SetTitle("p_{T}, GeV/c");
   all_cut->GetYaxis()->SetTitle("dN/dp_{T}");
   all_cut->SetMarkerStyle(41);
   all_cut->SetLineColor(kOrange);

   all->Scale(1./events);
   all->GetXaxis()->SetTitle("p_{T}, GeV/c");
   all->GetYaxis()->SetTitle("dN/dp_{T}");
   all->SetMarkerStyle(51);
   all->SetLineColor(kRed);


   auto legend_sp = new TLegend(0.1,0.1,0.38,0.3);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
     legend_sp->AddEntry(both_cut,"spectrum Disp&CPV Rejection","lep");
     legend_sp->AddEntry(both,"spectrum Disp&CPV","lep");
     legend_sp->AddEntry(all_cut,"spectrum All Rejection","lep");
     legend_sp->AddEntry(all,"spectrum All","lep");
 // TF1 *f6 = new TF1("f6","[0]*TMath::Exp(-[1]*x)", 6, 16);
 // etaspectrum->Fit("f6", "RMS");
  both->Draw("p e");
  both_cut->Draw("p e same");
  all_cut->Draw("p e same");
  all->Draw("p e same");
  legend_sp->Draw("p e same");

  c4->SaveAs("spec.root");
}
