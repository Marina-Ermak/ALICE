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

using namespace std;

Double_t CB1(Double_t * x, Double_t * par){
    Double_t kMean=0.134;
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

void IntFit(TCanvas *c, TH1D* _hf, TH1D* _hf_cb, char* name, TH1D* signal, TH1D* signal_bg, TH1D* signal_cb, TH1D* signal_bg_cb, double dpt, int i, TH1D * spectrum, TH1D * spectrum_cb, TH1D* position, TH1D* positionCB, TH1D* width, TH1D* widthCB);
void draw_spec(TH1D * pi0spectrum, TH1D * pi0spectrumCB, TH1D* pi0position, TH1D* pi0positionCB, TH1D* pi0width, TH1D* pi0widthCB, TString name, Double_t events);
TH1D* fill_hist(TH2D* _h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i);
void beauty_draw(TH1D* ( &h )[32][2], int i, TString name, TString hist_name);


void mod_pi0_LHC() {

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("all ","all ",1200,800);
  TCanvas *c2 = new TCanvas("Disp&CPV ","Disp&CPV",1200,800);

  c1->Divide(2);
  c2->Divide(2);

  TFile *_file0 = TFile::Open("LHC16_17.root");//watch the name!
  THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEtaPHI7;1");
  TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");
  TH2D *_h_all_CB1=(TH2D*)_l->FindObject("InvMass_All_E300");
  TH2D *_h_both=(TH2D*)_l->FindObject("InvMass_Both_E300");
  TH2D *_h_both_CB1=(TH2D*)_l->FindObject("InvMass_Both_E300");


  TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
  Double_t events=(Double_t)_ev->GetBinContent(2);
  cout<< "!!! Number of events - " << events<< "!!!"<<endl;

  Int_t bin1=0;
  Int_t bin2=0;
  double pTbins[33];// {0.5, 0.7, 0.9, 1.2, 1.5, 2, 2.5, 3.2, 4., 5., 7., 10.}
   for (int i=0; i<=32; i++)
       {
           pTbins[i]=3.5+i*1.;
           cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
       }

  TString hist_name;
  Int_t bins = 16;

  TH1D* signal_bg_all = new TH1D("signal/background bg all "," all", bins, pTbins);
  TH1D* signal_bg_both = new TH1D("signal/background bg Disp&CPV "," Disp&CPV", bins, pTbins);
  TH1D* signal_all = new TH1D("signal/background all "," all", bins, pTbins);
  TH1D* signal_both = new TH1D("signal/background all "," all", bins, pTbins);
  TH1D* signal_bg_all_cb = new TH1D("signal/background bg all cb"," all cb", bins, pTbins);
  TH1D* signal_bg_both_cb = new TH1D("signal/background bg Disp&CPV cb"," Disp&CPV cb", bins, pTbins);
  TH1D* signal_all_cb = new TH1D("signal/background all cb"," all cb", bins, pTbins);
  TH1D* signal_both_cb = new TH1D("signal/background all cb"," Disp&CPV cb", bins, pTbins);

  TH1D* pi0_spectrum_all = new TH1D("pi0_spectrum_all ","pi0 spectrum all ", bins, pTbins);
  TH1D* pi0_spectrum_both = new TH1D("pi0_spectrum_both ","pi0 spectrum both ", bins, pTbins);
  TH1D* pi0_spectrum_all_cb = new TH1D("pi0_spectrum_all cb","pi0 spectrum all cb", bins, pTbins);
  TH1D* pi0_spectrum_both_cb = new TH1D("pi0_spectrum_both cb ","pi0 spectrum both cb", bins, pTbins);

  TH1D* position = new TH1D("pi0position","pi0 position ", bins, pTbins);//11
  TH1D* width = new TH1D("pi0width", "pi0 width",bins, pTbins);//watch the number of bins!
  TH1D* position_cb = new TH1D("pi0position_cb","pi0 position_cb", bins, pTbins);//11
  TH1D* width_cb = new TH1D("pi0width_cb", "pi0 width cb", bins, pTbins);//watch the number of bins!

  TH1D* position_both = new TH1D("pi0position Disp&CPV","pi0 position Disp&CPV", bins, pTbins);//11
  TH1D* width_both = new TH1D("pi0width Disp&CPV", "pi0 width Disp&CPV", bins, pTbins);//watch the number of bins!
  TH1D* position_both_cb = new TH1D("pi0position_cb Disp&CPV","pi0 position_cb Disp&CPV", bins, pTbins);//11
  TH1D* width_both_cb = new TH1D("pi0width_cb Disp&CPV", "pi0 width cb Disp&CPV", bins, pTbins);//watch the number of bins!

  TH1D *_hf_all;
  TH1D *_hf_both;
  TH1D *_hf_all_cb;
  TH1D *_hf_both_cb;

TH1D * h[32][2];

  for(int i=1; i<bins; i++)
    {
    //  cout<<"work"<<i<<endl;
      bin1 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i]);
      bin2 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i+1]);
cout<<"!!!!!!!!!!!!!!!!!!"<<bin1<<" "<<bin2<<endl;

          hist_name.Form("%.2f < p_{T} < %.2f GeV/c",pTbins[i],pTbins[i+1]);
          cout<<"Hist_name - "<< hist_name<<endl;

          _hf_all=fill_hist(_h_all, bin1, bin2, "all", hist_name, i);
          _hf_both=fill_hist(_h_both, bin1, bin2, "both", hist_name, i);
          _hf_all_cb=fill_hist(_h_all, bin1, bin2, "all_cb", hist_name, i);
          _hf_both_cb=fill_hist(_h_both, bin1, bin2, "both_cb", hist_name, i);

          _hf_all->SetLineColor(kBlue);
          _hf_both->SetLineColor(kRed);
          _hf_all_cb->SetLineColor(kBlack);
          _hf_both_cb->SetLineColor(kGreen);

          _hf_all->SetMarkerColor(kBlue);
          _hf_both->SetMarkerColor(kRed);
          _hf_all_cb->SetMarkerColor(kBlack);
          _hf_both_cb->SetMarkerColor(kGreen);

          double dpt = pTbins[i+1]-pTbins[i];

          c1->cd(i+1);
          IntFit(c1,_hf_all, _hf_all_cb, "all", signal_all, signal_bg_all, signal_all_cb, signal_bg_all_cb, dpt, i, pi0_spectrum_all, pi0_spectrum_all_cb, position, position_cb, width, width_cb);
          c2->cd(i+1);
          IntFit(c2,_hf_both, _hf_both_cb, "Disp&CPV", signal_both, signal_bg_both, signal_both_cb, signal_bg_both_cb, dpt, i, pi0_spectrum_both, pi0_spectrum_both_cb, position_both, position_both_cb, width_both, width_both_cb);
          h[i][0]=_hf_all;
          h[i][1]=_hf_all_cb;
          beauty_draw(h, i, "All", hist_name);//!
          h[i][0]=_hf_both;
          h[i][1]=_hf_both_cb;
          beauty_draw(h, i, "Disp&CPV", hist_name);
}
draw_spec(pi0_spectrum_all, pi0_spectrum_all_cb, position, position_cb, width, width_cb, "all", events);
draw_spec(pi0_spectrum_both, pi0_spectrum_both_cb, position_both, position_both_cb, width_both, width_both_cb, "disp&cpv", events);
}

TH1D* fill_hist(TH2D* _h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i)
{
  TH1D*_hf=(TH1D*)_h->ProjectionX(name+Form("_%d",i),bin1,bin2);
  _hf->Rebin(10);
  _hf->Sumw2();
  _hf->SetName(name+Form("%d",i));//all SetName-s are different
  _hf->GetXaxis()->SetRangeUser(0.08, 0.25);
  _hf->SetTitle("");
  _hf->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
  return _hf;
}

void IntFit(TCanvas *c, TH1D* _hf, TH1D* _hf_cb, char* name, TH1D* signal, TH1D* signal_bg, TH1D* signal_cb, TH1D* signal_bg_cb, double dpt, int i, TH1D * spectrum, TH1D * spectrum_cb, TH1D* position, TH1D* positionCB, TH1D* width, TH1D* widthCB)
{
    char * fname = Form("f1_%s_%d",name,i);
    char * fname2 = Form("N_%s_%d",name,i);
    TF1  *f1 = new TF1(fname,"gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.1, 0.2);
        f1->SetLineStyle(2);

        f1->SetParLimits(1, 0.13, 0.15);
        f1->SetParLimits(2, 1e-03, 20e-03);
  f1->SetParameters(1.*_hf->GetBinContent(66),0.135,0.006,1.,3.,2.*_hf->GetBinContent(51)) ;

    _hf->Fit(fname, "r0");

  TFitResultPtr r = _hf->Fit(fname, "RMS+");//check if f1 works
  Double_t p[3],ep[3];
  Double_t p2[8],ep2[8];
  TMatrixDSym covrfit = r->GetCovarianceMatrix();
  TMatrixTSym<double> dataBG;
  covrfit.GetSub(0,2,dataBG);///
  TF1 *Number= new TF1(fname2,"gaus", 0.1, 0.2);

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

  position->SetBinContent(i+1, p[1]);
  position->SetBinError(i+1, ep[1]);

  width->SetBinContent(i+1, abs(p[2]));
  width->SetBinError(i+1, abs(ep[2]));

  TF1  *f_CB1 = new TF1("f_CB1",CB1, 0.08, 0.2,8);
    f_CB1->SetParameters(f1->GetParameter(0),f1->GetParameter(1),f1->GetParameter(2),1.,3.,2.*_hf_cb->GetBinContent(51),-0.001,0.0001) ;


  f_CB1->SetLineColor(kGreen);

  TF1 *Number_CB1 = new TF1("number_pi0_both_CB1", CB1S, 0.1, 0.25, 5);
  TFitResultPtr rCB = _hf_cb->Fit("f_CB1", "RMS+");
  TMatrixDSym covrfit1 = rCB->GetCovarianceMatrix();
  TMatrixTSym<double> dataBG1;
  covrfit1.GetSub(0,4,dataBG1);

  for(int j=0; j<8; j++)
  {
   p2[j]=f_CB1->GetParameter(j);
   ep2[j]=f_CB1->GetParError(j);
   Number_CB1->SetParameter(j,p2[j]);
   Number_CB1->SetParError(j,ep2[j]);
  }

  double num_cb =  Number_CB1->Integral(mean-2.*sigma, mean+2.*sigma)/_hf_cb->GetXaxis()->GetBinWidth(0);
  double err_cb =  Number_CB1->IntegralError(mean-2.*sigma, mean+2.*sigma,rCB->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/_hf_cb->GetXaxis()->GetBinWidth(0);

  cout<<"Integral both CB - "<< num_cb  <<endl;
  cout<<"Integral both CB error - "<< err_cb <<endl;

  double signal_v = Number->Integral(mean-2.*sigma, mean+2.*sigma)/_hf->GetXaxis()->GetBinWidth(0);
  double signal_v_error = Number->IntegralError(mean-2.*sigma, mean+2.*sigma,r->GetParams(), dataBG.GetMatrixArray())/_hf->GetXaxis()->GetBinWidth(0);
  cout<<"Signal: "<<signal_v<<endl;
  cout<<"Signal_error: "<<signal_v_error<<endl;

  double signal_v_cb = Number_CB1->Integral(mean-2.*sigma, mean+2.*sigma)/_hf_cb->GetXaxis()->GetBinWidth(0);
  double signal_v_error_cb = Number_CB1->IntegralError(mean-2.*sigma, mean+2.*sigma,rCB->GetParams(), dataBG1.GetMatrixArray())/_hf_cb->GetXaxis()->GetBinWidth(0);
  cout<<"Signal: "<<signal_v<<endl;
  cout<<"Signal_error: "<<signal_v_error<<endl;

  signal_cb->SetBinContent(i+1, signal_v_cb/dpt);//только сигнал
  signal_cb->SetBinError(i+1, signal_v_error_cb/dpt);

  signal_bg_cb->SetBinContent(i+1, (num_cb - signal_v_cb)/dpt);///
  signal_bg_cb->SetBinError(i+1, sqrt((err_cb)*(err_cb) + (signal_v_error_cb)*(signal_v_error_cb))/dpt);

  spectrum->SetBinContent(i+1, signal_v/dpt);
  spectrum->SetBinError(i+1, signal_v_error/dpt);

  spectrum_cb->SetBinContent(i+1, signal_v/dpt);
  spectrum_cb->SetBinError(i+1, signal_v_error/dpt);

  positionCB->SetBinContent(i+1, p2[1]);
  positionCB->SetBinError(i+1, ep2[1]);

  widthCB->SetBinContent(i+1, abs(p2[2]));
  widthCB->SetBinError(i+1, abs(ep2[2]));
}

void beauty_draw(TH1D* ( &h )[32][2], int i, TString name, TString hist_name)
{
  TCanvas * c1DP = new TCanvas("Sginal1","Signal1",600,800) ;
  gPad->SetBorderMode(0);
  gPad->SetBorderSize(2);
  gPad->SetFrameBorderMode(0);
  gPad->SetLeftMargin(0.18);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(0.03);
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
  box->SetMaximum(145000);


  h[i][1]->SetLineColor(kRed) ;
  h[i][1]->SetMarkerStyle(21);
  h[i][1]->SetMarkerColor(kBlack);

  h[i][1]->SetLineStyle(2);
  h[i][0]->SetLineColor(kGreen) ;

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
    box->Draw() ;
  h[i][1]->Draw() ;
  h[i][0]->Draw("same") ;

  TLegend * leg = new TLegend(0.67,0.66,0.95,0.74) ;
  leg->AddEntry(h[i][0], "CB function","l") ;
  leg->AddEntry(h[i][1], "gaus function","l") ;
  leg->SetLineColor(0) ;
  leg->Draw() ;
  TLatex t1 = TLatex();
  t1.SetNDC();
  t1.DrawLatex(0.56, 0.8, "#splitline{ALICE LHC}{work in progress}");//12.01.2021
  t1.SetTextSize(0.018) ;
  t1.Draw() ;
  TLatex t3 = TLatex();
  t3.SetNDC();
  t3.DrawLatex(0.48, 0.87, hist_name);
  t3.SetTextSize(0.018) ;
  t3.Draw() ;
  c1DP->SaveAs("pi0/LHC/"+name+Form("_%d",i)+".pdf");
}

void draw_spec(TH1D * pi0spectrum, TH1D * pi0spectrumCB, TH1D* pi0position, TH1D* pi0positionCB, TH1D* pi0width, TH1D* pi0widthCB, TString name, Double_t events){
//  TCanvas *c2 = new TCanvas("c2_results","The Fit Canvas",1200,800);
    TCanvas *c2 = new TCanvas(pi0spectrum->GetName(),pi0spectrum->GetName(),1200,800);

  c2->Divide(1,3);
  c2->cd(1)->SetLogy();
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
  pi0position->SetMinimum(0.135);
  pi0position->SetMaximum(0.137);

  TF1  *f4  = new TF1("f4","[0]", 0.5, 20);
  pi0position->Fit("f4", "RMS");
  pi0position->Draw("p e");
  pi0positionCB->SetLineColor(kGreen);
  pi0positionCB->SetMarkerStyle(22);
  pi0positionCB->SetMarkerColor(kGreen);

  pi0positionCB->Draw("p e same");

  c2->cd(3);
  //TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*TMath::Power(x, [5])", 5, 20);
  //TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*x+[5]*x*x", 0.5, 16);
  TF1 *f5 = new TF1("f5","TMath::Exp(-[0]*x)*([1]+[2]*x+[3]*x*x)", 5, 20);
  f5->SetParameter(0,0.08);
  f5->SetParameter(1,0.008);
  pi0width->GetXaxis()->SetTitle("p_{T}, GeV/c");
  pi0width->GetYaxis()->SetTitle("Width, Gev");
  pi0width->SetMinimum(0.004);
  pi0width->SetMaximum(0.007);

  pi0width->Fit("f5", "RMS");
  pi0width->SetMarkerStyle(21);
  pi0widthCB->SetLineColor(kGreen);
  pi0widthCB->SetMarkerColor(kGreen);
  pi0widthCB->SetMarkerStyle(22);

  pi0width->Draw("p e");
  TLegend * leg = new TLegend(0.67,0.66,0.95,0.74) ;
  leg->AddEntry(pi0width, "gaus", "lep");
  leg->AddEntry(pi0widthCB, "gaus", "lep");
  pi0widthCB->Draw("p e same");

  c2->SaveAs("pi0/LHC/"+name+"pi0spectrum_CB1_LHC.pdf");
}
