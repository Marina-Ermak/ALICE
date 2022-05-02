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

void IntFit(TH1D* _hf, char* name, TH1D* signal, TH1D* signal_bg, double dpt, int i, TH1D * spectrum,
            TH1D* position, TH1D* width);
void draw_spec(TH1D * etaspectrum, TH1D * etaspectrum_both, TH1D * etaspectrum_true, TH1D* etaposition, TH1D* etaposition_both,
                           TH1D* etapsiotion_true, TH1D* etawidth, TH1D* etawidth_both, TH1D* etawidth_true, TH1D *
                           etaspectrum_R, TH1D * etaspectrum_both_R, TH1D * etaspectrum_true_R, TH1D* etaposition_R, TH1D* etaposition_both_R,
                           TH1D* etaposition_true_R,TH1D* etawidth_R, TH1D* etawidth_both_R, TH1D* etawidth_true_R,
                           TString name, Double_t events, TH1D* _hf_eff, Double_t events1);
TH1D* fill_hist(TH2D* _h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i);
void beauty_draw(TH1D* ( &h )[32][6], int i, TString name, TString hist_name);


void efficiency_eta() {
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0);

TFile *_file1 = TFile::Open("histos_new.root");//watch the name!
THashList *_l1 =(THashList*)_file1->FindObjectAny("PHOSEta;1");

TH1D* _h_eff = (TH1D*)_l1->FindObject("EtaGeneratedSpectra");
_h_eff->Sumw2();
_h_eff->GetXaxis()->SetLimits(0., 15.);
_h_eff->SetTitle("gen_eta");
_h_eff->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");

TFile *_file0 = TFile::Open("LHC18f5b_2.root");//watch the name!
THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEta;1");

//TH1D *_h_eff=(TH1D*)_l->FindObject("EtaGeneratedSpectra");

TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");
TH2D *_h_both=(TH2D*)_l->FindObject("InvMass_Both_E300");
TH2D *_h_true=(TH2D*)_l->FindObject("True_eta_Inv_mass");
TH2D *_h_all_R=(TH2D*)_l->FindObject("InvMass_All_Rejection_E300");
TH2D *_h_both_R=(TH2D*)_l->FindObject("InvMass_Both_Rejection_E300");
TH2D *_h_true_R=(TH2D*)_l->FindObject("True_eta_Inv_mass_Rejection");

Int_t bins=4;
double pTbins[25]={1.5,3.,4.5,6.,7.5};

for (int i=0; i<5; i++)
    {
        _h_eff->SetBinContent(i, _h_eff->GetBinContent(i)/_h_eff->GetBinWidth(i));
    }
/*
pTbins[0]=0.5;// {0.5, 0.7, 0.9, 1.2, 1.5, 2, 2.5, 3.2, 4., 5., 7., 10.}
 for (int i=1; i<11; i++)
     {
         if(i<9)
          pTbins[i]=pTbins[i-1]+1;
         else
          pTbins[i]=pTbins[i-1]+3;
         cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;
     }
*/

TH1D *_ev=(TH1D*)_l1->FindObject("hSelEvents");
   Double_t events1=(Double_t)_ev->GetBinContent(2);
    cout<< "!!! Number of events  - " << events1<< "!!!"<<endl;

    TH1D *_ev1=(TH1D*)_l->FindObject("hSelEvents");
       Double_t events=(Double_t)_ev1->GetBinContent(2);
        cout<< "!!! Number of events  - " << events<< "!!!"<<endl;

   Int_t bin1=0;
   Int_t bin2=0;
TString hist_name;

    TH1D* signal_bg_all = new TH1D("signal/background bg all "," all", bins, pTbins);
    TH1D* signal_bg_both = new TH1D("signal/background bg Disp&CPV "," Disp&CPV", bins, pTbins);
    TH1D* signal_bg_true = new TH1D("signal/background bg True"," True", bins, pTbins);
    TH1D* signal_all = new TH1D("signal all "," all", bins, pTbins);
    TH1D* signal_both = new TH1D("signal Disp&CPV "," Disp&CPV", bins, pTbins);
    TH1D* signal_true = new TH1D("signal True"," True", bins, pTbins);
    TH1D* signal_bg_all_R = new TH1D("signal/background bg all R "," all R", bins, pTbins); //REJECTION
    TH1D* signal_bg_both_R = new TH1D("signal/background bg Disp&CPV R "," Disp&CPV R", bins, pTbins);
    TH1D* signal_bg_true_R = new TH1D("signal/background bg True R"," True R", bins, pTbins);
    TH1D* signal_all_R = new TH1D("signal all R"," all R", bins, pTbins);
    TH1D* signal_both_R= new TH1D("signal both R"," both R", bins, pTbins);
    TH1D* signal_true_R = new TH1D("signal True R"," True R", bins, pTbins);

    TH1D* eta_spectrum_all = new TH1D("eta_spectrum_all ","eta spectrum all ", bins, pTbins);
    TH1D* eta_spectrum_both = new TH1D("eta_spectrum_both ","eta spectrum both ", bins, pTbins);
    TH1D* eta_spectrum_true = new TH1D("eta_spectrum_true ","eta spectrum true ", bins, pTbins);
    TH1D* eta_spectrum_all_R = new TH1D("eta_spectrum_all R","eta spectrum all R", bins, pTbins); //REJECTION
    TH1D* eta_spectrum_both_R = new TH1D("eta_spectrum_both R","eta spectrum both R", bins, pTbins);
    TH1D* eta_spectrum_true_R = new TH1D("eta_spectrum_true R","eta spectrum true R", bins, pTbins);

    TH1D* position = new TH1D("etaposition","eta position ", bins, pTbins);//11
    TH1D* width = new TH1D("etawidth", "eta width",bins, pTbins);//watch the number of bins!
    TH1D* position_both = new TH1D("etaposition Disp&CPV","eta position Disp&CPV", bins, pTbins);//11
    TH1D* width_both = new TH1D("etawidth Disp&CPV", "eta width Disp&CPV", bins, pTbins);//watch the number of bins!
    TH1D* position_true = new TH1D("etaposition True","eta position True", bins, pTbins);//11
    TH1D* width_true = new TH1D("etawidth True", "eta width True", bins, pTbins);
    //REJECTION
    TH1D* position_R = new TH1D("etaposition R","eta position R", bins, pTbins);//11
    TH1D* width_R = new TH1D("etawidth R", "eta width R",bins, pTbins);//watch the number of bins!
    TH1D* position_both_R = new TH1D("etaposition Disp&CPV R","eta position Disp&CPV R", bins, pTbins);//11
    TH1D* width_both_R = new TH1D("etawidth Disp&CPV R", "eta width Disp&CPV R", bins, pTbins);//watch the number of bins!
    TH1D* position_true_R = new TH1D("etaposition True R","eta position True R", bins, pTbins);//11
    TH1D* width_true_R = new TH1D("etawidth True R", "eta width True R", bins, pTbins);//watch the number of bins!

    TH1D *_hf_all;
    TH1D *_hf_both;
    TH1D *_hf_true;

    TH1D *_hf_all_R;
    TH1D *_hf_both_R;
    TH1D *_hf_true_R;

    TH1D * h[32][6];

for(int i=0; i<bins; i++)
    {

      bin1 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i]);
      bin2 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i+1]);

        hist_name.Form("%.2f < p_{T} < %.2f GeV/c",pTbins[i],pTbins[i+1]);
        cout<<"########################################################################################Hist_name - "<< hist_name<<endl;

        _hf_all=fill_hist(_h_all, bin1, bin2, "all", hist_name, i);
        _hf_both=fill_hist(_h_both, bin1, bin2, "both", hist_name, i);
        _hf_true=fill_hist(_h_true, bin1, bin2, "true", hist_name, i);
        _hf_all_R=fill_hist(_h_all_R, bin1, bin2, "all_R", hist_name, i);
        _hf_both_R=fill_hist(_h_both_R, bin1, bin2, "both_R", hist_name, i);
        _hf_true_R=fill_hist(_h_true_R, bin1, bin2, "true_R", hist_name, i);

        _hf_all->SetLineColor(kRed);
        _hf_both->SetLineColor(kBlue);
        _hf_true->SetLineColor(kGreen);
        _hf_all->SetMarkerColor(kRed);
        _hf_both->SetMarkerColor(kBlue);
        _hf_true->SetMarkerColor(kGreen);
        _hf_all_R->SetLineColor(kRed);
        _hf_both_R->SetLineColor(kBlue);
        _hf_true_R->SetLineColor(kGreen);
        _hf_all_R->SetMarkerColor(kRed);
        _hf_both_R->SetMarkerColor(kBlue);
        _hf_true_R->SetMarkerColor(kGreen);

        double dpt = pTbins[i+1]-pTbins[i];


        IntFit(_hf_all, "all", signal_all, signal_bg_all, dpt, i, eta_spectrum_all, position, width);
        IntFit(_hf_both, "both", signal_both, signal_bg_both, dpt, i, eta_spectrum_both, position_both, width_both);
        IntFit(_hf_true, "true", signal_true, signal_bg_true, dpt, i, eta_spectrum_true, position_true, width_true);
        IntFit(_hf_all_R, "all R", signal_all_R, signal_bg_all_R, dpt, i, eta_spectrum_all_R, position_R, width_R);
        IntFit(_hf_both_R, "both R", signal_both_R, signal_bg_both_R, dpt, i, eta_spectrum_both_R, position_both_R, width_both_R);
        IntFit(_hf_true_R, "true R", signal_true_R, signal_bg_true_R, dpt, i, eta_spectrum_true_R, position_true_R, width_true_R);


        h[i][0]=_hf_all;
        h[i][1]=_hf_both;
        h[i][2]=_hf_true;
        h[i][3]=_hf_all_R;
        h[i][4]=_hf_both_R;
        h[i][5]=_hf_true_R;


        beauty_draw(h, i, "All and Disp&CPV and True", hist_name);//!

    }
draw_spec(eta_spectrum_all, eta_spectrum_both, eta_spectrum_true, position, position_both, position_true,
          width, width_both, width_true, eta_spectrum_all_R, eta_spectrum_both_R, eta_spectrum_true_R,
          position_R, position_both_R, position_true_R,width, width_both_R, width_true_R,
          "all and disp&cpv and true and rejection", events, _h_eff, events1);

}

void draw_spec(TH1D * etaspectrum, TH1D * etaspectrum_both, TH1D * etaspectrum_true, TH1D* etaposition, TH1D* etaposition_both,
               TH1D* etapsiotion_true, TH1D* etawidth, TH1D* etawidth_both, TH1D* etawidth_true, TH1D *
               etaspectrum_R, TH1D * etaspectrum_both_R, TH1D * etaspectrum_true_R, TH1D* etaposition_R, TH1D* etaposition_both_R,
               TH1D* etaposition_true_R,TH1D* etawidth_R, TH1D* etawidth_both_R, TH1D* etawidth_true_R,
               TString name, Double_t events, TH1D* _hf_eff, Double_t events1)
{
//  TCanvas *c2 = new TCanvas("c2_results","The Fit Canvas",1200,800);

    TCanvas *c3 = new TCanvas(etaspectrum->GetName(),etaspectrum->GetName(),1200,800);

    gPad->SetLogy();

    etaspectrum->Scale(1./events);
    etaspectrum_both->Scale(1./events);
    etaspectrum_true->Scale(1./events);
    etaspectrum_R->Scale(1./events);
    etaspectrum_both_R->Scale(1./events);
    etaspectrum_true_R->Scale(1./events);

    c3->cd(1);
    etaspectrum->GetXaxis()->SetTitle("p_{T}, GeV/c");
    etaspectrum->GetXaxis()->SetRangeUser(1e-10, 1);
    etaspectrum->GetYaxis()->SetTitle("dN/dp_{T}");
    etaspectrum->SetMarkerStyle(21);
    //TF1 *f3 = new TF1("f3","[0]*TMath::Exp(-[1]*x)", 0.5, 10);
    TF1 *f3 = new TF1("f3","[0]*TMath::Exp(-[1]*x)+TMath::Power([2]*(1+[3]*x*x), 13)", 0.5, 10);
    //etaspectrum->Fit("f3", "RMS");
    etaspectrum->Draw("p e1 ");
    etaspectrum_both->SetLineColor(kGreen);
    etaspectrum_both->SetMarkerColor(kGreen);
    etaspectrum_both->SetMarkerStyle(20);
    etaspectrum_both->Draw("p e1 same");

    etaspectrum_true->SetLineColor(kRed);
    etaspectrum_true->SetMarkerColor(kRed);
    etaspectrum_true->SetMarkerStyle(22);
    etaspectrum_true->Draw("p e1 same");

    TLegend * leg = new TLegend(0.7,0.6,0.9,0.9) ;
    leg->AddEntry(etaspectrum, "all", "lep");
    leg->AddEntry(etaspectrum_both, "Disp&CPV", "lep");
    leg->AddEntry(etaspectrum_true, "True", "lep");
    leg->Draw("p e same");

      TCanvas *c4 = new TCanvas("etaspectrum_R","etaspectrum_R",1200,800);
        c3->SaveAs("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/for_efficiency/eta/etaspectrum.pdf");

    gPad->SetLogy();
    etaspectrum_R->SetLineColor(kRed);
    etaspectrum_R->SetMarkerColor(kRed);
    etaspectrum_R->SetMarkerStyle(25);
      etaspectrum_R->GetXaxis()->SetRangeUser(0.00000001, 0.01);
    etaspectrum_R->Draw("p e1");


    etaspectrum_both_R->SetLineColor(kGreen);
    etaspectrum_both_R->SetMarkerColor(kGreen);
    etaspectrum_both_R->SetMarkerStyle(24);
    etaspectrum_both_R->Draw("p e1 same");

    etaspectrum_true_R->SetLineColor(kRed);
    etaspectrum_true_R->SetMarkerColor(kRed);
    etaspectrum_true_R->SetMarkerStyle(26);
    etaspectrum_true_R->Draw("p e1 same");

    TLegend * leg3 = new TLegend(0.7,0.6,0.9,0.9) ;
    leg3->AddEntry(etaspectrum_R, "al_Rl", "lep");
    leg3->AddEntry(etaspectrum_both_R, "Disp&CPV _R", "lep");
    leg3->AddEntry(etaspectrum_true_R, "True_R", "lep");
    leg3->Draw("p e same");

    c4->SaveAs("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/for_efficiency/eta/etaspectrum_rejection.pdf");
/*
    c3->cd(2);
    etaposition->GetXaxis()->SetTitle("p_{T}, GeV/c");
    etaposition->GetYaxis()->SetTitle("Peak position,9,pTbins, Gev");
    etaposition->SetMarkerStyle(21);
    etaposition->SetMinimum(0.5);
    etaposition->SetMaximum(0.7);

    TF1  *f4  = new TF1("f4","[0]", 0.5, 10);
  //  etaposition->SetMinimum(0.132);
  //  etaposition->SetMaximum(0.142);
    etaposition->Fit("f4", "RMS");
    etaposition->Draw("p e");
    etaposition_both->SetLineColor(kGreen);
    etaposition_both->SetMarkerStyle(22);
    etaposition_both->SetMarkerColor(kGreen);

    etaposition_both->Draw("p e same");

    c3->cd(3);
    //TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*TMath::Power(x, [5])", 5, 20);
    //TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*x+[5]*x*x", 0.5, 16);
    TF1 *f5 = new TF1("f5","TMath::Exp(-[0]*x)*([1]+[2]*x+[3]*x*x)", 5, 10);
//  f5->SetParameter(0,0.08);
//  f5->SetParameter(1,0.008);
    etawidth->GetXaxis()->SetTitle("p_{T}, GeV/c");
    etawidth->GetYaxis()->SetTitle("Width, Gev");
    etawidth->SetMinimum(0.);
    etawidth->SetMaximum(0.025);

    etawidth->Fit("f5", "RMS");
    etawidth->SetMarkerStyle(21);
    etawidth_both->SetLineColor(kGreen);
    etawidth_both->SetMarkerColor(kGreen);
    etawidth_both->SetMarkerStyle(22);

    etawidth->Draw("p e");
    TLegend * leg = new TLegend(0.67,0.66,0.95,0.74) ;
    leg->AddEntry(etawidth, "gaus", "lep");
    leg->AddEntry(etawidth_both, "gaus", "lep");
    etawidth_both->Draw("p e same");
    cout<<"work!"<<endl;
    c3->SaveAs("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/for_efficiency/eta/etaspectrum.pdf");



    //c2->Divide(1,3);
    c2->cd(1)->SetLogy();



    etaspectrum_R->Scale(1./events);
    etaspectrum_both_R->Scale(1./events);
    etaspectrum_true_R->Scale(1./events);
*/
    TCanvas *c2 = new TCanvas("Efficiency_eta_MC","Efficiency_eta_MC",1200,800);
    gPad->SetLogy();
    _hf_eff->Scale(1./events1);

    etaspectrum->GetXaxis()->SetLimits(0., 15.);
    etaspectrum_both->GetXaxis()->SetLimits(0., 15.);
    etaspectrum_true->GetXaxis()->SetLimits(0., 15.);
    etaspectrum_R->GetXaxis()->SetLimits(0., 15.);
    etaspectrum_both_R->GetXaxis()->SetLimits(0., 15.);
    etaspectrum_true_R->GetXaxis()->SetLimits(0., 15.);

    etaspectrum->Divide(_hf_eff);
    etaspectrum_both->Divide(_hf_eff);
    etaspectrum_true->Divide(_hf_eff);
    etaspectrum_R->Divide(_hf_eff);
    etaspectrum_both_R->Divide(_hf_eff);
    etaspectrum_true_R->Divide(_hf_eff);

    etaspectrum->GetXaxis()->SetTitle("p_{T}, GeV/c");
    etaspectrum->GetYaxis()->SetTitle("dN/dp_{T}");
    etaspectrum->SetLineColor(kRed);
    etaspectrum->SetMarkerColor(kRed);
    etaspectrum->SetMarkerStyle(21);

    etaspectrum_both->SetLineColor(kBlue);
    etaspectrum_both->SetMarkerColor(kBlue);
    etaspectrum_both->SetMarkerStyle(20);
    cout<<"here!"<<endl;
    etaspectrum_true->SetLineColor(kGreen);
    etaspectrum_true->SetMarkerColor(kGreen);
    etaspectrum_true->SetMarkerStyle(22);

    etaspectrum_R->SetLineColor(kRed);
    etaspectrum_R->SetMarkerColor(kRed);
    etaspectrum_R->SetMarkerStyle(25);

    etaspectrum_both_R->SetLineColor(kBlue);
    etaspectrum_both_R->SetMarkerColor(kBlue);
    etaspectrum_both_R->SetMarkerStyle(24);

    etaspectrum_true_R->SetLineColor(kGreen);
    etaspectrum_true_R->SetMarkerColor(kGreen);
    etaspectrum_true_R->SetMarkerStyle(26);

    TLegend * leg1 = new TLegend(0.7,0.6,0.9,0.9) ;
    leg1->AddEntry(etaspectrum, "all", "lep");
    leg1->AddEntry(etaspectrum_both, "Disp&CPV", "lep");
    leg1->AddEntry(etaspectrum_true, "True", "lep");
    leg1->AddEntry(etaspectrum_R, "al_Rl", "lep");
    leg1->AddEntry(etaspectrum_both_R, "Disp&CPV _R", "lep");
    leg1->AddEntry(etaspectrum_true_R, "True_R", "lep");

    etaspectrum->Draw("p e");
    etaspectrum_both->Draw("p e same");
    etaspectrum_true->Draw("p e same");
    etaspectrum_R->Draw("p e same");//!!!
    etaspectrum_both_R->Draw("p e same");
    etaspectrum_true_R->Draw("p e same");
    leg1->Draw("p e same");

    cout<<"here!"<<endl;

    c2->SaveAs("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/for_efficiency/eta/Efficiency_Eta_MC.pdf");

    cout<<"here!"<<endl;

 }



TH1D* fill_hist(TH2D* _h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i)
{
  TH1D*_hf=(TH1D*)_h->ProjectionX(name+Form("_%d",i),bin1,bin2);
  _hf->Rebin(10);//6 8
  _hf->Sumw2();
//  _hf->SetName(name+Form("%d",i));//all SetName-s are different
  _hf->GetXaxis()->SetRangeUser(0.4, 0.7);
  _hf->SetTitle(hist_name);
  _hf->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
  return _hf;
}



void IntFit(TH1D* _hf, char* name, TH1D* signal, TH1D* signal_bg, double dpt, int i,
            TH1D * spectrum, TH1D* position, TH1D* width)
{
    char * fname = Form("f1_%s_%d",name,i);
    char * fname2 = Form("N_%s_%d",name,i);
    TF1  *f1 = new TF1(fname,"gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.4, 0.7);
    f1->SetLineStyle(2);
  /*  if (name == "all")
      f1->SetParameters(2.5*_hf->GetBinContent(66),0.55,0.017,-2e5/(i+1),3e4/(i+1),2.5*_hf->GetBinContent(51)/(i+1)) ;
    if (name == "all R" or name == "both" or name == "true" or name == "true R")
      f1->SetParameters(1.*_hf->GetBinContent(66),0.5456,0.015,1.,3.,2.*_hf->GetBinContent(51)) ;
    if (name == "both R")
      f1->SetParameters(30.,0.5456,0.015,-2e3/(i+1),3e2/(i+1),2.5*_hf->GetBinContent(51)/(i+1)) ;//both_cut or just both?
*/
  //f1->SetParameters(30.,0.5456,0.015,0,0,2.*_hf->GetBinContent(51)) ;//both_cut or just both?
  f1->SetParLimits(1, 0.49, 0.56);
  f1->SetParLimits(0, 0,2*_hf->GetBinContent(56));
  f1->SetParameters(1.*_hf->GetBinContent(56),0.546,0.015,0,0,2.*_hf->GetBinContent(51)/(i+1));
/*
  if (i==0){
    f1->SetParLimits(1, 0.5, 0.56);
    f1->SetParLimits(0, 0,2*_hf->GetBinContent(56));
    f1->SetParameters(1.*_hf->GetBinContent(56),0.546,0.015,0,0,2.*_hf->GetBinContent(51)/(i+1));
    }
  if (i==1 and name=="all"){
      f1->FixParameter(1, 0.538);
      f1->FixParameter(2, 0.018);
    //  f1->FixParameters(10, 0.538, 0.018, -73., 18., 440.);
      }
  if (i==1 and name=="all R"){
          f1->FixParameter(1, 0.54);
          f1->FixParameter(2, 0.0234);
        //  f1->SetParameters(10, 0.538, 0.018, -73., 18., 440.);
          }
  if (i==2){
    f1->SetParameters(3.,0.546,0.015,0,0,2.*_hf->GetBinContent(51)/(i+1));//-1e4/(i+1)/(i+1),1e3/(i+1)/(i+1),1.7e4/(i+1)/(i+1));//-2e5/(i+1),3e4/(i+1),2.5*_hf->GetBinContent(51)/(i+1)) ;
  }

  if (i==2 and name=="both"){
          f1->FixParameter(0, 4.36);
          f1->FixParameter(1, 0.55);
          f1->FixParameter(2, 0.0135);
          }
if (i==2 and name=="all"){
  f1->FixParameter(0, 4.4775);
  f1->FixParameter(1, 0.53621);
  f1->FixParameter(2, 0.00685464);
  f1->FixParameter(3, -16.8183);
  f1->FixParameter(4, 4.47154);
  f1->FixParameter(5, 48.2514);
}

if (i==2 and name=="all R"){
  f1->FixParameter(0, 1.9);
  f1->FixParameter(1, 0.550);
  f1->FixParameter(2, 0.034);
  f1->FixParameter(3, -6.93);
  f1->FixParameter(4, 1.3);
  f1->FixParameter(5, 58.44);
}


if (i==3 and name=="both"){
  f1->FixParameter(0, 2.32);
  f1->FixParameter(1, 0.530);
  f1->FixParameter(2, 0.005);
  f1->FixParameter(3, -5.72);
  f1->FixParameter(4, 1.32);
  f1->FixParameter(5, 27.);
}

if (i==3 and name=="both R"){
  f1->FixParameter(0, 2);
  f1->FixParameter(1, 0.47);
  f1->FixParameter(2, 0.003);
  f1->FixParameter(3, 2.47);
  f1->FixParameter(4, 0.99);
  f1->FixParameter(5, 5.98e-5);
}
*/


    _hf->Fit(fname, "0");

    TFitResultPtr r = _hf->Fit(fname, "RMS+");//check if f1 works
    Double_t p[3],ep[3];
    TMatrixDSym covrfit = r->GetCovarianceMatrix();
    TMatrixTSym<double> dataBG;
    covrfit.GetSub(0,2,dataBG);///
    TF1 *Number = new TF1(fname2,"gaus", 0.4, 0.7);

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

    double signal_v = Number->Integral(mean-2.*sigma, mean+2.*sigma)/_hf->GetXaxis()->GetBinWidth(0);
    double signal_v_error = Number->IntegralError(mean-2.*sigma, mean+2.*sigma,r->GetParams(), dataBG.GetMatrixArray())/_hf->GetXaxis()->GetBinWidth(0);
    cout<<"Signal: "<<signal_v<<endl;
    cout<<"Signal_error: "<<signal_v_error<<endl;

    spectrum->SetBinContent(i+1, signal_v/dpt);
    spectrum->SetBinError(i+1, signal_v_error/dpt);

}



void beauty_draw(TH1D* ( &h )[32][6], int i, TString name, TString hist_name)
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
  TH1D * box = new TH1D("box","",200,0.401,0.699) ;
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
  box->SetMaximum(500);

  h[i][0]->SetMarkerStyle(24) ;
  h[i][0]->SetMarkerColor(kPink) ;
  h[i][0]->SetLineColor(kPink) ;
  h[i][0]->SetLineWidth(2) ;

  h[i][2]->SetMarkerStyle(24) ;
  h[i][2]->SetMarkerColor(kMagenta) ;
  h[i][2]->SetLineColor(kMagenta) ;
  h[i][2]->SetLineWidth(2) ;

  h[i][1]->Draw() ;//&&
  h[i][0]->Draw("same") ;
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
/*  TLatex t2 = TLatex();
  t2.SetNDC();
  t2.DrawLatex(0.227, 0.825, "12.01.2021");
  t2.SetTextSize(0.018) ;
  t2.Draw() ;
*/  TLatex t3 = TLatex();
  t3.SetNDC();
  t3.DrawLatex(0.229, 0.775, hist_name);
  t3.SetTextSize(0.018) ;
  t3.Draw() ;
  c1DP->SaveAs("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/for_efficiency/eta/"+name+Form("_%d",i)+".pdf");
}
