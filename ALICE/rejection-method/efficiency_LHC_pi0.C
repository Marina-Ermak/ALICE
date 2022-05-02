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
void draw_spec(TH1D * pi0spectrum, TH1D * pi0spectrum_both, TH1D* pi0position,
              TH1D* pi0position_both, TH1D* pi0width, TH1D* pi0width_both, TString name,
              Double_t events, TH1D* _hf_eff, Double_t events1,
            TH1D * pi0spectrum_lhc, TH1D * pi0spectrum_both_lhc, Double_t events2);
TH1D* fill_hist(TH2D* _h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i);
void beauty_draw(TH1D* ( &h )[32][4], int i, TString name, TString hist_name);


void efficiency_LHC_pi0() {
 gStyle->SetOptFit(0);
 gStyle->SetOptStat(0);

TFile *_file1 = TFile::Open("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/for_efficiency/pi0/histos.root");//watch the name!
THashList *_l1 =(THashList*)_file1->FindObjectAny("PHOSEta;1");

TH1D* _h_eff = (TH1D*)_l1->FindObject("Pi0GeneratedSpectra");
_h_eff->Sumw2();
_h_eff->GetXaxis()->SetLimits(0., 20.);
_h_eff->SetTitle("gen");
_h_eff->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");

TH1D* _h_eff_mod = (TH1D*) _h_eff->Clone();

for (int i=0; i<28; i++)
    {
        _h_eff_mod->SetBinContent(i, _h_eff_mod->GetBinContent(i)/_h_eff->GetBinWidth(i));
    }

TFile *_file0 = TFile::Open("histos_MC_09102020.root");//watch the name!
THashList *_l =(THashList*)_file0->FindObjectAny("PHOSEta;1");//watch the name!
TH2D *_h_all=(TH2D*)_l->FindObject("InvMass_All_E300");
TH2D *_h_both=(TH2D*)_l->FindObject("InvMass_Both_E300");

TFile *_file2 = TFile::Open("LHC16_17.root");//watch the name!
THashList *_l2 =(THashList*)_file2->FindObjectAny("PHOSEtaMB;1");//watch the name!
TH2D *_h_all_lhc=(TH2D*)_l2->FindObject("InvMass_All_E300");
TH2D *_h_both_lhc=(TH2D*)_l2->FindObject("InvMass_Both_E300");


Int_t bins=27;
double pTbins[28]={};
pTbins[0]=0.5;// {0.5, 0.7, 0.9, 1.2, 1.5, 2, 2.5, 3.2, 4., 5., 7., 10.}
 for (int i=1; i<28; i++)
     {
         if(i<18)
          pTbins[i]=pTbins[i-1]+0.5;
         else
          pTbins[i]=pTbins[i-1]+1;
         cout<<" pTBins ["<<i<<"] = "<<pTbins[i]<<endl;

     }


TH1D *_ev1=(TH1D*)_l1->FindObject("hSelEvents");
   Double_t events1=(Double_t)_ev1->GetBinContent(2);
    cout<< "!!! Number of events  - " << events1<< "!!!"<<endl;// events1 -- histos.root

    TH1D *_ev=(TH1D*)_l->FindObject("hSelEvents");
       Double_t events=(Double_t)_ev->GetBinContent(2);
        cout<< "!!! Number of events  - " << events<< "!!!"<<endl;//events -- MC

TH1D *_ev2=(TH1D*)_l2->FindObject("hSelEvents");
 Double_t events2=(Double_t)_ev2->GetBinContent(2);
cout<< "!!! Number of events  - " << events2<< "!!!"<<endl;//events2 -- LHC

   Int_t bin1=0;
   Int_t bin2=0;
   Int_t bin1_lhc=0;
   Int_t bin2_lhc=0;
   TString hist_name;



    TH1D* signal_bg_all = new TH1D("signal/background bg all","sg/bg all", bins, pTbins);
    TH1D* signal_bg_both = new TH1D("signal/background bg Disp&CPV","sg/bg Disp&CPV", bins, pTbins);
    TH1D* signal_all = new TH1D("signal all ","sg all", bins, pTbins);
    TH1D* signal_both = new TH1D("signal Disp&CPV","sg Disp&CPV", bins, pTbins);

    TH1D* pi0_spectrum_all = new TH1D("pi0_spectrum_all","pi0 spectrum all", bins, pTbins);
    TH1D* pi0_spectrum_both = new TH1D("pi0_spectrum_Disp&CPV","pi0 spectrum Disp&CPV", bins, pTbins);

    TH1D* position = new TH1D("pi0position all","pi0 position all", bins, pTbins);//11
    TH1D* width = new TH1D("pi0width all", "pi0 width all",bins, pTbins);//watch the number of bins!
    TH1D* position_both = new TH1D("pi0position Disp&CPV","pi0 position Disp&CPV", bins, pTbins);//11
    TH1D* width_both = new TH1D("pi0width Disp&CPV", "pi0 width Disp&CPV", bins, pTbins);//watch the number of bins!


    TH1D *_hf_all;
    TH1D *_hf_both;

    TH1D* signal_bg_all_lhc = new TH1D("signal/background bg all _lhc","sg/bg all _lhc", bins, pTbins);
    TH1D* signal_bg_both_lhc = new TH1D("signal/background bg Disp&CPV _lhc","sg/bg Disp&CPV _lhc", bins, pTbins);
    TH1D* signal_all_lhc = new TH1D("signal all_lhc ","sg all_lhc", bins, pTbins);
    TH1D* signal_both_lhc = new TH1D("signal Disp&CPV_lhc","sg Disp&CPV_lhc", bins, pTbins);

    TH1D* pi0_spectrum_all_lhc = new TH1D("pi0_spectrum_all_lhc","pi0 spectrum all_lhc", bins, pTbins);
    TH1D* pi0_spectrum_both_lhc = new TH1D("pi0_spectrum_Disp&CPV_lhc","pi0 spectrum Disp&CPV_lhc", bins, pTbins);

    TH1D* position_lhc = new TH1D("pi0position all_lhc","pi0 position all_lhc", bins, pTbins);//11
    TH1D* width_lhc = new TH1D("pi0width all_lhc", "pi0 width all_lhc",bins, pTbins);//watch the number of bins!
    TH1D* position_both_lhc = new TH1D("pi0position Disp&CPV_lhc","pi0 position Disp&CPV_lhc", bins, pTbins);//11
    TH1D* width_both_lhc = new TH1D("pi0width Disp&CPV_lhc", "pi0 width Disp&CPV_lhc", bins, pTbins);//watch the number of bins!

    TH1D *_hf_all_lhc;
    TH1D *_hf_both_lhc;

    TH1D * h[32][4];

for(int i=0; i<bins-1; i++)
    {

      bin1 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i]);
      bin2 = (Int_t) _h_all->GetYaxis()->FindBin(pTbins[i+1]);

      bin1_lhc = (Int_t) _h_all_lhc->GetYaxis()->FindBin(pTbins[i]);
      bin2_lhc = (Int_t) _h_all_lhc->GetYaxis()->FindBin(pTbins[i+1]);

          hist_name.Form("%.2f < p_{T} < %.2f GeV/c",pTbins[i],pTbins[i+1]);
          cout<<"Hist_name - "<< hist_name<<endl;

        _hf_all=fill_hist(_h_all, bin1, bin2, "all", hist_name, i);
        _hf_both=fill_hist(_h_both, bin1, bin2, "both", hist_name, i);

        _hf_all->SetLineColor(kRed);
        _hf_both->SetLineColor(kBlue);
        _hf_all->SetMarkerColor(kRed);
        _hf_both->SetMarkerColor(kBlue);

        _hf_all_lhc=fill_hist(_h_all_lhc, bin1_lhc, bin2_lhc, "all_lhc", hist_name, i);
        _hf_both_lhc=fill_hist(_h_both_lhc, bin1_lhc, bin2_lhc, "both_lhc", hist_name, i);

        _hf_all_lhc->SetLineColor(kBlack);
        _hf_both_lhc->SetLineColor(kMagenta);
        _hf_all_lhc->SetMarkerColor(kBlack);
        _hf_both_lhc->SetMarkerColor(kMagenta);

        double dpt = pTbins[i+1]-pTbins[i];

        IntFit(_hf_all, "all", signal_all, signal_bg_all, dpt, i, pi0_spectrum_all, position, width);
        IntFit(_hf_both, "both", signal_both, signal_bg_both, dpt, i, pi0_spectrum_both, position_both, width_both);

        IntFit(_hf_all_lhc, "all_lhc", signal_all_lhc, signal_bg_all_lhc, dpt, i, pi0_spectrum_all_lhc, position_lhc, width_lhc);
        IntFit(_hf_both_lhc, "both_lhc", signal_both_lhc, signal_bg_both_lhc, dpt, i, pi0_spectrum_both_lhc, position_both_lhc, width_both_lhc);

      //  h[i][0]=_hf_all;
      //  h[i][1]=_hf_both;

        h[i][0]=_hf_all_lhc;
        h[i][1]=_hf_both_lhc;


        beauty_draw(h, i, "All and Disp&CPV", hist_name);//!

    }
draw_spec(signal_all, signal_both, position, position_both, width, width_both, "all and disp&cpv", events, _h_eff_mod, events1,
pi0_spectrum_all_lhc, pi0_spectrum_both_lhc, events2);

}

void draw_spec(TH1D * pi0spectrum, TH1D * pi0spectrum_both, TH1D* pi0position,
  TH1D* pi0position_both, TH1D* pi0width, TH1D* pi0width_both, TString name,
  Double_t events, TH1D* _hf_eff, Double_t events1,
TH1D * pi0spectrum_lhc, TH1D * pi0spectrum_both_lhc, Double_t events2)
{
//  TCanvas *c2 = new TCanvas("c2_results","The Fit Canvas",1200,800);

    TCanvas *c3 = new TCanvas("Spectrum pi0 all  both","Spectrum pi0 all both",1200,800);
    c3->Divide(1,3);

    c3->cd(1)->SetLogy();

    pi0spectrum->Scale(1./events);
    pi0spectrum_both->Scale(1./events);
  //  pi0spectrum_both_lhc->Scale(1./events2);

    TH1D* cp_spec = (TH1D*) pi0spectrum->Clone();
    TH1D* cp_spec_both = (TH1D*) pi0spectrum_both->Clone();

    pi0spectrum->GetXaxis()->SetTitle("p_{T}, GeV/c");
    pi0spectrum->GetYaxis()->SetTitle("dN/dp_{T}");
    pi0spectrum->SetMarkerStyle(21);
    pi0spectrum->Draw("p e");

    pi0spectrum_both->SetLineColor(kGreen);
    pi0spectrum_both->SetMarkerColor(kGreen);
    pi0spectrum_both->SetMarkerStyle(22);
    pi0spectrum_both->Draw("same p e");

    c3->cd(2);
    pi0position->GetXaxis()->SetTitle("p_{T}, GeV/c");
    pi0position->GetYaxis()->SetTitle("Peak position, pTbins, Gev");
    pi0position->SetMarkerStyle(21);
    pi0position->SetMinimum(0.135);
    pi0position->SetMaximum(0.137);

  //  TF1  *f4  = new TF1("f4","[0]", 0.5, 10);
    pi0position->SetMinimum(0.132);
    pi0position->SetMaximum(0.142);
    pi0position->Fit("f4", "RMS");
    pi0position->Draw("p e");
    pi0position_both->SetLineColor(kGreen);
    pi0position_both->SetMarkerStyle(22);
    pi0position_both->SetMarkerColor(kGreen);
    pi0position_both->Draw("p e same");

    c3->cd(3);
    //TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*TMath::Power(x, [5])", 5, 20);
    //TF1 *f5 = new TF1("f5","[2]*TMath::Exp(-[0]*x+[1])+[3]+[4]*x+[5]*x*x", 0.5, 16);

    pi0width->GetXaxis()->SetTitle("p_{T}, GeV/c");
    pi0width->GetYaxis()->SetTitle("Width, Gev");
    pi0width->SetMinimum(0.);
    pi0width->SetMaximum(0.025);

  //  pi0width->Fit("f5", "RMS");
    pi0width->SetMarkerStyle(21);
    pi0width_both->SetLineColor(kRed);
    pi0width_both->SetMarkerColor(kRed);
    pi0width_both->SetMarkerStyle(22);

    pi0width->Draw("p e");

    TLegend * leg = new TLegend(0.67,0.66,0.95,0.74);
    leg->AddEntry(pi0width, "gaus", "lep");
  //  leg->AddEntry(pi0width_both, "gaus", "lep");
  //  pi0width_both->Draw("p e same");
    //leg->Draw("p e same");
    cout<<"work!"<<endl;
    c3->SaveAs("pi0spectrum_BOTH_MC.pdf");

    TCanvas *c2 = new TCanvas("Efficiency_pi0_MC","Efficiency_pi0_MC",1200,800);

    //c2->Divide(1,3);
    c2->cd(1)->SetLogy();

    _hf_eff->Scale(1./events1);
    _hf_eff->GetXaxis()->SetLimits(0., 20.);

    pi0spectrum->Divide(_hf_eff);
    pi0spectrum_both->Divide(_hf_eff);

    pi0spectrum->GetXaxis()->SetTitle("p_{T}, GeV/c");
    pi0spectrum->GetYaxis()->SetTitle("dN/dp_{T}");
    pi0spectrum->SetLineColor(kRed);
    pi0spectrum->SetMarkerColor(kRed);
    pi0spectrum->GetXaxis()->SetLimits(0., 20.);
    pi0spectrum->SetMarkerStyle(21);

    pi0spectrum_both->SetLineColor(kBlue);
    pi0spectrum_both->SetMarkerColor(kBlue);
    pi0spectrum_both->GetXaxis()->SetLimits(0., 20.);
    pi0spectrum_both->SetMarkerStyle(20);

    TF1 *f_pi0 = new TF1("f_pi0","[0]+[1]*log(x)", 0.5, 19.);
 //  TF1 *f_pi0 = new TF1("f_pi0","[0]+[1]*x+[2]*pow(x, 2)+[3]*pow(x, 3)", 0.5, 19.);
    pi0spectrum->Fit("f_pi0", "RMS+");
    TF1 *f_pi0_both = new TF1("f_pi0_both","[0]+[1]*log(x)", 0.5, 19.);
    TLegend * leg1 = new TLegend(0.5,0.36,0.8,0.3) ;
    pi0spectrum_both->Fit("f_pi0_both", "RMS+");
    leg1->AddEntry(pi0spectrum, "all", "lep");
    leg1->AddEntry(pi0spectrum_both, "Disp&CPV", "lep");

    pi0spectrum->Draw("p e");
    pi0spectrum_both->Draw("p e same");
    leg1->Draw("p e same");

    c2->SaveAs("Efficiency_pi0_MC.pdf");

    TCanvas *c = new TCanvas("Corrected_spectra_pi0_MC","Corrected_spectra_pi0_MC",1200,800);

    c->SetLogy();

    pi0spectrum_lhc->Scale(1./events2);
    pi0spectrum_both_lhc->Scale(1./events2);

    pi0spectrum_lhc->GetXaxis()->SetTitle("p_{T}, GeV/c");
    pi0spectrum_lhc->GetYaxis()->SetTitle("dN/dp_{T}");
    pi0spectrum_lhc->GetYaxis()->SetRangeUser(1e-10, 1.);

    pi0spectrum_lhc->Divide(f_pi0);
    pi0spectrum_both_lhc->Divide(f_pi0_both);

    pi0spectrum_lhc->GetXaxis()->SetTitle("p_{T}, GeV/c");
    pi0spectrum_lhc->GetYaxis()->SetTitle("dN/dp_{T}");
    pi0spectrum_lhc->SetLineColor(kRed);
    pi0spectrum_lhc->SetMarkerColor(kRed);
    pi0spectrum_lhc->GetXaxis()->SetLimits(0., 20.);
    pi0spectrum_lhc->SetMarkerStyle(21);

    pi0spectrum_both_lhc->SetLineColor(kBlue);
    pi0spectrum_both_lhc->SetMarkerColor(kBlue);
    pi0spectrum_both_lhc->GetXaxis()->SetLimits(0., 20.);
    pi0spectrum_both_lhc->SetMarkerStyle(20);

    _hf_eff->SetLineColor(kMagenta);
    _hf_eff->SetMarkerColor(kMagenta);
    _hf_eff->GetXaxis()->SetLimits(0., 20.);
    _hf_eff->SetMarkerStyle(22);

    TLegend * leg2 = new TLegend(0.5,0.5,0.9,0.9) ;
    leg2->AddEntry(pi0spectrum_lhc, "all LHC", "lep");
    leg2->AddEntry(pi0spectrum_both_lhc, "Disp&CPV LHC", "lep");
    leg2->AddEntry(_hf_eff, "Generated", "lep");

    pi0spectrum_lhc->Draw("p e");
    pi0spectrum_both_lhc->Draw("p e same");
    _hf_eff->Draw("p e same");
    leg2->Draw("p e same");

    c->SaveAs("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/for_efficiency/pi0/corrected_spectra_lhc.pdf");
 }



TH1D* fill_hist(TH2D* _h, Int_t bin1, Int_t bin2, TString name, TString hist_name, int i)
{
    TH1D*_hf=(TH1D*)_h->ProjectionX(name+Form("_%d",i),bin1,bin2);
    _hf->Rebin(10);
    _hf->Sumw2();
    _hf->SetName(name+Form("%d",i));//all SetName-s are different
    cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA:   "<<name+Form("%d",i)<<endl;
    _hf->GetXaxis()->SetRangeUser(0.08, 0.25);
    _hf->SetTitle(name+Form("%d",i));
    _hf->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV/c^{2})");
    return _hf;
}



void IntFit(TH1D* _hf, char* name, TH1D* signal, TH1D* signal_bg, double dpt, int i,
            TH1D * spectrum, TH1D* position, TH1D* width)
{
    char * fname = Form("f1_%s_%d", name,i);
    char * fname2 = Form("N1_%s_%d", name,i);
    cout<<"fname "<< fname <<endl;
    cout<<"fname2 "<< fname2 <<endl;
    TF1  *f1 = new TF1(fname,"gaus+[3]*(x-[1])+[4]+[5]*(x-[1])*(x-[1])", 0.1, 0.2);
    f1->SetLineStyle(2);

    f1->SetParLimits(1, 0.13, 0.15);
    f1->SetParLimits(2, 1e-03, 20e-03);
    f1->SetParameters(1.*_hf->GetBinContent(66),0.135,0.006,1.,3.,2.*_hf->GetBinContent(51)) ;

    _hf->Fit(fname, "r0");

    TFitResultPtr r = _hf->Fit(fname, "RMS+");//check if f1 works
    Double_t p[3],ep[3];
  //  Double_t p2[8],ep2[8];
    TMatrixDSym covrfit = r->GetCovarianceMatrix();
    TMatrixTSym<double> dataBG;
    covrfit.GetSub(0,2,dataBG);///
    TF1 *Number = new TF1(fname2,"gaus", 0.1, 0.2);

    for(int j=0; j<3; j++)
    {
        p[j]=f1->GetParameter(j);
        ep[j]=f1->GetParError(j);
        Number->SetParameter(j,p[j]);//to be corrected later
    }
    double mean = f1->GetParameter(1);
    double sigma = f1->GetParameter(2);

    double num = f1->Integral(mean-2.*sigma, mean+2.*sigma)/_hf->GetXaxis()->GetBinWidth(0);
    double num_err = f1->IntegralError(mean-2.*sigma, mean+2.*sigma,
      r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/_hf->GetXaxis()->GetBinWidth(0);

    cout<<"Integral: "<< num  <<endl;
    cout<<"Integral error: "<< num_err <<endl;

    position->SetBinContent(i+1, p[1]);
    position->SetBinError(i+1, ep[1]);

    width->SetBinContent(i+1, abs(p[2]));
    width->SetBinError(i+1, abs(ep[2]));

    double signal_v = Number->Integral(mean-2.*sigma,
      mean+2.*sigma)/_hf->GetXaxis()->GetBinWidth(0);
    double signal_v_error = Number->IntegralError(mean-2.*sigma,
      mean+2.*sigma,r->GetParams(), dataBG.GetMatrixArray())/_hf->GetXaxis()->GetBinWidth(0);

    cout<<"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB Signal: "<<signal_v<<endl;
    cout<<"BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBSignal_error: "<<signal_v_error<<endl;

    signal->SetBinContent(i+1, signal_v/ dpt);//только сигнал
    signal->SetBinError(i+1, signal_v_error/ dpt);

    signal_bg->SetBinContent(i+1, num/ dpt - signal_v/ dpt);///
    signal_bg->SetBinError(i+1, sqrt((num_err/ dpt)*(num_err/ dpt) + (signal_v_error/ dpt)*(signal_v_error/ dpt)));///!!!

    spectrum->SetBinContent(i+1, signal_v/ dpt);
    spectrum->SetBinError(i+1, signal_v_error/ dpt);
}



void beauty_draw(TH1D* ( &h )[32][4], int i, TString name, TString hist_name)
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


  //  h[i][1]->SetLineColor(kRed) ;
    h[i][1]->SetMarkerStyle(21);
  //  h[i][1]->SetMarkerColor(kBlack);

    h[i][1]->SetLineStyle(2);
  //  h[i][0]->SetLineColor(kGreen) ;
/*
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

  */  box->Draw() ;
    h[i][1]->Draw() ;
    h[i][0]->Draw("same") ;

    TLegend * leg = new TLegend(0.67,0.66,0.95,0.74) ;
    leg->AddEntry(h[i][0], "all function","l") ;
    leg->AddEntry(h[i][1], "both function","l") ;
    leg->SetLineColor(0) ;
    leg->Draw() ;
    TLatex t1 = TLatex();
    t1.SetNDC();
    t1.DrawLatex(0.56, 0.8, "#splitline{ALICE MC}{work in progress}");//12.01.2021
    t1.SetTextSize(0.018) ;
    t1.Draw() ;
    TLatex t3 = TLatex();
    t3.SetNDC();
    t3.DrawLatex(0.48, 0.87, hist_name);
    t3.SetTextSize(0.018) ;
    t3.Draw() ;
    c1DP->SaveAs("pi0/model/"+name+Form("mag_%d",i)+".pdf");
}
