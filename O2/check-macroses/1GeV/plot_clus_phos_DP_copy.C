/// \file plot_dig_phos.C
/// \brief Simple macro to plot PHOS digits per event

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <sstream>
#include <iostream>

#include "TROOT.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "DataFormatsPHOS/Cluster.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsPHOS/MCLabel.h"
#include "PHOSBase/Geometry.h"
#endif

void IntFit(TH1D* _h);


void plot_clus_phos_DP_copy(int ievent = 0, TString inputfile = "phosclusters.root")
{
  using namespace std;
  TFile *fKine = new TFile("o2sim_Kine.root");
  TTree *tKine = (TTree *)fKine->Get("o2sim");
  vector<o2::MCTrack> *mcArr = nullptr;
  tKine->SetBranchAddress("MCTrack", &mcArr);

  TFile *file1 = TFile::Open(inputfile.Data());
  TTree *cluTree = (TTree *)file1->Get("o2sim");
  o2::dataformats::MCTruthContainer<o2::phos::MCLabel> *mClustersMCArray = nullptr;
  cluTree->SetBranchAddress("PHOSClusterTrueMC", &mClustersMCArray);
  vector<o2::phos::Cluster> *mClustersArray = nullptr;
  cluTree->SetBranchAddress("PHOSCluster", &mClustersArray);
  vector<o2::phos::TriggerRecord> *mClustersTR = nullptr;
  cluTree->SetBranchAddress("PHOSClusterTrigRec", &mClustersTR);

//  o2::dataformats::MCTruthContainer<o2::phos::MCLabel> *clu_energy = nullptr;
//  cluTree->SetBranchAddress("mTruthArray", &clu_energy);


  if (!mClustersArray)
  {
    cout << "PHOS clusters not in tree. Exiting ..." << endl;
    return;
  }

  cout<<"events:"<<endl;
  cout << cluTree->GetEntries() << endl;//вывод всего, что есть в дереве
  cluTree->GetEvent(0);//выделяет память
  cout << "Read " << mClustersTR->size() << " events" << endl;
  cout << "Lables " << mClustersMCArray->getNElements() << endl;

  TH1D *energyMod_MC[4];
  TH1D *energyMod_Rec[4];
  TH1D *zMod[4];
  TH1D *xMod[4];
  TH1D *res_energy[4];


  TCanvas *canv[4];

  for (int m = 0; m < 4; m++)
  {
    energyMod_MC[m] = new TH1D(Form("energyModMC%d", m + 1), Form("energy MC Mod%d", m + 1), 100, 0.6, 1.1);
    energyMod_Rec[m] = new TH1D(Form("energyMod%d", m + 1), Form("energy Rec Mod%d", m + 1), 100, 0.6, 1.1);
    zMod[m] = new TH1D(Form("zMod%d", m + 1), Form("Z Mod%d", m + 1), 1000, 0, 20);
    xMod[m] = new TH1D(Form("xMod%d", m + 1), Form("X Mod%d", m + 1), 1000, 0, 20);
    canv[m] = new TCanvas(Form("clustersInMod%d", m+1), Form("PHOS clusters in module %d", m+1), 800,600);

    res_energy[m] = new TH1D(Form("res_energy_Mod%d", m + 1), Form("res_energy Mod%d", m + 1), 1000, -1, 1);
  }

  o2::phos::Geometry *geom = new o2::phos::Geometry("PHOS");
  int iev = 0;// index of event
  int count;
  TH2D *Cluster_per_event = new TH2D(Form("Cluster per event"), Form("Cluster per event"), 1000, 0, 1100, 1000, 0, 4);

  for (auto it = mClustersTR->begin(); it != mClustersTR->end(); it++)
  {
    count++;
    cout<<"iterator - "<<*it<<endl;
    cout << " --- " << endl;
    int firstClusterInEvent = it->getFirstEntry();
    int mLastClusterInEvent = firstClusterInEvent + it->getNumberOfObjects();
    tKine->GetEntry(iev++);//заполнение листа в дереве: запись данных в массив переданных внутрь по ссылке, переводим указатель mcTrack
    const auto &mcTrack = (*mcArr)[0];//загрузили инфу
    Printf(
      "Particle: pdg = %d, E = %f, x = %f, y = %f, z = %f",
      mcTrack.GetPdgCode(), mcTrack.GetEnergy(),
      460. * mcTrack.GetStartVertexMomentumX() / mcTrack.GetEnergy(),
      460. * mcTrack.GetStartVertexMomentumY() / mcTrack.GetEnergy(),
      460. * mcTrack.GetStartVertexMomentumZ() / mcTrack.GetEnergy()
    );

    Cluster_per_event->Fill(count, mLastClusterInEvent-firstClusterInEvent);

    for (int i = firstClusterInEvent; i < mLastClusterInEvent; i++) //бегаем по кластерам в событии
    {
      auto &clu = mClustersArray->at(i);
      float x, z;
      int mod;
      clu.getLocalPosition(x, z);//meh
      mod = clu.module();//модуль PHOS
      short absId;
      TVector3 helper;
      geom->local2Global(mod, x, z, helper);
      geom->relPosToAbsId(mod, x, z, absId);//преобразование в абсолютные координаты
      char relid[3];
      geom->absToRelNumbering(absId, relid);//преобразование в относительные координаты
      cout << "mod: " << mod << " x: " << x << " z: " << z << " E = " << clu.getEnergy() << endl;
      //MC info
      //Get list of labels corresponding to this cluster
      gsl::span<const o2::phos::MCLabel> spDigList = mClustersMCArray->getLabels(i);
      cout << "Primaries " << spDigList.size() << endl;
      for (auto &l : spDigList)//выводим свойства кластера (лейблы)
      {
        cout <<
        "Label: empty=" << l.isEmpty() <<
        " noise=" << l.isNoise() <<
        " valid=" << l.isValid() <<
        " fake=" << l.isFake() <<
        " prim=" << l.getTrackID() <<
        " deposited energy " << l.getEdep() << endl;
      }

      float energy = clu.getEnergy();
      float mOccCut = 0.;
      float x1, z1, x1_loc, z1_loc;
      clu.getLocalPosition(x1_loc, z1_loc);
      x1 = 460. * helper.X() / helper.Mag();
      z1 = 460. * helper.Z() / helper.Mag();
      cout<<" --- "<<endl;
      cout<<"Cluster: x1: "<<x1<<" z1: "<<z1<<" x1_loc: "<<x1_loc<<" z1_loc: "<<z1_loc<<endl;

      float energy_MC=mcTrack.GetEnergy();
      float x2 = 460. * mcTrack.GetStartVertexMomentumX() / mcTrack.GetEnergy();
      float z2 = 460. * mcTrack.GetStartVertexMomentumZ() / mcTrack.GetEnergy();

      cout<<"Monte_Carlo x2: "<<x2<<" z2: "<<z2<<endl;
      cout<<" --- "<<endl;
      if (0<= mod && mod <4 && energy > mOccCut)
      {
          energyMod_Rec[mod]->Fill(energy);
          energyMod_MC[mod]->Fill(energy_MC);
          res_energy[mod]->Fill((energy-energy_MC)/energy_MC);
          if(TMath::Abs(energy-energy_MC)/energy_MC<=0.1){
            xMod[mod]->Fill(TMath::Abs(x1-x2));
            zMod[mod]->Fill(TMath::Abs(z1-z2));
          }
      }
    }

  }

  TString name;
  for (int m = 0; m < 4; m++)
  {
    canv[m]->Divide(2,2);
    canv[m]->cd(1);
    energyMod_Rec[m]->SetLineColor(kGreen);

    energyMod_Rec[m]->Draw("colz");//e
    energyMod_MC[m]->Draw("same");//e
    name.Form("1GeV_%d.pdf", m);
    if(m==2)
    {
      IntFit(energyMod_Rec[m]);
    }
    canv[m]->cd(2);
    res_energy[m]->Draw("colz");
    canv[m]->cd(3);
    xMod[m]->Draw("colz");

    canv[m]->cd(4);
    zMod[m]->Draw("colz");

    canv[m]->SaveAs(name);
  }
  cout<<"count = "<<count<<endl;
  TCanvas *cluster_per_event = new TCanvas("Cluster per event","Cluster per event",600,400);
//  TH1D* Cluster_per_event1=Cluster_per_event->ProjectionX();
  cluster_per_event->cd(0);
  Cluster_per_event->SetMarkerStyle(21);
  //Cluster_per_event->GetXaxis()->SetTitle("events");
  //Cluster_per_event->GetYaxis()->SetTitle("clusters");
  Cluster_per_event->Draw();
  Cluster_per_event->SaveAs("Cluster_per_event_1GeV.pdf");
}


void IntFit(TH1D* _h)
{
  TF1  *f1 = new TF1("f1","gaus", 0.6, 1.1);
      f1->SetLineStyle(2);
//  TF1 *Number= new TF1("f2","gaus", 0.35, 0.55);

     f1->SetParLimits(1, 0.85, 1.);
     f1->SetParLimits(2, 0.04, 0.06);

    _h->Fit("f1", "r0");

  TFitResultPtr r = _h->Fit("f1", "RMS+");//check if f1 works
  Double_t p[3],ep[3];
  Double_t p2[8],ep2[8];
//  Double_t Number;
  TMatrixDSym covrfit = r->GetCovarianceMatrix();
  TMatrixTSym<double> dataBG;
  covrfit.GetSub(0,2,dataBG);///
  //TF1 *Number= new TF1(fname2,"gaus", 0.1, 0.2);

  for(int j=0; j<3; j++)
  {
      p[j]=f1->GetParameter(j);
      ep[j]=f1->GetParError(j);
    //  Number->SetParameter(j,p[j]);//to be corrected later
  }
  double mean = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);

  double num = f1->Integral(mean-2.*sigma, mean+2.*sigma)/_h->GetXaxis()->GetBinWidth(0);
  double num_err = f1->IntegralError(mean-2.*sigma, mean+2.*sigma,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/_h->GetXaxis()->GetBinWidth(0);

  cout<<"#############################################Integral: "<< num  <<endl;
  cout<<"#############################################Integral error: "<< num_err <<endl;

}
