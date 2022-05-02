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
  TProfile *h_NE_max =  new TProfile("h_NE_max","N_{clusters}(R)", 200, 0.0, 100., 0, 10);//????
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
  TH1D *Radius[4];

//  TCanvas *canv[4];
/*
  for (int m = 0; m < 4; m++)
  {

    Radius[m] = new TH1D(Form("Radius_Mod%d", m + 1), Form("Radius Mod%d", m + 1), 1000, 0, 200);
    canv[m] = new TCanvas(Form("clustersInMod%d", m+1), Form("PHOS clusters in module %d", m+1), 800,600);
  }
  */

  o2::phos::Geometry *geom = new o2::phos::Geometry("PHOS");
  int iev = 0;// index of event
  int count;
  TH2D *Cluster_per_event = new TH2D(Form("Cluster per event"), Form("Cluster per event"), 1000, 0, 1100, 1000, 0, 4);
  TH1D *h_energy = new TH2D("Energy", 1000, 0, 11.);

  for (auto it = mClustersTR->begin(); it != mClustersTR->end(); it++)
  {
    count++;
    cout<<"iterator - "<<*it<<endl;
    cout << " --- " << endl;
    int firstClusterInEvent = it->getFirstEntry();
    int mLastClusterInEvent = firstClusterInEvent + it->getNumberOfObjects();
    tKine->GetEntry(iev++);//заполнение листа в дереве: запись данных в массив переданных внутрь по ссылке, переводим указатель mcTrack
    const auto &mcTrack1 = (*mcArr)[0];//загрузили инфу
  /*  Printf(
      "Particle1: pdg = %d, E = %f, x = %f, y = %f, z = %f",
      mcTrack1.GetPdgCode(), mcTrack1.GetEnergy(),
      460. * mcTrack1.GetStartVertexMomentumX() / mcTrack1.GetEnergy(),
      460. * mcTrack1.GetStartVertexMomentumY() / mcTrack1.GetEnergy(),
      460. * mcTrack1.GetStartVertexMomentumZ() / mcTrack1.GetEnergy()
    );
*/

    const auto &mcTrack2 = (*mcArr)[1];//загрузили инфу
    /*Printf(
      "Particle1: pdg = %d, E = %f, x = %f, y = %f, z = %f",
      mcTrack2.GetPdgCode(), mcTrack2.GetEnergy(),
      460. * mcTrack2.GetStartVertexMomentumX() / mcTrack2.GetEnergy(),
      460. * mcTrack2.GetStartVertexMomentumY() / mcTrack2.GetEnergy(),
      460. * mcTrack2.GetStartVertexMomentumZ() / mcTrack2.GetEnergy()
    );
    */

    Double_t energy_MC1=mcTrack1.GetEnergy();
    Double_t x1 = 460. * mcTrack1.GetStartVertexMomentumX() / mcTrack1.GetEnergy();
    Double_t z1 = 460. * mcTrack1.GetStartVertexMomentumZ() / mcTrack1.GetEnergy();

    Double_t energy_MC2=mcTrack2.GetEnergy();
    Double_t x2 = 460. * mcTrack2.GetStartVertexMomentumX() / mcTrack2.GetEnergy();
    Double_t z2 = 460. * mcTrack2.GetStartVertexMomentumZ() / mcTrack2.GetEnergy();



    Cluster_per_event->Fill(count, mLastClusterInEvent-firstClusterInEvent);
    cout<<"NClus"<<mLastClusterInEvent-firstClusterInEvent<<endl;

    for (int i = firstClusterInEvent; i < mLastClusterInEvent; i++) //бегаем по кластерам в событии
    {
      auto &clu = mClustersArray->at(i);
      float x, z;
      int mod;
    //  Double_t NE_max = clu.getNExMax();
      Double_t energy = clu.getEnergy();
      clu.getLocalPosition(x, z);//meh
      mod = clu.module();//модуль PHOS
      short absId;
      TVector3 helper;
      geom->local2Global(mod, x, z, helper);
      geom->relPosToAbsId(mod, x, z, absId);//преобразование в абсолютные координаты
      char relid[3];
      geom->absToRelNumbering(absId, relid);//преобразование в относительные координаты
      cout << "mod: " << mod << " x: " << x << " z: " << z << " E = " << clu.getEnergy() << endl;
      h_energy->Fill(clu.getEnergy());
      //MC info
      //Get list of labels corresponding to this cluster
      gsl::span<const o2::phos::MCLabel> spDigList = mClustersMCArray->getLabels(i);
      cout << "Primaries " << spDigList.size() << endl;
      for (auto &l : spDigList)//выводим свойства кластера (лейблы)
      {
      /*  cout <<
        "Label: empty=" << l.isEmpty() <<
        " noise=" << l.isNoise() <<
        " valid=" << l.isValid() <<
        " fake=" << l.isFake() <<
        " prim=" << l.getTrackID() <<
        " deposited energy " << l.getEdep() << endl;
        */
      }

      cout<<"Monte_Carlo x2: "<<x2<<" z2: "<<z2<<endl;
      cout<<" --- "<<endl;
      if (0<= mod && mod <4)//0<= mod && mod <4
      {
        Double_t NE_max = clu.getNExMax();
        cout<<"############################### NE_max "<<NE_max<<endl;
      //  Radius[mod]->Fill(sqrt(pow((x1-x2), 2)+sqrt(pow((z1-z2), 2))));
        h_NE_max->Fill(sqrt(pow((x1-x2), 2)+pow((z1-z2), 2)), NE_max);
        // prof->Fill(R,locmax)
      }

    }

  }
/*
  cout<<"work1!"<<endl;
  TString name;
  for (int m = 0; m < 4; m++)
  {
    canv[m]->cd();
    cout<<"work2!"<<endl;
    Radius[m]->Draw("colz");//e
    cout<<"work3!"<<endl;
    name.Form("0.5GeV_%d.pdf", m);
    canv[m]->SaveAs(name);
  }
  cout<<"count = "<<count<<endl;
  */
//  TCanvas *cluster_per_event = new TCanvas("Cluster per event","Cluster per event",600,400);
//  TH1D* Cluster_per_event1=Cluster_per_event->ProjectionX();
//  cluster_per_event->cd(0);
//  Cluster_per_event->SetMarkerStyle(21);
  //Cluster_per_event->GetXaxis()->SetTitle("events");
  //Cluster_per_event->GetYaxis()->SetTitle("clusters");
//  Cluster_per_event->Draw();
//  Cluster_per_event->SaveAs("Cluster_per_event_0.5GeV.pdf");

  TCanvas *tEnergy = new TCanvas("Energy","Energy",800,600);
  tEnergy->cd(0);

  h_energy->Draw("pe");

  TCanvas *Nclu = new TCanvas("Number of clusters","Number of clusters",800,600);

  Nclu->cd(0);

  h_NE_max->SetLineColor(kBlack);
  h_NE_max->SetMarkerStyle(23);
  h_NE_max->SetMarkerColor(kBlack);

  h_NE_max->GetXaxis()->SetRangeUser(0., 10.);
  h_NE_max->GetYaxis()->SetRangeUser(0, 4);

  h_NE_max->GetXaxis()->SetTitle("R, cm");
  h_NE_max->GetYaxis()->SetTitle("N_{clusters}");

  h_NE_max->Draw("pe");
  Nclu->SaveAs("Nclu.pdf");
}


void IntFit(TH1D* _h)
{
  TF1  *f1 = new TF1("f1","gaus", 0.35, 0.55);
      f1->SetLineStyle(2);
//  TF1 *Number= new TF1("f2","gaus", 0.35, 0.55);

     f1->SetParLimits(1, 0.41, 0.48);
     f1->SetParLimits(2, 0.01, 0.025);

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

//  cout<<"#############################################Integral: "<< num  <<endl;
//  cout<<"#############################################Integral error: "<< num_err <<endl;

}
