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
#include "DataFormatsPHOS/Cell.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DataFormatsPHOS/MCLabel.h"
#include "PHOSBase/Geometry.h"
#endif

using namespace std;

void plot_cells_phos(int ievent = -1, TString inputfile = "phoscells.root")
{

    TFile* file1 = TFile::Open(inputfile.Data());
    TTree* cluTree = (TTree*)gFile->Get("o2sim");
    // o2::dataformats::MCTruthContainer<o2::phos::MCLabel>* mCellsMCArray = nullptr;
    // cluTree->SetBranchAddress("PHOSCellTrueMC", &mCellsMCArray);
    std::vector<o2::phos::Cell>* mCellArray = nullptr;
    cluTree->SetBranchAddress("PHOSCell", &mCellArray);
    std::vector<o2::phos::TriggerRecord>* mCellTR = nullptr;
    cluTree->SetBranchAddress("PHOSCellTrigRec", &mCellTR);


    if (!mCellArray) {
        cout << "PHOS cells not in tree. Exiting ..." << endl;
        return;
    }

    TH2D* vMod[4];
    TH1D* eMod[4];
    TCanvas* c[4];
    for(int m=0;m<4;m++){
        vMod[m]=new TH2D(Form("hMod%d", m+1), Form("hMod%d", m+1), 64, 0., 64., 56, 0., 56.);
        eMod[m]=new TH1D(Form("eMod%d", m+1), Form("eMod%d", m+1), 1000, 0, 1000.);
        c[m] = new TCanvas(Form("cellsInMod%d", m+1), Form("PHOS cells in module %d", m+1), 800,600);
    }

    cout<<cluTree->GetEntries()<<endl;
    for(int iev=0; iev<cluTree->GetEntries(); iev++){
    cluTree->GetEvent(iev);
    cout << "Read "<<mCellTR->size() << " events" << endl ;
    // cout << "Lables " << mCellsMCArray->getNElements() << endl;


     for(auto it = mCellTR->begin(); it != mCellTR->end(); it++) {

      int firstCellInEvent = it->getFirstEntry();
      int mLastCellInEvent = firstCellInEvent + it->getNumberOfObjects();
      cout << "event:" << it->getBCData() << " from " << firstCellInEvent << " to " << mLastCellInEvent << endl ;
    	for(int i=firstCellInEvent;i<mLastCellInEvent; i++ ){
    	  auto &clu = mCellArray->at(i)	;
        short absId = clu.getAbsId() ;
        char relid[3];
        o2::phos::Geometry::absToRelNumbering(absId, relid);//// Converts the absolute numbering into the following array
  //  relid[0] = PHOS Module number 1:module
  //  relid[1] = Row number inside a PHOS module (Phi coordinate)
  //  relid[2] = Column number inside a PHOS module (Z coordinate)
        int mod=relid[0]-1;
        float e = clu.getEnergy();
        if(clu.getLowGain())e*=16;
        float t = clu.getTime();
        float mOccCut=0.;

          if(mod>=0&&mod<4){
            if (e > mOccCut) {
                //                mHist2D[kOccupancyM1 + mod]->Fill(relid[1] - 0.5, relid[2] - 0.5);
                vMod[mod]->Fill(relid[1] - 0.5, relid[2] - 0.5);
                eMod[mod]->Fill(e);
                //cout<<"eneregy - "<<e<<endl;
            }
          }
        }
     }
    }

    for(int m=0;m<4;m++){
        c[m]->Divide(2,1);
        c[m]->cd(1);
        vMod[m]->Draw("colz");
        c[m]->cd(2);
        eMod[m]->SetMarkerStyle(21);
        eMod[m]->GetXaxis()->SetTitle("Energy");
        eMod[m]->Draw("e");
        TString number;
        number=Form("Mod%d", m + 1);
        c[m]->SaveAs("cells_5GeV_mod"+number+".pdf");
    }

}
