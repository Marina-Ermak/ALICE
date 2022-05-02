/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Analysis task for resonanse measurement in PHOS
// Authors: Dmitri Peresunko
#include <iostream>
#include <cstdio>
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "THashList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisPHOSEta.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTriggerAnalysis.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"

// const Double_t massPi0=0.1349770;

ClassImp(AliAnalysisPHOSEta) ;

//________________________________________________________________________
AliAnalysisPHOSEta::AliAnalysisPHOSEta(const char *name)
: AliAnalysisTaskSE(name),
fOutputContainer(0),
fGamma(0x0),
fStack(0x0),
isMC(0)

//fTriggerAnalysis(new AliTriggerAnalysis),

{
    // Output slots #0 write into a TH1 container
    DefineOutput(1,THashList::Class());
}

AliAnalysisPHOSEta:: AliAnalysisPHOSEta(const AliAnalysisPHOSEta& rh):
AliAnalysisTaskSE(rh.GetName()),
fOutputContainer(0x0),
fGamma(0x0),
fStack(0x0),
isMC(0)
{
    if(fOutputContainer)
        delete fOutputContainer;
    fOutputContainer = new THashList();
}

//________________________________________________________________________
AliAnalysisPHOSEta& AliAnalysisPHOSEta::operator=(const AliAnalysisPHOSEta& rh){
    // assignment operator

    this->~AliAnalysisPHOSEta();
    new(this) AliAnalysisPHOSEta(rh);
    return *this;
}
//________________________________________________________________________
AliAnalysisPHOSEta::~AliAnalysisPHOSEta()
{
    //Destructor
    if(fOutputContainer)
    {
        delete fOutputContainer;
        fOutputContainer=0x0 ;
    }
    //No need to delete histograms in array fhHistos[]!
    //They are deleted as content of fOutputContainer
}
//________________________________________________________________________
void AliAnalysisPHOSEta::UserCreateOutputObjects()
{
    // Create histograms
    // Called once

    // AOD histograms
    if(fOutputContainer != NULL)
    {
        delete fOutputContainer;
    };
    fOutputContainer = new THashList();
    fOutputContainer->SetOwner(kTRUE);
/*
    Int_t nPt=73;
    Double_t ptBins[74]={0.,0.1,0.2,0.3,0.4, 0.5,0.6,0.7,0.8,0.9, 1.0,1.1,1.2,1.3,1.4, 1.5,1.6,1.7,1.8,1.9, 2.0,2.2,2.4,2.5,2.6,
        2.8,3.,3.2,3.4,3.6, 3.8,4.0,4.5,4.8,5., 5.5,5.6,6.0,6.4,6.5, 7.0,7.2,7.5,8.0,8.5, 9.0,9.5,10.,11.,12., 13.,14.,15.,16.,17., 18.,19.,20.,22.,24., 25.,26.,28.,30.,32.,
        35.,40.,45.,50.,55., 60.,65.,70.,80.};
*/
Int_t nPt=4;
double ptBins[5]={1.5,3.,4.5,6.,7.5};
/*ptBins[0]=0.5;// {0.5, 0.7, 0.9, 1.2, 1.5, 2, 2.5, 3.2, 4., 5., 7., 10.}
 for (int i=1; i<11; i++)
     {
         if(i<9)
          ptBins[i]=ptBins[i-1]+1;
         else
          ptBins[i]=ptBins[i-1]+3;
         cout<<" pTBins ["<<i<<"] = "<<ptBins[i]<<endl;

     }

    Int_t nPt2=80;
    Double_t ptBins2[81]={0.,0.1,0.2,0.3,0.4, 0.5,0.6,0.7,0.8,0.9, 1.0,1.1,1.2,1.3,1.4, 1.5,1.6,1.7,1.8,1.9, 2.0,2.2,2.4,2.5,2.6,
        2.8,3.,3.2,3.4,3.6, 3.8,4.0,4.5,4.8,5., 5.5,5.6,6.0,6.4,6.5, 7.0,7.2,7.5,8.0,8.5, 9.0,9.5,10.,11.,12., 13.,14.,15.,16.,17., 18.,19.,20.,22.,24., 25.,26.,28.,30.,32.,
        35.,40.,45.,50.,55., 60.,65.,70.,80.,90.,100.,110.,120.,130.,140.,150.};

    fOutputContainer->Add(new TH2F("InvMass_All_Rejection_E100","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_All_Rejection_E300","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_Both_Rejection_E100","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_Both_Rejection_E300","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;

    //   fOutputContainer->Add(new TH2F("1_sigma","Invariant mass", 400,0., 1.,50,0,5.)) ;
    // fOutputContainer->Add(new TH2F("2_sigma","Invariant mass", 400,0., 1.,50,0,5.));
    fOutputContainer->Add(new TH2F("InvMass_All_E100","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_All_E300","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_Disp_E100","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_Disp_E300","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_CPV_E100","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_CPV_E300","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_Both_E100","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
    fOutputContainer->Add(new TH2F("InvMass_Both_E300","Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;

    fOutputContainer->Add(new TH2F("True_eta_Inv_mass","Invariant mass", 1500,0., 1.5,nPt,ptBins));
    fOutputContainer->Add(new TH2F("True_pi0_Inv_mass","Invariant mass", 1500,0., 1.5,nPt,ptBins));
    fOutputContainer->Add(new TH2F("True_etaprime_Inv_mass","Invariant mass", 1500,0., 1.5,nPt,ptBins));
    fOutputContainer->Add(new TH2F("True_eta_Inv_mass_Rejection","Invariant mass", 1500,0., 1.5,nPt,ptBins));

    fOutputContainer->Add(new TH1F("EtaGeneratedSpectra", "EtaSpectraGenerated", nPt,ptBins));
    fOutputContainer->Add(new TH1F("EtaGeneratedSpectra_PHOS", "EtaSpectraGenerated PHOS", nPt,ptBins));
    fOutputContainer->Add(new TH1F("Pi0GeneratedSpectra", "Pi0SpectraGenerated", nPt,ptBins));
    fOutputContainer->Add(new TH1F("Pi0GeneratedSpectra_PHOS", "Pi0SpectraGenerated PHOS", nPt, ptBins));


    for(int i=0;i<5;i++){
        fOutputContainer->Add(new TH2F(Form("InvMass_All_Rejection_E100_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
        fOutputContainer->Add(new TH2F(Form("InvMass_All_Rejection_E300_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
        fOutputContainer->Add(new TH2F(Form("InvMass_Both_Rejection_E100_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
        fOutputContainer->Add(new TH2F(Form("InvMass_Both_Rejection_E300_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;


        fOutputContainer->Add(new TH2F(Form("InvMass_All_E100_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
        fOutputContainer->Add(new TH2F(Form("InvMass_All_E300_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
        fOutputContainer->Add(new TH2F(Form("InvMass_Both_E100_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;
        fOutputContainer->Add(new TH2F(Form("InvMass_Both_E300_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins)) ;

        fOutputContainer->Add(new TH2F(Form("True_eta_Inv_mass_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins));
        fOutputContainer->Add(new TH2F(Form("True_pi0_Inv_mass_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins));
        fOutputContainer->Add(new TH2F(Form("True_eta_Inv_mass_Rejection_M%d",i),"Invariant mass", 1500,0., 1.5,nPt,ptBins));
    }

//    fOutputContainer->Add(new TH3F("Alpha","Phi12(Alpha)", 50, 0, TMath::Pi(), 50,0,1., nPt,ptBins));
//    fOutputContainer->Add(new TH3F("Alpha_sel_clu","Phi12(Alpha) of pi0", 50, 0, TMath::Pi(), 50,0,1., nPt,ptBins));

    fOutputContainer->Add(new TH1F("hSelEvents","Events selected",10,0.,10.));
    fOutputContainer->Add(new TH1F("hCellMultEvent", "PHOS cell multiplicity per event",2000,0,2000));
    fOutputContainer->Add(new TH1F("hClusterMult", "CaloCluster multiplicity", 100,0,100));
    fOutputContainer->Add(new TH1F("hPHOSClusterMult","PHOS cluster multiplicity",100,0,100));
    fOutputContainer->Add(new TH1F("hCellEnergy", "Cell energy", 500,0.,50.));
    fOutputContainer->Add(new TH1F("hClusterEnergy", "Cluster energy", 500,0.,50.));
    fOutputContainer->Add(new TH2F("hClusterEvsN", "Cluster energy vs digit multiplicity",     500,0.,50.,40,0.,40.));
    fOutputContainer->Add(new TH1F("hCellMultClu","Cell multiplicity per cluster",200,0,200));
    fOutputContainer->Add(new TH1F("hModule","Module events",5,0.,5.));
    fOutputContainer->Add(new TH2F("hClusterTOFvsE", "Cluster time vs energy", 500,-250.e-9,250.e-9,40,0.,40.));
    for(Int_t module=1; module<5; module++){
        fOutputContainer->Add(new TH2F(Form("hModule%d",module), Form("Cluster occupancy in module %d",module), 64,0.,64,56,0.,56.));
    }

    fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
    fOutputContainer->Add(new TH1F("hNvertexTracks",   "N of primary tracks from the primary vertex",150,0.,150.));
    fOutputContainer->Add(new TH1F("hT0TOF","T0 time (s)",2000,-1.0e-7,1.0e-7));
    fOutputContainer->Add(new TH1F("hTrackMult" ,"Charged track multiplicity",150,0.,150.));

    Int_t n=100;
    Double_t bins[100];
    for(int i=0;i<n;i++) bins[i]=10.*(double)i/(double)n;

    fOutputContainer->Add(new TH3F("hDispE", "Cluster dispersion vs energy", n-1,bins,n-1,bins,nPt2, ptBins2));

    PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisPHOSEta::UserExec(Option_t *)
{
    // Main loop, called for each event
    // Analyze AOD

    FillHistogram("hSelEvents",1) ;
    fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fEvent)
    {
        Printf("ERROR: Could not retrieve event");
        return;
    }
    FillHistogram("hSelEvents",2);

    //	AliAODInputHandler *esdH = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    // Checks if we have a primary vertex
    // Get primary vertices form AOD
    /*
     *	if      (fEvent->GetPrimaryVertexTracks()->GetNContributors()>0)
     *		fEventVtxExist    = kTRUE;
     *	else if (fEvent->GetPrimaryVertexSPD()   ->GetNContributors()>0)
     *		fEventVtxExist    = kTRUE;
     */

    if(!isMC){
        TString trigClasses = ((AliAODHeader*)fEvent->GetHeader())->GetFiredTriggerClasses();

        Bool_t isMB = (fInputHandler->IsEventSelected() & AliVEvent::kINT7)  ;
        Bool_t isPHI7 = (fInputHandler->IsEventSelected() & AliVEvent::kPHI7);


        if((fIsMB && !isMB) || (!fIsMB && !isPHI7)){
                PostData(1, fOutputContainer);
                return;
        }
    }

    FillHistogram("hSelEvents",3) ;

    const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();

    Double_t vtx5[3];
    vtx5[0] = esdVertex5->GetX();
    vtx5[1] = esdVertex5->GetY();
    vtx5[2] = esdVertex5->GetZ();
    TVector3 vtx(vtx5) ;

    FillHistogram("hNvertexTracks",esdVertex5->GetNContributors());
    FillHistogram("hZvertex"      ,esdVertex5->GetZ());
    if (TMath::Abs(esdVertex5->GetZ()) > 10. )
        return ;

    FillHistogram("hSelEvents",4) ;

    //	if (fEvent->IsPileupFromSPD())
    //          return ;

    //	FillHistogram("hSelEvents",5) ;

    // Fill event statistics for different selection criteria

    //Vtx class z-bin
    Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
    if(zvtx<0)zvtx=0 ;
    if(zvtx>9)zvtx=9 ;

    SelectGamma() ;

    // const Double_t massPi0=0.1349770; лучше потом раскомментить

    //Fill gamma-gamma and select pi0s
    Int_t nGamma = 0;
    if(fGamma) nGamma = fGamma->GetEntriesFast() ;


    for(Int_t i=0; i<nGamma; i++){
        AliCaloPhoton * pv1 = (AliCaloPhoton*)fGamma->At(i) ;

        FillHistogram("hClusterEnergy",pv1->E());
    }
    //fOutputContainer->Add(new TH2F("Alpha","Phi12(Alpha)", 400, 0, 1, 50,0,5.));
    // Post output data.
    PostData(1, fOutputContainer);
}

//________________________________________________________________________
void AliAnalysisPHOSEta::Terminate(Option_t *)
{
}
//_____________________________________________________________________________
void AliAnalysisPHOSEta::FillHistogram(const char * key,Double_t x)const{
    //FillHistogram
    TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
    if(hist)
        hist->Fill(x) ;
    else
        AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________fOutputContainer->Add(new TH2F("Alpha","Phi12(Alpha)", 400, 0, 1, 50,0,5.));


void AliAnalysisPHOSEta::FillHistogram(const char * key,Double_t x,Double_t y)const{
    //FillHistogram
    TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
    if(th1)
        th1->Fill(x, y) ;
    else
        AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}


//_____________________________________________________________________________
void AliAnalysisPHOSEta::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
    //Fills 1D histograms with key
    TObject * obj = fOutputContainer->FindObject(key);

    TH2 * th2 = dynamic_cast<TH2*> (obj);
    if(th2) {
        th2->Fill(x, y, z) ;
        return;
    }

    TH3 * th3 = dynamic_cast<TH3*> (obj);
    if(th3) {
        th3->Fill(x, y, z) ;
        return;
    }

    AliError(Form("can not findi histogram (of instance TH2) <%s> ",key)) ;
}



//_____________________________________________________________________________
double pi0_width(Double_t x)
{
    return 0.0006358*TMath::Exp(-0.1327*x+0.3698)-0.004659+0.01161*TMath::Power(x, -0.2168);
    // return 159600*TMath::Exp(-1.176*x)
    //  return  0.006434*TMath::Exp(-1.096*x)+0.005637;
    //    return 0.01729*TMath::Exp(-0.8375*x)+0.001515*x+(-2.538*TMath::Power(10, -5))*x*x*x ;
}

void AliAnalysisPHOSEta::SelectGamma(){
    Int_t inPHOS=0/*, iPi0Merged=0*/;
    const Double_t massPi0=0.1349770;
    if(fGamma)
        fGamma->Clear() ;
    else
        fGamma = new TClonesArray("AliCaloPhoton",100) ;

    const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();

    Double_t vtx5[3] ={esdVertex5->GetX(),esdVertex5->GetY(),esdVertex5->GetZ()};
    //const Double_t rcut = 1;

    Int_t multClust = fEvent->GetNumberOfCaloClusters();
    if(isMC){
        fStack = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());//
    }

    if (isMC && fStack){

      Int_t multMC = fStack->GetEntriesFast();
      for (Int_t iMC = 0; iMC < multMC; iMC++){
        AliAODMCParticle *particle = (AliAODMCParticle*)fStack->At(iMC);
        if (particle->GetPdgCode() == 221 && TMath::Abs(particle->Y()) < 0.5){


            if (abs(particle->Y()) < 0.5)
              FillHistogram("EtaGeneratedSpectra", particle->Pt());

            if (TMath::Abs(particle->Phi()) < 0.611 && TMath::Abs(particle->Eta()) < 0.125){
              FillHistogram("EtaGeneratedSpectra_PHOS", particle->Pt());
            }

        }

        if (particle->GetPdgCode() == 111 && TMath::Abs(particle->Y()) < 0.5){


            if (abs(particle->Y()) < 0.5)
              FillHistogram("Pi0GeneratedSpectra", particle->Pt());

            if (TMath::Abs(particle->Phi()) < 0.611 && TMath::Abs(particle->Eta()) < 0.125){
              FillHistogram("Pi0GeneratedSpectra_PHOS", particle->Pt());
            }

        }
        }
      }

    Int_t multBin = 0;
    if(multClust<5) multBin = 0;
    if(multClust>=5 && multClust<10) multBin = 1;
    if(multClust>=10 && multClust<20) multBin = 2;
    if(multClust>=20 && multClust<30) multBin = 3;
    if(multClust>=30) multBin = 4;


    for (Int_t i=0; i<multClust; i++)
    {
        AliAODCaloCluster *clu = fEvent->GetCaloCluster(i);
        if (clu->GetType() !=AliVCluster::kPHOSNeutral)
            continue;


        if(clu->E()<0.1)
            continue;

        if(clu->E()>2. && clu->GetNCells()<3)
            continue ;

        if(clu->E()>2. && clu->GetM02()<0.1)
            continue ;

        FillHistogram("hClusterTOFvsE", clu->GetTOF(),clu->E());

        if((clu->GetTOF()>30.e-9) || (clu->GetTOF() <-30.e-9) )
            continue;

        FillHistogram("hDispE", clu->GetM02(), clu->GetM20(), clu->E());

        TLorentzVector pv1;
        clu->GetMomentum(pv1, vtx5);
        if(inPHOS >= fGamma->GetSize())
        {
            fGamma->Expand(inPHOS+50);
        }

        TLorentzVector momentum ;
        clu->GetMomentum(momentum, vtx5);

/*        Double_t cluE = clu->E()*1.028*(1.- 0.018*TMath::Exp(-clu->E()*0.5)); //non-linearity

        Double_t sc = cluE/clu->E();
        AliCaloPhoton *p = new ((*fGamma)[inPHOS++]) AliCaloPhoton(sc*pv1.Px(),sc*pv1.Py(),sc*pv1.Pz(),sc*clu->E() );
*/

        AliCaloPhoton * p = new((*fGamma)[inPHOS++]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E());

        //    p->SetDispBit(clu->Chi2()<2.5*2.5);
        p->SetDispBit(clu->GetDispersion()<2.5*2.5) ;

        //    p->SetTOFBit((clu->GetTOF()>-50.e-9) && (clu->GetTOF() <50.e-9) );
        p->SetCPVBit(clu->GetEmcCpvDistance()>2.5);
        p->SetCPV2Bit(clu->GetEmcCpvDistance()>1.);
        //p->SetPrimary(0) ; //no matched partner yet
        int num_it=0;//!
        if(isMC){
            Int_t primLabel=clu->GetLabelAt(num_it);//!!
            while(primLabel>-1) //previously if
            {
                AliAODMCParticle *prim = (AliAODMCParticle*)fStack->At(primLabel);
                Int_t iparent = primLabel;
                //        printf ("Code of the daughter particle %d\n", prim->GetPdgCode());
                AliAODMCParticle *parent = prim;
                iparent = parent->GetMother();
                if(iparent>-1){
                    parent = (AliAODMCParticle*)fStack->At(iparent);
                    if(parent->GetPdgCode()==221 || parent->GetPdgCode()==111 || parent->GetPdgCode()==331){
                        p->SetPrimaryAtVertex(iparent);
                        break;//!
                    }
                    else{
                        p->SetPrimaryAtVertex(primLabel);
                    }

                    //        printf ("Code of the oldest parent particle: check %d\n", parent->GetPdgCode());
                }

                p->SetPrimary(primLabel);
                primLabel=clu->GetLabelAt(num_it++);
            }
        }
    }


    for (int i=0; i<inPHOS-1; i++){
        //     AliAODCaloCluster *clu1 = fEvent->GetCaloCluster(i);
        AliCaloPhoton *pv1 = static_cast<AliCaloPhoton*>(fGamma->At(i));
        //AliCaloPhoton *pv1 = (AliCaloPhoton*)(fGamma->At(i));
        //            TLorentzVector pv1;
        //            clu1->GetMomentum(pv1, vtx5);

        for(int j=i+1; j<inPHOS; j++){

            //     AliAODCaloCluster *clu2 = fEvent->GetCaloCluster(j);
            //             TLorentzVector pv2;
            //             clu2->GetMomentum(pv2, vtx5);
            AliCaloPhoton *pv2 = static_cast<AliCaloPhoton*>(fGamma->At(j));

            Double_t alpha=TMath::Abs(pv1->E()-pv2->E())/(pv1->E()+pv2->E());
            Int_t pos1 = pv1->GetPrimaryAtVertex();
            Int_t pos2 = pv2->GetPrimaryAtVertex();
            //	 AliAODMCParticle *phot1 = (AliAODMCParticle*)(fStack->At(pos1));
            //	 AliAODMCParticle *phot2 = (AliAODMCParticle*)(fStack->At(pos2));
            //	 Int_t Pdg1 = phot1->GetPdgCode();
            //	 Int_t Pdg2 = phot2->GetPdgCode();
            if(isMC){
                if (pos1 == pos2 /* && (*pv1+*pv2).M()>0.1*/) {
                    AliAODMCParticle *prim = (AliAODMCParticle*)fStack->At(pos1);
                    if (prim->GetPdgCode()==221)
                    {
                        FillHistogram("True_eta_Inv_mass", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                        FillHistogram(Form("True_eta_Inv_mass_M%d",multBin), (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                    }
                    else if (prim->GetPdgCode()==111)
                    {
                        FillHistogram("True_pi0_Inv_mass", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                        FillHistogram(Form("True_pi0_Inv_mass_M%d",multBin), (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                    }
                    else if (prim->GetPdgCode()==331)
                    {
                        FillHistogram("True_etaprime_Inv_mass", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                        FillHistogram(Form("True_etaprime_Inv_mass_M%d",multBin), (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                    }
                    //		 printf ("If-condition");
                    FillHistogram("Alpha_sel_clu",TMath::ACos(1-massPi0*massPi0/(2*(pv1->E())*(pv2->E()))), alpha, pv1->E());
                    FillHistogram("Alpha_sel_clu",TMath::ACos(1-massPi0*massPi0/(2*(pv1->E())*(pv2->E()))), alpha, pv2->E());
                }
            }

            FillHistogram("Alpha",TMath::ACos(1-massPi0*massPi0/(2*(pv1->E())*(pv2->E()))), alpha, pv1->E());
            FillHistogram("Alpha",TMath::ACos(1-massPi0*massPi0/(2*(pv1->E())*(pv2->E()))), alpha, pv2->E());

            FillHistogram("InvMass_All_E100", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
            FillHistogram(Form("InvMass_All_E100_M%d",multBin), (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
            if(pv1->E()>0.3 && pv2->E()>0.3){
                FillHistogram("InvMass_All_E300", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                FillHistogram(Form("InvMass_All_E300_M%d",multBin), (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
            }
            if(pv1->IsDispOK() && pv2->IsDispOK()){
                FillHistogram("InvMass_Disp_E100", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                if(pv1->E()>0.3 && pv2->E()>0.3)
                    FillHistogram("InvMass_Disp_E300", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                if(pv1->IsCPVOK() && pv2->IsCPVOK()){
                    FillHistogram("InvMass_Both_E100", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                    FillHistogram(Form("InvMass_Both_E100_M%d",multBin), (*pv1+*pv2).M(), (*pv1+*pv2).Pt());

                    if(pv1->E()>0.3 && pv2->E()>0.3){
                        FillHistogram("InvMass_Both_E300", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                        FillHistogram(Form("InvMass_Both_E300_M%d",multBin), (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                    }

                }
            }

            if(pv1->IsCPVOK() && pv2->IsCPVOK()){
                FillHistogram("InvMass_CPV_E100", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());
                if(pv1->E()>0.3 && pv2->E()>0.3)
                    FillHistogram("InvMass_CPV_E300", (*pv1+*pv2).M(), (*pv1+*pv2).Pt());

            }

        }
    }
    //   TF1 *pi0_width = new TF1("pi0_width","[3]*TMath::Exp(-[0]*x)+[1]*x+[2]*x*x*x", 0.5, 5);
    int pointed_clusters[inPHOS*inPHOS] ;
    int k=0 ;
    Double_t M_pt[2];
    for (int i=0; i<inPHOS-1; i++)
        for(int j=i+1; j<inPHOS; j++)
        {
            AliCaloPhoton *photon1 = static_cast<AliCaloPhoton*>(fGamma->At(i)) ;
            AliCaloPhoton *photon2 = static_cast<AliCaloPhoton*>(fGamma->At(j)) ;
            M_pt[0]=(*photon1+*photon2).M() ;//mass
            M_pt[1]=(*photon1+*photon2).Pt() ;//momentum
            if (M_pt[0]>(massPi0-2*pi0_width(M_pt[1])) && M_pt[0]<(massPi0+2*pi0_width((M_pt[1]))))
            {
                //                            printf("PI0_WIDTH - %lf \n", pi0_width(M_pt[1])) ;
                pointed_clusters[k]=i ;
                k++ ;
                pointed_clusters[k]=j ;
                k++ ;
            }
        }
        int flag=0 ;
        for (int i=0; i<inPHOS-1; i++)
            for (int j=i+1; j<inPHOS; j++)
            {
                AliCaloPhoton *photon1_ = static_cast<AliCaloPhoton*>(fGamma->At(i)) ;
                AliCaloPhoton *photon2_ = static_cast<AliCaloPhoton*>(fGamma->At(j)) ;
                for (int l=0; l<k; l++)
                    if(i==pointed_clusters[l] || j==pointed_clusters[l])
                    {
                        flag=1 ;
                    }

                    if (flag!=1)
                    {


                        if(isMC){
                            Int_t pos1 = photon1_->GetPrimaryAtVertex();
                            Int_t pos2 = photon2_->GetPrimaryAtVertex();
                            if (pos1 == pos2 /* && (*pv1+*pv2).M()>0.1*/) {
                                AliAODMCParticle *prim = (AliAODMCParticle*)fStack->At(pos1);
                                if (prim->GetPdgCode()==221)
                                {
                                    FillHistogram("True_eta_Inv_mass_Rejection", (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt());
                                    FillHistogram(Form("True_eta_Inv_mass_Rejection_M%d",multBin), (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt());
                                }
                            }
                        }
                        if(photon1_->IsCPVOK() && photon2_->IsCPVOK() && photon1_->IsDispOK() && photon2_->IsDispOK()){
                            FillHistogram("InvMass_Both_Rejection_E100", (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                            FillHistogram(Form("InvMass_Both_Rejection_E100_M%d",multBin), (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                            if(photon1_->E()>0.3 && photon2_->E()>0.3){
                                FillHistogram("InvMass_Both_Rejection_E300", (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                                FillHistogram(Form("InvMass_Both_Rejection_E300_M%d",multBin), (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                            }
                        }
                        FillHistogram("InvMass_All_Rejection_E100", (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                        FillHistogram(Form("InvMass_All_Rejection_E100_M%d",multBin), (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                        if(photon1_->E()>0.3 && photon2_->E()>0.3){
                            FillHistogram("InvMass_All_Rejection_E300", (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                            FillHistogram(Form("InvMass_All_Rejection_E300_M%d",multBin), (*photon1_+*photon2_).M(), (*photon1_+*photon2_).Pt()) ;
                        }
                    }
                    flag=0 ;
            }
}
//_____________________________________________________________________________
Double_t AliAnalysisPHOSEta::Pi0Mass (Double_t /*pt*/){
    return 0.137 ;
}
//_____________________________________________________________________________
Double_t AliAnalysisPHOSEta::Pi0Width (Double_t /*pt*/){
    return 0.012;   //2sigma
}
//_____________________________________________________________________________
Double_t AliAnalysisPHOSEta::EtaMass (Double_t /*pt*/){
    return 0.555 ;
}
//_____________________________________________________________________________
Double_t AliAnalysisPHOSEta::EtaWidth (Double_t /*pt*/){
    return 0.030;   //2sigma
}
//______________________________________________________
Double_t AliAnalysisPHOSEta::PionDispCut(Double_t m02, Double_t m20, Double_t E){
    //Returns ditance to pi0 peak center in sigmas

    //No Disp cut for soft energies
    if(E<25.)
        return 999;

    //Parameterization using single pi0 simulation
    Double_t longMpi = 1.857398e+00 + 1.208331e+01*TMath::Exp(-4.977723e-02*E) ;
    Double_t longSpi = 3.820707e-01 + 1.000542e+00*TMath::Exp(-3.877147e-02*E) ;
    Double_t shortMpi = 1.152118e+00 - 4.076138e-01*TMath::Exp(-2.372902e-02*E) ;
    Double_t shortSpi = 1.517538e-01 + 9.382205e+00*TMath::Exp(-1.563037e-01*E) ;
    Double_t powerNpi = 2.055773e+00 + 9.616408e+03*TMath::Exp(-2.664167e-01*E) ;

    Double_t dx = (m02-longMpi)/longSpi ;
    Double_t dy = (m20-shortMpi)/shortSpi ;

    //we have non-gaussian power, so re-calculate in Gaussian sigmas
    return TMath::Sign(TMath::Sqrt(TMath::Power(TMath::Abs(dx),powerNpi)+TMath::Power(TMath::Abs(dy),powerNpi)),dx) ;


}
