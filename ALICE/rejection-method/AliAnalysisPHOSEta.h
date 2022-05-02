#ifndef AliAnalysisPHOSEta_cxx
#define AliAnalysisPHOSEta_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for extracting spectra of resonances in PHOS plus tracks
// Authors: Dmitri Peresunko
// 12-July-2018

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class AliESDtrackCuts;
class AliPIDResponse;
class AliAODTrack ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisPHOSEta : public AliAnalysisTaskSE {
public:
	AliAnalysisPHOSEta(const char *name = "AliAnalysisPHOSEta");
        AliAnalysisPHOSEta(const AliAnalysisPHOSEta&); 
	AliAnalysisPHOSEta& operator=(const AliAnalysisPHOSEta&); 
	virtual ~AliAnalysisPHOSEta();

	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
    void SetMC(Bool_t b){isMC=b;}
    void SetMB(Bool_t b){fIsMB=b;}
    
private:

    
	void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
// 	void FillHistogram(Int_t remi,Channels_t ch,PID_t pid,Double_t x) const ; //Fill 1D histogram witn name key
	
	void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
//	void FillHistogram(Int_t remi,Channels_t ch,PID_t pid,Double_t x, Double_t y, Double_t z) const ; //Fill 2D histogram witn name key
	void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
	void SelectGamma() ;
        Bool_t AcceptTrack(const AliAODTrack *t);
        Double_t Pi0Mass(Double_t pt) ;
        Double_t Pi0Width(Double_t pt) ;
        Double_t EtaMass(Double_t pt) ;
        Double_t EtaWidth(Double_t pt) ;
        Double_t PionDispCut(Double_t m02, Double_t m20, Double_t E) ;
       
 
private:
	THashList * fOutputContainer;  //! final histogram container

	
	AliAODEvent  * fEvent ;      //! Current event
	TClonesArray * fGamma;       //! List of selected photons
	TClonesArray *fStack;
    Bool_t isMC;
    Bool_t fIsMB;
    
	ClassDef(AliAnalysisPHOSEta, 1); 
};
#endif
