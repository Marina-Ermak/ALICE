#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH(./)
// #include <TAlienCollection.h>
// #include <TAlienResult.h>
#include "AliAODInputHandler.h"
#include <OADB/macros/AddTaskPhysicsSelection.C>
#include <OADB/macros/AddTaskCentrality.C>
#include <ANALYSIS/macros/AddTaskPIDResponse.C>
#include <PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C>
#include "AliAnalysisPHOSEta.h"
R__LOAD_LIBRARY(AliAnalysisPHOSEta_cxx)
#endif
void Tagging(const char* dataset="collection.xml")
{
  gErrorIgnoreLevel=2001 ;

  TChain* chain = new TChain("aodTree");


  // Connect to alien
//  TGrid::Connect("alien://");

//  cout << "Pi0Analysis: processing collection " << dataset << endl;
//  TGridCollection * collection = dynamic_cast<TGridCollection*>(TAlienCollection::Open(dataset));

//  TAlienResult* result = (TAlienResult *)collection->GetGridResult("",0 ,0);
//  TList* rawFileList = result->GetFileInfoList();

/*
  for (Int_t counter=0 ; counter < rawFileList->GetEntries() ; counter++) {
    TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(counter)) ;
    const char * rawFile = fi->GetCurrentUrl()->GetUrl() ;
    printf("Processing %s\n", rawFile) ;
    chain->Add(rawFile);
    printf("Chain: %lld entries.\n",chain->GetEntries());
  }
*/

char folder[1000];
  char kostil[1000];
  for(int i=1; i<183; i++)
  {
      memset(folder, 0, sizeof(folder));
      strcpy(folder,"/home/hellaweas/271868/00");
      cout<<"Folder after strcpy "<<folder<<endl;
      if (i<10)
      {
          strcat(folder,"0");
          sprintf(kostil,"%d",i);
            cout<<"kostil is "<<kostil<<endl;
            strcat(folder, kostil);
             cout<<"After kostil "<<folder<<endl;
             strcat(folder, "/AliAOD.root");
            cout<<"after second strcat Folder is "<<folder<<endl;
             chain->Add(folder);

      }
      else if(i<=99)
      {
          sprintf(kostil,"%d",i);

           cout<<"kostil is "<<kostil<<endl;
           strcat(folder, kostil);
           cout<<"After kostil "<<folder<<endl;
           strcat(folder, "/AliAOD.root");
           cout<<"after second strcat Folder is "<<folder<<endl;
           chain->Add(folder);



      }
      else{
          strcpy(folder,"/home/hellaweas/271868/0");
          sprintf(kostil,"%d",i);
          cout<<"kostil is "<<kostil<<endl;
          strcat(folder, kostil);
          cout<<"After kostil "<<folder<<endl;
          strcat(folder, "/AliAOD.root");
          cout<<"after second strcat Folder is "<<folder<<endl;
          chain->Add(folder);
     }
   }

      //chain->Add("/home/hellaweas/alice/AliPhysics/PWGGA/PHOSTasks/PHOS_Resonances/AliAOD.root");
      cout<<" Compiling Task..."<<endl;


//chain->Add("../data/0002/AliAOD.root");// old part

//     chain->Add("/home/prsnko/Tagging/pPb8TeV/MC/Code/AliAOD.root");

//    chain->Add("/home/prsnko/Calibr/Run2/CalibPi02/Alignment/Code/AliAOD.root");

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Spectrum");

  // ESD input handler
  AliAODInputHandler* esdH = new AliAODInputHandler();
//  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

  // Debug level
//  mgr->SetDebugLevel(2);

//  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C") ;
  AddTaskPhysicsSelection(false,true) ;


  AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","Run2NoNcellCut",1,true) ;
  AliPHOSTenderSupply* supply = tenderPHOS->GetPHOSTenderSupply() ;
  supply->SetNonlinearityVersion("Run2TuneMCNoNcell");
  supply->ForceUsingBadMap("BadMap.root");



  // Add my task
  AliAnalysisPHOSEta* task1 = new AliAnalysisPHOSEta("PHOSEta");
  task1->SetMC(kTRUE);///1
  mgr->AddTask(task1);

  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer();

//   Create containers for input/output
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("PHOSEta",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");


  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);


  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }

}
