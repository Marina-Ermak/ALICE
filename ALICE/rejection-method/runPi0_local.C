#include "AliAODInputHandler.h"
#include <cstring>
#include <iostream>
R__LOAD_LIBRARY($ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C)
R__LOAD_LIBRARY(AddTaskPHOSEta.C)// eta or resonances

void runPi0_local(){
    // Load common libraries
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSISaliceBase");
    gSystem->Load("libCORRFW");
    gSystem->Load("libOADB");
    //gSystem->Load("libCORRFW.so");
    // Use AliRoot includes to compile our task
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    
    AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
    
    //Create input handler
    AliAODInputHandler* esdH = new AliAODInputHandler();
    mgr->SetInputEventHandler(esdH);
    
 TChain* chain;
  chain = new TChain("aodTree");
   // chain->Add("~/data/AliAOD.1.root");
   // chain->Add("~/data/AliAOD.2.root");
   // chain->Add("~/data/AliAOD.3.root");
char folder[1000];
//strcpy(folder,"/home/hellaweas/271868/00");
char kostil[1000];
for(int i=1; i<3; i++)
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
    cout<<" Compiling Task..."<<endl;
    // compile our task
//    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_Resonances/AddTaskPHOSResonances.C");
    //Tender Supplies
//gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C");
  //    AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","Run2NoNcellCut",1,true) ;
//  AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","Run2TuneMCNoNcell",1,true) ;
//  AliPHOSTenderSupply* supply = tenderPHOS->GetPHOSTenderSupply() ;
 // supply->ForceUsingBadMap("BadMap.root");
 AliPHOSTenderTask *tenderPHOS = AddAODPHOSTender("PHOSTenderTask","PHOStender","Run2NoNcellCut",1,true) ;

    AliPHOSTenderSupply* supply = tenderPHOS->GetPHOSTenderSupply() ;
    supply->SetNonlinearityVersion("Run2TuneMCNoNcell");
    supply->ForceUsingBadMap("BadMap.root");
    // create our task
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    AliAnalysisPHOSResonances *taskpi0 = AddTaskPHOSEta();//Eta or resonances
    // Enable debug printouts
    mgr->SetDebugLevel(9);
    
    if (!mgr->InitAnalysis())
        return;
    
    mgr->PrintStatus();
    // Start analysis in grid.
    mgr->StartAnalysis("local", chain);

};
