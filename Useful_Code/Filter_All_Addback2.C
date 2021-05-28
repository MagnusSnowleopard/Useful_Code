#include <stdio.h>
#include <map>

#include </home/derosa/packages/GRUTinizer/include/GPeak.h>
#include </home/derosa/packages/GRUTinizer/include/GRootCommands.h>
#include </home/derosa/packages/GRUTinizer/include/TChannel.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TTree.h>

// This code displays all the fits for all the detectors (54 graphs)



void add_back_save(){
  
 TFile* f = new TFile("All_file_basic_addback.root","recreate"); 
 TChain* tree = new TChain("tree"); 

// tree->Add("/home/derosa/packages/simple_data_reduction/Data_Reduction/Data/basic_master.root ");
 tree->Add("/home/derosa/packages/simple_data_reduction/Data_Reduction/Data/basic1222-00.root ");

 TH2D* Addback_summary = (TH2D*) add_back_sum(tree);
   
 f->Write(); 

 f->Close(); 

return; 
}













