#include <TFile.h>
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "EVGCore/EventRecord.h"
#include "nusystematics/artless/response_helper.hh"
#include <stdio.h>

genie::NtpMCEventRecord * NtpMC_mcrec;
genie::EventRecord * mcrec;

std::vector< TArrayF* >   weightArrays;
std::vector< TArrayF* >   knotArrays;
std::vector< TTree* >     sysTrees;
std::vector <int> resplist;
TArrayF dVals;
TArrayF wVals;
 
// main loop function
void loop(TFile * filename, TFile * outf,  TTree * tree, TTree * gtree, std::string fhicl_filename )
{
  // DUNE reweight getter
  int realparam = 5;
  nusyst::response_helper rh(fhicl_filename);
  mcrec = NULL;
  NtpMC_mcrec = NULL;
  gtree->SetBranchAddress( "genie_record", &NtpMC_mcrec );
  gtree->GetEntry(0);
  mcrec = NtpMC_mcrec->event;
  systtools::event_unit_response_w_cv_t response = rh.GetEventVariationAndCVResponse(*mcrec);
  //DUNE method: Get list of variations, and make CAF branch for each one
  outf->cd();
  std::vector<unsigned int> parIds = rh.GetParameters();
  for( unsigned int i = 0; i < parIds.size(); ++i ) {
    resplist.push_back(i);
    systtools::SystParamHeader head = rh.GetHeader(parIds[i]);
    std::cout << "Parameter " << i << "with name " << head.prettyName.c_str() << " has the ID: " << parIds[i] << std::endl;
    //printf( "Adding reweight branch %u for %s with %lu shifts\n", parIds[i], head.prettyName.c_str(), head.paramVariations.size());
    
    char ttempname[100];
    string tempname; 
    
    sprintf(ttempname,"%s",head.prettyName.c_str());
    tempname = ttempname;
    //std::cout << "Parameter with name: " << tempname << " has CV value of " << head.centralParamValue << std::endl; 
    
    TTree* sysWeightsTree  = new TTree((tempname+"_tree").c_str(), ("Contains dial values and associated weights for the " + tempname + " systematic parameter.").c_str());
    
    sysWeightsTree->Branch((tempname+"_weights").c_str(),&wVals);
    sysWeightsTree->Branch((tempname+"_knots").c_str(),&dVals);
    sysTrees.push_back(sysWeightsTree); 
} 


  int k = 0;
  // Main event loop
  int N = tree->GetEntries();
  TTree * looptree;
  for ( int ii = 0; ii < N; ++ii ) {
	if (ii==0) {std::cout << "Starting Event loop"<< std::endl;}

    tree->GetEntry(ii);
    if( ii % 1 == 0 ) printf( "Event %d of %d...\n", ii, N );
    //if( ii % 100000 == 0 ) {std::cout << "Event " << ii << " of " << N << std::endl;}

    // get GENIE event record
    gtree->GetEntry(ii);
    mcrec = NtpMC_mcrec->event;
    genie::Interaction * in = mcrec->Summary();

    //Add DUNErw weights to the CAF
    systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*mcrec);
    int j = 0;
    std::vector <int> noresp;
    noresp = resplist;
    for( systtools::event_unit_response_w_cv_t::iterator it = resp.begin(); it != resp.end(); ++it ) {  
      char lloopname[100];
      string loopname;
      string treename;
      systtools::SystParamHeader head = rh.GetHeader((*it).pid);
      sprintf(lloopname,"%s",head.prettyName.c_str());
      //std::cout << "Param name = " << head.prettyName.c_str() << " and it.PID = " << (*it).pid << std::endl;
      loopname = lloopname;
      noresp.erase(std::remove(noresp.begin(), noresp.end(), (*it).pid), noresp.end());
      int resp_size = (*it).responses.size();
      looptree = sysTrees[(*it).pid];
      wVals.Set(resp_size);
      dVals.Set(resp_size);
      j = j+1;
      k = k+1;
        for( int i = 0; i < resp_size; ++i ) {
          wVals.SetAt((*it).responses[i]*(*it).CV_response,i);
          dVals.SetAt(head.paramVariations[i],i);
          }
      looptree->Fill();
      
     }
  for(int nr : noresp) {
         
      looptree = sysTrees[nr];
      systtools::SystParamHeader head = rh.GetHeader(nr);
      int resp_size = head.paramVariations.size();
      wVals.Set(resp_size);
      dVals.Set(resp_size);
      j = j+1;
      k = k+1;
      for( int i = 0; i < resp_size; ++i ) {
	    wVals.SetAt(1.0,i);
        dVals.SetAt(head.paramVariations[i],i);
      }
      looptree->Fill();
	  wVals.Reset();
	  dVals.Reset();
    }
	delete mcrec;
  }
 }

int main( int argc, char const *argv[] ) 
{

  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

// get command line options
  std::string infile;
  std::string outfile;
  std::string fhicl_filename;
  std::string outfile2;
  std::string param_num; 
 
  const char* path  = "/vols/dune/ljw20/DUNE_2023_ND_CAFs/"; 
  //const char* path  = "/vols/dune/ljw20/ND_FHC_CAFs_split/"; 
  //const char* outpath  = "/vols/t2k/users/ljw20/data/DUNE_2021/DUNE_2021_new_xsec_weights/"; 
  const char* outpath  = "/vols/dune/ljw20/DUNE_2023_ND_xsec_weights_split/"; 
  //const char* outpath  = "testWeights/."; 
 
  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--infile") ) {
      infile = argv[i+1];
      i += 2;
    } //else if( argv[i] == std::string("--gfile") ) {
      //gfile = argv[i+1];
      //i += 2;} 
      else if( argv[i] == std::string("--param") ) {
      param_num = argv[i+1];
      i += 2; }
      else if( argv[i] == std::string("--fhicl") ) {
      fhicl_filename = argv[i+1];
      i += 2;} 
    else i += 1; // look for next thing
  }

  outfile2 = infile.substr(0, infile.find_last_of('.')) + "_weight_" + param_num + ".root";
  char filepath[200];
  char outfilepath[200];
  sprintf(filepath, "%s/%s", path, infile.c_str());
  sprintf(outfilepath, "%s/%s", outpath, outfile2.c_str());
  std::cout << "Output File: " << outfilepath <<std::endl;
  
  TFile * out = new TFile(outfilepath, "RECREATE" );
  TFile * tf = new TFile(filepath);
  
  TTree * tree = (TTree*) tf->Get( "caf" ); 
  TTree * gtree = (TTree*) tf->Get( "genieEvt" );

  loop(tf, out, tree, gtree, fhicl_filename );

  std::cout << "Writing output ROOT File" << std::endl;
  out->Write();
  out->Close();
}
