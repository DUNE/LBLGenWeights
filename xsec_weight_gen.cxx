#include <TFile.h>
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "EVGCore/EventRecord.h"
#include "nusystematics/artless/response_helper.hh"
#include <stdio.h>

genie::EventRecord * mcrec;

// main loop function
void loop(TFile * filename, TFile * outf,  TTree * tree, TTree * gtree, std::string fhicl_filename )
{
  // DUNE reweight getter
  nusyst::response_helper rh(fhicl_filename);
  
  //DUNE method: Get list of variations, and make CAF branch for each one
  outf->cd();
  std::vector<unsigned int> parIds = rh.GetParameters();
  for( unsigned int i = 0; i < parIds.size(); ++i ) {
    systtools::SystParamHeader head = rh.GetHeader(parIds[i]);
    printf( "Adding reweight branch %u for %s with %lu shifts\n", parIds[i], head.prettyName.c_str(), head.paramVariations );
    char ttempname[100];
    string tempname; 
    sprintf(ttempname,"%s",head.prettyName.c_str());
    tempname = ttempname;
    TTree* sysWeightsTree  = new TTree((tempname+"_tree").c_str(), ("Contains dial values and associated weights for the " + tempname + " systematic parameter.").c_str());
} 


  mcrec = NULL;
  gtree->SetBranchAddress( "genie_record", &mcrec );
  int k = 0;
  // Main event loop
  int N = tree->GetEntries();
  TTree * looptree;
  for ( int ii = 0; ii < N; ++ii ) {
   // if (ii > 9) {
     // break; }
    tree->GetEntry(ii);
    if( ii % 100 == 0 ) printf( "Event %d of %d...\n", ii, N );

    // get GENIE event record
    gtree->GetEntry(ii);
    genie::Interaction * in = mcrec->Summary();

    //Add DUNErw weights to the CAF
    systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*mcrec);
    int j = 0;
    for( systtools::event_unit_response_w_cv_t::iterator it = resp.begin(); it != resp.end(); ++it ) {  
      char lloopname[100];
      string loopname;
      string treename;
      systtools::SystParamHeader head = rh.GetHeader(parIds[j]);
      sprintf(lloopname,"%s",head.prettyName.c_str());
      loopname = lloopname;
      if (loopname != "MaCCQE") {
        continue; }
      treename = loopname + "_tree"; 
      looptree = (TTree*) outf->Get(treename.c_str());
      TArrayF dVals = TArrayF((*it).responses.size());
      TArrayF wVals = TArrayF((*it).responses.size());
      looptree->Branch((loopname+"_weights").c_str(),&wVals);
      looptree->Branch((loopname+"_knots").c_str(),&dVals);
      j = j+1;
      k = k+1;
      for( int i = 0; i < (*it).responses.size(); ++i ) {
        dVals.SetAt(i-3,i);
	wVals.SetAt((*it).responses[i],i);
        }
      
      looptree->Fill();
      }
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
 
  const char* path  = "/vols/t2k/users/ljw20/software/m3_dune/inputs/mtuples/GenieSampleCuts"; 
  const char* outpath  = "/vols/t2k/users/ljw20/software/m3_dune/inputs/mtuples/WeightFiles"; 
  //std::string outfolder = "./weights";
 
  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--infile") ) {
      infile = argv[i+1];
      i += 2;
    } //else if( argv[i] == std::string("--gfile") ) {
      //gfile = argv[i+1];
      //i += 2;} 
      //else if( argv[i] == std::string("--outfile") ) {
      //outfile = argv[i+1];
      //i += 2; }
      else if( argv[i] == std::string("--fhicl") ) {
      fhicl_filename = argv[i+1];
      i += 2;} 
    else i += 1; // look for next thing
  }

  outfile2 = infile.substr(0, infile.find_last_of('.')) + "_weight.root";
  char filepath[200];
  char outfilepath[200];
  sprintf(filepath, "%s/%s", path, infile.c_str());
  sprintf(outfilepath, "%s/%s", outpath, outfile2.c_str());
  std::cout <<"Output file: " << outfilepath <<std::endl;
  TFile * out = new TFile(outfilepath, "CREATE" );
  TFile * tf = new TFile(filepath);
  
  TTree * tree = (TTree*) tf->Get( "caf" );
   
  TTree * gtree = (TTree*) tf->Get( "genieEvt" );

  loop(tf, out, tree, gtree, fhicl_filename );

  std::cout << "Writing output ROOT File" << std::endl;
  out->Write();
  out->Close();
}
