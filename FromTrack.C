#define FromTrack_cxx
#include "FromTrack.h"
#include <TH2.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;




void FromTrack::Loop()
{
//   In a ROOT session, you can do:
//      root> .L FromTrack.C
//      root> FromTrack t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    
   
   if (fChain == 0) return;
   TFile *fOut = new TFile("trackingResults.root", "RECREATE");
   fOut->cd();
   TH1F *h_TrackingPlane_Track_x[16];
   TH1F *h_TrackingPlane_Track_y[16];
   TH1F *h_TrackingPlane_vtx_z[16];
   TH2F *h_TrackingPlane_track_xy[16];
   TH2F *h_TrackingPlane_vtxz_vtxx[16];
   TH2F *h_TrackingPlane_vtxz_vtxy[16];
   TH2F *h_TrackingPlane_vtxx_vtxy[16];
   
   for(int i = 0; i < 16; ++i){
       
//       stringstream ss;
//       ss << "h_TrackingPlane_Track_x_Stave"+std::to_string(i);
      //cout << ss.str() << endl;
      TString trackPlaneX         = "h_TrackingPlane_Track_x_Stave"+std::to_string(i);
      TString trackPlaneY         = "h_TrackingPlane_Track_y_Stave"+std::to_string(i);
      TString trackPlaneVtxZ      = "h_TrackingPlane_vtx_z_Stave"+std::to_string(i);
      TString trackPlaneTrackXY   = "h_TrackingPlane_track_xy_Stave"+std::to_string(i);
      TString trackPlaneVtxzVtxx  = "h_TrackingPlane_vtxz_vtxx_Stave"+std::to_string(i);
      TString trackPlaneVtxzVtxy  = "h_TrackingPlane_vtxz_vtxy_Stave"+std::to_string(i);
      TString trackPlaneVtxxVtxy  = "h_TrackingPlane_vtxx_vtxy_Stave"+std::to_string(i);
       
      h_TrackingPlane_Track_x[i]   = new TH1F(trackPlaneX, trackPlaneX, 1300, -650.0, 650.0);
      h_TrackingPlane_Track_y[i]   = new TH1F(trackPlaneY, trackPlaneY, 50, -25.0, 25.0);
      h_TrackingPlane_vtx_z[i]     = new TH1F(trackPlaneVtxZ, trackPlaneVtxZ, 30000, -10000.0, 20000.0);
      
      h_TrackingPlane_track_xy[i]  = new TH2F(trackPlaneTrackXY, trackPlaneTrackXY, 1300, -650.0, 650.0, 50, -25.0, 25.0);
      h_TrackingPlane_vtxz_vtxx[i] = new TH2F(trackPlaneVtxzVtxx, trackPlaneVtxzVtxx, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
      h_TrackingPlane_vtxz_vtxy[i] = new TH2F(trackPlaneVtxzVtxy,trackPlaneVtxzVtxy, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
      h_TrackingPlane_vtxx_vtxy[i] = new TH2F(trackPlaneVtxxVtxy, trackPlaneVtxxVtxy, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
      
   }
   
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%1000==0)std::cout << "processed: " << jentry << std::endl;
//       std::cout << "eventid: " << eventid << std::endl;
//       std::cout << "detid size: " << detid->size() << ", pdg size: " << pdg->size() << " energy size: " << E->size() << std::endl;
      for(size_t j=0; j < detid->size(); ++j){
          if(detid->at(j) < 1000 || detid->at(j) > 1015)continue;
          int det = detid->at(j) - 1000;
          cout << detid->at(j) << endl;
          h_TrackingPlane_Track_x[det]->Fill(x->at(j),weight);
          h_TrackingPlane_Track_y[det]->Fill(y->at(j),weight);
          h_TrackingPlane_vtx_z[det]->Fill(vtxz->at(j), weight);
          h_TrackingPlane_track_xy[det]->Fill(x->at(j), y->at(j), weight);
          h_TrackingPlane_vtxz_vtxx[det]->Fill(vtxz->at(j), vtxx->at(j), weight);
          h_TrackingPlane_vtxz_vtxy[det]->Fill(vtxz->at(j), vtxy->at(j), weight);
          h_TrackingPlane_vtxx_vtxy[det]->Fill(vtxx->at(j), vtxy->at(j), weight);
      }

      // if (Cut(ientry) < 0) continue;
   }
   
   fOut->Write();
   fOut->Close();
}
