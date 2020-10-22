//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  1 18:02:00 2020 by ROOT version 6.22/00
// from TTree Tracks/Tracks hitting volumes marked for track interception
// found on file: luxe_hics_signal_165gev_3031nm_cv7emstd_1umstep_tv4_hv2_1_t0.root
//////////////////////////////////////////////////////////

#ifndef FromTrack_h
#define FromTrack_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class FromTrack {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventid;
   vector<int>     *trackid;
   vector<int>     *detid;
   vector<int>     *pdg;
   vector<int>     *physproc;
   vector<double>  *E;
   vector<double>  *x;
   vector<double>  *y;
   vector<double>  *z;
   vector<double>  *t;
   vector<double>  *vtxx;
   vector<double>  *vtxy;
   vector<double>  *vtxz;
   vector<double>  *px;
   vector<double>  *py;
   vector<double>  *pz;
   vector<double>  *theta;
   vector<double>  *phi;
   vector<double>  *xlocal;
   vector<double>  *ylocal;
   vector<double>  *zlocal;
   Double_t        weight;
   vector<int>     *ptrackid;
   vector<int>     *nsecondary;

   // List of branches
   TBranch        *b_eventid;   //!
   TBranch        *b_trackid;   //!
   TBranch        *b_detid;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_physproc;   //!
   TBranch        *b_E;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_t;   //!
   TBranch        *b_vtxx;   //!
   TBranch        *b_vtxy;   //!
   TBranch        *b_vtxz;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_xlocal;   //!
   TBranch        *b_ylocal;   //!
   TBranch        *b_zlocal;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_ptrackid;   //!
   TBranch        *b_nsecondary;   //!

   FromTrack(TTree *tree=0);
   virtual ~FromTrack();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FromTrack_cxx
FromTrack::FromTrack(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("luxe_hics_signal_165gev_3031nm_cv7emstd_1umstep_tv4_hv2_1_t0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("luxe_hics_signal_165gev_3031nm_cv7emstd_1umstep_tv4_hv2_1_t0.root");
      }
      f->GetObject("Tracks",tree);

   }
   Init(tree);
}

FromTrack::~FromTrack()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FromTrack::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FromTrack::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FromTrack::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trackid = 0;
   detid = 0;
   pdg = 0;
   physproc = 0;
   E = 0;
   x = 0;
   y = 0;
   z = 0;
   t = 0;
   vtxx = 0;
   vtxy = 0;
   vtxz = 0;
   px = 0;
   py = 0;
   pz = 0;
   theta = 0;
   phi = 0;
   xlocal = 0;
   ylocal = 0;
   zlocal = 0;
   ptrackid = 0;
   nsecondary = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventid", &eventid, &b_eventid);
   fChain->SetBranchAddress("trackid", &trackid, &b_trackid);
   fChain->SetBranchAddress("detid", &detid, &b_detid);
   fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
   fChain->SetBranchAddress("physproc", &physproc, &b_physproc);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
   fChain->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
   fChain->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("xlocal", &xlocal, &b_xlocal);
   fChain->SetBranchAddress("ylocal", &ylocal, &b_ylocal);
   fChain->SetBranchAddress("zlocal", &zlocal, &b_zlocal);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("ptrackid", &ptrackid, &b_ptrackid);
   fChain->SetBranchAddress("nsecondary", &nsecondary, &b_nsecondary);
   Notify();
}

Bool_t FromTrack::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FromTrack::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FromTrack::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FromTrack_cxx
