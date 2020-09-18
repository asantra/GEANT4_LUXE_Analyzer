//
//

#ifndef __RUN_PROC_HITS_TREE__

void process_hits_tree_draw_v2(const char *fnlist = 0, const char *commentstr = 0)
{
   gROOT->ProcessLineSync("#define __RUN_PROC_HITS_TREE__ 1");
   gROOT->ProcessLineSync(".L MHists.C+");
   gROOT->ProcessLine("#include \"process_hits_tree_draw_v2.C\"");
   gROOT->ProcessLine(fnlist ? Form("process_hits_tree_draw(\"%s\")", fnlist) : "process_hits_tree_draw(0)");
   gROOT->ProcessLine("#undef __RUN_PROC_HITS_TREE__");
}

#else

//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <sstream>
//#include <string>
//#include <vector>
//#include <algorithm>
//#include <stdexcept>

//#include "TChain.h"

//#include "MHists.h"

int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist);
void CreateHistograms(MHists *mh);


int process_hits_tree_draw(const char *fnlist = 0, const char *commentstr = 0)
{

  int debugl = 0; //1;

  if (!fnlist) {
    std::cout << "Usage: root -l process_hits_tree_draw_v2.C\'(\"file_with_list_of_mc_files\")\'\n";
    return -1;
  }
  
  std::string fnamelist(fnlist);
  std::vector<std::string>  flist;
  ProcessList(fnamelist, flist);  
  if (debugl) {
    std::cout << "The following files will be processed:\n";
    std::for_each(flist.begin(), flist.end(), [](const std::string ss) {std::cout << ss << std::endl;});
  }

  MHists *lhist = new MHists();
  CreateHistograms(lhist);
   
  TChain *hitstree = new TChain("Hits");
  std::for_each(flist.begin(), flist.end(), [hitstree](const std::string ss) {hitstree->Add(ss.c_str());} );

  double e_dep, ev_weight;
  int det_id, ev_id, layer_id, hit_id, cell_x, cell_y;
  std::vector<int> track_list;
  std::vector<double> track_x, track_y, track_z;

  std::vector<int> *track_list_p = &track_list;
  std::vector<double> *track_x_p = &track_x; 
  std::vector<double> *track_y_p = &track_y;
  std::vector<double> *track_z_p = &track_z;

  hitstree->SetBranchAddress("eventid", &ev_id);
  hitstree->SetBranchAddress("detid", &det_id);
  hitstree->SetBranchAddress("layerid", &layer_id);
  hitstree->SetBranchAddress("hitid", &hit_id);
  hitstree->SetBranchAddress("cellx", &cell_x);
  hitstree->SetBranchAddress("celly", &cell_y);
  hitstree->SetBranchAddress("edep",  &e_dep);
  hitstree->SetBranchAddress("weight", &ev_weight);
  hitstree->SetBranchAddress("track_list", &track_list_p);
  hitstree->SetBranchAddress("trackx", &track_x_p);
  hitstree->SetBranchAddress("tracky", &track_y_p);
  hitstree->SetBranchAddress("trackz", &track_z_p);
  
  const double ipmagcutx = 165.0; 
  const double ipmagcuty = 54.0; 
  const double zproj = 2748.0;

  Long64_t nevproc = hitstree->GetEntries();
  for (Long64_t ii = 0; ii < nevproc; ++ii) {
    hitstree->GetEntry(ii);

    if (!(ii % 10000000)) { std::cout << "Event " << ii << std::endl; }

    if (det_id >= 1000 && det_id <= 1008) {
      const int nlayers = 16;  
      int det = det_id - 1000;
      int detlayer = det * nlayers + layer_id;
      lhist->FillHistW("tracking_planes_hits_x", detlayer, cell_x, ev_weight);
      lhist->FillHistW("tracking_planes_hits_y", detlayer, cell_x, ev_weight);
      lhist->FillHistW("tracking_planes_hits_xy", detlayer, cell_x, cell_y, ev_weight);
      lhist->FillHistW("tracking_planes_hits_xy_edep", detlayer, cell_x, cell_y, e_dep*ev_weight);
    }

    if (det_id >= 2000 && det_id <= 2001) {
      const int nlayers = 21;  
      int det = det_id - 2000;
      int detlayer = det * nlayers + layer_id;
      lhist->FillHistW("ecal_hits_x", detlayer, cell_x, ev_weight);
      lhist->FillHistW("ecal_hits_edep", detlayer, e_dep, ev_weight);
      lhist->FillHistW("ecal_hits_xy", detlayer, cell_x, cell_y, ev_weight);
      lhist->FillHistW("ecal_hits_xy_edep", detlayer, cell_x, cell_y, e_dep*ev_weight);
      int ntrck = track_list.size();
      lhist->FillHistW("ecal_hits_ntrck", detlayer, ntrck, ev_weight);
      lhist->FillHistW("ecal_hits_xy_ntrck", detlayer, cell_x, cell_y, ntrck*ev_weight);
    }

    if (det_id >= 3000 && det_id <= 3001) {
      int det = det_id - 3000;
      lhist->FillHistW("lyso_hits_x", det, cell_x, ev_weight);
      lhist->FillHistW("lyso_hits_xy", det, cell_x, cell_y, ev_weight);
      lhist->FillHistW("lyso_hits_xy_edep", det, cell_x, cell_y, e_dep*ev_weight);
      lhist->FillHistW("lyso_hits_edep", det, e_dep, ev_weight);
    }

    if (det_id >= 4000 && det_id <= 4007) {
      int det = det_id - 4000;
      lhist->FillHistW("gammamon_hit_edep", det, e_dep, ev_weight);
    }

  } // tree loop
    
 //*********************To draw histos
  lhist->DrawHist1D_BT15("tracking_planes_hits_x");

  lhist->DrawHist1D_BT15("ecal_hits_edep");
  lhist->DrawHist2D_BT15("ecal_hits_xy", "ecal_hits_xy", "x (mm)", "y (mm)", "N", "colz");

//   lhist->DrawHist1D_BT15("lyso_hits_edep");

  std::string suffix("_hits");
  std::string foutname = fnamelist.substr(fnamelist.find_last_of("/")+1);
  foutname = foutname.substr(0, foutname.find_last_of("."));
  foutname += suffix + std::string(".root");
  lhist->SaveHists(foutname);

//  delete hitstree;
  return 0;  
}



int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist)
{
  std::fstream  fdata;
  fdata.open(fnamelist, std::ios::in);
  if (!fdata.is_open()) {
    throw std::runtime_error(std::string("Error reding data from the file ") + fnamelist);
  }
  
  unsigned long lid = 0;
  while (!fdata.eof()) {
    std::string  ffname;
    double fweight;
    fdata >> ffname;
    if (!fdata.fail()) { 
//       std::cout << "File name " << ffname << " is read from the list file" << std::endl;
      flist.push_back(ffname);
    }
    else if (fdata.eof()) { break; }
    else {
      std::cout << "ProcessList(..)  :  Error reading data from the file " << fnamelist 
                << ",  line: " << lid << ". Exit." << std::endl;
      fdata.close();          
      return -2;
    }
    ++lid;
  }
  
  fdata.close();

  return 0;
}  
  


void CreateHistograms(MHists *mh) 
{
  std::cout << "Creating histograms\n";
  int nlayers = 16;
  int nsensor = 9;
  int ndet = nsensor * nlayers;
  int npixx = 1024;
  int npixy = 512;
  mh->AddHists("tracking_planes_hits_x", ndet, npixx, 0.0, npixx);
  mh->AddHists("tracking_planes_hits_y", ndet, npixy, 0.0, npixy);
  mh->AddHists("tracking_planes_hits_xy", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
  mh->AddHists("tracking_planes_hits_xy_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);

  nlayers = 21;
  nsensor = 2;
  ndet = nsensor * nlayers;
  npixx = 110;
  npixy = 11;
  mh->AddHists("ecal_hits_x", ndet, npixx, 0.0, npixx);
  mh->AddHists("ecal_hits_edep", ndet, 10000, 0.0, 10.0);
  mh->AddHists("ecal_hits_xy", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
  mh->AddHists("ecal_hits_ntrck", ndet, 1000, 0.0, 1000);
  mh->AddHists("ecal_hits_xy_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
  mh->AddHists("ecal_hits_xy_ntrck", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
 
  ndet = 2;
  npixx = 150;
  npixy = 25;
  mh->AddHists("lyso_hits_x", ndet, npixx, 0.0, npixx);
  mh->AddHists("lyso_hits_xy", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
  mh->AddHists("lyso_hits_xy_edep", ndet, npixx, 0.0, npixx, npixy, 0.0, npixy);
  mh->AddHists("lyso_hits_edep", ndet, 10000, 0.0, 10.0);

  ndet = 8;
  mh->AddHists("gammamon_hit_edep", ndet, 10000, 0.0, 10000.0);

}


#endif

