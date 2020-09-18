//
//


#ifndef __RUN_PROC_TREE__

void process_track_tree_draw_v2(const char *fnlist = 0, const char *commentstr = 0)
{
   gROOT->ProcessLineSync("#define __RUN_PROC_TREE__ 1");
   gROOT->ProcessLineSync(".L MHists.C+");
   gROOT->ProcessLine("#include \"process_track_tree_draw_v2.C\"");
   gROOT->ProcessLine(fnlist ? Form("run_process_track_tree_draw(\"%s\")", fnlist) : "run_process_track_tree_draw(0)");
   gROOT->ProcessLine("#undef __RUN_PROC_TREE__");
}

#else

// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <algorithm>
// #include <stdexcept>
// #include <vector>
// #include <tuple>

// #include "TChain.h"

// #include "MHists.h"

int ProcessList(const std::string &fnamelist, std::vector<std::string> &flist);
void CreateHistograms(MHists *mh);

typedef std::tuple<int,int,int, double,double,double, double,double,double,double, double,double,double, double> TrckType;
//track_id, det_id, pdg, xx, yy, zz, E, pxx, pyy, pzz, vtx_x, vtx_y, vtx_z, ev_weight

void TestTracks(const int evid, std::vector<TrckType> evtracks, MHists *mhist);
void DumpTrack(const int evid, const TrckType &trk);

int run_process_track_tree_draw(const char *fnlist = 0, const char *commentstr = 0)
{
  int debugl = 0; //1;

  if (!fnlist) {
    std::cout << "Usage: root -l process_track_tree_draw_v2.C\'(\"file_with_list_of_mc_files\")\'\n";
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
   
  TChain *track_tree = new TChain("Tracks");
  std::for_each(flist.begin(), flist.end(), [track_tree](const std::string ss) {track_tree->Add(ss.c_str());} );

  double eneg, xx, yy, zz, tt, pxx, pyy, pzz, vtx_x, vtx_y, vtx_z, ev_weight;
  int det_id, pdg, ev_id, track_id;
  track_tree->SetBranchAddress("x", &xx);
  track_tree->SetBranchAddress("y", &yy);
  track_tree->SetBranchAddress("z", &zz);
  track_tree->SetBranchAddress("t", &tt);
  track_tree->SetBranchAddress("px", &pxx);
  track_tree->SetBranchAddress("py", &pyy);
  track_tree->SetBranchAddress("pz", &pzz);
  track_tree->SetBranchAddress("E",  &eneg);
  track_tree->SetBranchAddress("vtxx", &vtx_x);
  track_tree->SetBranchAddress("vtxy", &vtx_y);
  track_tree->SetBranchAddress("vtxz", &vtx_z);
  track_tree->SetBranchAddress("detid", &det_id);
  track_tree->SetBranchAddress("pdg", &pdg);
  track_tree->SetBranchAddress("eventid", &ev_id);
  track_tree->SetBranchAddress("trackid", &track_id);
  track_tree->SetBranchAddress("weight", &ev_weight);

  const double ipmagcutx = 165.0; 
  const double ipmagcuty = 54.0; 
  const double zproj = 2748.0;

  std::vector<TrckType> event_tracks;
  int pevent = -1;
  int nevproc = track_tree->GetEntries();
  std::cout << "Total number of events: " << nevproc << std::endl;
  
  for (Long64_t ii = 0; ii < nevproc; ++ii) {
    track_tree->GetEntry(ii);
//     if (eneg < 1.0) continue;

    double xproj = pxx/pzz*(zproj-zz) + xx;
    double yproj = pyy/pzz*(zproj-zz) + yy;

    if (!(ii % 10000000)) { std::cout << "Event " << ii << std::endl; }

    if (det_id >= 1000 && det_id <= 1015) {
      int det = det_id - 1000;
      lhist->FillHistW("tracking_planes_track_x", det, xx, ev_weight);
      lhist->FillHistW("tracking_planes_track_y", det, yy, ev_weight);
      lhist->FillHistW("tracking_planes_vtx_z", det, vtx_z, ev_weight);
      lhist->FillHistW("tracking_planes_track_xy", det, xx, yy, ev_weight);
      lhist->FillHistW("tracking_planes_vtxz_vtxx", det, vtx_z, vtx_x, ev_weight);
      lhist->FillHistW("tracking_planes_vtxz_vtxy", det, vtx_z, vtx_y, ev_weight);
      lhist->FillHistW("tracking_planes_vtxx_vtxy", det, vtx_x, vtx_y, ev_weight);

      if (pdg == 11)  {  
        lhist->FillHistW("tracking_planes_track_e_electrons", det, eneg, ev_weight); 
        lhist->FillHistW("tracking_planes_vtx_z_electrons", det, vtx_z, ev_weight);
        lhist->FillHistW("tracking_planes_track_xy_electrons", det, xx, yy, ev_weight);
        lhist->FillHistW("tracking_planes_vtxz_vtxx_electrons", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("tracking_planes_vtxz_vtxy_electrons", det, vtx_z, vtx_y, ev_weight);
        lhist->FillHistW("tracking_planes_vtxx_vtxy_electrons", det, vtx_x, vtx_y, ev_weight);
        if (vtx_z < 3000.0 && fabs(vtx_x) < 165.0 && fabs(vtx_y) < ipmagcutx) {     
          lhist->FillHistW("tracking_planes_track_xy_electrons_vtxcut", det, xx, yy, ev_weight);
          lhist->FillHistW("tracking_planes_track_x_electrons_vtxcut", det, xx, ev_weight);
        }
        if (fabs(xproj) < ipmagcutx && fabs(yproj) < ipmagcuty) {     
          lhist->FillHistW("tracking_planes_track_xy_electrons_pcut", det, xx, yy, ev_weight);
          lhist->FillHistW("tracking_planes_track_x_electrons_pcut", det, xx, ev_weight);
        }
      }

      if (pdg == -11) {  
        lhist->FillHistW("tracking_planes_track_e_positrons", det, eneg, ev_weight); 
        lhist->FillHistW("tracking_planes_vtx_z_positrons", det, vtx_z, ev_weight);
        lhist->FillHistW("tracking_planes_track_xy_positrons", det, xx, yy, ev_weight);
        lhist->FillHistW("tracking_planes_vtxz_vtxx_positrons", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("tracking_planes_vtxz_vtxy_positrons", det, vtx_z, vtx_y, ev_weight);
        lhist->FillHistW("tracking_planes_vtxx_vtxy_positrons", det, vtx_x, vtx_y, ev_weight);
        if (vtx_z < 3000.0 && fabs(vtx_x) < 165.0 && fabs(vtx_y) < ipmagcutx) {     
          lhist->FillHistW("tracking_planes_track_xy_positrons_vtxcut", det, xx, yy, ev_weight);
          lhist->FillHistW("tracking_planes_track_x_positrons_vtxcut", det, xx, ev_weight);
        }
        if (fabs(xproj) < ipmagcutx && fabs(yproj) < ipmagcuty) {     
          lhist->FillHistW("tracking_planes_track_xy_positrons_pcut", det, xx, yy, ev_weight);
          lhist->FillHistW("tracking_planes_track_x_positrons_pcut", det, xx, ev_weight);
        }
      }

      if (pdg == 22)  {  
        lhist->FillHistW("tracking_planes_track_e_gamma", det, eneg, ev_weight); 
        lhist->FillHistW("tracking_planes_vtx_z_gamma", det, vtx_z, ev_weight);
        lhist->FillHistW("tracking_planes_track_xy_gamma", det, xx, yy, ev_weight);
        lhist->FillHistW("tracking_planes_vtxz_vtxx_gamma", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("tracking_planes_vtxz_vtxy_gamma", det, vtx_z, vtx_y, ev_weight);
        lhist->FillHistW("tracking_planes_vtxx_vtxy_gamma", det, vtx_x, vtx_y, ev_weight);
        if (vtx_z < 3000.0 && fabs(vtx_x) < 165.0 && fabs(vtx_y) < ipmagcutx) {     
          lhist->FillHistW("tracking_planes_track_xy_gamma_vtxcut", det, xx, yy, ev_weight);
        }
        if (fabs(xproj) < ipmagcutx && fabs(yproj) < ipmagcuty) {     
          lhist->FillHistW("tracking_planes_track_xy_gamma_pcut", det, xx, yy, ev_weight);
        }
      }
    }

    double ecal_ytop = 26.5;
    double ecal_ybot = -26.5;
    double ecal_zfrnt = 4260.0;
    double ecal_zback = 4345.0;
    double ecal_xleft = 578.5;
    double ecal_xrght = -578.5;
    double ecal_xinleft = 30.5;
    double ecal_xinrght = -30.5;
    if (det_id >= 2000 && det_id <= 2001) {
      int det = det_id - 2000;
      lhist->FillHistW("ecal_track_x", det, xx, ev_weight);
      lhist->FillHistW("ecal_track_y", det, yy, ev_weight);
      lhist->FillHistW("ecal_track_vtx_z", det, vtx_z, ev_weight);
      lhist->FillHistW("ecal_track_xy", det, xx, yy, ev_weight);
      lhist->FillHistW("ecal_track_vtxz_vtxx", det, vtx_z, vtx_x, ev_weight);
      lhist->FillHistW("ecal_track_vtxz_vtxy", det, vtx_z, vtx_y, ev_weight);
      lhist->FillHistW("ecal_track_vtxx_vtxy", det, vtx_x, vtx_y, ev_weight);
      if (zz > ecal_zfrnt && zz < ecal_zback) {
        if (yy > ecal_ytop) {
          lhist->FillHistW("ecal_track_vtxz_vtxx_top", det, vtx_z, vtx_x, ev_weight);
          lhist->FillHistW("ecal_track_vtxz_vtxy_top", det, vtx_z, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_vtxx_vtxy_top", det, vtx_x, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_zx_top", det, zz, xx, ev_weight);
          lhist->FillHistW("ecal_track_e_top", det, eneg, ev_weight);
        }
        if (yy < ecal_ybot) {
          lhist->FillHistW("ecal_track_vtxz_vtxx_bottom", det, vtx_z, vtx_x, ev_weight);
          lhist->FillHistW("ecal_track_vtxz_vtxy_bottom", det, vtx_z, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_vtxx_vtxy_bottom", det, vtx_x, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_zx_bottom", det, zz, xx, ev_weight);
          lhist->FillHistW("ecal_track_e_bottom", det, eneg, ev_weight);
        }
        if (xx > ecal_xleft) {
          lhist->FillHistW("ecal_track_vtxz_vtxx_left", det, vtx_z, vtx_x, ev_weight);
          lhist->FillHistW("ecal_track_vtxz_vtxy_left", det, vtx_z, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_vtxx_vtxy_left", det, vtx_x, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_zy_left", det, zz, yy, ev_weight);
          lhist->FillHistW("ecal_track_e_left", det, eneg, ev_weight);
        }
        if (xx < ecal_xrght) {
          lhist->FillHistW("ecal_track_vtxz_vtxx_right", det, vtx_z, vtx_x, ev_weight);
          lhist->FillHistW("ecal_track_vtxz_vtxy_right", det, vtx_z, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_vtxx_vtxy_right", det, vtx_x, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_zy_right", det, zz, yy, ev_weight);
          lhist->FillHistW("ecal_track_e_right", det, eneg, ev_weight);
        }
        if (xx > ecal_xinrght && xx < ecal_xinleft) {
          lhist->FillHistW("ecal_track_vtxz_vtxx_inner", det, vtx_z, vtx_x, ev_weight);
          lhist->FillHistW("ecal_track_vtxz_vtxy_inner", det, vtx_z, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_vtxx_vtxy_inner", det, vtx_x, vtx_y, ev_weight);
          lhist->FillHistW("ecal_track_zy_inner", det, zz, yy, ev_weight);
          lhist->FillHistW("ecal_track_e_inner", det, eneg, ev_weight);
        }
        lhist->FillHistW("ecal_track_vtxx_vtxy_zcut", det, vtx_x, vtx_y, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxx_zcut", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxy_zcut", det, vtx_z, vtx_y, ev_weight);
        lhist->FillHistW("ecal_track_pdg_zcut", det, pdg, ev_weight);
        lhist->FillHistW("ecal_track_e_zcut", det, eneg, ev_weight);
      }
      
      if (zz < ecal_zfrnt) {
        lhist->FillHistW("ecal_track_vtxx_vtxy_front", det, vtx_x, vtx_y, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxx_front", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxy_front", det, vtx_z, vtx_y, ev_weight);
        lhist->FillHistW("ecal_track_xy_front", det, xx, yy, ev_weight);
        lhist->FillHistW("ecal_track_pdg_front", det, pdg, ev_weight);
        lhist->FillHistW("ecal_track_e_front", det, eneg, ev_weight);
      }
      if (zz > ecal_zback) {
        lhist->FillHistW("ecal_track_vtxx_vtxy_back", det, vtx_x, vtx_y, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxx_back", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxy_back", det, vtx_z, vtx_y, ev_weight);
        lhist->FillHistW("ecal_track_xy_back", det, xx, yy, ev_weight);
        lhist->FillHistW("ecal_track_pdg_back", det, pdg, ev_weight);
        lhist->FillHistW("ecal_track_e_back", det, eneg, ev_weight);
      }

      if (pdg == 11)  {  
        lhist->FillHistW("ecal_track_e_electrons", det, eneg, ev_weight); 
        lhist->FillHistW("ecal_track_vtx_z_electrons", det, vtx_z, ev_weight);
        lhist->FillHistW("ecal_track_xy_electrons", det, xx, yy, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxx_electrons", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxy_electrons", det, vtx_z, vtx_y, ev_weight);
      }

      if (pdg == -11) {  
        lhist->FillHistW("ecal_track_e_positrons", det, eneg, ev_weight); 
        lhist->FillHistW("ecal_track_vtx_z_positrons", det, vtx_z, ev_weight);
        lhist->FillHistW("ecal_track_xy_positrons", det, xx, yy, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxx_positrons", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxy_positrons", det, vtx_z, vtx_y, ev_weight);
      }

      if (pdg == 22)  {  
        lhist->FillHistW("ecal_track_e_gamma", det, eneg, ev_weight); 
        lhist->FillHistW("ecal_track_vtx_z_gamma", det, vtx_z, ev_weight);
        lhist->FillHistW("ecal_track_xy_gamma", det, xx, yy, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxx_gamma", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("ecal_track_vtxz_vtxy_gamma", det, vtx_z, vtx_y, ev_weight);
      }
    }

    if (det_id >= 3000 && det_id <= 3001) {
      int det = det_id - 3000;
      lhist->FillHistW("lyso_track_x", det, xx, ev_weight);
      lhist->FillHistW("lyso_track_y", det, yy, ev_weight);
      lhist->FillHistW("lyso_track_vtx_z", det, vtx_z, ev_weight);
      lhist->FillHistW("lyso_track_xy", det, xx, yy, ev_weight);
      lhist->FillHistW("lyso_track_vtxz_vtxx", det, vtx_z, vtx_x, ev_weight);
      lhist->FillHistW("lyso_track_vtxz_vtxy", det, vtx_z, vtx_y, ev_weight);
      if (pdg == 11)  {  
        lhist->FillHistW("lyso_track_e_electrons", det, eneg, ev_weight); 
      lhist->FillHistW("lyso_track_vtx_z_electrons", det, vtx_z, ev_weight);
      lhist->FillHistW("lyso_track_xy_electrons", det, xx, yy, ev_weight);
      lhist->FillHistW("lyso_track_vtxz_vtxx_electrons", det, vtx_z, vtx_x, ev_weight);
      lhist->FillHistW("lyso_track_vtxz_vtxy_electrons", det, vtx_z, vtx_y, ev_weight);
      }

      if (pdg == -11) {  
        lhist->FillHistW("lyso_track_e_positrons", det, eneg, ev_weight); 
      lhist->FillHistW("lyso_track_vtx_z_positrons", det, vtx_z, ev_weight);
      lhist->FillHistW("lyso_track_xy_positrons", det, xx, yy, ev_weight);
      lhist->FillHistW("lyso_track_vtxz_vtxx_positrons", det, vtx_z, vtx_x, ev_weight);
      lhist->FillHistW("lyso_track_vtxz_vtxy_positrons", det, vtx_z, vtx_y, ev_weight);
      }

      if (pdg == 22)  {  
        lhist->FillHistW("lyso_track_e_gamma", det, eneg, ev_weight); 
        lhist->FillHistW("lyso_track_vtx_z_gamma", det, vtx_z, ev_weight);
        lhist->FillHistW("lyso_track_xy_gamma", det, xx, yy, ev_weight);
        lhist->FillHistW("lyso_track_vtxz_vtxx_gamma", det, vtx_z, vtx_x, ev_weight);
        lhist->FillHistW("lyso_track_vtxz_vtxy_gamma", det, vtx_z, vtx_y, ev_weight);
      }
    }

    if (det_id >= 4000 && det_id <= 4007) {
      int det = det_id - 4000;
      lhist->FillHistW("gammamon_track_z", det, zz, ev_weight);
      lhist->FillHistW("gammamon_track_vtx_z", det, vtx_z, ev_weight);
      lhist->FillHistW("gammamon_track_xy", det, xx, yy, ev_weight);
      lhist->FillHistW("gammamon_track_vtxz_vtxx", det, vtx_z, vtx_x, ev_weight);
      lhist->FillHistW("gammamon_track_vtxz_vtxy", det, vtx_z, vtx_y, ev_weight);
      if (pdg == 11)  {  
        lhist->FillHistW("gammamon_track_e_electrons", det, eneg, ev_weight); 
      }

      if (pdg == -11) {  
        lhist->FillHistW("gammamon_track_e_positrons", det, eneg, ev_weight); 
      }

      if (pdg == 22)  {  
        lhist->FillHistW("gammamon_track_e_gamma", det, eneg, ev_weight); 
      }
    }
   
    // Accumulate tracks for one event 
    if (pevent == ev_id) {
      event_tracks.push_back(std::make_tuple(track_id, det_id, pdg, xx, yy, zz,eneg, pxx, pyy, pzz, 
                                              vtx_x, vtx_y, vtx_z, ev_weight));
    } else if (pevent != ev_id) {
      TestTracks(pevent, event_tracks, lhist);
      event_tracks.clear();
      event_tracks.push_back(std::make_tuple(track_id, det_id, pdg, xx, yy, zz, eneg, pxx, pyy, pzz, 
                                              vtx_x, vtx_y, vtx_z, ev_weight));
      pevent = ev_id;
    } else { //(pevent > ev_id)
      std::cout << "New event ID is smaller than previous!\n";
      delete track_tree;
      delete lhist;
      return -1;
    }

  } // tree loop
    
 //*********************To draw histos
//  lhist->DrawHist1D_BT15("tracking_planes_track_x");
//  lhist->DrawHist1D_BT15("tracking_planes_track_y");
//  lhist->DrawHist1D_BT15("tracking_planes_vtx_z");

//  lhist->DrawHist1D_BT15("tracking_planes_track_e_electrons");
//  lhist->DrawHist1D_BT15("tracking_planes_track_e_positrons");
//  lhist->DrawHist1D_BT15("tracking_planes_track_e_gamma");
 
   lhist->DrawHist2D_BT15("tracking_planes_track_xy", "tracking_planes_track_xy", "x (mm)", "y (mm)", "N", "colz");
   lhist->DrawHist2D_BT15("tracking_planes_track_xy_charged_tracks", "tracking_planes_track_xy_charged_tracks", "x (mm)", "y (mm)", "N", "colz");
//   lhist->DrawHist2D_BT15("tracking_planes_track_xy_electrons");
//   lhist->DrawHist2D_BT15("tracking_planes_track_xy_positrons");
//   lhist->DrawHist2D_BT15("tracking_planes_track_xy_gamma");
//   lhist->DrawHist2D_BT15("tracking_planes_vtxz_vtxx");
//   lhist->DrawHist2D_BT15("tracking_planes_vtxz_vtxy");
//   lhist->DrawHist2D_BT15("tracking_planes_vtxx_vtxy");

  std::string suffix("");
  std::string foutname = fnamelist.substr(fnamelist.find_last_of("/")+1);
  foutname = foutname.substr(0, foutname.find_last_of("."));
  foutname += suffix + std::string(".root");
  lhist->SaveHists(foutname);

//  delete track_tree;
  return 0;  
}



void TestTracks(const int evid, std::vector<TrckType> evtracks, MHists *mhist)
{
//   DumpTrack(evid, trk);
  std::stable_sort(evtracks.begin(), evtracks.end(), [](TrckType x, TrckType y) {if (std::get<0>(x) < std::get<0>(y)) 
                                                                                 return true; else return false;} );
  std::map<int, int> track_count;
  int ltid = -100;
  int tcount = 0;
  for (const auto &trk : evtracks) {
//     std::cout << "Sorted Track id: " << std::get<0>(trk) << " det id: " << std::get<1>(trk) << std::endl;
    int detid = std::get<1>(trk);
//     if (detid == 2000 || detid == 2001) DumpTrack(evid, trk);
    if (detid < 1000 || detid > 1015) continue;
    int trackid = std::get<0>(trk);
    if (trackid == ltid) {
      ++tcount;  
    } else {
      if (ltid > 0) {
        track_count[ltid] = tcount; 
        mhist->FillHistW("tracking_planes_n_track_cross", 0, tcount, std::get<13>(trk));  
      }
      ltid = trackid;
      tcount = 1;
    }
  }
  
  const double ipmagcutx = 165.0; 
  const double ipmagcuty = 54.0; 
  const double zproj = 2748.0;
  
  int ntrckmax = 3;
  for (const auto &tid : track_count) {
    if (tid.second >= ntrckmax) {
      for (const auto &trk : evtracks) {
        int detid = std::get<1>(trk);
        if (detid < 1000 || detid > 1015) continue;
        int trackid = std::get<0>(trk);
        if (trackid == tid.first) {
          int pdg = std::get<2>(trk);
          int det = detid-1000;
          if (pdg==11 || pdg==-11) {
            double xx = std::get<3>(trk);
            double yy = std::get<4>(trk);
            double zz = std::get<5>(trk);
            double pxx = std::get<7>(trk);
            double pyy = std::get<8>(trk);
            double pzz = std::get<9>(trk);
            double ev_weight = std::get<13>(trk);
            mhist->FillHistW("tracking_planes_track_xy_charged_tracks", det, xx, yy, ev_weight);
            
            double xproj = pxx/pzz*(zproj-zz) + xx;
            double yproj = pyy/pzz*(zproj-zz) + yy;
            if (fabs(xproj) < ipmagcutx && fabs(yproj) < ipmagcuty) {
              mhist->FillHistW("tracking_planes_track_xy_charged_tracks_cut", det, xx, yy, ev_weight);
            }
          }
        }
      }
    }
  }
  
}


void DumpTrack(const int evid, const TrckType &trk)
{
 std::cout << "Event: " << evid << "  ";
   std::cout << "Track (id / det_id / pdg / E / xyz / p_xyz / vtx_xyz/ w /:" 
             << "  " << std::get<0>(trk) 
             << "  " << std::get<1>(trk)
             << "  " << std::get<2>(trk)
             << "  " << std::get<6>(trk)
             << "  " << std::get<3>(trk)
             << "  " << std::get<4>(trk)
             << "  " << std::get<5>(trk)
             << "  " << std::get<7>(trk)
             << "  " << std::get<8>(trk)
             << "  " << std::get<9>(trk)
             << "  " << std::get<10>(trk)
             << "  " << std::get<11>(trk)
             << "  " << std::get<12>(trk)
             << "  " << std::get<13>(trk) << std::endl;
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
  int ndet = 16;
  mh->AddHists("tracking_planes_track_x", ndet, 1300, -650.0, 650.0);
  mh->AddHists("tracking_planes_track_x_electrons_vtxcut", ndet, 1300, -650.0, 650.0);
  mh->AddHists("tracking_planes_track_x_positrons_vtxcut", ndet, 1300, -650.0, 650.0);
  mh->AddHists("tracking_planes_track_x_electrons_pcut", ndet, 1300, -650.0, 650.0);
  mh->AddHists("tracking_planes_track_x_positrons_pcut", ndet, 1300, -650.0, 650.0);
  mh->AddHists("tracking_planes_track_y", ndet, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_electrons", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_positrons", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_gamma", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_electrons_vtxcut", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_positrons_vtxcut", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_gamma_vtxcut", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_electrons_pcut", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_positrons_pcut", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_gamma_pcut", ndet, 1300, -650.0, 650.0, 50, -25.0, 25.0);
  mh->AddHists("tracking_planes_vtx_z", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("tracking_planes_vtx_z_electrons", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("tracking_planes_vtx_z_gamma", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("tracking_planes_vtx_z_positrons", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("tracking_planes_vtxz_vtxx", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxz_vtxy", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxx_vtxy", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxz_vtxx_electrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxz_vtxy_electrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxx_vtxy_electrons", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxz_vtxx_gamma", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxz_vtxy_gamma", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxx_vtxy_gamma", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxz_vtxx_positrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxz_vtxy_positrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("tracking_planes_vtxx_vtxy_positrons", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);

  mh->AddHists("tracking_planes_track_e_electrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("tracking_planes_track_e_positrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("tracking_planes_track_e_gamma", ndet, 2000, 0.0, 20.0);

  mh->AddHists("tracking_planes_track_xy_charged_tracks", ndet, 650, -650.0, 650.0, 25, -25.0, 25.0);
  mh->AddHists("tracking_planes_track_xy_charged_tracks_cut", ndet, 650, -650.0, 650.0, 25, -25.0, 25.0);

  mh->AddHists("tracking_planes_n_track_cross", 1, 20, 0.0, 20.0);
  
  ndet = 2;
  mh->AddHists("ecal_track_x", ndet, 1300, -650.0, 650.0);
  mh->AddHists("ecal_track_y", ndet, 100, -50.0, 50.0);
  mh->AddHists("ecal_track_xy", ndet, 1300, -650.0, 650.0, 100, -50.0, 50.0);
  mh->AddHists("ecal_track_xy_electrons", ndet, 1300, -650.0, 650.0, 100, -50.0, 50.0);
  mh->AddHists("ecal_track_xy_gamma", ndet, 1300, -650.0, 650.0, 100, -50.0, 50.0);
  mh->AddHists("ecal_track_xy_positrons", ndet, 1300, -650.0, 650.0, 100, -50.0, 50.0);
  mh->AddHists("ecal_track_vtx_z", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("ecal_track_vtx_z_electrons", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("ecal_track_vtx_z_gamma", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("ecal_track_vtx_z_positrons", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("ecal_track_vtxz_vtxx", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxy", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxx_vtxy", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxx_vtxy_zcut", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxx_zcut", ndet, 2000, 3000.0, 5000.0, 2000, -2000.0, 2000.0);
  mh->AddHists("ecal_track_vtxz_vtxy_zcut", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxx_electrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxy_electrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxx_vtxy_electrons", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxx_gamma", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxy_gamma", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxx_vtxy_gamma", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxx_positrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxz_vtxy_positrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_vtxx_vtxy_positrons", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("ecal_track_e_electrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("ecal_track_e_positrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("ecal_track_e_gamma", ndet, 2000, 0.0, 20.0);
  mh->AddHists("ecal_track_pdg_zcut", ndet, 200, -100.0, 100.0);
  mh->AddHists("ecal_track_e_zcut", ndet, 2000, 0.0, 20.0);

  int nvtxz = 14000;
  int nvtxx = 200;
  int nvtxy = 200;
  double vtxz_min = -8000.0;
  double vtxz_max =  6000.0;
  double vtxx_min = -100.0;
  double vtxx_max =  100.0;
  double vtxy_min = -100.0;
  double vtxy_max =  100.0;
  mh->AddHists("ecal_track_vtxz_vtxx_top", ndet, nvtxz, vtxz_min, vtxz_max, nvtxx, vtxx_min, vtxx_max);
  mh->AddHists("ecal_track_vtxz_vtxy_top", ndet, nvtxz, vtxz_min, vtxz_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_vtxx_vtxy_top", ndet, nvtxx, vtxx_min, vtxx_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_zx_top", ndet, 100, 4200.0, 4400.0, 600, -600.0, 600.0);
  mh->AddHists("ecal_track_e_top", ndet, 2000, 0.0, 20.0);

  mh->AddHists("ecal_track_vtxz_vtxx_bottom", ndet, nvtxz, vtxz_min, vtxz_max, nvtxx, vtxx_min, vtxx_max);
  mh->AddHists("ecal_track_vtxz_vtxy_bottom", ndet, nvtxz, vtxz_min, vtxz_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_vtxx_vtxy_bottom", ndet, nvtxx, vtxx_min, vtxx_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_zx_bottom", ndet, 100, 4200.0, 4400.0, 600, -600.0, 600.0);
  mh->AddHists("ecal_track_e_bottom", ndet, 2000, 0.0, 20.0);

  mh->AddHists("ecal_track_vtxz_vtxx_left", ndet, nvtxz, vtxz_min, vtxz_max, nvtxx, vtxx_min, vtxx_max);
  mh->AddHists("ecal_track_vtxz_vtxy_left", ndet, nvtxz, vtxz_min, vtxz_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_vtxx_vtxy_left", ndet, nvtxx, vtxx_min, vtxx_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_zy_left", ndet, 100, 4200.0, 4400.0, 35, -35.0, 35.0);
  mh->AddHists("ecal_track_e_left", ndet, 2000, 0.0, 20.0);

  mh->AddHists("ecal_track_vtxz_vtxx_right", ndet, nvtxz, vtxz_min, vtxz_max, nvtxx, vtxx_min, vtxx_max);
  mh->AddHists("ecal_track_vtxz_vtxy_right", ndet, nvtxz, vtxz_min, vtxz_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_vtxx_vtxy_right", ndet, nvtxx, vtxx_min, vtxx_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_zy_right", ndet, 100, 4200.0, 4400.0, 35, -35.0, 35.0);
  mh->AddHists("ecal_track_e_right", ndet, 2000, 0.0, 20.0);

  mh->AddHists("ecal_track_vtxz_vtxx_inner", ndet, nvtxz, vtxz_min, vtxz_max, nvtxx, vtxx_min, vtxx_max);
  mh->AddHists("ecal_track_vtxz_vtxy_inner", ndet, nvtxz, vtxz_min, vtxz_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_vtxx_vtxy_inner", ndet, nvtxx, vtxx_min, vtxx_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_zy_inner", ndet, 100, 4200.0, 4400.0, 35, -35.0, 35.0);
  mh->AddHists("ecal_track_e_inner", ndet, 2000, 0.0, 20.0);

  mh->AddHists("ecal_track_vtxz_vtxx_front", ndet, nvtxz, vtxz_min, vtxz_max, nvtxx, vtxx_min, vtxx_max);
  mh->AddHists("ecal_track_vtxz_vtxy_front", ndet, nvtxz, vtxz_min, vtxz_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_vtxx_vtxy_front", ndet, nvtxx, vtxx_min, vtxx_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_xy_front", ndet, 600, -600.0, 600.0, 35, -35.0, 35.0);
  mh->AddHists("ecal_track_pdg_front", ndet, 200, -100.0, 100.0);
  mh->AddHists("ecal_track_e_front", ndet, 2000, 0.0, 20.0);

  mh->AddHists("ecal_track_vtxz_vtxx_back", ndet, nvtxz, vtxz_min, vtxz_max, nvtxx, vtxx_min, vtxx_max);
  mh->AddHists("ecal_track_vtxz_vtxy_back", ndet, nvtxz, vtxz_min, vtxz_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_vtxx_vtxy_back", ndet, nvtxx, vtxx_min, vtxx_max, nvtxy, vtxy_min, vtxy_max);
  mh->AddHists("ecal_track_xy_back", ndet, 600, -600.0, 600.0, 35, -35.0, 35.0);
  mh->AddHists("ecal_track_pdg_back", ndet, 200, -100.0, 100.0);
  mh->AddHists("ecal_track_e_back", ndet, 2000, 0.0, 20.0);

  
  ndet = 2;
  mh->AddHists("lyso_track_x", ndet, 800, -400.0, 400.0);
  mh->AddHists("lyso_track_y", ndet, 100, -50.0, 50.0);
  mh->AddHists("lyso_track_xy", ndet, 800, -400.0, 400.0, 100, -50.0, 50.0);
  mh->AddHists("lyso_track_xy_electrons", ndet, 800, -400.0, 400.0, 100, -50.0, 50.0);
  mh->AddHists("lyso_track_xy_gamma", ndet, 800, -400.0, 400.0, 100, -50.0, 50.0);
  mh->AddHists("lyso_track_xy_positrons", ndet, 800, -400.0, 400.0, 100, -50.0, 50.0);
  mh->AddHists("lyso_track_vtx_z", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("lyso_track_vtx_z_electrons", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("lyso_track_vtx_z_gamma", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("lyso_track_vtx_z_positrons", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("lyso_track_vtxz_vtxx", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxz_vtxy", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxx_vtxy", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxz_vtxx_electrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxz_vtxy_electrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxx_vtxy_electrons", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxz_vtxx_gamma", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxz_vtxy_gamma", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxx_vtxy_gamma", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxz_vtxx_positrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxz_vtxy_positrons", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_vtxx_vtxy_positrons", ndet, 300, -3000.0, 3000.0, 300, -3000.0, 3000.0);
  mh->AddHists("lyso_track_e_electrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("lyso_track_e_positrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("lyso_track_e_gamma", ndet, 2000, 0.0, 20.0);

  ndet = 8;
  mh->AddHists("gammamon_track_z", ndet, 2000, 13000.0, 15000.0);
  mh->AddHists("gammamon_track_xy", ndet, 400, -200.0, 200.0, 400, -200.0, 200.0);
  mh->AddHists("gammamon_track_vtx_z", ndet, 30000, -10000.0, 20000.0);
  mh->AddHists("gammamon_track_vtxz_vtxx", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("gammamon_track_vtxz_vtxy", ndet, 3000, -10000.0, 20000.0, 300, -3000.0, 3000.0);
  mh->AddHists("gammamon_track_e_electrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("gammamon_track_e_positrons", ndet, 2000, 0.0, 20.0);
  mh->AddHists("gammamon_track_e_gamma", ndet, 2000, 0.0, 20.0);

}

#endif

