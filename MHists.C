
// This class allows to create and manage histograms of similar type.
// E.g. create 128 identical histograms to study signal of each channel of one APV chip.
// There is MHist_test_v1(...) function for test and demostration.



// #if !defined(__CINT__) || defined(__MAKECINT__)
#if !defined(__CLING__) || defined(__ROOTCLING__)  // for root 6

 #include <iostream>
 #include <sstream>
 #include <string>
 #include <vector>
 #include <map>
 #include <algorithm> 
 #include <iterator>
 #include "stdlib.h"
 #include "TH1D.h"
 #include "TH2D.h"
 #include "TH3D.h"
 #include "TCanvas.h"
 #include "TLegend.h"
 #include "TFile.h"
 #include "TClass.h"
 #include "TLatex.h"
 #include "TPaveStats.h"
 #include "TStyle.h"
 #include "TPad.h"
 
 #include "TCollection.h"
 #include "TList.h"
 #include "TF1.h"

 #include "TROOT.h"
 #include "TMath.h"

 #include "MHists.h"
 
#endif




MHists::MHists(): debugl(0)
{
  const Int_t nclrs = 11;
  Int_t clrs[nclrs] = {2, 4, 3, 6, 12, 30, 36, 38, 39, 40, 46};
  for (Int_t ii = 0; ii < nclrs; ++ii) clrarr.push_back(clrs[ii]);
}


MHists::~MHists()
{
  for (std::map< std::string, std::vector <TH1*> >::iterator hmitr = fhist_map.begin(); hmitr != fhist_map.end(); ++hmitr) {
    for (std::vector <TH1*>::iterator hvitr = hmitr->second.begin(); hvitr != hmitr->second.end(); ++hvitr) {
      if (*hvitr) delete *hvitr;
    }
    hmitr->second.clear();
  }
  fhist_map.clear();

  for (std::map< std::string, std::vector <TLegend*> >::iterator hmitr = flgnd_map.begin(); hmitr != flgnd_map.end(); ++hmitr) {
    for (std::vector <TLegend*>::iterator hvitr = hmitr->second.begin(); hvitr != hmitr->second.end(); ++hvitr) {
      if (*hvitr) delete *hvitr;
    }
    hmitr->second.clear();
  }
  flgnd_map.clear();
}



Int_t MHists::SaveHists (const std::string &fname)
{
  if (fname.empty()) {
    std::cout << "Warning: MHists::SaveHists: File name is empty. Return.\n";
    return 0;
  }
  TFile *pfile = TFile::Open(fname.c_str(), "RECREATE");
  if (!pfile) {
    std::cout << "Warning: MHists::SaveHists: Can not open file " << fname << "! Return!\n";
    return 0;
  }
  TH1D  *hinfo = new TH1D("fhinfo", "fhinfo", fhist_map.size(), 0, fhist_map.size());
  Int_t nh = 0;
  for (std::map< std::string, std::vector <TH1*> >::iterator hmitr = fhist_map.begin(); hmitr != fhist_map.end(); ++hmitr) {
    hinfo->Fill(hmitr->first.c_str(), hmitr->second.size());
    for (std::vector <TH1*>::iterator hvitr = hmitr->second.begin(); hvitr != hmitr->second.end(); ++hvitr) {
      if (*hvitr) (*hvitr)->Write();
    }
    nh+=hmitr->second.size();
  }
  hinfo->Write();
  pfile->Close();
  delete pfile;
  return nh;
}



Int_t MHists::SaveHists (const std::vector<std::string> hmapids, const std::string &fname)
{
  if (fname.empty()) {
    std::cout << "Warning: MHists::SaveHists: File name is empty. Return.\n";
    return 0;
  }
  TFile *pfile = TFile::Open(fname.c_str(), "RECREATE");
  if (!pfile) {
    std::cout << "Warning: MHists::SaveHists: Can not open file " << fname << "! Return!\n";
    return 0;
  }
  TH1D  *hinfo = new TH1D("fhinfo", "fhinfo", hmapids.size(), 0, hmapids.size());
  Int_t nh = 0;
  for (std::vector<std::string>::const_iterator mnitr = hmapids.begin(); mnitr != hmapids.end(); ++mnitr) {
    std::map< std::string, std::vector <TH1*> >::iterator hmitr = fhist_map.find(*mnitr);
    if (hmitr == fhist_map.end()) continue;
    hinfo->Fill(hmitr->first.c_str(), hmitr->second.size());
    for (std::vector <TH1*>::iterator hvitr = hmitr->second.begin(); hvitr != hmitr->second.end(); ++hvitr) {
      if (*hvitr) (*hvitr)->Write();
    }
    nh+=hmitr->second.size();
  }
  hinfo->Write();
  pfile->Close();
  delete pfile;
  return nh;
}



Int_t MHists::LoadFromFile(const std::string &fname)
{
  std::string infoh_name("fhinfo");
  if (fname.empty()) {
    std::cout << "Warning: MHists::LoadFromFile: File name is empty. Return.\n";
    return 0;
  }
  TFile *pfile = TFile::Open(fname.c_str(), "READ");
  if (!pfile) {
    std::cout << "Warning: MHists::LoadFromFile: Can not open file " << fname << "! Return!\n";
    return 0;
  }
  TH1D  *hinfo = dynamic_cast<TH1D*>(pfile->Get(infoh_name.c_str()));
  if (!hinfo) {
    std::cout << "Warning: MHists::LoadFromFile: Can not read " << infoh_name << " historgam from the file " << fname << "! Return!\n";
    pfile->Close();
    delete pfile;
    return 0;
  }
  Int_t nhread = 0;
  Int_t map_size = hinfo->GetNbinsX();
  for (Int_t binid = 0; binid < map_size; ++binid) {
    std::string mapid = hinfo->GetXaxis()->GetBinLabel(binid+1);
    Int_t n_hists = hinfo->GetBinContent(binid+1);
    if (!AddHists(mapid, n_hists)) continue;
    for (Int_t hid = 0; hid < n_hists; ++hid) {
      std::stringstream hh_name("");
      hh_name << mapid << "_" << hid;
      TH1 *hh = dynamic_cast<TH1*>(pfile->Get(hh_name.str().c_str()));
      if (!hh) {
        std::cout << "Warning: MHists::LoadFromFile: Can not read " << hh_name.str() << " historgam from the file " << fname << "! Continue!\n";
        continue;
      }
      hh->SetDirectory(0);
      AddHists(mapid, hh, hid);
      ++nhread;
    }
  }
  pfile->Close();
  delete pfile;
  return nhread;
}



Int_t MHists::AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1)
{ // 1D
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH1D(hname.str().c_str(), hname.str().c_str(), n1, a1, b1);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
    if (debugl > 1) std::cout << "Histogram " << hname.str() << "  " << " has been created.\n";
  }
  return nhists - hfail;
}



Int_t MHists::AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1,  
                                                                     Int_t n2, Double_t a2, Double_t b2)
{ // 2D
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH2D(hname.str().c_str(), hname.str().c_str(), n1, a1, b1, n2, a2, b2);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
  }
  return nhists - hfail;
}



Int_t MHists::AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1,  
                                                                     Int_t n2, Double_t a2, Double_t b2,
                                                                     Int_t n3, Double_t a3, Double_t b3)
{ // 3D
  if (nhists < 0) return 0;
  int hfail = 0;
  for (int ii = 0; ii < nhists; ++ii) {
    std::stringstream  hname("");
    hname << hmapid << "_" << ii;
    TH1 *hh = new TH3D(hname.str().c_str(), hname.str().c_str(), n1, a1, b1, n2, a2, b2, n3, a3, b3);
    if (hh) {
      hh->SetDirectory(0);
      hh->Sumw2();
      fhist_map[hmapid].push_back(hh);
    } else  { ++hfail; std::cerr << "Error:  MHists::AddHists: Problem to create histogram: " << hname.str() << ". Continue.\n"; }
  }
  return nhists - hfail;
}



Int_t MHists::AddHists(const std::string hmapid, const Int_t nhists)
{
  if (nhists < 0) return 0;
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    std::cout << "Warning:  MHists::AddHists: " << hmapid << " already exists. Nothing was done.\n";
    return 0;
  }
  fhist_map[hmapid].reserve(nhists);
  fhist_map[hmapid].assign(nhists, (TH1D*)0);
  return nhists;
}



Int_t MHists::AddHists(const std::string hmapid, TH1 *hist, const Int_t posid)
{
  if (posid < 0 || !hist) return 0;
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    if ( static_cast<unsigned int>(posid) < fhist_map[hmapid].size() ) {
      fhist_map[hmapid][posid] = hist;
      return 1;
    }
  } else {
    std::cout << "Warning:  MHists::AddHists: Cannnot assign histogram, collection " << hmapid << " does not exist. Nothing was done.\n";
  }
  return 0;
}




TH1* MHists::GetHist(const std::string hmapid, Int_t histid)
{
  std::map< std::string, std::vector <TH1*> >::iterator  hitr = fhist_map.find(hmapid);
  if ( hitr != fhist_map.end() ) {
    if (histid >= 0 && hitr->second.size() > static_cast<unsigned int>(histid))
      return hitr->second[histid];
  } 
  return 0;  
}


// TH3* MHists::GetHist(const std::string hmapid, Int_t histid, string hClass)
// {
//   hClass = "TH3";
//   std::map< std::string, std::vector <TH3*> >::iterator  hitr = fhist_map.find(hmapid);
//   if ( hitr != fhist_map.end() ) {
//     if (histid >= 0 && hitr->second.size() > static_cast<unsigned int>(histid))
//       return hitr->second[histid];
//   } 
//   return 0;  
// }



void MHists::FillHist(const std::string hmapid, Int_t histid, Double_t val) 
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH1") != std::string::npos) hh->Fill(val);
    else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 1D. Continue.\n";
  } else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}

void MHists::SetBinHist(const std::string hmapid, Int_t histid, Int_t binVal, Double_t val) 
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH1") != std::string::npos) hh->SetBinContent(binVal, val);
    else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 1D. Continue.\n";
  } else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}

void MHists::FillHist(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2)
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH2") != std::string::npos) hh->Fill(val1, val2);
    else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 2D. Continue.\n";
  } else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}



void MHists::FillHistW(const std::string hmapid, Int_t histid, Double_t val, Double_t weight) 
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH1") != std::string::npos) hh->Fill(val, weight);
    else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 1D. Continue.\n";
  } else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}


void MHists::FillHistW(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2, Double_t weight)
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH2") != std::string::npos) dynamic_cast<TH2*>(hh)->Fill(val1, val2, weight);
    else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 2D. Continue.\n";
  } else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}


void MHists::FillHistW(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2, Double_t val3, Double_t weight)
{ 
  TH1 *hh = GetHist(hmapid, histid);
  if (hh) { 
    std::string cln(hh->ClassName());
    if (cln.find("TH3") != std::string::npos) dynamic_cast<TH3*>(hh)->Fill(val1, val2, val3, weight);
    else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] is not 3D. Continue.\n";
  } else  std::cerr << "Error:  MHists::FillHist: Histogram: " << hmapid << "[" << histid << "] not found. Continue.\n";
}


void MHists::CanvasDesign(Int_t nhist, Int_t &rows, Int_t &col)
{
  if (nhist <= 0) rows = col = 0;
  else {
    rows = static_cast<Int_t>(sqrt(nhist));
    col = static_cast<Int_t>(nhist/rows);
    if (col*rows < nhist) ++col;
  }
  if (debugl > 5) std::cout << "MHists::CanvasDesign:  N_Hists: " << nhist << "   N_Rows: " << rows  << "   N_Clns: " << col << std::endl; 
}



TCanvas *MHists::DrawHist2D_BT15(const std::string hmapid, const char *canvas_name, const char *x_name, const char *y_name, 
                                                                        const char *z_name, const char *drawopt) 
{
  std::vector<std::string> titles;
  if (canvas_name) titles.push_back(canvas_name);
  else titles.push_back("");
  if (x_name) titles.push_back(x_name);
  else titles.push_back("");
  if (y_name) titles.push_back(y_name);
  else titles.push_back("");
  if (z_name) titles.push_back(z_name);
  else titles.push_back("");
  if (drawopt) titles.push_back(drawopt);
  return DrawHist2D_BT15(hmapid, titles);
}



TCanvas *MHists::DrawHist2D_BT15(const std::string hmapid, const std::vector<std::string> &titles)
{
//  titles: canvas, x_axis, y_axis, z_axis 
  TH2D  *hh;
  Int_t nrow, ncol, hii; 
  std::string c_title;

  if (fhist_map.find(hmapid) == fhist_map.end()) {
    std::cout << "Error MHists::DrawHist2D_BT15:  Key string " << hmapid << " is wrong\n";
    return 0;
  }

  if ( titles.size() > 0  && !titles[0].empty()) c_title = titles[0];
  else c_title = hmapid;
  TCanvas *c5 = new TCanvas(c_title.c_str(), c_title.c_str(), 1800, 600);
  CanvasDesign(fhist_map[hmapid].size(), nrow, ncol);
  
  c5->Divide(ncol, nrow);
  Double_t zmax, zmin;
  zmax = zmin = 0.0;
  hii = 0;
  for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr, ++hii) {
    if (fhist_map[hmapid].size()==7) { int npos = (hii&1) ? (hii>>1)+1 : (hii>>1)+5; c5->cd(npos); }
    else c5->cd(hii+1);
    hh = dynamic_cast<TH2D*>(*itr);
    if (!hh) continue;
    if (titles.size() > 1) hh->GetXaxis()->SetTitle(titles[1].c_str());
    if (titles.size() > 2) hh->GetYaxis()->SetTitle(titles[2].c_str());
    if (titles.size() > 3) hh->GetZaxis()->SetTitle(titles[3].c_str());
    if (titles.size() > 4) hh->Draw(titles[4].c_str());
    else hh->Draw("lego2");
    
    if (zmax < hh->GetBinContent(hh->GetMaximumBin()) ) zmax = hh->GetBinContent(hh->GetMaximumBin());
    if (zmin > hh->GetBinContent(hh->GetMinimumBin()) ) zmin = hh->GetBinContent(hh->GetMinimumBin());
  }
  
  for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
    hh = dynamic_cast<TH2D*>(*itr);
    if (!hh) continue;
    hh->GetZaxis()->SetRangeUser(zmin, zmax);
  }
  return c5;
}


TCanvas *MHists::DrawHist1D_BT15(const std::string hmapid, const char *canvas_name, const char *x_name, const char *y_name, const char *dopt) 
{
  std::vector<std::string> titles;
  if (canvas_name) titles.push_back(canvas_name);
  else titles.push_back("");
  if (x_name) titles.push_back(x_name);
  else titles.push_back("");
  if (y_name) titles.push_back(y_name);
  else titles.push_back("");
  if (dopt) titles.push_back(dopt);
  return DrawHist1D_BT15(hmapid, titles);
}


TCanvas *MHists::DrawHist1D_BT15(const std::string hmapid, const std::vector<std::string> &titles)
{ //  titles: canvas, x_axis, y_axis, draw_opt 
  // Creats a canvas and divide it in pads and draws 1D histograms one in each pad.
  TLegend   *lgnd;
  TH1D      *hh;
  TList     *hfnlist;
  Int_t nrow, ncol, hii, nclrs; 
  std::string c_title, draw_opt;

  if (fhist_map.find(hmapid) == fhist_map.end()) {
    std::cout << "Error MHists::DrawHist1D_BT15:  Key string " << hmapid << " is wrong\n";
    return 0;
  }

  flgnd_map[hmapid].clear();
  if ( titles.size() > 0  && !titles[0].empty()) c_title = titles[0];
  else c_title = hmapid;
  TCanvas   *c4 = new TCanvas(c_title.c_str(), c_title.c_str(), 900, 600);
  CanvasDesign(fhist_map[hmapid].size(), nrow, ncol);

  c4->Divide(ncol, nrow);
  nclrs = clrarr.size();
  if (titles.size() > 3  && !titles[3].empty() ) draw_opt = titles[3];
  else draw_opt = "hist";
  
  hii = 0;
  for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr, ++hii) {
    if (fhist_map[hmapid].size()==7) { int npos = (hii&1) ? (hii>>1)+1 : (hii>>1)+5; c4->cd(npos); }
    else c4->cd(hii+1);
    hh = dynamic_cast<TH1D*>(*itr);
    if (!hh) continue;
    lgnd = new TLegend(0.75, 0.55, 0.88, 0.88);
    hh->SetStats(1);
    hh->SetMarkerStyle((hii%nclrs)+22);
    hh->SetLineColor(clrarr[hii%nclrs]); 
    hh->SetMarkerColor(clrarr[hii%nclrs]);
    hh->SetMarkerSize(1.5);
    hh->SetLineWidth(2);
    if (titles.size() > 1) hh->GetXaxis()->SetTitle(titles[1].c_str());
    if (titles.size() > 2) hh->GetYaxis()->SetTitle(titles[2].c_str());
    hh->Draw(draw_opt.c_str());
    lgnd->AddEntry(hh, hh->GetName(), (draw_opt.find_first_of("pP") == std::string::npos) ? "L":"LP");
    if ( (hfnlist = hh->GetListOfFunctions()) && (hfnlist->GetSize() > 0) ) {
      TIter next(hfnlist);
      TF1 *fn;
      while ( (fn = dynamic_cast<TF1*>(next())) ) {
        fn->Draw("same");
        lgnd->AddEntry(fn, fn->GetName(), "l");
      }
    }
//hh->Print(/*"ALL"*/);
//  htempl.front()->GetYaxis()->SetRangeUser(0.0, 1.0);
//  lgnd->SetTextSize(0.0275);
    lgnd->SetBorderSize(0);
    lgnd->SetFillColor(0);
    lgnd->SetTextFont(42);
    lgnd->SetLineColor(1);
    lgnd->SetFillStyle(1001);
    lgnd->Draw();
    flgnd_map[hmapid].push_back(lgnd);
  }
  return c4;
}



TCanvas *MHists::DrawHist1D(const std::string hmapid, const std::vector<std::string> &titles)
{ //  titles: canvas, x_axis, y_axis, draw_opt. 
  // Draws 1D histograms in current, all in one.
  TLegend   *lgnd;
  TH1D      *hh;
  Int_t hii, nclrs; 
  std::string c_title, draw_opt;
  TCanvas   *c4 = 0;

  if (fhist_map.find(hmapid) == fhist_map.end()) {
    std::cout << "Error MHists::DrawHist1D:  Key string " << hmapid << " is wrong\n";
    return 0;
  }
  
  if ( titles.size() > 0  && !titles[0].empty()) c_title = titles[0];
  else c_title = hmapid;
  TObject *cobj = gROOT->GetListOfCanvases()->FindObject(c_title.c_str());
  if (cobj /*&&  std::string(cobj->ClassName()).find("TCanvas") !=  std::string::npos */ ) {
    c4 = dynamic_cast<TCanvas*>(cobj);
    if (debugl > 5) std::cout <<  "Found canvas " << c_title << std::endl;
  } else {
    c4 = new TCanvas(c_title.c_str(), c_title.c_str(), 900, 600);
  }

//   flgnd_map[hmapid].clear();
  lgnd = new TLegend(0.75, 0.55, 0.88, 0.88);

  nclrs = clrarr.size();
  if (titles.size() > 3  && !titles[3].empty() ) draw_opt = titles[3];
  else draw_opt = "hist";
  
  hii = 0;
  for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr, ++hii) {
    hh = dynamic_cast<TH1D*>(*itr);
    if (!hh) continue;
    hh->SetStats(0);
    hh->SetMarkerStyle((hii%nclrs)+22);
    hh->SetLineColor(clrarr[hii%nclrs]); 
    hh->SetMarkerColor(clrarr[hii%nclrs]);
    hh->SetMarkerSize(1.5);
    hh->SetLineWidth(2);
    if (!hii) {
      hh->Draw(draw_opt.c_str());
      if (titles.size() > 1 && !titles[1].empty() ) hh->GetXaxis()->SetTitle(titles[1].c_str());
      if (titles.size() > 2 && !titles[2].empty()) hh->GetYaxis()->SetTitle(titles[2].c_str());
    } else hh->Draw(Form("%s same", draw_opt.c_str()));
    lgnd->AddEntry(hh, hh->GetName(), "LP");
//     hh->Print(/*"ALL"*/);
  }
//   htempl.front()->GetYaxis()->SetRangeUser(0.0, 1.0);
//   lgnd->SetTextSize(0.0275);
  lgnd->SetBorderSize(0);
  lgnd->SetFillColor(0);
  lgnd->SetTextFont(42);
  lgnd->SetLineColor(1);
  lgnd->SetFillStyle(1001);
  lgnd->Draw();
  flgnd_map[hmapid].push_back(lgnd);
  return c4;
}



TCanvas *MHists::DrawHist1D(const std::vector<std::string> hmapids, const std::vector<std::string> &titles)
{
//  titles: canvas, x_axis, y_axis, draw_opt 
  Int_t nrow, ncol, hii; 
  std::vector<std::string> loc_titles = titles;
  
  if (hmapids.size()==0) return 0;
  
  std::string c_title("");
  if ( loc_titles.size() > 0 && !loc_titles[0].empty() ) c_title = loc_titles[0];
  else {
    for (std::vector <std::string>::const_iterator sitr = hmapids.begin(); sitr != hmapids.end(); ++sitr) {
      if ( c_title.empty() ) c_title += *sitr;
      else c_title += ", " + *sitr; 
    }
    if ( loc_titles.size() > 0) loc_titles[0] = c_title; else loc_titles.push_back(c_title);
  }
  
  TCanvas   *c4 = new TCanvas(c_title.c_str(), c_title.c_str(), 900, 600);
  CanvasDesign(hmapids.size(), nrow, ncol);
  c4->Divide(ncol, nrow);

  hii = 0;
  for (std::vector <std::string>::const_iterator sitr = hmapids.begin(); sitr != hmapids.end(); ++sitr, ++hii) {
    c4->cd(hii+1);
    DrawHist1D(*sitr, loc_titles);
  }
  return c4;
}



TCanvas *MHists::DrawHist1DCmp(const std::vector<std::string> hmapids, const char *canvas_name, const char *x_name, 
                                                                       const char *y_name, const char *dopt)
{
  std::vector<std::string> titles;
  if (canvas_name) titles.push_back(canvas_name);
  else titles.push_back("");
  if (x_name) titles.push_back(x_name);
  else titles.push_back("");
  if (y_name) titles.push_back(y_name);
  else titles.push_back("");
  if (dopt) titles.push_back(dopt);
  return DrawHist1DCmp(hmapids, titles);
}



TCanvas *MHists::DrawHist1DCmp(const std::vector<std::string> hmapids, const std::vector<std::string> &titles)
{
  Int_t nclrs, ncol, nrow;
  TLegend   *lgnd;
  std::string draw_opt;
  if (hmapids.size() <= 0) return 0;
  unsigned int hii, maxmap = 0;
  for (std::vector<std::string>::const_iterator itr = hmapids.begin(); itr != hmapids.end(); ++itr) {
    if (fhist_map.find(*itr) != fhist_map.end()) {
      if (debugl > 1) std::cout << "Size of the vector assignes to " << *itr << ":  " << fhist_map[*itr].size() << std::endl;
      if (fhist_map[*itr].size() > maxmap) maxmap = fhist_map[*itr].size();
    }
  }
  
  std::string c_title("");
  if ( titles.size() > 0 && !titles[0].empty() ) c_title = titles[0];
  else {
    for (std::vector <std::string>::const_iterator sitr = hmapids.begin(); sitr != hmapids.end(); ++sitr) {
      if ( c_title.empty() ) c_title += *sitr;
      else c_title += ", " + *sitr; 
    }
  }

  TCanvas   *c4 = new TCanvas(c_title.c_str(), c_title.c_str(), 900, 600);
  CanvasDesign(maxmap, nrow, ncol);
  c4->Divide(ncol, nrow);

  nclrs = clrarr.size();
  if (titles.size() > 3  && !titles[3].empty() ) draw_opt = titles[3];
  else draw_opt = "hist";

  hii = 0;
  for (hii = 0; hii < maxmap; ++hii) {
    if (maxmap==7) { int npos = (hii&1) ? (hii>>1)+1 : (hii>>1)+5; c4->cd(npos); }
    else c4->cd(hii+1);
    Int_t hdrwn = 0;
    for (std::vector <std::string>::const_iterator sitr = hmapids.begin(); sitr != hmapids.end(); ++sitr) {
      if (fhist_map.find(*sitr) == fhist_map.end()) {
        std::cout << "Warning MHists::DrawHist1DCmp:  Key string " << *sitr << " not found. Continue.\n";
        continue;
      }
      if (fhist_map[*sitr].size() == 0) continue;
      TH1D *hh = dynamic_cast<TH1D*>(fhist_map[*sitr][hii % fhist_map[*sitr].size()]);
      if (!hh) continue;
      hh->SetStats(0);
      hh->SetMarkerStyle((hdrwn%nclrs)+22);
      hh->SetLineColor(clrarr[hdrwn%nclrs]); 
      hh->SetMarkerColor(clrarr[hdrwn%nclrs]);
//       hh->SetMarkerSize(1.5);
      hh->SetLineWidth(2);
      if (!hdrwn) {
        hh->Draw(draw_opt.c_str());
        if (titles.size() > 1 && !titles[1].empty() ) hh->GetXaxis()->SetTitle(titles[1].c_str());
        if (titles.size() > 2 && !titles[2].empty()) hh->GetYaxis()->SetTitle(titles[2].c_str());
        lgnd = new TLegend(0.75, 0.55, 0.88, 0.88);
      } else hh->Draw(Form("%s same", draw_opt.c_str()));
      lgnd->AddEntry(hh, hh->GetName(), "LP");
      ++hdrwn;
    }
    lgnd->SetBorderSize(0);
    lgnd->SetFillColor(0);
    lgnd->SetTextFont(42);
    lgnd->SetLineColor(1);
    lgnd->SetFillStyle(1001);
    lgnd->Draw();
    flgnd_map[hmapids.front()].push_back(lgnd);
  }
  return c4;
}



void MHists::SetRange(const std::string hmapid, const Double_t val1, const Double_t val2, const Int_t axis)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      switch (axis) {
        case 0: hh->GetXaxis()->SetRangeUser(val1, val2); break;
        case 1: hh->GetYaxis()->SetRangeUser(val1, val2); break;
        case 2: hh->GetZaxis()->SetRangeUser(val1, val2); break;
      }
    }
  }
}


void MHists::SetLineColor(const std::string hmapid, const Int_t clr)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      hh->SetLineColor(clr);
      hh->SetMarkerColor(clr);
    }
  }
}


void MHists::SetLineWidth(const std::string hmapid, const Int_t wdth)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      hh->SetLineWidth(wdth);
    }
  }
}


void MHists::SetMarkerSize(const std::string hmapid, const Double_t msize)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1D *hh = dynamic_cast<TH1D*>(*itr);
      if (!hh) continue;
      hh->SetMarkerSize(msize);
    }
  }
}


void MHists::Scale(const std::string hmapid, const Double_t x)
{
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1 *hh = dynamic_cast<TH1*>(*itr);
      if (!hh) continue;
      hh->Scale(x);
    }
  }
}


double MHists::Integral(const std::string hmapid)
{ 
  double integral = 0;
  if (fhist_map.find(hmapid) != fhist_map.end()) {
    for (std::vector <TH1*>::iterator itr = fhist_map[hmapid].begin(); itr != fhist_map[hmapid].end(); ++itr) {
      TH1 *hh = dynamic_cast<TH1*>(*itr);
      if (!hh) continue;
      integral = hh->Integral();
    }
  }
  return integral;
}


int MHists::SetLog (TCanvas *cc, const int axis)
{
  TPad *pp;
  if (!cc) return 0;
  int nn = 1;
  while ( (pp = dynamic_cast<TPad*>(cc->GetPad(nn))) ) {
    switch (axis) {
      case 0 : pp->SetLogx(); break;
      case 1 : pp->SetLogy(); break;
      case 2 : pp->SetLogz(); break;
    }
    ++nn;
  }
  return nn;
}


Int_t MHist_test_v1(const Int_t n_hist = 4)
{
  if (n_hist <= 0 ) { 
    std::cout << "Number of histograms is " << n_hist << ". Nothing to do.\n"; 
    return 0; 
  }
  MHists *prj_hists = new MHists();
  prj_hists->SetDebug(0);
  prj_hists->AddHists("gaus_test", n_hist, 100, -10.0, 10.0);
  prj_hists->AddHists("gaus_test2", n_hist, 100, -10.0, 10.0);
  
  for (Int_t ii = 0; ii < n_hist; ++ii) {
    TH1D *hh = dynamic_cast<TH1D*>(prj_hists->GetHist("gaus_test", ii));
    if (hh) hh->FillRandom("gaus", 10000);
    else std::cout << "MHists::GetHist returned 0\n";
  }
  
  for (Int_t ii = 0; ii < n_hist; ++ii) {
    TH1D *hh = dynamic_cast<TH1D*>(prj_hists->GetHist("gaus_test2", ii));
    if (hh) hh->FillRandom("gaus", 10000);
    else std::cout << "MHists::GetHist returned 0\n";
  }
  
  std::string hhtitle[] = {"test_gaus", "data", "dN", "hist"};
  std::vector<std::string> hhstr;
  for (Int_t ii = 0; ii < 4; ++ii) hhstr.push_back(hhtitle[ii]);

  prj_hists->DrawHist1D_BT15("gaus_test", hhstr);
  
  std::vector< std::string > v_str, v_str1;
  v_str.push_back("gaus_test");
  v_str.push_back("gaus_test2");
  prj_hists->DrawHist1D(v_str, v_str1);
  hhstr[0] = "test_gaus2";
  hhstr[3] = "e1p";
  prj_hists->DrawHist1D("gaus_test2", hhstr);
  
//   prj_hists->SaveHists("test_gaus_12.root");
  
  return 0;
}
  

  


