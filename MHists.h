
// This class allows to create and manage histograms of similar type.
// E.g. create 128 identical histograms to study signal of each channel of one APV chip.

#if !defined MHISTS
#define MHISTS

// #if !defined(__CINT__) || defined(__MAKECINT__)
#if !defined(__CLING__) || defined(__ROOTCLING__)  // for root 6

 #include <string>
 #include <vector>
 #include <map>

 #include <Rtypes.h>
 
#endif

class TH1;
class TLegend;
class TCanvas;

class MHists
{
public:  
   MHists();
  virtual ~MHists();
  Int_t AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1); // 1D
  Int_t AddHists(const std::string hmapid, const Int_t nhists, Int_t n1, Double_t a1, Double_t b1,  // 2D
                                                                   Int_t n2, Double_t a2, Double_t b2);
  Int_t AddHists(const std::string hmapid, const Int_t nhists);
  Int_t AddHists(const std::string hmapid, TH1 *hist, const Int_t posid);

  TH1* GetHist(const std::string hmapid, const Int_t histid);
  std::vector<TH1*>& GetHist(const std::string hmapid) { return fhist_map[hmapid]; }
  
  void FillHist(const std::string hmapid, Int_t histid, Double_t val);
  void FillHist(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2);

  void FillHistW(const std::string hmapid, Int_t histid, Double_t val, Double_t weight);
  void FillHistW(const std::string hmapid, Int_t histid, Double_t val1, Double_t val2, Double_t weight);
  
  TCanvas *DrawHist2D_BT15(const std::string hmapid, const std::vector<std::string> &titles);
  TCanvas *DrawHist2D_BT15(const std::string hmapid, const char *canvas_name = 0, const char *x_name = 0, const char *y_name = 0, 
                                                     const char *z_name = 0, const char *drawopt = 0); 
  
  TCanvas *DrawHist1D_BT15(const std::string hmapid, const std::vector<std::string> &titles);
  TCanvas *DrawHist1D_BT15(const std::string hmapid, const char *canvas_name = 0, const char *x_name = 0, const char *y_name = 0, const char *dopt = 0); 

  TCanvas *DrawHist1D(const std::vector<std::string> hmapids, const std::vector<std::string> &titles);  
  TCanvas *DrawHist1D(const std::string hmapid, const std::vector<std::string> &titles);

  TCanvas *DrawHist1DCmp(const std::vector<std::string> hmapids, const std::vector<std::string> &titles);  
  TCanvas *DrawHist1DCmp(const std::vector<std::string> hmapids, const char *canvas_name = 0, const char *x_name = 0, const char *y_name = 0, const char *dopt = 0); 
  
  void SetRange(const std::string hmapid, const Double_t val1, const Double_t val2, const Int_t axis = 0);
  void SetLineColor(const std::string hmapid, const Int_t clr = 0);
  void SetLineWidth(const std::string hmapid, const Int_t wdth = 1);
  void SetMarkerSize(const std::string hmapid, const Double_t msize = 1.0);
  void Scale(const std::string hmapid, const Double_t x);
  static int SetLog (TCanvas *cc, const int axis = 1);

  Int_t SaveHists (const std::string &fname);
  Int_t SaveHists (const std::vector<std::string> hmapids, const std::string &fname);
  Int_t LoadFromFile(const std::string &fname);

  void CanvasDesign(Int_t nhist, Int_t &rows, Int_t &col);
  void SetDebug(const Int_t dbgl = 0) { debugl = dbgl; }
  int GetSize(const std::string hmapid) const {try { return fhist_map.at(hmapid).size(); } catch (const std::out_of_range& rr) { return -1; }} 

protected:

  std::map< std::string, std::vector <TH1*> > fhist_map;
  
  std::map< std::string, std::vector<TLegend*> > flgnd_map;
  std::vector<Int_t> clrarr;
  
  Int_t debugl;

};

#endif





