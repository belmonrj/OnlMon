#ifndef CEMC_CEMCMONDRAW_H
#define CEMC_CEMCMONDRAW_H

#include <onlmon/OnlMonDraw.h>

#include <string>  // for allocator, string
#include <TH2.h>
#include <TStyle.h>

class OnlMonDB;
class TCanvas;
class TGraphErrors;
class TPad;

class CemcMonDraw : public OnlMonDraw
{
 public:
  CemcMonDraw(const std::string &name = "CEMCMON"); // same name as server!

  ~CemcMonDraw() override {}

  int Init() override;
  int Draw(const std::string &what = "ALL") override;
  int MakePS(const std::string &what = "ALL") override;
  int MakeHtml(const std::string &what = "ALL") override;
  void setSave(int s) {save = s;}

 protected:
  int MakeCanvas(const std::string &name);
  int DrawFirst(const std::string &what = "ALL");
  int DrawSecond(const std::string &what = "ALL");
  int DrawThird(const std::string &what = "ALL");
  int DrawFourth(const std::string &what = "ALL");
  int DrawHistory(const std::string &what = "ALL");
  int FindHotTower(TPad *warn,TH2D* );
  std::vector<int>  getBadPackets(TH1 *hist, int what, float cutoff);


  int TimeOffsetTicks = -1;
  int save = 0;
  TCanvas *TC[9] = {nullptr};
  TPad *transparent[9] = {nullptr};
  TPad *Pad[18] = {nullptr};
  TPad *warning[18] = {nullptr};
  TGraphErrors *gr[2] = {nullptr};
  OnlMonDB *dbvars = nullptr;
  TStyle* cemcStyle = nullptr;
};

#endif /* CEMC_CEMCMONDRAW_H */