#include "CommonFuncs.C"

#include <onlmon/sepd/SepdMonDraw.h>

#include <onlmon/OnlMonClient.h>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libonlsepdmon_client.so)

void sepdDrawInit(const int online = 0)
{
  OnlMonClient *cl = OnlMonClient::instance();
  // register histos we want with monitor name
  // first histogram uses the TH1->GetName() as key
  cl->registerHisto("h_ADC_corr", "SEPDMON_0");
  cl->registerHisto("h_hits_corr", "SEPDMON_0");
  cl->registerHisto("h_ADC_all_channel", "SEPDMON_0");
  cl->registerHisto("h_hits_all_channel", "SEPDMON_0");
  cl->registerHisto("h_event", "SEPDMON_0");
  cl->registerHisto("h1_waveform_twrAvg", "SEPDMON_0");
  cl->registerHisto("h1_waveform_time", "SEPDMON_0");
  cl->registerHisto("h1_waveform_pedestal", "SEPDMON_0");
  cl->registerHisto("h2_sepd_waveform", "SEPDMON_0");
  cl->registerHisto("h1_packet_chans", "SEPDMON_0");
  cl->registerHisto("h1_packet_length", "SEPDMON_0");
  cl->registerHisto("h1_packet_number", "SEPDMON_0");
  cl->registerHisto("h1_packet_event", "SEPDMON_0");

  for ( int i = 0; i < 768; ++i )
    {
      cl->registerHisto(Form("h_ADC_channel_%d",i), "SEPDMON_0");
    }

  // for local host, just call ohcalDrawInit(2)
  CreateSubsysHostlist("sepd_hosts.list", online);

  // says I know they are all on the same node
  cl->requestHistoBySubSystem("SEPDMON_0", 1);
  OnlMonDraw *sepdmon = new SepdMonDraw("SEPDMONDRAW");  // create Drawing Object
  cl->registerDrawer(sepdmon);                           // register with client framework
}

void sepdDraw(const char *what = "ALL")
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->requestHistoBySubSystem("SEPDMON_0", 1);  // update histos
  cl->Draw("SEPDMONDRAW", what);                // Draw Histos of registered Drawers
}

void sepdSavePlot()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->SavePlot("SEPDMONDRAW");                  // Save Plots
  return;
}

void sepdHtml()
{
  OnlMonClient *cl = OnlMonClient::instance();  // get pointer to framewrk
  cl->MakeHtml("SEPDMONDRAW");                  // Create html output
  return;
}
