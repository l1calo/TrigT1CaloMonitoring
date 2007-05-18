// ********************************************************************
//
// NAME:        MonHelpers.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: collection of useful functions that are used by several
//              classes
//
// ********************************************************************

#include <sstream>
#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>

#include "TrigT1CaloMonitoring/MonHelper.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"

/*---------------------------------------------------------*/
HistoBooker::HistoBooker(ManagedMonitorToolBase::MonGroup* MonitoringGroup, 
			 MsgStream::MsgStream* log, std::string DataType) 
/*---------------------------------------------------------*/
{
  m_MonGroup = MonitoringGroup;
  m_DataType = DataType;
  m_log = log;
}

/*---------------------------------------------------------*/
HistoBooker::~HistoBooker()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
TH1F* HistoBooker::book1F(std::string HistoName, std::string HistoTitle, 
	     int NoBins, double xmin, double xmax, std::string xAxisTitle, std::string yAxisTitle)
/*---------------------------------------------------------*/
{
  //set Histo values
  TH1F* hist = new TH1F(TString(m_DataType + "_" + HistoName), TString(m_DataType + ": " + HistoTitle), NoBins, xmin, xmax);
  hist -> GetXaxis() -> SetTitle(TString(xAxisTitle));
  hist -> GetYaxis() -> SetTitle(TString(yAxisTitle));
 
 //book histo to corresponding mongroup
  StatusCode regtest = m_MonGroup->regHist (hist);

  if (regtest != StatusCode::SUCCESS) 
    {
      *m_log << MSG::WARNING << "Could not register histogram : " 
	     << hist->GetName()  << endreq;
    }

  if (regtest = StatusCode::SUCCESS) 
    {
      *m_log << MSG::VERBOSE << "Registered histogram : " 
	     << hist->GetName()  << endreq;
    }
  
  
  return hist;
}


/*---------------------------------------------------------*/
TH2F* HistoBooker::book2F(std::string HistoName, std::string HistoTitle, 
	     int xBins, double xmin, double xmax, int yBins, double ymin, double ymax, 
	     std::string xAxisTitle, std::string yAxisTitle)
/*---------------------------------------------------------*/
{
  //set Histo values
  TH2F* hist = new TH2F(TString(m_DataType + "_" + HistoName), TString(m_DataType + ": " + HistoTitle), xBins, xmin, xmax, yBins, ymin, ymax);
  hist -> GetXaxis() -> SetTitle(TString(xAxisTitle));
  hist -> GetYaxis() -> SetTitle(TString(yAxisTitle));
  hist-> SetOption ("colz");

  //book histo to corresponding mongroup
  StatusCode regtest = m_MonGroup->regHist (hist);

  if (regtest != StatusCode::SUCCESS) 
    {
      *m_log << MSG::WARNING << "Could not register histogram : " 
	     << hist->GetName()  << endreq;
    }

  if (regtest = StatusCode::SUCCESS) 
    {
      *m_log << MSG::VERBOSE << "Registered histogram : " 
	     << hist->GetName()  << endreq;
    }
 
  return hist;
}

/*---------------------------------------------------------*/
Helper::Helper() 
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
Helper::~Helper()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
std::string Helper::Binary(unsigned int Hits, int NumberBits) 
/*---------------------------------------------------------*/
{
  std::string temp = "";
  int temp2=0;

  for (int i=0; i<NumberBits; i++) //change hits into binary form
    //like it is in the bytestream
    {
      temp2= ((Hits >> i) &0x1 );
      //note: the bit-word is ordered right to left:
      //right: lowest thresholds, left: highest thresholds
      //e.g. JEM Hits in Main direction: 
      // thr7|thr6|thr5|thr5|thr4|thr3|thr2|thr0
      if (temp2==0) temp = "0" + temp;
      if (temp2==1) temp = "1" + temp;
    }

  return temp;
}


/*---------------------------------------------------------*/
int Helper::Multiplicity(std::string BinaryHitMap, 
			 int ThresholdNumber, int BitsPerThresh) 
/*---------------------------------------------------------*/
{
  std::string temp = "";
  int mult = 0;

  temp = BinaryHitMap;

  temp.assign(temp,temp.length()-(ThresholdNumber * BitsPerThresh)-BitsPerThresh, BitsPerThresh);
  //assigns to temp only the interesting bits
  //note: the bit-word is ordered right to left:
  //right: lowest thresholds, left: highest thresholds
  //e.g. JEM Hits in Main direction: 
  // thr7|thr6|thr5|thr5|thr4|thr3|thr2|thr0

  for (int i=0; i<BitsPerThresh; i++) // converts the binary
    // display into an integer number
    {
      if ((temp.substr()[i])=='0') mult= mult + int(pow(2,(BitsPerThresh-1-i))) * 0;
      if ((temp.substr()[i])=='1') mult= mult + int(pow(2,(BitsPerThresh-1-i))) * 1;
    }

  return mult;
}
