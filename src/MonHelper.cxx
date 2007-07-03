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
  std::string name,title;
  if (m_DataType=="")
    {
      name=HistoName;
      title=HistoTitle;
    }
  else
    {
      name=m_DataType + "_" + HistoName;
      title=m_DataType + ": " + HistoTitle;
    }

  TH1F* hist = new TH1F(TString(name), TString(title), NoBins, xmin, xmax);
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
  std::string name,title;
  if (m_DataType=="")
    {
      name=HistoName;
      title=HistoTitle;
    }
  else
    {
      name=m_DataType + "_" + HistoName;
      title=m_DataType + ": " + HistoTitle;
    }

  TH2F* hist = new TH2F(TString(name), TString(title), xBins, xmin, xmax, yBins, ymin, ymax);
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
  NoTTEtaBins=67;
  NoTTPhiBins=65;
  TTEtaBins = new double[NoTTEtaBins];
  TTPhiBins = new double[NoTTPhiBins];

  TTEtaBins[0]=-4.9;
  TTEtaBins[1]=-4.475;
  TTEtaBins[2]=-4.050;
  TTEtaBins[3]=-3.625;
  TTEtaBins[4]=-3.2;
  
  BinContent=-3.1;
  for (int i=5; i<8; i++)
    {
      TTEtaBins[i]=BinContent;
      BinContent=BinContent+0.2;
    }
  for (int i=8; i<58; i++)
    {
      TTEtaBins[i]=BinContent;
      BinContent=BinContent+0.1;
    }
  for (int i=58; i<62; i++)
    {
      TTEtaBins[i]=BinContent;
      BinContent=BinContent+0.2;
    }

  TTEtaBins[62]=3.2;
  TTEtaBins[63]=3.625;
  TTEtaBins[64]=4.050;
  TTEtaBins[65]=4.475;
  TTEtaBins[66]=4.9;


  
  BinContent=0;
  for (int i=0; i<65; i++)
    {
      TTPhiBins[i]=BinContent;
      BinContent=BinContent+2*TMath::Pi()/64;
    } 


  NoJEEtaBins=33;
  NoJEPhiBins=33;
  JEEtaBins = new double[NoJEEtaBins];
  JEPhiBins = new double[NoJEPhiBins];

  BinContent=-2.4;

  JEEtaBins[0]=-4.9;
  JEEtaBins[1]=-3.2;
  JEEtaBins[2]=-2.9;
  JEEtaBins[3]=-2.7;
  
  JEEtaBins[29]=2.7;
  JEEtaBins[30]=2.9;
  JEEtaBins[31]=3.2;
  JEEtaBins[32]=4.9;

  for (int i=4; i<29; i++)
    {
      JEEtaBins[i]=BinContent;
      BinContent=BinContent+0.2;
    }
  
  BinContent=0;
  for (int i=0; i<33; i++)
    {
      JEPhiBins[i]=BinContent;
      BinContent=BinContent+2*TMath::Pi()/32;
    } 
}

/*---------------------------------------------------------*/
Helper::~Helper()
/*---------------------------------------------------------*/
{
}

/*---------------------------------------------------------*/
double* Helper::TTEtaBinning() 
/*---------------------------------------------------------*/
{
  return TTEtaBins;
}

/*---------------------------------------------------------*/
double* Helper::TTPhiBinning() 
/*---------------------------------------------------------*/
{
  return TTPhiBins;
}

/*---------------------------------------------------------*/
double* Helper::JEEtaBinning() 
/*---------------------------------------------------------*/
{
  return JEEtaBins;
}

/*---------------------------------------------------------*/
double* Helper::JEPhiBinning() 
/*---------------------------------------------------------*/
{
  return JEPhiBins;
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

/*---------------------------------------------------------*/
void Helper::FillHitsHisto(TH1F* Histo, std::string HitWord, int minThresh, int NoThresh, int ThreshOffset, 
			    int BitsPerThresh, MsgStream::MsgStream* log)
/*---------------------------------------------------------*/
{
  int min=(minThresh+ThreshOffset);
  int max=(NoThresh+ThreshOffset);

  for (int i=min; i<max; i++)
    {
      int JetMult = Multiplicity(HitWord,i,BitsPerThresh);

      *log<<MSG::VERBOSE<<"Thresh.No: "<<(i-ThreshOffset)<<" Multiplicity: "<< JetMult<<endreq;
      Histo->Fill((i-ThreshOffset),JetMult);
    }
  return;
}

