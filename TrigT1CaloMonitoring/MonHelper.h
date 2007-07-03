// ********************************************************************
//
// NAME:        MonHelper.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: collection of useful functions that are used by several
//              classes
//
// ********************************************************************

#ifndef MonHelper_H
#define MonHelper_H

#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/MsgStream.h"

#include "TH1.h"
#include "TH2.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "AthenaMonitoring/AthenaMonManager.h"

class TH1F;
class TH2F;

class HistoBooker 
{
 public:

     HistoBooker(ManagedMonitorToolBase::MonGroup* MonitoringGroup, MsgStream::MsgStream* log, std::string DataType); 
     virtual ~HistoBooker();

     //books TH1F-Histos to a MonGroup
     TH1F* book1F(std::string HistoName, std::string HistoTitle, 
		  int NoBins, double xmin, double xmax, std::string xAxisTitle, std::string yAxisTitle = "#");

     //books TH2F-Histos to a MonGroup
     TH2F* book2F(std::string HistoName, std::string HistoTitle, 
		  int xBins, double xmin, double xmax, int yBins, double ymin, double ymax, 
		  std::string xAxisTitle, std::string yAxisTitle="#");
 protected:
     ManagedMonitorToolBase::MonGroup* m_MonGroup;
     std::string m_DataType;
     MsgStream::MsgStream* m_log;
};

class Helper 
{
 public:

     Helper(); 
     virtual ~Helper();

     //returns array with eta and phi binning for TT and JE
     double* TTEtaBinning();
     double* TTPhiBinning();
     double* JEEtaBinning();
     double* JEPhiBinning();

     // converts the integer Hits value into a string with lenth = NumberBits  
     std::string Binary(unsigned int Hits, int NumberBits);

     // gives back the multiplicity a certain threshold (ThresNo), where BitsPerThresh
     // is the number of bits that are reserved for each Threshold
     int Multiplicity(std::string BinaryHitMap,int ThreshNo, int BitsPerThresh);

     //fills Hit Histos
     void FillHitsHisto(TH1F* Histo, std::string HitWord, int minThresh, int NoThresh, int ThreshOffset, 
			 int BitsPerThresh, MsgStream::MsgStream* log);
			  
 protected:
     int NoTTEtaBins;
     int NoTTPhiBins;
     double *TTEtaBins;
     double *TTPhiBins;
     int NoJEEtaBins;
     int NoJEPhiBins;
     double *JEEtaBins;
     double *JEPhiBins;

     double BinContent;
};


#endif
