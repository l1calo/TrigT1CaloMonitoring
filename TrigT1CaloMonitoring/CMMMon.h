// ********************************************************************
//
// NAME:        CMMMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on CMM level
//
// ********************************************************************

#ifndef CMMMon_H
#define CMMMon_H

#include <string>

#include "GaudiKernel/ToolHandle.h"

#include "DataModel/DataVector.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class TH1F_LW;
class TH2F_LW;
class TH2I_LW;
class TProfile2D_LW;

class StatusCode;

namespace LVL1 {
  class CMMJetHits;
  class CMMEtSums;
}

class TrigT1CaloMonErrorTool;
class TrigT1CaloLWHistogramTool;

class CMMMon : public ManagedMonitorToolBase
{
public:

   typedef DataVector<LVL1::CMMJetHits> CMMJetHitsCollection;
   typedef DataVector<LVL1::CMMEtSums>  CMMEtSumsCollection;
  
   CMMMon( const std::string & type, const std::string & name,
	                             const IInterface* parent ); 

   virtual ~CMMMon();

   virtual StatusCode initialize();
   virtual StatusCode bookHistograms( bool isNewEventsBlock,
                                      bool isNewLumiBlock, bool isNewRun );
   virtual StatusCode fillHistograms();
   virtual StatusCode procHistograms( bool isEndOfEventsBlock,
                                      bool isEndOfLumiBlock, bool isEndOfRun );

private:

   enum SummaryErrors { JetStatus, EnergyStatus, JetParity, EnergyParity,
                        RoIParity, NumberOfSummaryBins };
   enum RoIParityErrors { ExParity, EyParity, EtParity, JetEtParity,
                          NumberOfRoIParityBins };

   // Tool to retrieve bytestream errors
   ToolHandle<TrigT1CaloMonErrorTool>    m_errorTool;
   ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;

   /** location of data */
   std::string m_CMMJetHitsLocation;
   std::string m_CMMEtSumsLocation;
   std::string m_CMMRoILocation;

   std::string m_PathInRootFile;   
   std::string m_ErrorPathInRootFile;
   bool m_histBooked;

   /** Histos */   
   // CMM Jet Hits
   TH1F_LW* m_h_CMMJetHits_MainJets;
   TH1F_LW* m_h_CMMJetHits_FwdJetsRight;
   TH1F_LW* m_h_CMMJetHits_FwdJetsLeft;
   TH1F_LW* m_h_CMMJetHits_EtMap;
   // JEM Hits
   TH1F_LW* m_h_CMMJetHits_JEM_MainHits;
   TH1F_LW* m_h_CMMJetHits_JEM_FwdHitsRight;
   TH1F_LW* m_h_CMMJetHits_JEM_FwdHitsLeft;
   TH1F_LW* m_h_CMMJetHits_JEM_Crate0ParityError;
   TH1F_LW* m_h_CMMJetHits_JEM_Crate1ParityError;

   // CMM Et Sums
   TH1F_LW* m_h_CMMEtSums_Ex;
   TH1F_LW* m_h_CMMEtSums_Ey;
   TH1F_LW* m_h_CMMEtSums_Et;
   TH1F_LW* m_h_CMMEtSums_MissingEtMap;
   TH1F_LW* m_h_CMMEtSums_SumEtMap;
   TH1F_LW* m_h_CMMEtSums_MissingEtSigMap;
   TProfile2D_LW* m_h_CMMEtSums_Overflow;
   // JEM Et Sums
   TH1F_LW*  m_h_CMMEtSums_JEM_Ex;
   TH1F_LW*  m_h_CMMEtSums_JEM_Ey;
   TH1F_LW*  m_h_CMMEtSums_JEM_Et; 

   // CMM RoI
   TH1F_LW* m_h_CMMRoI_JetEtHits;
   TH1F_LW* m_h_CMMRoI_SumEtHits;
   TH1F_LW* m_h_CMMRoI_MissingEtHits;
   TH1F_LW* m_h_CMMRoI_MissingEtSigHits;

   TH1F_LW* m_h_CMMRoI_Ex;
   TH1F_LW* m_h_CMMRoI_Ey;
   TH1F_LW* m_h_CMMRoI_Et;

   //errors
   TH2F_LW* m_h_CMMJet_error;
   TH2F_LW* m_h_CMMEnergy_error;
   TH2F_LW* m_h_CMMJet_parity;
   TH2F_LW* m_h_CMMEnergy_parity;
   TH1F_LW* m_h_CMMRoI_error;
   TH1F_LW* m_h_CMM_ErrorSummary;
   TH1F_LW* m_h_TriggeredSlice;	  
   TH2I_LW* m_h_CMM_Events;
};


#endif
