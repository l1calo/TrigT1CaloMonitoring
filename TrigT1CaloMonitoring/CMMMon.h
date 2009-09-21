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

#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TH1.h"
#include "TH2.h"

#include "DataModel/DataVector.h"
#include "TrigT1CaloEvent/CMMJetHits.h"
#include "TrigT1CaloEvent/CMMEtSums.h"


#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"


class CMMMon : public ManagedMonitorToolBase
{
public:
/*         typedef DataVector<const LVL1::CMMJetHits*> pCMMJetHitsCollection; */
        typedef DataVector<LVL1::CMMJetHits> CMMJetHitsCollection;
        typedef DataVector<const LVL1::CMMJetHits> cCMMJetHitsCollection;
        typedef DataVector<LVL1::CMMEtSums> CMMEtSumsCollection;
  

	CMMMon( const std::string & type, const std::string & name,
	                 const IInterface* parent ); 

	virtual ~CMMMon();

	virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
	virtual StatusCode fillHistograms();
	virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

protected:

   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   /** location of data */
   std::string m_CMMJetHitsLocation;
   std::string m_CMMEtSumsLocation;
   std::string m_CMMRoILocation;

   std::string m_DataType;   
   std::string m_PathInRootFile;   
   std::string m_ErrorPathInRootFile;
   int m_NoEvents;
   int m_MaxEnergyRange;
   bool m_Offline;


  /** Histos */   
   // CMM Jet Hits
   TH1F* m_h_CMMJetHits_MainJets;
   TH1F* m_h_CMMJetHits_FwdJetsRight;
   TH1F* m_h_CMMJetHits_FwdJetsLeft;
   TH1F* m_h_CMMJetHits_EtMap;
   // JEM Hits
   TH1F* m_h_CMMJetHits_JEM_MainHits;
   TH1F* m_h_CMMJetHits_JEM_FwdHitsRight;
   TH1F* m_h_CMMJetHits_JEM_FwdHitsLeft;
   TH1F* m_h_CMMJetHits_JEM_Crate0ParityError;
   TH1F* m_h_CMMJetHits_JEM_Crate1ParityError;

   // CMM Et Sums
   TH1F* m_h_CMMEtSums_Ex;
   TH1F* m_h_CMMEtSums_Ey;
   TH1F* m_h_CMMEtSums_Et;
   TH1F* m_h_CMMEtSums_MissingEtMap;
   TH1F* m_h_CMMEtSums_SumEtMap;
   // JEM Et Sums
   TH1F*  m_h_CMMEtSums_JEM_Ex;
   TH1F*  m_h_CMMEtSums_JEM_Ey;
   TH1F*  m_h_CMMEtSums_JEM_Et; 

   // CMM RoI
   TH1F* m_h_CMMRoI_JetEtHits;
   TH1F* m_h_CMMRoI_SumEtHits;
   TH1F* m_h_CMMRoI_MissingEtHits;

   TH1F* m_h_CMMRoI_Ex;
   TH1F* m_h_CMMRoI_Ey;
   TH1F* m_h_CMMRoI_Et;

   //errors
   TH2F* m_h_CMMJet_error;
   TH2F* m_h_CMMEnergy_error;
   TH1F* m_h_CMMRoI_error;
   TH1F* m_h_CMM_ErrorSummary;
   TH1F* m_h_TriggeredSlice;	  
};


#endif
