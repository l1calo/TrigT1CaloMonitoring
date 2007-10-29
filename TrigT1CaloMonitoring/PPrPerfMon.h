// ********************************************************************
//
// NAME:     PPrPerfMon.h
// PACKAGE:  TrigT1CaloMonitoring
//
// AUTHOR:   Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//	     
//
// ********************************************************************
#ifndef PPRPERFMON_H
#define PPRPERFMON_H

#include <map>
#include "AthenaMonitoring/AthenaMonManager.h"
#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "GaudiKernel/StatusCode.h"
#include "CLHEP/Units/SystemOfUnits.h"


#include "TH1.h"
#include "TH2.h"

class PPrPerfMon: public ManagedMonitorToolBase
{

 public:
  
  PPrPerfMon(const std::string & type, const std::string & name,
		  const IInterface* parent);
    

  virtual ~PPrPerfMon();

  virtual StatusCode bookHistograms( bool isNewEventsBlock, bool isNewLumiBlock, bool isNewRun );
  virtual StatusCode fillHistograms();
  virtual StatusCode procHistograms( bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun );

private:
  std::string m_DataType;
  std::string m_PathInRootFile;
      
protected:
   /** a handle on Store Gate for access to the Event Store */
   StoreGateSvc* m_storeGate;

   // StoreGate service
   StoreGateSvc* m_detStore;
   
};

#endif
