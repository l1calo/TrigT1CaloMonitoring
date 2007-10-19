#include "TrigT1CaloMonitoring/TrigT1CaloCpmMonTool.h"
#include "TrigT1CaloMonitoring/PPrMon.h"
#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/CMMMon.h"
#include "TrigT1CaloMonitoring/JEPTransPerfMon.h"
#include "TrigT1CaloMonitoring/PPrPerfMon.h"


#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_TOOL_FACTORY(TrigT1CaloCpmMonTool )
DECLARE_TOOL_FACTORY(PPrMon)
DECLARE_TOOL_FACTORY(JEMMon)
DECLARE_TOOL_FACTORY(CMMMon)
DECLARE_TOOL_FACTORY(JEPTransPerfMon)
DECLARE_TOOL_FACTORY(PPrPerfMon)

DECLARE_FACTORY_ENTRIES(TrigT1CaloMonitoring) {
  DECLARE_ALGTOOL(TrigT1CaloCpmMonTool );
  DECLARE_ALGTOOL(PPrMon);
  DECLARE_ALGTOOL(JEMMon);
  DECLARE_ALGTOOL(CMMMon);
  DECLARE_ALGTOOL(JEPTransPerfMon);
  DECLARE_ALGTOOL(PPrPerfMon);

}

