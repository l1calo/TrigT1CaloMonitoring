#include "TrigT1CaloMonitoring/TrigT1CaloMonTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloBSMonTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloCpmMonTool.h"
#include "TrigT1CaloMonitoring/JetElementMon.h"
#include "TrigT1CaloMonitoring/TriggerTowerMon.h"
#include "TrigT1CaloMonitoring/JEMMon.h"
#include "TrigT1CaloMonitoring/CMMMon.h"
#include "TrigT1CaloMonitoring/SimBSMon.h"

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_TOOL_FACTORY(TrigT1CaloMonTool )
DECLARE_TOOL_FACTORY(TrigT1CaloBSMonTool )
DECLARE_TOOL_FACTORY(TrigT1CaloCpmMonTool )
DECLARE_TOOL_FACTORY(TriggerTowerMon);
DECLARE_TOOL_FACTORY(JetElementMon);
DECLARE_TOOL_FACTORY(JEMMon);
DECLARE_TOOL_FACTORY(CMMMon);
DECLARE_TOOL_FACTORY(SimBSMon);

DECLARE_FACTORY_ENTRIES(TrigT1CaloMonitoring) {
  DECLARE_ALGTOOL(TrigT1CaloMonTool )
  DECLARE_ALGTOOL(TrigT1CaloBSMonTool )
  DECLARE_ALGTOOL(TrigT1CaloCpmMonTool )
  DECLARE_ALGTOOL(TriggerTowerMon);
  DECLARE_ALGTOOL(JetElementMon);
  DECLARE_ALGTOOL(JEMMon);
  DECLARE_ALGTOOL(CMMMon);
  DECLARE_ALGTOOL(SimBSMon);

}

