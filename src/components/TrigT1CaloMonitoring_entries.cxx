#include "TrigT1CaloMonitoring/TrigT1CaloMonTool.h"
#include "TrigT1CaloMonitoring/TrigT1CaloBSMonTool.h"
#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_TOOL_FACTORY(TrigT1CaloMonTool )
DECLARE_TOOL_FACTORY(TrigT1CaloBSMonTool )

DECLARE_FACTORY_ENTRIES(TrigT1CaloMonitoring) {
  DECLARE_ALGTOOL(TrigT1CaloMonTool )
  DECLARE_ALGTOOL(TrigT1CaloBSMonTool )

}

