####################################################################################################################
##################################### Calo Monitoring ##############################################################
####################################################################################################################

# *************************
# Application configuration
# *************************

#theApp.setup( MONTECARLO )
theApp.Dlls   += [ "AthenaMonitoring" ]
theApp.Dlls   += [ "TrigT1CaloMonitoring" ]
theApp.TopAlg += [ "AthenaMonManager/PrimaryManager" ]

# ************************
# Monitoring configuration
# ************************

## Setup the output file(s):
THistSvc = Algorithm( "THistSvc" )

## The string "TestMon" in the argument below is the 'FileKey'
## used by Athena to access the output file internally
THistSvc.Output = ["CaloMon DATAFILE='CaloMon.root' OPT='RECREATE'"]

## AthenaMonManager is the Algorithm that manages many classes inheriting
## from ManagedMonitorToolBase
monMan = Algorithm( "PrimaryManager" )

## Add all the ManagedMonitorToolBase objects
monMan.AthenaMonTools += [ "TrigT1CaloMonTool/L1CaloTool" ]

## get a handle on the ToolSvc
#from AthenaCommon.AppMgr import ToolSvc as toolSvc
ToolSvc = Algorithm( "ToolSvc" )

ToolSvc.L1CaloTool.TriggerTowerLocation = "LVL1TriggerTowers"
ToolSvc.L1CaloTool.JetElementLocation = "LVL1JetElements"
ToolSvc.L1CaloTool.OutputLevel = DEBUG


## FileKey must match that given to THistSvc
monMan.FileKey = "CaloMon"

## Set global monitoring parameters: see the AthenaMonManager class
## in the Control/AthenaMonitoring package
monMan.ManualDataTypeSetup = True
monMan.DataType            = "collisions"
monMan.Environment         = "tier0"
monMan.ManualRunLBSetup    = True
monMan.Run                 = 1
monMan.LumiBlock           = 1

####################################################################################################################
###################################### Ende ########################################################################
####################################################################################################################
