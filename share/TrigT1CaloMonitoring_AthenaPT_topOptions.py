######################################################################
# TopOptions to execute RecExCommissioning in PT
#               (and provide data for online event display)
#
# 21/08/2007  Jiri.Masik@cern.ch
######################################################################

#default flags
include("RecExCommission/RecExCommissionFlags_jobOptions.py")

ATLASCosmicFlags.EvtMax = 10

#crucial
ATLASCosmicFlags.doSim=False
ATLASCosmicFlags.JiveXML=False
ATLASCosmicFlags.doOnline=False

#optional switches
ATLASCosmicFlags.doMonitoring=False
ATLASCosmicFlags.doLVL1Calo=True

ATLASCosmicFlags.doInDet=False
ATLASCosmicFlags.doMuons=False
ATLASCosmicFlags.doLAr=False
ATLASCosmicFlags.doTile=False
ATLASCosmicFlags.doCaloTrkMuId=False
ATLASCosmicFlags.doLVL1CaloCPJEPMonitoring=False
ATLASCosmicFlags.doHLT=False
ATLASCosmicFlags.doCTPMon=False

include( "RecExCommission/RecExCommission_topOptions.py" )

#Settuing up the monitoring infrastructure
from AthenaMonitoring.AthenaMonitoringConf import AthenaMonManager
topSequence += AthenaMonManager( "ManagedAthenaMonPilot" )
ManagedAthenaMonPilot = topSequence.ManagedAthenaMonPilot
ManagedAthenaMonPilot.FileKey             = "GLOBAL"
#FIXME: The FileKey is not propagated to the following MonManagers.
#FIXME: We have to make sure that the subdetector fragements use the same FileKey
#FIXME: alternavely, we could add all the FileKeys with more  THistSvc.Output +=[..] statements.
ManagedAthenaMonPilot.ManualDataTypeSetup = True
ManagedAthenaMonPilot.DataType            = "cosmics"
##  if doOnline:
##      ManagedAthenaMonPilot.Environment     = "online"
##  else:
ManagedAthenaMonPilot.Environment         = "tier0"
ManagedAthenaMonPilot.ManualRunLBSetup    = False
#ManagedAthenaMonPilot.Run                = 1
#ManagedAthenaMonPilot.LumiBlock          = 1

include ( "TrigT1CaloMonitoring/TrigT1CaloMonitoring_forAthenaPT.py")
toolSvc.L1TriggerTowerMonTool.DistPerChannel = False


#-----------------------------------------------------------------------------
# Monitoring output
#-----------------------------------------------------------------------------
theApp.Dlls += [ "GaudiSvc" ]
THistSvc = Algorithm( "THistSvc" )
#    if os.path.exists(RootHistOutputFilename):
#        os.remove(RootHistOutputFilename)
#THistSvc.Output += [MonitorOutput+" DATAFILE='"+ATLASCosmicFlags.MonOutputFile+"'OPT='NEW'"]
#THistSvc.Output += ["SHIFT DATAFILE='"+ATLASCosmicFlags.MonOutputFile+"'OPT='NEW'"]
THistSvc.Output += ["GLOBAL DATAFILE='"+ATLASCosmicFlags.MonOutputFile+"'OPT='RECREATE'"]
