####################################################################################################################
########################################## L1Calo Simulation #######################################################
####################################################################################################################
# Import the configurable algorithms for TrigT1Calo
from AthenaCommon.GlobalFlags  import globalflags

from TrigT1Calo.TrigT1CaloConf import LVL1__TriggerTowerMaker
from TrigT1Calo.TrigT1CaloConf import LVL1__CPMTowerMaker
from TrigT1Calo.TrigT1CaloConf import LVL1__JetElementMaker
from TrigT1Calo.TrigT1CaloConf import LVL1__EmTauTrigger
from TrigT1Calo.TrigT1CaloConf import LVL1__JetTrigger
from TrigT1Calo.TrigT1CaloConf import LVL1__EnergyTrigger
from TrigT1Calo.TrigT1CaloConf import LVL1__CPCMMMaker
from TrigT1Calo.TrigT1CaloConf import LVL1__JEPCMMMaker
from TrigT1Calo.TrigT1CaloConf import LVL1__ROD
from TrigT1Calo.TrigT1CaloConf import LVL1__Tester
from TrigT1Calo.TrigT1CaloConf import LVL1__DumpTriggerObjects

# Get the algorithm sequence
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

job += LVL1__CPMTowerMaker( 'CPMTowerMaker' )
job += LVL1__JetElementMaker( 'JetElementMaker' )
job += LVL1__EmTauTrigger( 'EmTauTrigger' )
job += LVL1__JetTrigger( 'JetTrigger' )
job += LVL1__EnergyTrigger( 'EnergyTrigger' )
job += LVL1__ROD( 'ROD' )
job += LVL1__CPCMMMaker( 'CPCMMMaker' )
job += LVL1__JEPCMMMaker( 'JEPCMMMaker' )


#*************************************************************
# Set input and output locations  
#*************************************************************
# CPMTowerMaker
job.CPMTowerMaker.TriggerTowerLocation = "TriggerTowers"
job.CPMTowerMaker.CPMTowerLocation = "Sim_CPMTowers"

# JetElementMaker
job.JetElementMaker.TriggerTowerLocation ="TriggerTowers"  
job.JetElementMaker.JetElementLocation ="Sim_JetElements"
#job.JetElementMaker.OutputLevel = VERBOSE

# EmTauTrigger
job.EmTauTrigger.TriggerTowerLocation = "TriggerTowers"
job.EmTauTrigger.CPMHitsLocation = "Sim_CPMHits"

# JetTrigger
job.JetTrigger.JetElementLocation="JetElements"
job.JetTrigger.JEMHitsLocation = "Sim_JEMHits"

# EnergyTrigger
job.EnergyTrigger.JetElementLocation="JetElements"
job.EnergyTrigger.JEMEtSumsLocation="Sim_JEMEtSums"

# JEP CMMs
job.JEPCMMMaker.JetElementLocation="JetElements"
job.JEPCMMMaker.JEMHitsLocation = "JEMHits"
job.JEPCMMMaker.JEMEtSumsLocation = "JEMEtSums"
job.JEPCMMMaker.CMMJetHitsLocation = "Sim_CMMJetHits"
job.JEPCMMMaker.CMMEtSumsLocation = "Sim_CMMEtSums"
job.JEPCMMMaker.JEMRoILocation = "JEMRoIs"
job.JEPCMMMaker.CMMRoILocation = "Sim_CMMRoIs"
#job.JEPCMMMaker.OutputLevel = VERBOSE

# CP CMMs
job.CPCMMMaker.CPMTowerLocation = "PMTowers"
job.CPCMMMaker.CPMHitsLocation = "CPMHits"
job.CPCMMMaker.CMMCPHitsLocation = "Sim_CMMCPHits"
job.CPCMMMaker.CPMRoILocation = "Sim_CPMRoIs"
