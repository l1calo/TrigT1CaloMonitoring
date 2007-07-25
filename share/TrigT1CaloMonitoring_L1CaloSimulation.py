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
#from TrigT1Calo.TrigT1CaloConf import LVL1__Tester
#from TrigT1Calo.TrigT1CaloConf import LVL1__DumpTriggerObjects

# Get the algorithm sequence
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

job += LVL1__CPMTowerMaker( 'CPMTowerMaker' )
job += LVL1__JetElementMaker( 'JetElementMaker' )
job += LVL1__EmTauTrigger( 'EmTauTrigger' )
job += LVL1__EnergyTrigger( 'EnergyTrigger' )
job += LVL1__JetTrigger( 'JetTrigger' )
job += LVL1__JEPCMMMaker( 'JEPCMMMaker' )
job += LVL1__CPCMMMaker( 'CPCMMMaker' )
job += LVL1__ROD( 'ROD' )


#*************************************************************
# Set input and output locations  
#*************************************************************
# CPMTowerMaker
job.CPMTowerMaker.TriggerTowerLocation = "TriggerTowers"
job.CPMTowerMaker.CPMTowerLocation = "Sim_CPMTowers"
#job.CPMTowerMaker.OutputLevel = DEBUG

# JetElementMaker
job.JetElementMaker.TriggerTowerLocation ="TriggerTowers"  
job.JetElementMaker.JetElementLocation ="Sim_JetElements"
#job.JetElementMaker.OutputLevel = DEBUG

# EmTauTrigger
job.EmTauTrigger.TriggerTowerLocation = "TriggerTowers"
job.EmTauTrigger.CPMHitsLocation = "Sim_CPMHits"
job.EmTauTrigger.EmTauROILocation = "Sim_EmTauRoIs"
#job.EmTauTrigger.LVL1ConfigSvc = ""
#job.EmTauTrigger.OutputLevel = DEBUG

# JetTrigger
job.JetTrigger.JetElementLocation="Sim_JetElements"
job.JetTrigger.JEMHitsLocation = "Sim_JEMHits"
job.JetTrigger.JetROIOutputLocation= "Sim_JetRoIs"
#job.JetTrigger.LVL1ConfigSvc = ""
#job.JetTrigger.NumberOfThreshSets = ""
#job.JetTrigger.DefaultClusterThresholds = ""
#job.JetTrigger.DefaultMultiplicities = ""
#job.JetTrigger.DefaultAlgorithms = ""
#job.JetTrigger.JetElementThreshold = ""
#job.JetTrigger.OutputLevel=DEBUG

# EnergyTrigger
job.EnergyTrigger.JetElementLocation="Sim_JetElements"
job.EnergyTrigger.JEMEtSumsLocation="Sim_JEMEtSums"
job.EnergyTrigger.EnergyRoILocation = "Sim_EnergyRoIs"
#job.EnergyTrigger.EnergyCTPLocation = ""
#job.EnergyTrigger.LVL1ConfigSvc = ""
#job.EnergyTrigger.EtSumJEThresh = ""
#job.EnergyTrigger.EtMissJEThresh = ""
#job.EnergyTrigger.ExEyRanges = ""
#job.EnergyTrigger.OutputLevel=DEBUG

# JEP CMMs
job.JEPCMMMaker.JetElementLocation="Sim_JetElements"
job.JEPCMMMaker.JEMHitsLocation = "Sim_JEMHits"
job.JEPCMMMaker.JEMEtSumsLocation = "Sim_JEMEtSums"
job.JEPCMMMaker.CMMJetHitsLocation = "Sim_CMMJetHits"
job.JEPCMMMaker.CMMEtSumsLocation = "Sim_CMMEtSums"
job.JEPCMMMaker.JEMRoILocation = "Sim_JEMRoIs"
job.JEPCMMMaker.CMMRoILocation = "Sim_CMMRoIs"
job.JEPCMMMaker.JetEtRoILocation = "Sim_JetEtRoIs"
job.JEPCMMMaker.EnergyRoILocation = "Sim_EnergyRoIs"
job.JEPCMMMaker.EtMapsLocation = "Sim_CMMEtSumsMAPS"
job.JEPCMMMaker.JetRoILocation = "Sim_JetRoIs"
job.JEPCMMMaker.JEPBSCollectionLocation = "Sim_JEPBS"
job.JEPCMMMaker.JEPRoIBSCollectionLocation = "Sim_JEPRoIBS"
job.JEPCMMMaker.OutputLevel = DEBUG

# CP CMMs
job.CPCMMMaker.CPMTowerLocation = "Sim_CPMTowers"
job.CPCMMMaker.CPMHitsLocation = "Sim_CPMHits"
job.CPCMMMaker.CMMCPHitsLocation = "Sim_CMMCPHits"
job.CPCMMMaker.CPMRoILocation = "Sim_CPMRoIs"
job.CPCMMMaker.EmTauROILocation = "Sim_EmTauRoIs"
job.CPCMMMaker.CPBSCollectionLocation = "Sim_CPBS"
#job.CPCMMMaker.OutputLevel = DEBUG

#job.ROD.LVL1ConfigSvc = ""
job.ROD.EmTauRoILocation = "Sim_EmTauRoIs"
#job.ROD.EmTauCTPLocation = "Sim_EmTauCTP"
#job.ROD.EmTauSlinkLocation = "Sim_EmTauSlink"
#job.ROD.JEPSlinkLocation = "Sim_JEPSlink"
#job.ROD.EnergySlinkLocation = "Sim_EnergySlink"
job.ROD.JetRoILocation = "Sim_JetRoIs"
#job.ROD.JetCTPLocation = "Sim_JetCTP"
job.ROD.JetEtRoILocation = "Sim_JetEtRoIs"
job.ROD.EnergyRoILocation = "Sim_EnergyRoIs"
#job.ROD.OutputLevel = DEBUG
