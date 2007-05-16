####################################################################################################################
##################################### L1CaloSimulation #############################################################
####################################################################################################################

 
# TriggerTowerMaker parameters
TriggerTowerMaker = Algorithm( "TriggerTowerMaker" )
#TriggerTowerMaker.OutputLevel = DEBUG

# JetElementMaker
JetElementMaker = Algorithm( "JetElementMaker" )
#JetElementMaker.OutputLevel = DEBUG

# JetTrigger
JetTrigger = Algorithm( "JetTrigger" )
JetTrigger.NumberOfThreshSets = 1
JetTrigger.JetElementThreshold = 0
#JetTrigger.OutputLevel = DEBUG

# EnergyTrigger
EnergyTrigger = Algorithm( "EnergyTrigger" )
EnergyTrigger.EtSumJEThresh = 1
EnergyTrigger.EtMissJEThresh = 0
#EnergyTrigger.OutputLevel = DEBUG

# JEP CMMs
JEPCMMMaker = Algorithm( "JEPCMMMaker" )
#JEPCMMMaker.OutputLevel = DEBUG

# ROD
ROD = Algorithm( "ROD" )
#ROD.OutputLevel = DEBUG
