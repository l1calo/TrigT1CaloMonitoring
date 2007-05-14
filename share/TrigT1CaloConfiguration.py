####################################################################################################################
##################################### L1CaloSimulation #############################################################
####################################################################################################################

 
# TriggerTowerMaker parameters
#TriggerTowerMaker = Algorithm( "TriggerTowerMaker" )

# JetElementMaker
JetElementMaker = Algorithm( "JetElementMaker" )

# JetTrigger
JetTrigger = Algorithm( "JetTrigger" )
JetTrigger.NumberOfThreshSets = 1
JetTrigger.JetElementThreshold = 0
#JetTrigger.OutputLevel = VERBOSE

# EnergyTrigger
EnergyTrigger = Algorithm( "EnergyTrigger" )
EnergyTrigger.EtSumJEThresh = 1
EnergyTrigger.EtMissJEThresh = 0
#EnergyTrigger.OutputLevel = VERBOSE

# JEP CMMs
JEPCMMMaker = Algorithm( "JEPCMMMaker" )

# ROD
ROD = Algorithm( "ROD" )
