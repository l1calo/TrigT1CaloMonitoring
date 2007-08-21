######################################################################
# TopOptions to execute RecExCommissioning in PT
#               (and provide data for online event display)
#
# 21/08/2007  Jiri.Masik@cern.ch
######################################################################

#default flags
include("RecExCommission/RecExCommissionFlags_jobOptions.py")

#crucial
ATLASCosmicFlags.doSim=False
ATLASCosmicFlags.JiveXML=False
ATLASCosmicFlags.doOnline=False

#optional switches
ATLASCosmicFlags.doMonitoring=True
ATLASCosmicFlags.doLVL1Calo=True
ATLASCosmicFlags.doInDet=False
ATLASCosmicFlags.doMuons=False
ATLASCosmicFlags.doLAr=False
ATLASCosmicFlags.doTile=False
ATLASCosmicFlags.doCaloTrkMuId=False

include( "RecExCommission/RecExCommission_topOptions.py" )

include ( "TrigT1CaloMonitoring/TrigT1CaloMonitoring_forAthenaPT.py")
