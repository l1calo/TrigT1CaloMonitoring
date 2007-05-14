####################################################################################################################
##################################### Bytestream Optionen ##########################################################
####################################################################################################################

#StreamBS = Algorithm( "StreamBS" )
#StreamBS.ItemList    = [ ]

include( "TriggerRelease/TrigWriteBS_topOptions.py" )
 
include ( "TrigT1CaloByteStream/WritePpmBS_jobOptions.py" )
include ( "TrigT1CaloByteStream/WriteJepBS_jobOptions.py" )
include ( "TrigT1CaloByteStream/WriteCpBS_jobOptions.py" )

# Algs
theApp.Dlls += [ "TrigT1CaloByteStream" ]

#theApp.TopAlg += [ "PpmTester" ]
#Algorithm( "PpmTester" ).TriggerTowerLocation ="LVL1TriggerTowers"

#theApp.TopAlg += [ "JemTester" ]
#Algorithm( "JemTester" ).JetElementLocation ="LVL1JetElements"

#theApp.TopAlg += [ "JepContainerMaker" ]
#Algorithm( "JepContainerMaker" ).JetElementLocation ="LVL1JetElements"

#ByteStreamCnvSvc = Service( "ByteStreamCnvSvc" )
#ByteStreamCnvSvc.JepByteStreamTool.DataFormat = 1
#ByteStreamCnvSvc.PpmByteStreamTool.DataFormat = 1

#StreamBS = Algorithm( "StreamBS" )
#StreamBS.ItemList   += [ "1321678566#*" ]
#StreamBS.ItemList   += [ "6207#*" ]

#ByteStreamCnvSvc.PpmByteStreamTool.OutputLevel = VERBOSE
