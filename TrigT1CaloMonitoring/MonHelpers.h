// ********************************************************************
//
// NAME:        MonHelpers.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: collection of useful functions that are used by several
//              classes
//
// ********************************************************************

#ifndef MonHelpers_H
#define MonHelpers_H



std::string Binary(unsigned int Hits, int NumberBits);
// converts the integer value "Hits" into a string with lenth = NumberBits

int Multiplicity(std::string BinaryHitMap,int ThreshNo, int BitsPerThresh);
// gives back the multiplicity a certain threshold (ThresNo), where BitsPerThresh
// is the number of bits that are reserved for each Threshold


#endif
