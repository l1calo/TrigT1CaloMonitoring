// ********************************************************************
//
// NAME:        MonHelpers.cxx
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: collection of useful functions that are used by several
//              classes
//
// ********************************************************************

#include <sstream>


#include <algorithm>
#include <math.h>
#include <functional>
#include <iostream>


//_______________________________ Binary ___________________________________________
std::string Binary(unsigned int Hits, int NumberBits) 
{
  std::string temp = "";
  int temp2=0;

  for (int i=0; i<NumberBits; i++) //change hits into binary form
    {
      temp2= ((Hits >> i) &0x1 );
      temp=temp+static_cast<char>(temp2+48);
    }

  return temp;
}


//_______________________________ Multiplicity ___________________________________________
int Multiplicity(std::string BinaryHitMap, int ThresholdNumber, int BitsPerThresh) 
{
  std::string temp = "";
  int mult = 0;

  temp = BinaryHitMap;

  temp.assign(temp,(ThresholdNumber * BitsPerThresh), BitsPerThresh);
  // assigns to temp only the interesting bits

  for (int i=0; i<BitsPerThresh; i++) // converts the binary
    // display into an integer number
    {
      mult= mult + int(pow(2,i)) * (static_cast<int>(temp.substr()[i])-48);
    }

  return mult;
}
