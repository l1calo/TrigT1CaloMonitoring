// ********************************************************************
//
// NAME:        JEMMon.h
// PACKAGE:     TrigT1CaloMonitoring  
//
// AUTHOR:      Johanna Fleckner (Johanna.Fleckner@uni-mainz.de)
//           
// DESCRIPTION: Monitoring of the JEP on JEM level
//
// ********************************************************************

#ifndef JEMMon_H
#define JEMMon_H

#include <string>
#include <vector>

#include "GaudiKernel/ToolHandle.h"

#include "DataModel/DataVector.h"

#include "AthenaMonitoring/ManagedMonitorToolBase.h"

class TH1F_LW;
class TH2F_LW;
class TH2I_LW;

class StatusCode;

class TrigT1CaloMonErrorTool;
class TrigT1CaloLWHistogramTool;

namespace LVL1 {
  class JEMHits;
  class JEMEtSums;
  class JetElement;
  class JEMRoI;
}

/** Monitoring of the JEP on JEM level.
 *
 *  Produces histograms of JEM input data, output data and hardware errors.
 *
 *  <b>ROOT Histogram Directories:</b>
 *
 *  <table>
 *  <tr><th> Directory                       </th><th> Contents                             </th></tr>
 *  <tr><td> @c L1Calo/JEM/Input             </td><td> JetElement energy, eta and phi distributions <br>
 *                                                     Eta/phi hitmaps for each slice       <br>
 *                                                     Eta/phi hitmaps weighted by energy   <br>
 *                                                     Triggered slice number               </td></tr>
 *  <tr><td> @c L1Calo/JEM/Output/RoI        </td><td> Hit multiplicity per threshold       <br>
 *                                                     Eta/phi hits passing each threshold  <br>
 *                                                     Eta.phi of RoI saturation            </td></tr>
 *  <tr><td> @c L1Calo/JEM/Output/Thresholds </td><td> Hit multiplicity per threshold       <br>
 *                                                     Hits per JEM by threshold            </td></tr>
 *  <tr><td> @c L1Calo/JEM/Output/EnergySums </td><td> Ex, Ey, Et sum distributions         </td></tr>
 *  <tr><td> @c L1Calo/JEM/Errors/Hardware   </td><td> Error summary                        <br>
 *                                                     Substatus errors                     <br>
 *                                                     Parity and LinkDown errors           <br>
 *                                                     Error event numbers                  </td></tr>
 *  </table>
 *
 *  <!--
 *  <b>Notes on Particular Histograms:</b>
 *
 *  <table>
 *  <tr><th> Histogram                       </th><th> Comment                        </th></tr>
 *  <tr><td>                                 </td><td>                                </td></tr>
 *  </table>
 *  -->
 *
 *  <b>Custom Merges Used (Tier0):</b>
 *
 *  <table>
 *  <tr><th> Merge                                    </th><th> Used For                 </th></tr>
 *  <tr><td> @ref MergesUsedsection "@c eventSample " </td><td> Error event number plots </td></tr>
 *  </table>
 *
 *  <b>StoreGate Containers Used:</b>
 *
 *  <table>
 *  <tr><th> Container                 </th><th> Comment                                  </th></tr>
 *  <tr><td> @c DataVector
 *           @c <LVL1::JetElement>     </td><td> Jet element data                         </td></tr>
 *  <tr><td> @c DataVector
 *           @c <LVL1::JEMRoI>         </td><td> Jet RoI data                             </td></tr>
 *  <tr><td> @c DataVector
 *           @c <LVL1::JEMHits>        </td><td> Jet hits data                            </td></tr>
 *  <tr><td> @c DataVector
 *           @c <LVL1::JEMEtSums>      </td><td> Energy sums data                         </td></tr>
 *  <tr><td> @c std::vector<int>       <br>
 *           @c "L1CaloJEMErrorVector" </td><td> Output.
 *                                               Error summary bits for global histograms </td></tr>
 *  </table>
 *
 *  <b>Tools Used:</b>
 *
 *  <table>
 *  <tr><th> Tool                         </th><th> Description          </th></tr>
 *  <tr><td> @c TrigT1CaloMonErrorTool    </td><td> @copydoc m_errorTool </td></tr>
 *  <tr><td> @c TrigT1CaloLWHistogramTool </td><td> @copydoc m_histTool  </td></tr>
 *  </table>
 *
 *  <b>JobOption Properties:</b>
 *
 *  <table>
 *  <tr><th> Property               </th><th> Description                    </th></tr>
 *  <tr><td> @c JetElementLocation  </td><td> @copydoc m_JetElementLocation  </td></tr>
 *  <tr><td> @c JEMHitsLocation     </td><td> @copydoc m_JEMHitsLocation     </td></tr>
 *  <tr><td> @c JEMEtSumsLocation   </td><td> @copydoc m_JEMEtSumsLocation   </td></tr>
 *  <tr><td> @c JEMRoILocation      </td><td> @copydoc m_JEMRoILocation      </td></tr>
 *  <tr><td> @c NumberOfSlices      </td><td> @copydoc m_SliceNo             </td></tr>
 *  <tr><td> @c MaxEnergyRange      </td><td> @copydoc m_MaxEnergyRange      </td></tr>
 *  <tr><td> @c PathInRootFile      </td><td> @copydoc m_PathInRootFile      </td></tr>
 *  <tr><td> @c ErrorPathInRootFile </td><td> @copydoc m_ErrorPathInRootFile </td></tr>
 *  </table>
 *
 *  <b>Related Documentation:</b>
 *
 *  <a href="http://hepwww.rl.ac.uk/Atlas-L1/Modules/ROD/ROD-spec-version1_2_2.pdf">
 *  ATLAS Level-1 Calorimeter Trigger - Read-out Driver</a>
 *
 *  @authors Johanna Fleckner, Andrea Neusiedl, Peter Faulkner
 *
 */

class JEMMon : public ManagedMonitorToolBase
{
public:

   JEMMon( const std::string & type, const std::string & name,
   	                             const IInterface* parent ); 

   virtual ~JEMMon();

   virtual StatusCode initialize();
   virtual StatusCode bookHistogramsRecurrent();
   virtual StatusCode fillHistograms();
   virtual StatusCode procHistograms();

private:

   /// Error summary bins
   enum SummaryErrors { EMParity, HadParity, EMLink, HadLink, JEMStatus,
                        RoIParity, NumberOfSummaryBins };

   typedef DataVector<LVL1::JetElement> JECollection;
   typedef DataVector<LVL1::JEMHits> JEMHitsCollection;
   typedef DataVector<LVL1::JEMEtSums> JEMEtSumsCollection;
   typedef DataVector<LVL1::JEMRoI> JemRoiCollection;

   /// Tool to retrieve bytestream errors
   ToolHandle<TrigT1CaloMonErrorTool>    m_errorTool;
   /// Histogram helper tool
   ToolHandle<TrigT1CaloLWHistogramTool> m_histTool;

   /** location of data */
   /// JetElement collection StoreGate key
   std::string m_JetElementLocation;
   /// JEMHits collection StoreGate key
   std::string m_JEMHitsLocation;
   /// JEMEtSums collection StoreGate key
   std::string m_JEMEtSumsLocation;   
   /// JEMRoI collection StoreGate key
   std::string m_JEMRoILocation;

   /// Maximum number of slices
   int m_SliceNo;
   /// Maximum jet element energy
   int m_MaxEnergyRange;
   /// Histograms booked flag
   bool m_histBooked;

   /// Directory in ROOT
   std::string m_PathInRootFile;   
   /// Error directory in ROOT
   std::string m_ErrorPathInRootFile;

   /** Histos */
   TH1F_LW* m_h_jem_em_1d_jetEl_Eta;              ///< EM TowerSum distribution per eta  --  JEM input
   TH1F_LW* m_h_jem_had_1d_jetEl_Eta;             ///< HAD TowerSum distribution per eta  --  JEM input
   TH1F_LW* m_h_jem_em_1d_jetEl_Phi;              ///< EM TowerSum distribution per phi  --  JEM input
   TH1F_LW* m_h_jem_had_1d_jetEl_Phi;             ///< HAD TowerSum distribution per phi  --  JEM input
   TH1F_LW* m_h_jem_em_1d_jetEl_Energy;           ///< TowerSum EM energy distribution  --  JEM input
   TH1F_LW* m_h_jem_had_1d_jetEl_Energy;          ///< TowerSum HAD energy distribution  --  JEM input

   // HitMaps
   TH2F_LW* m_h_jem_em_2d_etaPhi_jetEl_HitMapWeighted;            ///< eta-phi map of EM TowerSum weighted with energy  --  JEM input
   TH2F_LW* m_h_jem_had_2d_etaPhi_jetEl_HitMapWeighted;           ///< eta-phi map of HAD TowerSum weighted with energy  --  JEM input
   std::vector<TH2F_LW*> m_v_jem_em_2d_etaPhi_jetEl_HitMapSlice;  ///< eta-phi map of EM TowerSum for each Timeslice  --  JEM input
   std::vector<TH2F_LW*> m_v_jem_had_2d_etaPhi_jetEl_HitMapSlice; ///< eta-phi map of HAD TowerSum for each Timeslice  --  JEM input

   // error maps
   TH2F_LW* m_h_jem_2d_Status;                    ///< Error reports from JEM SubStatus Word
   TH2F_LW* m_h_jem_em_2d_etaPhi_jetEl_Parity;    ///< Jet Element EM Parity Errors
   TH2F_LW* m_h_jem_had_2d_etaPhi_jetEl_Parity;   ///< Jet Element Had Parity Errors
   TH2F_LW* m_h_jem_em_2d_etaPhi_jetEl_LinkDown;  ///< Jet Element EM Link Down Errors
   TH2F_LW* m_h_jem_had_2d_etaPhi_jetEl_LinkDown; ///< Jet Element Had Link Down Errors
   
   // number of triggered slice
   TH1F_LW* m_h_jem_1d_TriggeredSlice;            ///< Number of the Triggered Slice for JE

   // JEM Hits
   TH1F_LW* m_h_jem_1d_thresh_MainHits;           ///< Main Jet Hit Multiplicity per Threshold  --  JEM DAQ
   TH1F_LW* m_h_jem_1d_thresh_FwdHitsRight;       ///< Fwd Right Jet Hit Multiplicity per Threshold  --  JEM DAQ
   TH1F_LW* m_h_jem_1d_thresh_FwdHitsLeft;        ///< Fwd Left Jet Hit Multiplicity per Threshold  --  JEM DAQ
   TH2F_LW* m_h_jem_2d_thresh_HitsPerJem;         ///< HitMap of Hits per JEM
   // JEM Et Sums
   TH1F_LW*  m_h_jem_1d_energy_SubSumsEx;         ///< JEM Ex  --  JEM DAQ
   TH1F_LW*  m_h_jem_1d_energy_SubSumsEy;         ///< JEM Ey  --  JEM DAQ
   TH1F_LW*  m_h_jem_1d_energy_SubSumsEt;         ///< JEM Et  --  JEM DAQ

   // JEM RoI
   TH1F_LW* m_h_jem_1d_roi_MainHits;              ///< Main Jet Hit Multiplicity per Threshold  --  JEM RoI
   TH1F_LW* m_h_jem_1d_roi_FwdHitsRight;          ///< Forward Right Jet Hit Multiplicity per Threshold  --  JEM RoI
   TH1F_LW* m_h_jem_1d_roi_FwdHitsLeft;           ///< Forward Left Jet Hit Multiplicity per Threshold  --  JEM RoI
  
   std::vector<TH2F_LW*> m_v_jem_2d_etaPhi_roi_MainThresh; ///< eta-phi Map of Main Hits passing each Threshold   --  JEM RoI
   std::vector<TH2F_LW*> m_v_jem_2d_etaPhi_roi_FwdThresh;  ///< eta-phi Map of Fwd Hits passing each Threshold   --  JEM RoI
 
   // errors and saturation
   TH2F_LW* m_h_jem_2d_etaPhi_roi_Parity;         ///< JEM RoI Parity Errors
   TH2F_LW* m_h_jem_2d_etaPhi_roi_Saturation;     ///< JEM RoI Saturation
   
   // Error Summary
   TH1F_LW* m_h_jem_1d_ErrorSummary;              ///< Summary of JEM Data Errors
   TH2I_LW* m_h_jem_2d_ErrorEventNumbers;         ///< JEM Error Event Numbers
};


#endif
