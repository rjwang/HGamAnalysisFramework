#ifndef HGAMMAINCLUDES_H
#define HGAMMAINCLUDES_H

// STL includes
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include <time.h>


// Root includes
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2F.h>

#include <EventLoop/StatusCode.h>
#include "EventLoop/OutputStream.h"
#include "CPAnalysisExamples/errorcheck.h"
#include "xAODRootAccess/tools/Message.h"

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#include "PATInterfaces/SystematicRegistry.h"

#ifndef __CINT__
// xAOD includes
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/AuxInfoBase.h"

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODEgamma/PhotonAuxContainer.h"
#include "xAODEgamma/PhotonxAODHelpers.h"
#include "xAODEgamma/ElectronxAODHelpers.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/ElectronAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/MuonAuxContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODBTagging/BTagging.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODMissingET/MissingETContainer.h"

// Tools
#include "AsgTools/AnaToolHandle.h"

#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h"
#include "ElectronPhotonSelectorTools/AsgPhotonIsEMSelector.h"
#include "ElectronEfficiencyCorrection/AsgElectronEfficiencyCorrectionTool.h"
#include "PhotonEfficiencyCorrection/AsgPhotonEfficiencyCorrectionTool.h"

#include "MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h"
#include "MuonEfficiencyCorrections/MuonTriggerScaleFactors.h"
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"

#include "JetCalibTools/JetCalibrationTool.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "xAODBTaggingEfficiency/BTaggingEfficiencyTool.h"

#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "PileupReweighting/IPileupReweightingTool.h"

#include "TrigDecisionTool/TrigDecisionTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "xAODTrigEgamma/TrigPhotonContainer.h"
#include "xAODTrigEgamma/TrigElectronContainer.h"
#include "TrigMuonMatching/TrigMuonMatching.h"

#include "PathResolver/PathResolver.h"

#include "PhotonVertexSelection/PhotonVertexHelpers.h"

#include "FourMomUtils/xAODP4Helpers.h"

#endif

// constants and typedefs

//! \brief Hgamma namespace
namespace HG {
  
  /// Helper macro to check xAOD::TReturnCode return values
  /// See Atilla's slide 8-9 at: https://indico.cern.ch/event/362819/
#define EL_CHECK( COMMENT, EXP )                        \
  do {                                                  \
    if ( ! EXP.isSuccess() ) {                          \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),  \
             #EXP );                                    \
      return EL::StatusCode::FAILURE;                   \
    }                                                   \
  } while( false );

#define CP_CHECK( COMMENT, EXP )                        \
  do {                                                  \
    if ( EXP != CP::SystematicCode::Ok ) {              \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),  \
             #EXP );                                    \
      return CP::SystematicCode::Unsupported;           \
    }                                                   \
  } while( false );

#define CC_CHECK( COMMENT, EXP )                   \
  do {                                                  \
    if ( EXP != CP::CorrectionCode::Ok ) {              \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),  \
             #EXP );                                    \
    }                                                   \
  } while( false );

#define HG_CHECK( COMMENT, EXP )                        \
  do {                                                  \
    if ( ! EXP.isSuccess() ) {                          \
      Error( COMMENT,                                   \
             XAOD_MESSAGE("\n  Failed to execute %s"),  \
             #EXP );                                    \
    }                                                   \
  } while( false );

  template <typename T>
  struct Identity {
    typedef T type;
  };

  // 1*GeV = 1000*MeV
  static const double GeV(1000), invGeV(1.0/GeV);
 
#ifndef __MAKECINT__
  typedef std::pair< xAOD::PhotonContainer*, xAOD::ShallowAuxContainer*> PhotonShallowCopies;
#endif

  typedef std::vector<CP::SystematicSet> SystematicList;

  namespace Iso {
    enum IsolationType {
      LooseTrackOnly,
      Loose,
      Gradient,
      GradientLoose,
      FixedCutTight,
      FixedCutTightTrackOnly,
      FixedCutLoose,
      FixedCutTightCaloOnly,
      FixedCutLooseCaloOnly,
      UserDefined,
      Undefined
    }; // enum IsolationType
  } // namespace Iso

  const unsigned int PhotonLoosePrime5 = egammaPID::PhotonLoose |
    0x1 << egammaPID::ClusterStripsDeltaEmax2_Photon |
    0x1 << egammaPID::ClusterStripsEratio_Photon;

  const unsigned int PhotonLoosePrime4 = PhotonLoosePrime5 |
    0x1 << egammaPID::ClusterStripsWtot_Photon;

  const unsigned int PhotonLoosePrime3 = PhotonLoosePrime4 |
    0x1 << egammaPID::ClusterStripsDEmaxs1_Photon;

  const unsigned int PhotonLoosePrime2 = PhotonLoosePrime3 |
    0x1 << egammaPID::ClusterStripsDeltaE_Photon;

}

#ifndef __CINT__
// Hgamma includes
#include "HGamAnalysisFramework/HgammaUtils.h"
#include "HGamAnalysisFramework/Config.h"
#include "HGamAnalysisFramework/TruthUtils.h"
#endif


#endif
