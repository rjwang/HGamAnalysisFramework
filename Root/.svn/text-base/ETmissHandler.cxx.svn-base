#include "HGamAnalysisFramework/ETmissHandler.h"

#include "HGamAnalysisFramework/HgammaUtils.h"

// MET EDM
#include "xAODMissingET/MissingET.h"
#include "xAODMissingET/MissingETContainer.h"

// MET Aux and association map
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETAssociationMap.h"

#include "xAODBase/IParticleHelpers.h"

namespace HG {

  /*! \brief Class that
   *  \author Luis March
   */

  //______________________________________________________________________________
  ETmissHandler::ETmissHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : HgammaHandler(name, event, store)
  , m_metMakerTool(nullptr)
  , m_metSysTool(nullptr)
  { }

  //______________________________________________________________________________
  ETmissHandler::~ETmissHandler()
  {
    SafeDelete(m_metMakerTool);
    SafeDelete(m_metSysTool);
  }

  //______________________________________________________________________________
  EL::StatusCode ETmissHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    m_isAFII = config.getBool("IsAFII", false);

    // Read in configuration information
    m_containerName = config.getStr(m_name+".ContainerName");
    m_truthName     = config.getStr(m_name+".TruthContainerName", "MET_Truth");

    m_assocMapName  = config.getStr(m_name+".METAssociactionMapName", "METAssoc_AntiKt4EMTopo");
    m_coreName      = config.getStr(m_name+".METCoreName", "METAssoc_AntiKt4EMTopo");

    m_metCst        = config.getStrV(m_name+".METCST");
    m_metTypes      = config.getStrV(m_name+".METTypes");
    
    // METMaker tool
    m_metMakerTool = new met::METMaker("METMaker");   // Defining m_metMakerTool as re-builder of ETmiss tool
    
    CP_CHECK(m_name, m_metMakerTool->setProperty("CustomJetJvtCut", config.getNum("JetHandler.Selection.JVT", 0.941) ));
    if (m_metMakerTool->initialize().isFailure())
      fatal("Failed to initialize METMakerTool");

    // METSystematics tool
    m_metSysTool = new met::METSystematicsTool("METSystematicsTool");
    CP_CHECK(m_name, m_metSysTool->setProperty("ConfigJetTrkFile", "JetTrackSyst.config"                                                   ));
    CP_CHECK(m_name, m_metSysTool->setProperty("JetColl"         , config.getStr("JetHandler.JetContainerName", "AntiKt4EMTopoJets").Data()));
    // AtlFast-II (Fast Simulation)
    if (m_isAFII) {
      CP_CHECK(m_name, m_metSysTool->setProperty("ConfigPrefix", "METUtilities/data15_13TeV/Dec15v2"                                       ));
      CP_CHECK(m_name, m_metSysTool->setProperty("ConfigSoftTrkFile", "TrackSoftTerms_afii.config"                                         ));
    }
    // Full Simulation
    else {
      CP_CHECK(m_name, m_metSysTool->setProperty("ConfigPrefix", "METUtilities/data15_13TeV/Dec15v1"                                       ));
      CP_CHECK(m_name, m_metSysTool->setProperty("ConfigSoftTrkFile", "TrackSoftTerms.config"                                              ));
    }
    if (m_metSysTool->initialize().isFailure())
      fatal("Failed to initialize METSystematicTool");

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer ETmissHandler::getCorrectedContainer()
  {
    xAOD::MissingETContainer dummyContainer;
    return dummyContainer;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer ETmissHandler::getCorrectedContainer(const xAOD::PhotonContainer   *photons  ,
                                                                const xAOD::JetContainer      *jets     ,
                                                                const xAOD::ElectronContainer *electrons,
                                                                const xAOD::MuonContainer     *muons    )
  {
    // Get Shallow copy from TEvent/TStore
    bool calib = false;
    // second argument false --> make empty raw contianer (not shallow copy of xAOD), since it's rebuilt below
    xAOD::MissingETContainer shallowContainer = getShallowContainer(calib, false);
    if (calib) {
      xAOD::MissingETContainer corrected(SG::VIEW_ELEMENTS);
      corrected.push_back(shallowContainer["CST"]);
      corrected.push_back(shallowContainer["TST"]);
      return corrected;
    }

    // Retrieve the MET association map: Needed for METMaker tool
    const xAOD::MissingETAssociationMap *metMap = nullptr;
    if (m_event->retrieve(metMap, m_assocMapName).isFailure())
      fatal("Unable to retrieve " + m_assocMapName + " from TEvent");

    // Retrieve the MET core container: Needed for METMaker tool
    const xAOD::MissingETContainer *coreMet  = nullptr;
    if (m_event->retrieve(coreMet, m_coreName).isFailure())
      fatal("Unable to retrieve " + m_coreName + " from TEvent");

    // For MET rebuilding, need access to the actual container linked with the auxdata
    // which is already in TStore (HgammaHandler magic)
    TString shallowName = "Shallow" + m_containerName + m_sysName;
    xAOD::MissingETContainer *shallowMet  = nullptr;
    if (m_store->retrieve(shallowMet, shallowName.Data()).isFailure())
      fatal("Unable to retrieve ShallowMET from TEvent");

    // Reset the MET map before each building
    metMap->resetObjSelectionFlags();

    // Rebuild ETmiss RefGamma, RefEle, and Muons terms
    HG_CHECK(m_name, m_metMakerTool->rebuildMET("RefGamma", xAOD::Type::Photon  , shallowMet, photons  , metMap));
    HG_CHECK(m_name, m_metMakerTool->rebuildMET("RefEle"  , xAOD::Type::Electron, shallowMet, electrons, metMap));
    HG_CHECK(m_name, m_metMakerTool->rebuildMET("Muons"   , xAOD::Type::Muon    , shallowMet, muons    , metMap));

    // Jet and Soft Terms: CST (no JVT cut)
    HG_CHECK(m_name, m_metMakerTool->rebuildJetMET("RefJet", "SoftClus", "PVSoftTrk", shallowMet, jets, coreMet, metMap, false));

    // Apply possible systematic uncertainty shifts: CST
    if (m_isMC) {
      xAOD::MissingET *softClusMet = (*shallowMet)["SoftClus"];
      if (softClusMet == nullptr)
        fatal("Couldn't retrieve SoftClus from shallowMet, exiting!");
      CC_CHECK(m_name, m_metSysTool->applyCorrection(*softClusMet));
    }

    // Rebuild full MET: CST
    HG_CHECK(m_name, m_metMakerTool->buildMETSum("CST", shallowMet, MissingETBase::Source::LCTopo));

    // Corrected MET container to return: Saved only RefGamma, RefEle, Muons and CST
    // RefJet, SoftClus and PVSoftTrk not saved because we will only save them for TST definition
    xAOD::MissingETContainer corrected(SG::VIEW_ELEMENTS);
    corrected.push_back((*shallowMet)["RefGamma"]);
    corrected.push_back((*shallowMet)["RefEle"]);
    corrected.push_back((*shallowMet)["Muons"]);
    corrected.push_back((*shallowMet)["CST"]);

    // We clear RefJet, SoftClus and PVSoftTrk from shallowMet
    for (auto metCST: *shallowMet) {
      if (std::find(m_metCst.begin(), m_metCst.end(), metCST->name().c_str()) != m_metCst.end())
	metCST->clear();
    }
    
    // Jet and Soft Terms: TST (JVT cut)
    HG_CHECK(m_name, m_metMakerTool->rebuildJetMET("RefJet", "SoftClus", "PVSoftTrk", shallowMet, jets, coreMet, metMap, true));

    // Apply possible systematic uncertainty shifts: TST
    if (m_isMC) {
      xAOD::MissingET *softTrkMet = (*shallowMet)["PVSoftTrk"];
      if (softTrkMet == nullptr)
        fatal("Couldn't retrieve PVSoftTrk from shallowMet, exiting!");
      CC_CHECK(m_name, m_metSysTool->applyCorrection(*softTrkMet));
    }

    // Rebuild full MET: TST
    HG_CHECK(m_name, m_metMakerTool->buildMETSum("TST", shallowMet, MissingETBase::Source::Track ));

    // RefJet, SoftClus and PVSoftTrk saved into the MET container to return
    int check = 0;
    for (auto metCST: *shallowMet) {
      if (std::find(m_metCst.begin(), m_metCst.end(), metCST->name().c_str()) != m_metCst.end())
	{
	  ++check;
	  if ( check > 3 ) corrected.push_back(metCST);
	}
    }
    // TST also saved into the MET container to return
    corrected.push_back((*shallowMet)["TST"]);

    // Return MET container: RefGamma, RefEle, Muons             => Common to CST and TST 
    //                       CST                                 => Only CST MET saved (needed an extra shallowMet container for RefJet and SoftClus)
    //                       RefJet, PVSoftTrk, SoftClus and TST => From TST
    return corrected;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer ETmissHandler::applySelection(xAOD::MissingETContainer &container)
  {
    xAOD::MissingETContainer selected(SG::VIEW_ELEMENTS);
    for (auto met: container) {
      // Limit MET to the types specified in config (TST, CST, ...)
      if (std::find(m_metTypes.begin(), m_metTypes.end(), met->name().c_str()) == m_metTypes.end())
        continue;
      
      selected.push_back(met);
    }

    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode ETmissHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    bool isAffected = false;
    for (auto var: sys) {
      if (m_metSysTool->isAffectedBySystematic(var)) {
        isAffected = true;
        break;
      }
    }

    if (isAffected) {
      if (m_isMC) CP_CHECK(m_name, m_metSysTool->applySystematicVariation(sys));
    } else {
      if (m_isMC) CP_CHECK(m_name, m_metSysTool->applySystematicVariation(CP::SystematicSet()));
    }

    // For MET, always make new container since shifting hard objects indirectly affects result
    m_sysName = sys.name() == "" ? "" : "_"+sys.name();

    return CP::SystematicCode::Ok;
  }

}
