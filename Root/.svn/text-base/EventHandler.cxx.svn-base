#include "HGamAnalysisFramework/EventHandler.h"

#include "HGamAnalysisFramework/HGamVariables.h"
#include "HGamAnalysisFramework/JetHandler.h"

#include "PhotonVertexSelection/PhotonVertexHelpers.h"

#include "xAODEventShape/EventShape.h"
#include "GoodRunsLists/TGRLCollection.h"
#include "VertexPositionReweighting/VertexPositionReweightingTool.h"

#include "TrigEgammaMatchingTool/TrigEgammaMatchingTool.h"

namespace HG {

  SG::AuxElement::Decorator<float> EventHandler::trigSF("SF_trig");
  SG::AuxElement::Decorator<unsigned int> EventHandler::RandomRunNumber("RandomRunNumber");

  //______________________________________________________________________________
  EventHandler::EventHandler(xAOD::TEvent *event, xAOD::TStore *store)
  : m_event(event)
  , m_store(store)
  , m_grl(nullptr)
//  , m_pileupRW(nullptr)
  , m_vtxRW(nullptr)
  , m_configTool(nullptr)
  , m_trigDecTool(nullptr)
  , m_trigMuonMatchTool(nullptr)
  , m_trigElectronMatchTool(nullptr)
  , m_trigMuonScaleFactors(nullptr)
  , m_trigElectronScaleFactors(nullptr)
  , m_trigElectronMCEfficiency(nullptr)
  , m_checkGRL(false)
  , m_checkTile(false)
  , m_checkLAr(false)
  , m_checkCore(false)
  , m_checkBkg(false)
  , m_checkTrig(false)
  , m_checkSCT(false)
  { }

  //______________________________________________________________________________
  EventHandler::~EventHandler()
  {
    SafeDelete(m_grl);
//    SafeDelete(m_pileupRW);
    SafeDelete(m_vtxRW);
    SafeDelete(m_configTool);
    SafeDelete(m_trigDecTool);
    SafeDelete(m_trigMuonMatchTool);
    SafeDelete(m_trigMuonScaleFactors);
    SafeDelete(m_trigElectronMatchTool);
    SafeDelete(m_trigElectronScaleFactors);
    SafeDelete(m_trigElectronMCEfficiency);

    for (auto dec: m_trigDec) SafeDelete(dec.second);
    for (auto dec: m_trigAcc) SafeDelete(dec.second);
  }

  //______________________________________________________________________________
  EL::StatusCode EventHandler::initialize(Config &config)
  {
    const char *APP_NAME = "HG::EventHandler";

    // General options
    const xAOD::EventInfo *eventInfo = nullptr;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }

    m_isMC    = eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION);
    m_isMxAOD = config.getBool("IsMxAOD", false);
    m_jvt     = config.getNum("JetHandler.Selection.JVT", 0.941);

    m_is50ns = config.getBool("Is50ns", false);
    m_isAFII = config.getBool("IsAFII", false);
    m_mcType = config.getStr("MonteCarloType", "MC15a");

    if (m_isMC) {
      int mcid = eventInfo->mcChannelNumber();
      TString option = TString::Format("MonteCarloType.%d", mcid);
      if (config.isDefined(option))
        m_mcType = config.getStr(option);
    }

    // GRL selection
    std::vector<std::string> vecGRL;
    for (auto grl: config.getStrV("EventHandler.GRL"))
      vecGRL.push_back(PathResolverFindCalibFile(grl.Data()));

    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    CHECK(m_grl->setProperty("GoodRunsListVec", vecGRL));

    if (m_grl->initialize().isFailure())
      fatal("Failed to initialize GRL tool");

    // Pileup weighting
    std::vector<std::string> confFiles;
    std::vector<std::string> lcalcFiles;
    if (m_is50ns) {
      for (TString val: config.getStrV("EventHandler.PRW.ConfigFiles50ns"))
        confFiles.push_back(val.Data());

      for (TString val: config.getStrV("EventHandler.PRW.LumiCalcFiles50ns"))
        lcalcFiles.push_back(val.Data());
    } else {
      for (TString val: config.getStrV("EventHandler.PRW.ConfigFiles"+m_mcType))
        confFiles.push_back(val.Data());

      for (TString val: config.getStrV("EventHandler.PRW.LumiCalcFiles"))
        lcalcFiles.push_back(val.Data());
    }

    m_prwSF     = config.getNum("EventHandler.PRW.DataScaleFactor"    , 0.862069);
    int defChan = config.getNum("EventHandler.PRW.DefaultChannel"+m_mcType, 341000);

    double prwSFup     = config.getNum("EventHandler.PRW.DataScaleFactorUP"  , 0.917431);
    double prwSFdown   = config.getNum("EventHandler.PRW.DataScaleFactorDOWN", 0.813008);

    // m_pileupRW = new CP::PileupReweightingTool("Pileup");
    // CHECK(m_pileupRW->setProperty("ConfigFiles"    , confFiles ));
    // CHECK(m_pileupRW->setProperty("LumiCalcFiles"  , lcalcFiles));
    // CHECK(m_pileupRW->setProperty("DataScaleFactor", m_prwSF   ));
    // CHECK(m_pileupRW->setProperty("DataScaleFactorUP"  , prwSFup  ));
    // CHECK(m_pileupRW->setProperty("DataScaleFactorDOWN", prwSFdown));
    // CHECK(m_pileupRW->setProperty("DefaultChannel" , defChan   ));
    // if (m_pileupRW->initialize().isFailure())
    //   fatal("Failed to initialize PRW tool");

    m_pileupRW.make("CP::PileupReweightingTool/HGam");
    CHECK(m_pileupRW.setProperty("ConfigFiles"    , confFiles ));
    CHECK(m_pileupRW.setProperty("LumiCalcFiles"  , lcalcFiles));
    CHECK(m_pileupRW.setProperty("DataScaleFactor", m_prwSF   ));
    CHECK(m_pileupRW.setProperty("DataScaleFactorUP"  , prwSFup  ));
    CHECK(m_pileupRW.setProperty("DataScaleFactorDOWN", prwSFdown));
    CHECK(m_pileupRW.setProperty("DefaultChannel" , defChan   ));
    if (m_pileupRW.retrieve().isFailure())
      fatal("Failed to initialize PRW tool");

    // Vertex weighting
    m_vtxRW = new CP::VertexPositionReweightingTool("VertexPosition");
    CHECK(m_vtxRW->setProperty("DataMean" , -20.0));
    CHECK(m_vtxRW->setProperty("DataSigma",  42.0));

    if (m_vtxRW->initialize().isFailure())
      fatal("Failed to initialize vertex RW tool");

    if (config.getBool("EventHandler.CheckDuplicate" , true)) m_checkDuplic = true;
    if (config.getBool("EventHandler.CheckGRL"       , true)) m_checkGRL    = true;
    if (config.getBool("EventHandler.CheckTile"      , true)) m_checkTile   = true;
    if (config.getBool("EventHandler.CheckLAr"       , true)) m_checkLAr    = true;
    if (config.getBool("EventHandler.CheckCore"      , true)) m_checkCore   = true;
    if (config.getBool("EventHandler.CheckBackground", true)) m_checkBkg    = true;
    if (config.getBool("EventHandler.CheckVertex"    , true)) m_checkVertex = true;
    if (config.getBool("EventHandler.CheckTriggers"  , true)) m_checkTrig   = true;
    if (config.getBool("EventHandler.CheckSCT"       , true)) m_checkSCT    = true;

    m_truthPtclName = config.getStr("TruthParticles.ContainerName","TruthParticle");

    //if (not m_isMC) m_checkDalitz = false;
    if (m_isMC) m_checkGRL = false;

    m_requiredTriggers = config.getStrV("EventHandler.RequiredTriggers");

    // Trigger decision tool
    if (not m_isMxAOD) {
      m_configTool = new TrigConf::xAODConfigTool("xAODConfigTool");
      ToolHandle<TrigConf::ITrigConfigTool> configHandle(m_configTool);
      if (configHandle->initialize().isFailure()) {
        fatal("Failed to initialize trigger config handle");
      }

      m_trigDecTool = new Trig::TrigDecisionTool("TrigDecisionTool");
      CHECK(m_trigDecTool->setProperty("ConfigTool"     , configHandle  ));
      CHECK(m_trigDecTool->setProperty("TrigDecisionKey","xTrigDecision"));

      if (m_trigDecTool->initialize().isFailure()) {
        fatal("Failed to initialize trigger decision tool");
      }

      // Set up trigger matching map
      for (auto trig: m_requiredTriggers) {
        m_trigMatch[trig] = config.getStr("EventHandler.TriggerMatch."+trig, "");
      }

      // Trigger matching for muons
      m_trigMuonMatchTool = new Trig::TrigMuonMatching("TrigMuonMatching");
      CHECK(m_trigMuonMatchTool->setProperty("TriggerTool", ToolHandle<Trig::TrigDecisionTool>(m_trigDecTool)));
      if (m_trigMuonMatchTool->initialize().isFailure()) {
        fatal("Failed to initialize TrigMuonMatchingTool");
      }

      m_trigElectronMatchTool = new Trig::TrigEgammaMatchingTool("TrigEgammaMatchingTool");

      // Determine if triggers only apply to certain runs
      for (TString trig: m_requiredTriggers) {
        TString first = "", last = "";
        int index = -1;
        StrV ranges;
        if (config.isDefined("EventHandler.RunNumbers."+trig))
          ranges = config.getStrV("EventHandler.RunNumbers."+trig);
        for (TString range: ranges) {
          if (range.Contains("-")) {
            // If it's a range, figure out which runs within the range are in GRL
            index  = range.Index("-");
            first  = range(0, index);
            last = range(index+1, range.Length());
            for (int run = first.Atoi(); run <= last.Atoi(); ++run)
              if (m_grl->getGRLCollection()->HasRun(run))
                m_trigRunNumbers[trig].insert(run);
          } else {
            // If it's just a single run, push it back
            m_trigRunNumbers[trig].insert(range.Atoi());
          }
        }
      }
    }
    
    m_trigMuonScaleFactors = new CP::MuonTriggerScaleFactors("TrigMuonScaleFactors");
    CHECK(m_trigMuonScaleFactors->setProperty("MuonQuality", TString(config.getStr("MuonHandler.Selection.PID", "FRED"))             .Data()));
    CHECK(m_trigMuonScaleFactors->setProperty("Isolation",   TString(config.getStr("MuonHandler.Selection.Isolation", "IsoGradient")).Data()));


    if (! m_trigMuonScaleFactors->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonTriggerScaleFactors Tool. Exiting." );
    }
    

    // Electron trigger scale factor per electron
    m_trigElectronScaleFactors = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyTrigCorrectionTool");
    TString bcid = m_is50ns ? "50ns" : "25ns";
    std::string file_trigSF   = config.getStr("EventHandler.ScaleFactor.TrigCorrectionFileName"+bcid).Data();
    std::vector< std::string > correctionFilesTrigSF;
    correctionFilesTrigSF.push_back(file_trigSF);
    CHECK(m_trigElectronScaleFactors->setProperty("CorrectionFileNameList", correctionFilesTrigSF));
    CHECK(m_trigElectronScaleFactors->setProperty("ForceDataType", m_isAFII ? 3 : 1));

    if (m_trigElectronScaleFactors->initialize().isFailure()) {
      fatal("Failed to initialize AsgElectronEfficiencyTrigCorrectionTool");
    }


    m_trigElectronMCEfficiency = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyTrigMCEffCorrectionTool");
    std::string file_trigMCEff   = config.getStr("EventHandler.ScaleFactor.TrigMCEffCorrectionFileName"+bcid).Data();
    std::vector< std::string > correctionFilesTrigMCEff;
    correctionFilesTrigMCEff.push_back(file_trigMCEff);
    CHECK(m_trigElectronMCEfficiency->setProperty("CorrectionFileNameList", correctionFilesTrigMCEff));
    CHECK(m_trigElectronMCEfficiency->setProperty("ForceDataType", m_isAFII ? 3 : 1));

    if (m_trigElectronMCEfficiency->initialize().isFailure()) {
      fatal("Failed to initialize AsgElectronEfficiencyTrigMCEffCorrectionTool");
    }


    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  double EventHandler::mcWeight()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }

    if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) {
      double mcweight = 1.0;

      const std::vector<float> weights = eventInfo->mcEventWeights();
      if (weights.size() > 0) mcweight = weights[0];

      return mcweight;
    }

    return 1.0;
  }

  //______________________________________________________________________________
  double EventHandler::vertexWeight()
  {
    if (var::vertexWeight.exists())
      return var::vertexWeight();

    var::vertexWeight.setValue(1.0);

    if (not m_isMC)
      return var::vertexWeight();

    double vtxWeight = 0.0;
    switch (m_vtxRW->getWeight(vtxWeight)) {
      case CP::CorrectionCode::Error:
        Error("vertexWeight()", "m_vtxWeight->getWeight returned error, returning 1.0");
        var::vertexWeight.setValue(1.0);
        break;
      case CP::CorrectionCode::OutOfValidityRange:
        var::vertexWeight.setValue(1.0);
        Error("vertexWeight()", "m_vtxWeight->getWeight out of validity range, returning 1.0");
        break;
      case CP::CorrectionCode::Ok:
        var::vertexWeight.setValue(vtxWeight);
    }

    return var::vertexWeight();
  }

  //______________________________________________________________________________
  double EventHandler::pileupWeight()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }

    var::pileupWeight.setValue(1.0);

    if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION))
      var::pileupWeight.setValue(m_pileupRW->getCombinedWeight(*eventInfo));

    return var::pileupWeight();
  }


  //______________________________________________________________________________
  double EventHandler::triggerPrescaleWeight(TString triggerList, bool muDependent)
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }

    if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION))
      return 1.0;
    else
      return m_pileupRW->getDataWeight(*eventInfo,triggerList,muDependent);
  }

  //______________________________________________________________________________
  double EventHandler::triggerPrescale(TString trigger)
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }
    
    if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION))
      return 1.0;
    else{
      auto cg = m_trigDecTool->getChainGroup(std::string(trigger.Data()));
      return cg->getPrescale();
    }
  }
  
  //______________________________________________________________________________
  double EventHandler::integratedLumi()
  {
    return m_pileupRW->GetIntegratedLumi();
  }

  //______________________________________________________________________________
  bool EventHandler::pass()
  {
    if (var::isPassedBasic.exists())
      return var::isPassedBasic();
    
    if (m_checkDuplic && isDuplicate  ()) return false;
    if (m_checkTrig   && !passTriggers()) return false;
    return passDQ();
  }
  
  bool EventHandler::passDQ()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
      HG::fatal("Cannot access EventInfo");
    
    if (m_checkGRL    && !passGRL       (eventInfo)) return false;
    if (m_checkLAr    && !passLAr       (eventInfo)) return false;
    if (m_checkTile   && !passTile      (eventInfo)) return false;
    if (m_checkCore   && !passCore      (eventInfo)) return false;
    if (m_checkBkg    && !passBackground(eventInfo)) return false;
    if (m_checkSCT    && !passSCT       (eventInfo)) return false;
    if (m_checkVertex && !passVertex    (eventInfo)) return false;
    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::isDuplicate()
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
      HG::fatal("Cannot access EventInfo");
    
    static SG::AuxElement::Decorator<char> acc_isDuplicate("isDuplicate");
    if (acc_isDuplicate.isAvailable(*eventInfo))
      return acc_isDuplicate(*eventInfo);

    unsigned int runNumber = eventInfo->runNumber();
    if (m_isMC)  runNumber = eventInfo->mcChannelNumber();

    acc_isDuplicate(*eventInfo) = not m_eventNumberSet[runNumber].insert(eventInfo->eventNumber()).second;
    return acc_isDuplicate(*eventInfo);
  }

  //______________________________________________________________________________
  bool EventHandler::isDalitz()
  {
    if (m_isMC) {
      if (var::isDalitzEvent.exists())
        return var::isDalitzEvent();
    
      const xAOD::TruthParticleContainer *truthParticles = nullptr;
      if (m_event->retrieve(truthParticles, m_truthPtclName.Data()).isFailure())
        HG::fatal("Can't access TruthParticleContainer");

      return HG::isDalitz(truthParticles);
    }
    
    // By default (for data) return false
    return false;
  }

  //______________________________________________________________________________
  double EventHandler::triggerScaleFactor(xAOD::ElectronContainer *Electrons, xAOD::MuonContainer *Muons)
  {
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      HG::fatal("Cannot access EventInfo");
    }

    if ( Electrons  && Electrons->size() > 0 )
    { 
      double _effElectronTrigSF = 0.0;
      getElectronTriggerScaleFactor(*Electrons, _effElectronTrigSF);
      trigSF(*eventInfo) = _effElectronTrigSF;
    }
    else if ( Muons && Muons->size() > 0 )
    { 
      double _effMuonTrigSF = 0.0;
      m_trigMuonScaleFactors->getTriggerScaleFactor(*Muons, _effMuonTrigSF, "HLT_mu20_iloose_L1MU15_OR_HLT_mu50");
      trigSF(*eventInfo) = _effMuonTrigSF;
    }
    else
    {
      HG::fatal("Unrecognised particle type for trigger scale factor, returning 1.0");
      trigSF(*eventInfo) = 1.0;
    }

    return trigSF(*eventInfo);

  }

  //______________________________________________________________________________
  bool EventHandler::passTriggers()
  {
    // Require at least one trigger to be passed
    for (auto trig: m_requiredTriggers) {
      if (passTrigger(trig.Data())) return true;
    }
    return false;
  }

  //______________________________________________________________________________
  StrV EventHandler::getPassedTriggers()
  {
    StrV passedTrigs;
    for (auto trig: m_requiredTriggers) {
      if (passTrigger(trig))
        passedTrigs.push_back(trig);
    }
    return passedTrigs;
  }

  //______________________________________________________________________________
  StrV EventHandler::getRequiredTriggers()
  {
    return m_requiredTriggers;
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_SinglePhoton(const TString &trig,
                                                   const xAOD::Photon &ph)
  {
    std::string str = trig.Data();
    return m_trigElectronMatchTool->matchHLT(&ph, str);
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_DiPhoton(const TString &trig,
                                               const xAOD::Photon &photon1,
                                               const xAOD::Photon &photon2)
  {
    if (passTriggerMatch_SinglePhoton(trig, photon1) &&
        passTriggerMatch_SinglePhoton(trig, photon2))
      return true;

    // Return false if both photons weren't matched to a trigger object
    return false;
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_DiMuon(const TString &trig,
                                             const xAOD::Muon &muon1,
                                             const xAOD::Muon &muon2)
  {
    std::pair<bool, bool> isMatchedMuon1, isMatchedMuon2;
    isMatchedMuon1 = std::make_pair(false, false);
    isMatchedMuon2 = std::make_pair(false, false);

    m_trigMuonMatchTool->matchDimuon(&muon1, &muon2, trig.Data(), isMatchedMuon1, isMatchedMuon2);

    return isMatchedMuon1.first && isMatchedMuon2.first;
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_SingleMuon(const TString &trig,
                                                 const xAOD::Muon &muon)
  {
    return m_trigMuonMatchTool->match(&muon, trig.Data());
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_SingleElectron(const TString &trig,
                                                     const xAOD::Electron &el)
  {
    std::string str = trig.Data();
    return m_trigElectronMatchTool->matchHLT(&el, str);
  }

  //______________________________________________________________________________
  bool EventHandler::passTriggerMatch_DiElectron(const TString &trig,
                                                 const xAOD::Electron &el1,
                                                 const xAOD::Electron &el2)
  {
    return (passTriggerMatch_SingleElectron(trig, el1) &&
	    passTriggerMatch_SingleElectron(trig, el2));
  }

  //______________________________________________________________________________
  double EventHandler::selectedVertexZ()
  {
    if (var::selectedVertexZ.exists())
      return var::selectedVertexZ();

    if (!m_store->contains<ConstDataVector<xAOD::VertexContainer> >("HGamVertices"))
      return -999;

    ConstDataVector<xAOD::VertexContainer> *hgamvertices = nullptr;
    if (m_store->retrieve(hgamvertices, "HGamVertices").isFailure())
      return -999;

    if (hgamvertices && hgamvertices->size() > 0) {
      var::selectedVertexZ.setValue((*hgamvertices)[0]->z());
      return var::selectedVertexZ();
    }

    return -999;
  }

  //______________________________________________________________________________
  double EventHandler::hardestVertexZ()
  {
    if (var::hardestVertexZ.exists())
      return var::hardestVertexZ();

    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Cannot access PrimaryVertices");

    const xAOD::Vertex *hardest = xAOD::PVHelpers::getHardestVertex(vertices);

    if (hardest) {
      var::hardestVertexZ.setValue(hardest->z());
      return var::hardestVertexZ();
    }

    return -999;
  }

  //______________________________________________________________________________
  double EventHandler::eventShapeDensity()
  {
    if (var::eventShapeDensity.exists())
      return var::eventShapeDensity();

    const xAOD::EventShape *eventShape = nullptr;
    if (m_event->retrieve(eventShape, "Kt4EMTopoEventShape").isFailure())
      HG::fatal("Couldn't retrieve Kt4EMTopoEventShape from TEvent");

    double rho;
    eventShape->getDensity(xAOD::EventShape::Density, rho);

    var::eventShapeDensity.setValue(rho);

    return rho;
  }

  //______________________________________________________________________________
  double EventHandler::mu()
  {
    if (var::mu.exists())
      return var::mu();

    const xAOD::EventInfo *eventInfo = nullptr;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
      HG::fatal("Cannot access EventInfo in mu()");

    double mu = eventInfo->averageInteractionsPerCrossing();
    if (not m_isMC)
      mu = m_pileupRW->getCorrectedMu(*eventInfo)*m_prwSF;

    var::mu.setValue(mu);

    return mu;
  }

  //______________________________________________________________________________
  int EventHandler::runNumber()
  {
    const xAOD::EventInfo *eventInfo = nullptr;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
      HG::fatal("Cannot access EventInfo in runNumber()");

    int RunNumber = eventInfo->runNumber();
    if (m_isMC) {
      if (RandomRunNumber.isAvailable(*eventInfo)) {
        RunNumber = RandomRunNumber(*eventInfo);
      } else {
        RunNumber = m_pileupRW->getRandomRunNumber(*eventInfo);
        RandomRunNumber(*eventInfo) = RunNumber;
      }
    }

    return RunNumber;
  }

  //______________________________________________________________________________
  int EventHandler::numberOfPrimaryVertices()
  {
    if (var::numberOfPrimaryVertices.exists())
      return var::numberOfPrimaryVertices();

    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Cannot access PrimaryVertices");

    int NPV = 0;
    for (auto vertex: *vertices) {
      if (vertex->vertexType() == xAOD::VxType::PriVtx ||
          vertex->vertexType() == xAOD::VxType::PileUp)
        NPV++;
    }

    var::numberOfPrimaryVertices.setValue(NPV);

    return NPV;
  }


  //______________________________________________________________________________
  xAOD::Photon* EventHandler::getClosestHLTObject(const TString &trig, const xAOD::Photon &pho)
  {
#ifndef __DC14__
    if (m_trigElectronMatchTool->matchHLT(&pho, trig.Data()))
    {
      return (xAOD::Photon*) m_trigElectronMatchTool->closestHLTObject(&pho, trig.Data());
    }
#endif
    return NULL;
  }

  //______________________________________________________________________________
  xAOD::Electron* EventHandler::getClosestHLTObject(const TString &trig, const xAOD::Electron &el)
  {
#ifndef __DC14__
    if (m_trigElectronMatchTool->matchHLT(&el, trig.Data()))
    {
      return (xAOD::Electron*) m_trigElectronMatchTool->closestHLTObject(&el, trig.Data());
    }
#endif
    return NULL;
  }

  //______________________________________________________________________________
  xAOD::Muon* EventHandler::getClosestHLTObject(const TString &/*trig*/, const xAOD::Muon &/*mu*/)
  {
    // This function intentionally left blank.  No Muon tool method to retrieve closest HLT object
    return NULL;
  }

  //______________________________________________________________________________
  bool EventHandler::passGRL(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        !m_grl->passRunLB(*eventInfo))
      return false;

    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passTile(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error)
      return false;

    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passLAr(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error)
      return false;

    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passBackground(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->isEventFlagBitSet(xAOD::EventInfo::Background, 20))
      return false;

    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passCore(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18))
      return false;

    return true;
  }

  //______________________________________________________________________________
  bool EventHandler::passVertex(const xAOD::EventInfo * /*eventInfo*/)
  {
    // Retrieve PV collection from TEvent
    const xAOD::VertexContainer* vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure()) {
      HG::fatal("Couldn't retrieve PrimaryVertices, exiting.");
      return false;
    }

    for (auto vertex: *vertices)
      if (vertex->vertexType() == xAOD::VxType::VertexType::PriVtx ||
          vertex->vertexType() == xAOD::VxType::VertexType::PileUp)
        return true;

    return false;
  }


  //______________________________________________________________________________
  bool EventHandler::passSCT(const xAOD::EventInfo *eventInfo)
  {
    if (!eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION) &&
        eventInfo->errorState(xAOD::EventInfo::SCT) == xAOD::EventInfo::Error)
      return false;

    return true;
  }

  //______________________________________________________________________________
  void EventHandler::getElectronTriggerScaleFactor(xAOD::ElectronContainer &electrons, double &trigSF)
  {

    double _eleTrigSF = 1.;
    double _eleTrigMCEff = 1.;

    double _numSF = 1.;
    double _denomSF = 1.;

    for (auto ele : electrons )
    {
       _eleTrigSF = 0.;
       _eleTrigMCEff = 0.;

       double cl_eta = 10.;
       const xAOD::CaloCluster* cluster = ele->caloCluster();
       if(cluster)  cl_eta = cluster->eta();

       if (m_isMC && std::abs(cl_eta) < 2.47 && ele->pt() >= 15000.)
       {
          if(m_trigElectronScaleFactors->getEfficiencyScaleFactor     (*ele, _eleTrigSF   ) == CP::CorrectionCode::Error)
                 fatal("ElectronEfficiencyTrigCorrection returned CP::CorrectionCode::Error");

          if(m_trigElectronMCEfficiency->getEfficiencyScaleFactor     (*ele, _eleTrigMCEff   ) == CP::CorrectionCode::Error)
                 fatal("ElectronEfficiencyTrigMCEffCorrection returned CP::CorrectionCode::Error");
       }

       _numSF *=  (1. - (_eleTrigSF*_eleTrigMCEff));
       _denomSF *=  (1. - _eleTrigMCEff);

    }
    if ( (1.- _denomSF) > 0. )
       trigSF = (1.- _numSF) / (1.- _denomSF);

    return;
  }

  //______________________________________________________________________________
  EL::StatusCode EventHandler::writeVars(TString eventName)
  {
    if (eventName == "")
      eventName = HG::VarHandler::getInstance()->getEventInfoName();

    const xAOD::EventInfo *eventInfo = 0;
    if (m_store->retrieve(eventInfo, eventName.Data()).isFailure()) {
      if (m_event->retrieve(eventInfo, eventName.Data()).isFailure()) {
        HG::fatal("EventHandler::write() cannot access " + eventName + ". Exiting.");
      }
    }

    // create containers
    xAOD::EventInfo *output       = new xAOD::EventInfo();
    xAOD::AuxInfoBase * outputAux = new xAOD::AuxInfoBase();
    output->setStore(outputAux);

    *output = *eventInfo;

    // record event info (yes, can be done before setting the actual values)
    if (!m_event->record(output, eventName.Data())) { return EL::StatusCode::FAILURE; }

    eventName += "Aux.";
    if (!m_event->record(outputAux, eventName.Data())) { return EL::StatusCode::FAILURE; }

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  EL::StatusCode EventHandler::writeEventInfo()
  {
    if (m_event->copy("EventInfo").isFailure())
      Warning("EventHandler::writeEventInfo()", "Couldn't copy EventInfo to output.");

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  CP::SystematicCode EventHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    bool isAffected = false;
    for (auto var: sys) {
      if (  m_pileupRW->isAffectedBySystematic(var) || 
            m_trigMuonScaleFactors->isAffectedBySystematic(var) ||
            m_trigElectronScaleFactors->isAffectedBySystematic(var) ||
            m_trigElectronMCEfficiency->isAffectedBySystematic(var) 
         ) {
        isAffected = true;
        break;
      }
    }

    if (isAffected) {
      CP_CHECK("EventHandler", m_pileupRW->applySystematicVariation(sys));
      CP_CHECK("EventHandler", m_trigMuonScaleFactors->applySystematicVariation(sys));
      CP_CHECK("EventHandler", m_trigElectronScaleFactors->applySystematicVariation(sys));
      CP_CHECK("EventHandler", m_trigElectronMCEfficiency->applySystematicVariation(sys));
      m_sysName = sys.name() == "" ? "" : "_"+sys.name();
    } else {
      CP_CHECK("EventHandler", m_pileupRW->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK("EventHandler", m_trigMuonScaleFactors->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK("EventHandler", m_trigElectronScaleFactors->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK("EventHandler", m_trigElectronMCEfficiency->applySystematicVariation(CP::SystematicSet()));
      m_sysName = "";
    }
    
    return CP::SystematicCode::Ok;
  }

} // namespace HG
