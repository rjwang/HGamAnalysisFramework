#ifndef HGamAnalysisFramework_EventHandler_H
#define HGamAnalysisFramework_EventHandler_H

#include "HGamAnalysisFramework/HgammaIncludes.h"

namespace Trig {
  class TrigEgammaMatchingTool;
}

namespace CP {
  class VertexPositionReweightingTool;
}

namespace HG {

  class EventHandler {
  protected:
    xAOD::TEvent              *m_event;
    xAOD::TStore              *m_store;
    GoodRunsListSelectionTool *m_grl;
    // CP::PileupReweightingTool *m_pileupRW;
    asg::AnaToolHandle<CP::IPileupReweightingTool> m_pileupRW;
    CP::VertexPositionReweightingTool *m_vtxRW;
    TrigConf::xAODConfigTool  *m_configTool;
    Trig::TrigDecisionTool    *m_trigDecTool;
    Trig::TrigMuonMatching    *m_trigMuonMatchTool;
    Trig::TrigEgammaMatchingTool *m_trigElectronMatchTool;
    CP::MuonTriggerScaleFactors  *m_trigMuonScaleFactors;
    AsgElectronEfficiencyCorrectionTool  *m_trigElectronScaleFactors;
    AsgElectronEfficiencyCorrectionTool  *m_trigElectronMCEfficiency;
    std::map<TString, SG::AuxElement::Decorator<char>* >     m_trigDec;
    std::map<TString, SG::AuxElement::ConstAccessor<char>* > m_trigAcc;
    std::map<TString, std::set<int> > m_trigRunNumbers;

    bool                       m_isMC;
    bool                       m_isAFII;
    bool                       m_isMxAOD;
    bool                       m_is50ns;
    TString                    m_mcType;
    bool                       m_checkDuplic;
    bool                       m_checkGRL;
    bool                       m_checkTile;
    bool                       m_checkLAr;
    bool                       m_checkCore;
    bool                       m_checkBkg;
    bool                       m_checkVertex;
    bool                       m_checkTrig;
    bool                       m_checkSCT;
    StrV                       m_requiredTriggers;
    std::map<TString, TString> m_trigMatch;
    TString                    m_truthPtclName;
    float                      m_jvt;
    float                      m_prwSF;
    std::map<unsigned int, std::set<unsigned int> > m_eventNumberSet;

    TString       m_sysName;

  public:
    static SG::AuxElement::Decorator<unsigned int> RandomRunNumber;
    static SG::AuxElement::Decorator<float> trigSF;


  public:
    EventHandler(xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~EventHandler();

    virtual EL::StatusCode initialize(Config &config);

    // Apply systematic variation to appropriate tools
    virtual CP::SystematicCode applySystematicVariation(const CP::SystematicSet &sys);

    // Global event pass (all selections)
    bool pass();
    
    // DQ selection (GRL+LAr+Tile+Core) + vertex
    bool passDQ();
    
    // Individual event selections
    bool isDalitz();
    bool isDuplicate();
    bool passGRL(const xAOD::EventInfo *eventInfo);
    bool passTile(const xAOD::EventInfo *eventInfo);
    bool passLAr(const xAOD::EventInfo *eventInfo);
    bool passCore(const xAOD::EventInfo *eventInfo);
    bool passBackground(const xAOD::EventInfo *eventInfo);
    bool passVertex(const xAOD::EventInfo *eventInfo);
    bool passSCT(const xAOD::EventInfo *eventInfo);
    bool passTriggers();

    // Event weights
    double mcWeight();
    double pileupWeight();
    double vertexWeight();
    double triggerPrescaleWeight(TString triggerList = "HLT_g10_loose", bool muDependent = true);
    double triggerPrescale(TString triggerList = "HLT_g10_loose");

    // Helper functions
    int numberOfPrimaryVertices();
    double selectedVertexZ();
    double hardestVertexZ();
    double eventShapeDensity();
    double mu();
    int runNumber();
    double integratedLumi();

    inline bool passTrigger(const TString &trig = "");
    bool passTriggerMatch_SinglePhoton(const TString &trig, const xAOD::Photon &photon1);
    bool passTriggerMatch_DiPhoton(const TString &trig, const xAOD::Photon &photon1, const xAOD::Photon &photon2);
    bool passTriggerMatch_SingleMuon(const TString &trig, const xAOD::Muon &muon);
    bool passTriggerMatch_DiMuon(const TString &trig, const xAOD::Muon &muon1, const xAOD::Muon &muon2);
    bool passTriggerMatch_SingleElectron(const TString &trig, const xAOD::Electron &el);
    bool passTriggerMatch_DiElectron(const TString &trig, const xAOD::Electron &el1, const xAOD::Electron &el2);
    StrV getPassedTriggers();
    StrV getRequiredTriggers();
    double triggerScaleFactor(xAOD::ElectronContainer *Electrons, xAOD::MuonContainer *Muons);
    void getElectronTriggerScaleFactor(xAOD::ElectronContainer &Electrons, double &trigSF);


    xAOD::Photon*   getClosestHLTObject(const TString &trig, const xAOD::Photon &pho);
    xAOD::Electron* getClosestHLTObject(const TString &trig, const xAOD::Electron &el);
    xAOD::Muon*     getClosestHLTObject(const TString &trig, const xAOD::Muon &mu);

    template<typename T>
    void storeVar(SG::AuxElement::Decorator<T> &dec, T value);

    template<typename T>
    void storeTruthVar(const char *name, T value);

    template<typename T>
    void storeVar(const char *name, T value);

    template<typename T>
    T getTruthVar(const char *name);

    template<typename T>
    T getVar(const char *name);

    virtual EL::StatusCode writeVars(TString name = "");
    virtual EL::StatusCode writeEventInfo();

    // these are deprecated!!
    template<typename T>
    void storeVariable(const char *name, T value);

    virtual EL::StatusCode write() { return writeEventInfo(); }



  };

  //______________________________________________________________________________
  bool EventHandler::passTrigger(const TString &trig)
  { 
    // Define decorator/accessor if not already available
    if (m_trigDec.find(trig) == m_trigDec.end()) {
      m_trigAcc[trig] = new SG::AuxElement::ConstAccessor<char>(("passTrig_"+trig).Data());
      m_trigDec[trig] = new SG::AuxElement::Decorator    <char>(("passTrig_"+trig).Data());
    }

    const xAOD::EventInfo *eventInfo = nullptr;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure())
      HG::fatal("EventHandler::passTrigger() : Cannot access EventInfo");

    // If the decision tool is defined, use it
    if (m_trigDecTool) {
      // In case it's not checked, default to false
      (*m_trigDec[trig])(*eventInfo) = false;

      // Check if only certain runs should be used for this trigger
      if (m_trigRunNumbers.find(trig) == m_trigRunNumbers.end()) {
        (*m_trigDec[trig])(*eventInfo) = m_trigDecTool->isPassed(trig.Data());
      } else {
        if (m_trigRunNumbers.at(trig).count(runNumber()))
          (*m_trigDec[trig])(*eventInfo) = m_trigDecTool->isPassed(trig.Data());
      }

      return (*m_trigDec[trig])(*eventInfo);
    }

    // If there's no decision tool, check for decision in EventInfo
    if (m_trigAcc[trig]->isAvailable(*eventInfo))
      return (*m_trigAcc[trig])(*eventInfo);

    Error("EventHandler::passTrigger()", "Trigger '%s' could not be checked!", trig.Data());
    fatal("TrigDecisionTool not defined. This is either an MxAOD without this trigger stored, or initialize() went wrong. Exiting.");

    // Should never get here
    return false;
  }

}

#include "HGamAnalysisFramework/EventHandler.hpp"

#endif // HGamAnalysisFramework_EventHandler_H
