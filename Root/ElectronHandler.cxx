#include "HGamAnalysisFramework/ElectronHandler.h"

#include "IsolationSelection/IsolationSelectionTool.h"
#include "IsolationCorrections/IsolationCorrectionTool.h"
#include "ElectronPhotonShowerShapeFudgeTool/ElectronPhotonShowerShapeFudgeTool.h"

namespace HG {
  
  SG::AuxElement::Accessor<float> ElectronHandler::effIDSF("SF_IDeff");
  SG::AuxElement::Accessor<float> ElectronHandler::effRecoSF("SF_Recoeff");
  SG::AuxElement::Accessor<float> ElectronHandler::effIsoSF("SF_Isoeff");
  SG::AuxElement::Accessor<float> ElectronHandler::scaleFactor("scaleFactor");
  SG::AuxElement::Accessor<float> ElectronHandler::Ecalib_ratio("Ecalib_ratio");
  SG::AuxElement::Accessor<float> ElectronHandler::Ereso("Ereso");
  SG::AuxElement::Accessor<float> ElectronHandler::eta_s2("eta_s2");
  SG::AuxElement::Accessor<char>  ElectronHandler::passIPCut("passIPCut");
  SG::AuxElement::Accessor<char>  ElectronHandler::isTight("isTight");
  SG::AuxElement::Accessor<char>  ElectronHandler::isMedium("isMedium");
  SG::AuxElement::Accessor<char>  ElectronHandler::isLoose("isLoose");
  
  //______________________________________________________________________________
  ElectronHandler::ElectronHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store)
  : HgammaHandler(name, event, store)
  , m_electronCalibTool(nullptr)
  , m_isoCorrTool(nullptr)
  , m_electronIDSF(nullptr)
  , m_electronRecoSF(nullptr)
  , m_electronIsoSF(nullptr)
  , m_fudgeTool(nullptr)
  { }

  //______________________________________________________________________________
  ElectronHandler::~ElectronHandler()
  {
    SafeDelete(m_electronCalibTool);
    SafeDelete(m_isoCorrTool);
    SafeDelete(m_electronIDSF);
    SafeDelete(m_electronRecoSF);
    SafeDelete(m_electronIsoSF);
    SafeDelete(m_fudgeTool);

    for (auto sel: m_electronSelectors) SafeDelete(sel.second);
    m_electronSelectors.clear();
    
    for (auto dec: m_pidAcc) SafeDelete(dec.second);
    m_pidAcc.clear();

    for (auto iso: m_isoTools) SafeDelete(iso.second);
    m_isoTools.clear();
    
    for (auto dec: m_isoAcc) SafeDelete(dec.second);
    m_isoAcc.clear();
  }

  //______________________________________________________________________________
  EL::StatusCode ElectronHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    // General configs
    m_isAFII = config.getBool("IsAFII", false);
    m_is50ns = config.getBool("Is50ns", false);
    TString bcid = m_is50ns ? "50ns" : "25ns";

    // egamma energy scale (data) and extra smearing correction (MC)
    m_electronCalibTool = new CP::EgammaCalibrationAndSmearingTool("ElectronCalibrationAndSmearingTool");

    for (TString prop : {"ESModel", "decorrelationModel"}) {
      CP_CHECK(m_name, m_electronCalibTool->setProperty(prop.Data(), config.getStr(m_name+".Calibration."+prop).Data()));
    }
    
    if (m_isAFII)
      CP_CHECK(m_name, m_electronCalibTool->setProperty("useAFII", true));

    if (m_electronCalibTool->initialize().isFailure()) {
      fatal("Failed to initialize EgammaCalibrationAndSmearingTool");
    }

    //electron selection
    m_doPidCut   = config.getBool(m_name+".Selection.ApplyPIDCut", true);
    m_pidCuts    = config.getStrV(m_name+".Selection.PID");
    
    if (m_pidCuts.size() < 1) fatal("You must specify at least one PID criteria for electrons");

    // loop over PID selections
    if (m_doPidCut && m_pidCuts.size() < 1)
      fatal("Electron PID cut requested, but no working point supplied. Exiting!");
    
    for (size_t i = 0; i < m_pidCuts.size(); ++i) {
      TString pid = m_pidCuts[i];
      m_pidAcc[pid] = new SG::AuxElement::Accessor<char>(("is"+m_pidCuts[i]).Data());
      
      if (i == 0) m_defaultPid = pid;
      
      TString cfgFile(config.getStr(m_name+".Selection.ConfigFile."+m_pidCuts[i]+bcid));

      m_electronSelectors[pid] = new AsgElectronLikelihoodTool("ElectronLikelihoodTool");

      CP_CHECK(m_name, m_electronSelectors[pid]->setProperty("primaryVertexContainer",config.getStr("PrimaryVertices.ContainerName","PrimaryVertices").Data()));
      CP_CHECK(m_name, m_electronSelectors[pid]->setProperty("ConfigFile",cfgFile.Data()));
      
      if (!m_electronSelectors[pid]->initialize().isSuccess()) {
        fatal(TString::Format("Failed to initialize %sElectronLikelihoodTool",pid.Data()));
      }
    }
    
    // isolation tools
    m_doIsoCut = config.getBool(m_name+".Selection.ApplyIsoCut", true);
    m_isoCuts  = config.getStrV(m_name+".Selection.IsoCriteria");
    if (m_doIsoCut && m_isoCuts.size() < 1)
      fatal("Isolation cut requested, but no working point supplied. Exiting!");
    for (size_t i = 0; i < m_isoCuts.size(); ++i) {
      HG::Iso::IsolationType iso = getIsoType(m_isoCuts[i]);
      m_isoAcc[iso] = new SG::AuxElement::Accessor<char>(("isIso"+m_isoCuts[i]).Data());

      // first isolation in the list is the default one to apply
      if (i == 0) m_defaultIso = iso;

      m_isoTools[iso] = new CP::IsolationSelectionTool(m_isoCuts[i].Data());
      CP_CHECK(m_name, m_isoTools[iso]->setProperty("ElectronWP", m_isoCuts[i].Data()));

      if (m_isoTools[iso]->initialize().isFailure()) 
        fatal("Failed to initialize IsolationSelectionTool with WP"+m_isoCuts[i]);
      
      CP_CHECK(m_name, m_isoTools[iso]->setIParticleCutsFrom(xAOD::Type::Electron));
    }

    // Isolation correction tool
    m_isoCorrTool = new CP::IsolationCorrectionTool("IsolationCorrectionTool");
    CP_CHECK(m_name, m_isoCorrTool->setProperty("IsMC", isMC()));
    
    if (m_isoCorrTool->initialize().isFailure())
      fatal("Failed to initialize IsolationCorrectionTool");

    
    // efficiency scale factor tool
    m_electronIDSF = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyIDCorrectionTool");
    m_electronRecoSF = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyRecoCorrectionTool");
    m_electronIsoSF = new AsgElectronEfficiencyCorrectionTool("AsgElectronEfficiencyIsoCorrectionTool");
    
    std::string file_ID   = config.getStr(m_name+".ScaleFactor.IDCorrectionFileName"+bcid).Data();
    std::string file_Reco = config.getStr(m_name+".ScaleFactor.RecoCorrectionFileName"+bcid).Data();
    std::string file_Iso = config.getStr(m_name+".ScaleFactor.IsoCorrectionFileName"+bcid).Data();
    std::vector< std::string > correctionFilesID;
    std::vector< std::string > correctionFilesReco;
    std::vector< std::string > correctionFilesIso;
    correctionFilesID.push_back(file_ID);
    correctionFilesReco.push_back(file_Reco);
    correctionFilesIso.push_back(file_Iso);
    
    CP_CHECK(m_name, m_electronIDSF->setProperty("CorrectionFileNameList", correctionFilesID));
    CP_CHECK(m_name, m_electronIDSF->setProperty("ForceDataType", m_isAFII ? 3 : 1)); 
    CP_CHECK(m_name, m_electronRecoSF->setProperty("CorrectionFileNameList", correctionFilesReco));
    CP_CHECK(m_name, m_electronRecoSF->setProperty("ForceDataType", m_isAFII ? 3 : 1)); 
    CP_CHECK(m_name, m_electronIsoSF->setProperty("CorrectionFileNameList", correctionFilesIso));
    CP_CHECK(m_name, m_electronIsoSF->setProperty("ForceDataType", m_isAFII ? 3 : 1)); 


    if (m_electronIDSF->initialize().isFailure()) {
      fatal("Failed to initialize AsgElectronEfficiencyIDCorrectionTool");
    }

    if (m_electronRecoSF->initialize().isFailure()) {
      fatal("Failed to initialize AsgElectronEfficiencyRecoIDCorrectionTool");
    }

    if (m_electronIsoSF->initialize().isFailure()) {
      fatal("Failed to initialize AsgElectronEfficiencyIsoIDCorrectionTool");
    }


    // electron fudge tool - OFF by 
    if (m_isMC && not m_isAFII) {
      m_fudgeTool  = new ElectronPhotonShowerShapeFudgeTool("FudgeTool");
      m_fudgeSet = config.getNum(m_name+".Calibration.FFSet"  , 16);
      CP_CHECK(m_name, m_fudgeTool->setProperty("Preselection", m_fudgeSet));
      if (m_fudgeTool->initialize().isFailure()) {
        fatal("Failed to initialize FudgeTool");
      }
    }

    
    // Read in configuration information
    m_containerName = config.getStr(m_name+".ContainerName", "ElectronCollection");

    m_etaCut      = config.getNum(m_name+".Selection.MaxAbsEta"  , 2.47);
    m_ptCut       = config.getNum(m_name+".Selection.PtPreCutGeV", 25.0)*GeV;

    m_crackReject = config.getBool(m_name+".Selection.ApplyCrackRejection", true);
    m_barrelMax   = config.getNum(m_name+".Selection.BarrelMaxAbsEta"  , 1.37);
    m_endcapMin   = config.getNum(m_name+".Selection.EndcapMinAbsEta"  , 1.52);
    
    m_applyIPCuts   = config.getBool(m_name+".Selection.ApplyIPCuts", false);
    m_d0BySigd0Cut  = config.getNum(m_name+".Selection.d0BySigd0Max", 5.0);
    m_z0Cut         = config.getNum(m_name+".Selection.z0Max", 0.5);

    //OFF by default - not recommended 
    m_doFudge = config.getBool(m_name+".Calibration.DoFudgeFactor", false);
    if (m_isData) m_doFudge = false;
    

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  xAOD::ElectronContainer ElectronHandler::getCorrectedContainer()
  {
    // get the event info
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      fatal("Cannot access EventInfo");
    }
      
    bool calib = false;
    xAOD::ElectronContainer shallowContainer = getShallowContainer(calib);
    if (calib) return shallowContainer;

    //PID on not calibrated objects, SF on calibrated objects, isolation on calibrated objects
    for (auto electron: shallowContainer) {
      calibrateAndSmearElectron(electron, eventInfo, m_electronCalibTool);
      if (m_doFudge) applyFudgeFactor(electron, eventInfo);
      correctIsoLeakage(*electron);
      decoratePID(*electron);
      decorateIso(*electron);
      decorateIPCut(*electron);
      applyScaleFactor(electron, eventInfo);

      eta_s2(*electron) = electron->caloCluster()->etaBE(2);
    }
    
    // sort the electrons
    shallowContainer.sort(comparePt);
    
    return shallowContainer;
  }

  //______________________________________________________________________________
  void ElectronHandler::correctIsoLeakage(xAOD::Electron &electron)
  {
    // // Tool complains if you correct electrons outside of interpolation limits
    // double eta = electron.caloCluster()->eta();
    // if (eta < -2.5 || 2.5 < eta)
    //   return;

    // Tool complains if you correct electrons outside of interpolation limits
    double eta = electron.caloCluster()->etaBE(1);
    if (eta < -2.5 || 2.5 < eta)
      return;

    eta = electron.caloCluster()->etaBE(2);
    if (eta < -2.5 || 2.5 < eta)
      return;

    // Check that the leakage is correctly updated
    if (not m_isoCorrTool->applyCorrection(electron))
      fatal("Couldn't correct electron isolation leakage?");
  }
  
  //______________________________________________________________________________
  xAOD::ElectronContainer ElectronHandler::applySelection(xAOD::ElectronContainer &container)
  {
    xAOD::ElectronContainer selected(SG::VIEW_ELEMENTS);
    for (auto electron: container) {
      // pT and eta cuts
      if (!passPtEtaCuts(electron)) continue;

      // d0/z0 cuts
      if (m_applyIPCuts && !passIPCuts(electron)) continue;

      // PID LH identification
      if (m_doPidCut && !passPIDCut(electron)) continue;

      // isolation cuts
      if (m_doIsoCut && !passIsoCut(electron)) continue;

      selected.push_back(electron);
    }
    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode ElectronHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {
    setVertexCorrected(false);

    bool isAffected = false;
    for (auto var: sys) {
      if (m_electronCalibTool->isAffectedBySystematic(var) ||
          m_electronIDSF->isAffectedBySystematic(var)      ||
          m_electronRecoSF->isAffectedBySystematic(var)    ||
	  m_electronIsoSF->isAffectedBySystematic(var)     ||
          m_isoCorrTool->isAffectedBySystematic(var)       ) {
        isAffected = true;
        break;
      }
    }

    if (isAffected) {
      CP_CHECK(m_name, m_electronCalibTool->applySystematicVariation(sys));
      CP_CHECK(m_name, m_electronIDSF     ->applySystematicVariation(sys));
      CP_CHECK(m_name, m_electronRecoSF   ->applySystematicVariation(sys));
      CP_CHECK(m_name, m_electronIsoSF    ->applySystematicVariation(sys));
      CP_CHECK(m_name, m_isoCorrTool      ->applySystematicVariation(sys));
      m_sysName = sys.name() == "" ? "" : "_"+sys.name();
    } else {
      CP_CHECK(m_name, m_electronCalibTool->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_electronIDSF     ->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_electronRecoSF   ->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_electronIsoSF    ->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_isoCorrTool      ->applySystematicVariation(CP::SystematicSet()));
      m_sysName = "";
    }
    
    return CP::SystematicCode::Ok;
  }

  //______________________________________________________________________________
  void ElectronHandler::calibrateAndSmearElectron(xAOD::Electron *ele,
                                                  const xAOD::EventInfo *evtInfo,
                                                  CP::EgammaCalibrationAndSmearingTool *electronCalibTool) {
    double cl_eta = 10.;
    const xAOD::CaloCluster* cluster = ele->caloCluster();
    if(cluster) cl_eta = cluster->eta();
    
    // Apply smearing 
    if(std::abs(cl_eta) < 2.47 && ele->pt() >= 20000.){
      electronCalibTool->setRandomSeed(evtInfo->eventNumber()*100+ele->index());
      
      // Calibrate the electron
      double E_before = ele->e();
      CP::CorrectionCode cc = electronCalibTool->applyCorrection(*ele);
      if (cc==CP::CorrectionCode::Error)
	Error("calibratedElectron()","Error calibrating current electron");
      if (cc==CP::CorrectionCode::OutOfValidityRange)
	Warning("calibratedElectron()","Current electron has no valid calibration due to out-of-range");
      
      Ecalib_ratio(*ele) = ele->e()/E_before;
#ifndef __DC14__
      Ereso(*ele) = electronCalibTool->getResolution(*ele);
#endif
    }
  }

  //______________________________________________________________________________
  void ElectronHandler::applyFudgeFactor(xAOD::Electron *ele, const xAOD::EventInfo *evtInfo)
  {
    // Don't actually apply the fudging in AFII
    if (m_isAFII)
      return;
    int isMC = evtInfo->eventType(xAOD::EventInfo::IS_SIMULATION);
    if (isMC) {
      CP::CorrectionCode cc = m_fudgeTool->applyCorrection(*ele);
      if (cc==CP::CorrectionCode::Error)
        Error("applyFudgeFactorElectron()","Fudging returned error");
      if (cc==CP::CorrectionCode::OutOfValidityRange)
	Warning("applyFudgeFactorElectron()","Current electron has no valid fudging due to out-of-range");
    }
  }

  
  //______________________________________________________________________________
  void ElectronHandler::applyScaleFactor(xAOD::Electron *ele, const xAOD::EventInfo *evtInfo)
  {
    double _effIDSF = 1.0, _effRecoSF = 1.0, _effIsoSF = 1.0;
    
    double cl_eta = 10.;
    const xAOD::CaloCluster* cluster = ele->caloCluster();
    if(cluster)  cl_eta = cluster->eta();
    
    if (m_isMC && std::abs(cl_eta) < 2.47 && ele->pt() >= 15000.){
      if(m_electronIDSF->getEfficiencyScaleFactor     (*ele, _effIDSF   ) == CP::CorrectionCode::Error)
        fatal("ElectronEfficiencyIDCorrection returned CP::CorrectionCode::Error");
      if(m_electronRecoSF->getEfficiencyScaleFactor   (*ele, _effRecoSF   ) == CP::CorrectionCode::Error)
	fatal("ElectronEfficiencyRecoCorrection returned CP::CorrectionCode::Error");
      if(m_electronIsoSF->getEfficiencyScaleFactor   (*ele, _effIsoSF   ) == CP::CorrectionCode::Error)
	fatal("ElectronEfficiencyRecoCorrection returned CP::CorrectionCode::Error");
    }

    effIDSF(*ele) = _effIDSF;
    effRecoSF(*ele) = _effRecoSF;
    effIsoSF(*ele) = _effIsoSF;
    scaleFactor(*ele) = _effIDSF*_effRecoSF*_effIsoSF;
  }

  //______________________________________________________________________________
  bool ElectronHandler::passPtEtaCuts(const xAOD::Electron *ele) 
  {
    double abs_eta = fabs(eta_s2(*ele));
    if (abs_eta > m_etaCut) return false;
    if (m_crackReject && (abs_eta > m_barrelMax && abs_eta < m_endcapMin)) return false;
    
    if (ele->pt() < m_ptCut) return false;
  
    return true;
  }

  //______________________________________________________________________________
  void ElectronHandler::decorateIso(xAOD::Electron &ele)
  {
    for (auto dec: m_isoAcc) {
      if (m_isoTools[dec.first]->accept(ele))
        (*dec.second)(ele) = true;
      else
        (*dec.second)(ele) = false;
    }
  }

  //______________________________________________________________________________
  bool ElectronHandler::passIsoCut(const xAOD::Electron *ele, HG::Iso::IsolationType iso)
  {
    /// applies isolation cut specified in config file
    if (iso == HG::Iso::Undefined) {
      if (!m_isoAcc[m_defaultIso]->isAvailable(*ele))
        return true;
      return (*m_isoAcc[m_defaultIso])(*ele);
    }
    
    if (m_isoTools.find(iso) != m_isoTools.end()) {
      if (!m_isoAcc[iso]->isAvailable(*ele))
        return true;
      return (*m_isoAcc[iso])(*ele);
    }

    fatal("Electron isolation cut requested that wasn't specified in config file. Exiting.");
    return false;
  }

  //______________________________________________________________________________
  void ElectronHandler::decoratePID(xAOD::Electron &ele)
  {
    for (auto dec: m_pidAcc) {
      if (m_electronSelectors[dec.first]->accept(&ele))
        (*dec.second)(ele) = true;
      else
        (*dec.second)(ele) = false;
    }
  }
  
  //______________________________________________________________________________
  bool ElectronHandler::passPIDCut(const xAOD::Electron *ele, TString pid)
  {

    // applies PID cut specified in config file
    if (pid == "Default") {
      if (!m_pidAcc[m_defaultPid]->isAvailable(*ele))
        return true;
      return (*m_pidAcc[m_defaultPid])(*ele);
    }
    
    if (m_pidAcc.find(pid) != m_pidAcc.end()) {
      if (!m_pidAcc[pid]->isAvailable(*ele))
        return true;
      return (*m_pidAcc[pid])(*ele);
    }

    fatal(TString::Format("%s: PID cut requested that wasn't specified in config file. Exiting.",m_name.Data()));
    return false;
  }

  //______________________________________________________________________________
  bool ElectronHandler::passIPCuts(const xAOD::Electron *ele) 
  {
    if(passIPCut.isAvailable(*ele) && !passIPCut(*ele)) return false;
    return true;
  }

  //______________________________________________________________________________
  void ElectronHandler::decorateIPCut(xAOD::Electron &ele)
  {
    passIPCut(ele) = false;


    const xAOD::TrackParticle *tp = ele.trackParticle();
    if (tp==nullptr) return;
    
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) {
      fatal("Cannot access EventInfo");
    }

    double d0sig = xAOD::TrackingHelpers::d0significance(tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if ( fabs(d0sig) > m_d0BySigd0Cut) return;
    
    const xAOD::VertexContainer* vertexCont =0;
    if (m_event->retrieve(vertexCont,"PrimaryVertices").isFailure()) return;
    
    const xAOD::Vertex* pvx = xAOD::PVHelpers::getHardestVertex(vertexCont);
    if (pvx == nullptr) return;
    
    double z0 = tp->z0() + tp->vz() - pvx->z();
    z0 = z0 * sin(tp->theta());
    if ( fabs(z0) > m_z0Cut ) return;

    passIPCut(ele) = true;
  }
  
  //______________________________________________________________________________
  HG::Iso::IsolationType ElectronHandler::getIsoType(TString isoName) {
    if      (isoName == "LooseTrackOnly") return HG::Iso::LooseTrackOnly;
    else if (isoName == "Loose") return HG::Iso::Loose;
    else if (isoName == "Gradient") return HG::Iso::Gradient;
    else if (isoName == "GradientLoose") return HG::Iso::GradientLoose;
    else if (isoName == "FixedCutTight") return HG::Iso::FixedCutTight;
    else if (isoName == "FixedCutTightTrackOnly") return HG::Iso::FixedCutTightTrackOnly;
    else if (isoName == "FixedCutLoose") return HG::Iso::FixedCutLoose;
    else if (isoName == "UserDefined") return HG::Iso::UserDefined;
    else fatal("Isolation "+isoName+" read from: "+
               m_name+".Selection.IsolationCriteria is not Tight, Gradient, Loose, or UserDefined. Exiting.");
    return HG::Iso::Undefined;
  }
    
}
