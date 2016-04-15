#include "HGamAnalysisFramework/MuonHandler.h"

#ifdef __DC14__
#include "ElectronIsolationSelection/IsolationSelectionTool.h"
#else
#include "IsolationSelection/IsolationSelectionTool.h"
#endif

namespace HG {
  //______________________________________________________________________________
  SG::AuxElement::Accessor<float> MuonHandler::effSF("SF_eff");
  SG::AuxElement::Accessor<float> MuonHandler::effSFIso("SF_eff_iso");
  SG::AuxElement::Accessor<float> MuonHandler::effSFTTVA("SF_eff_TTVA");
  SG::AuxElement::Accessor<float> MuonHandler::scaleFactor("scaleFactor");
  SG::AuxElement::Accessor<char> MuonHandler::isAccepted("isAccepted");
  SG::AuxElement::Accessor<char> MuonHandler::passIPCut("passIPCut");

  //______________________________________________________________________________
  MuonHandler::MuonHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store) : HgammaHandler(name, event, store)
  { }

  //______________________________________________________________________________
  MuonHandler::~MuonHandler()
  {
      SafeDelete(m_muonEffScaleFactors);
      SafeDelete(m_muonEffScaleFactorsIso);
      SafeDelete(m_muonEffScaleFactorsTTVA);
      SafeDelete(m_muonCalibTool);
      SafeDelete(m_muonSelectTool);

      for (auto iso: m_isoTools) SafeDelete(iso.second);
      m_isoTools.clear();

      for (auto dec: m_isoDecorators) SafeDelete(dec.second);
      m_isoDecorators.clear();
  }

  //______________________________________________________________________________
  EL::StatusCode MuonHandler::initialize(Config &config)
  {
    HgammaHandler::initialize(config);

    // Selecting muons
    m_muonSelectTool = new CP::MuonSelectionTool("MuonSelectionTool");
    m_MaxEta = config.getNum(m_name+".Selection.MaxEta", 2.7);
    CP_CHECK(m_name, m_muonSelectTool->setProperty("MaxEta", m_MaxEta));

    m_pidCut = config.getStr(m_name+".Selection.PID", "Medium");
    if (m_pidCut == "Tight")
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::Tight)));
    }
    else if (m_pidCut == "Medium")
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::Medium)));
    }
    else if (m_pidCut == "Loose") 
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::Loose)));
    }
    else if (m_pidCut == "VeryLoose") 
    {
      CP_CHECK(m_name, m_muonSelectTool->setProperty("MuQuality", int(xAOD::Muon::VeryLoose)));
    }
    else
      fatal(TString::Format("Value: %s for key: %s.Selection.PID not vlid MuQuality. Exiting.", m_name.Data(), m_pidCut.Data()).Data());

    if (!m_muonSelectTool->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonSelection Tool. Exiting.");
    }
   

    // Efficiency scale factors
    m_muonEffScaleFactors = new CP::MuonEfficiencyScaleFactors("MuonEfficiencyScaleFactors");
    for (TString prop : {"WorkingPoint", "CalibrationRelease"}) 
    {
      CP_CHECK(m_name, m_muonEffScaleFactors->setProperty(prop.Data(), config.getStr(m_name+".Efficiency."+prop).Data()));
    }

    if (! m_muonEffScaleFactors->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonEfficiencyScaleFactors Tool. Exiting." );
    }


    // Efficiency scale factors - second instance for isolation SF
    m_muonEffScaleFactorsIso = new CP::MuonEfficiencyScaleFactors("MuonEfficiencyScaleFactorsIso");
    CP_CHECK(m_name, m_muonEffScaleFactorsIso->setProperty("WorkingPoint",       config.getStr(m_name+".Efficiency.WorkingPointIso").Data()));
    CP_CHECK(m_name, m_muonEffScaleFactorsIso->setProperty("CalibrationRelease", config.getStr(m_name+".Efficiency.CalibrationRelease").Data()));
    if (! m_muonEffScaleFactorsIso->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonEfficiencyScaleFactors Tool, Isolation instance. Exiting." );
    }


    // Efficiency scale factors - third instance for track-to-vertex-association.  Isn't their a better way of setting this up???
    m_muonEffScaleFactorsTTVA = new CP::MuonEfficiencyScaleFactors("MuonEfficiencyScaleFactorsTTVA");
    CP_CHECK(m_name, m_muonEffScaleFactorsTTVA->setProperty("WorkingPoint",       config.getStr(m_name+".Efficiency.WorkingPointTTVA").Data()));
    CP_CHECK(m_name, m_muonEffScaleFactorsTTVA->setProperty("CalibrationRelease", config.getStr(m_name+".Efficiency.CalibrationRelease").Data()));

    if (! m_muonEffScaleFactorsTTVA->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonEfficiencyScaleFactors Tool, Isolation instance. Exiting." );
    }



    // Calibrate and smear
    m_muonCalibTool = new CP::MuonCalibrationAndSmearingTool( "MuonCalibrationAndSmearingTool" );
    for (TString prop : {"Year", "Algo", "SmearingType"}) 
    {
      CP_CHECK(m_name, m_muonCalibTool->setProperty(prop.Data(), config.getStr(m_name+".Calibration."+prop).Data()));
    }
    
    if (! m_muonCalibTool->initialize().isSuccess() )
    {
      fatal("Failed to properly initialize the MuonCalibrationAndSmearingTool Tool. Exiting." );
    }

          
    // Isolation tools
    m_doIsoCut = config.getBool(m_name+".Selection.ApplyIsoCut", true);
    m_isoCuts  = config.getStrV(m_name+".Selection.IsoCriteria");
    if (m_doIsoCut && m_isoCuts.size() < 1)
      fatal("Isolation cut requested, but no working point supplied. Exiting!");
    for (size_t i = 0; i < m_isoCuts.size(); ++i) {
      HG::Iso::IsolationType iso = getIsoType(m_isoCuts[i]);
      m_isoDecorators[iso] = new SG::AuxElement::Accessor<char>(("isIso"+m_isoCuts[i]).Data());
      
      // first isolation in the list is the default one to apply
      if (i == 0) m_defaultIso = iso;
      
      m_isoTools[iso] = new CP::IsolationSelectionTool(m_isoCuts[i].Data());
#ifdef __DC14__
      CP_CHECK(m_name, m_isoTools[iso]->setProperty("WorkingPoint", m_isoCuts[i].Data()));
#else
      CP_CHECK(m_name, m_isoTools[iso]->setProperty("MuonWP", m_isoCuts[i].Data()));
#endif
      
      if (m_isoTools[iso]->initialize().isFailure()) {
        fatal("Failed to initialize IsolationSelectionTool with WP: "+m_isoCuts[i]);
      }
    }
    
    
    // Read in remaining configuration information
    m_containerName = config.getStr(m_name+".ContainerName", "Muons");
    m_ApplyPtCut    = config.getBool(m_name+".Selection.ApplyPtCut", true);
    m_PtCut         = config.getNum(m_name+".Selection.PtCutGeV",   6.0);

    m_ApplyIPCuts   = config.getBool(m_name+".Selection.ApplyIPCuts", false);
    m_d0BySigd0Cut  = config.getNum(m_name+".Selection.d0BySigd0Max", 3.0);
    m_z0Cut         = config.getNum(m_name+".Selection.z0Max",  0.5);

    
    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  xAOD::MuonContainer MuonHandler::getCorrectedContainer()
  {
    // get the event info
    // const xAOD::EventInfo *eventInfo = 0;
    //if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) fatal("Cannot access EventInfo");

    bool calib = false;
    xAOD::MuonContainer shallowContainer = getShallowContainer(calib);
    if (calib) return shallowContainer;

    for (auto muon: shallowContainer)
    {
      if (fabs(muon->eta()) <= m_MaxEta) calibrateAndSmearMuon(muon, m_muonCalibTool);
      applyScaleFactor(muon);
      
      // Add selection decorations
      isAccepted(*muon) = m_muonSelectTool->accept(muon);
      decorateIPCut(*muon);
      decorateIso(*muon);
    }
    

    // sort the muons
    shallowContainer.sort(comparePt);
    
    return shallowContainer;
  }

  //______________________________________________________________________________
  xAOD::MuonContainer MuonHandler::applySelection(xAOD::MuonContainer &container)
  {
    xAOD::MuonContainer selected(SG::VIEW_ELEMENTS);
    for (auto muon: container)
    {
      // Apply selection cuts
      if (!passSelection(muon)) continue;
      
      // require Isolation
      if (m_doIsoCut && !passIsoCut(muon)) continue;
      
      selected.push_back(muon);
    }
    return selected;
  }

  //______________________________________________________________________________
  CP::SystematicCode MuonHandler::applySystematicVariation(const CP::SystematicSet &sys)
  {

    setVertexCorrected(false);

    bool isAffected = false;
    for (auto var: sys) {
      if (m_muonEffScaleFactors->isAffectedBySystematic(var) ||
          m_muonEffScaleFactorsIso->isAffectedBySystematic(var)      ||
          m_muonEffScaleFactorsTTVA->isAffectedBySystematic(var)    ||
          m_muonCalibTool->isAffectedBySystematic(var)   ) {
        isAffected = true;
        break;
      }
    }


    if (isAffected) {
      CP_CHECK(m_name, m_muonEffScaleFactors->applySystematicVariation(sys));
      CP_CHECK(m_name, m_muonEffScaleFactorsIso->applySystematicVariation(sys));
      CP_CHECK(m_name, m_muonEffScaleFactorsTTVA->applySystematicVariation(sys));
      CP_CHECK(m_name, m_muonCalibTool->applySystematicVariation(sys));
      m_sysName = sys.name() == "" ? "" : "_"+sys.name();
    } else {
      CP_CHECK(m_name, m_muonEffScaleFactors->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_muonEffScaleFactorsIso->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_muonEffScaleFactorsTTVA->applySystematicVariation(CP::SystematicSet()));
      CP_CHECK(m_name, m_muonCalibTool->applySystematicVariation(CP::SystematicSet()));
      m_sysName = "";
   }

   return CP::SystematicCode::Ok;

  }

  
  //______________________________________________________________________________
  void MuonHandler::applyScaleFactor(xAOD::Muon *muon)
  {
    float EfficiencyScaleFactor = 1.0;
    CP::CorrectionCode cc = m_muonEffScaleFactors->applyEfficiencyScaleFactor( *muon );
    if (cc==CP::CorrectionCode::Error)
      Error("applyScaleFactor()","Error applying efficiency scale factor to current muon");

    if(m_muonEffScaleFactors->getEfficiencyScaleFactor(*muon, EfficiencyScaleFactor) == CP::CorrectionCode::Error)
      fatal("applyScaleFactor():  Error retrieving applied scale factor from muon");

    effSF(*muon) = EfficiencyScaleFactor;
    scaleFactor(*muon) = EfficiencyScaleFactor;

    // Apply Iso SF
    float EfficiencyScaleFactorIso = 1.0;
    cc = m_muonEffScaleFactorsIso->applyEfficiencyScaleFactor( *muon );
    if (cc==CP::CorrectionCode::Error)
      Error("applyScaleFactor()","Error applying iso efficiency scale factor to current muon");
    
    if(m_muonEffScaleFactorsIso->getEfficiencyScaleFactor(*muon, EfficiencyScaleFactorIso) == CP::CorrectionCode::Error)
      fatal("applyScaleFactor():  Error retrieving applied iso scale factor from muon");
    effSFIso(*muon)    = EfficiencyScaleFactorIso;
    scaleFactor(*muon) *= EfficiencyScaleFactorIso;


    // Apply TTVA SF
    float EfficiencyScaleFactorTTVA = 1.0;
    cc = m_muonEffScaleFactorsTTVA->applyEfficiencyScaleFactor( *muon );
    if (cc==CP::CorrectionCode::Error)
      Error("applyScaleFactor()","Error applying TTVA efficiency scale factor to current muon");
    
    if(m_muonEffScaleFactorsTTVA->getEfficiencyScaleFactor(*muon, EfficiencyScaleFactorTTVA) == CP::CorrectionCode::Error)
      fatal("applyScaleFactor():  Error retrieving applied TTVA scale factor from muon");
    effSFTTVA(*muon)    = EfficiencyScaleFactorTTVA;
    scaleFactor(*muon) *= EfficiencyScaleFactorTTVA;
  }
  

  
  //______________________________________________________________________________
  void MuonHandler::calibrateAndSmearMuon(xAOD::Muon *muon, CP::MuonCalibrationAndSmearingTool *muonCalibTool)
  {
    // applies calibration to a muon

    // Apply smearing (?)     <<<<<< Does not apply to muons??
    //photonCalibTool->setRandomSeed(evtInfo->eventNumber()*100+gam->index());
    
    // Calibrate the muon
    double E_before = muon->e();
    CP::CorrectionCode cc = muonCalibTool->applyCorrection( *muon );
    if (cc==CP::CorrectionCode::Error)
      Error("calibratedAndSmearMuon()","Error calibrating current muon");
    if (cc==CP::CorrectionCode::OutOfValidityRange)
      Warning("calibratedAndSmearMuon()","Current muon has no valid calibration due to out-of-range");
    
    // decorate the muon with the calibration factor
    muon->auxdata< float >( "Ecalib_ratio" ) = muon->e()/E_before;
  }
    
  
  //______________________________________________________________________________
  //bool MuonHandler::passPtEtaCuts(const xAOD::Muon *muon) {
    /// applies kinematic preselection cut
    // eta cuts
  //  if (fabs(muon->eta()) > m_etaCut) return false;
  
    // pt cuts
  //  if (muon->pt() < m_ptCut) return false;
  //  return true;
  //}
  
  //______________________________________________________________________________
  bool MuonHandler::passSelection(const xAOD::Muon *muon)
  {
    // Muon selector performs the following cuts:
    //if (!m_muonSelectTool->getQuality(muon) <= xAOD::Muon::Medium) return false;
    //if (!m_muonSelectTool->passedIDCuts(muon)) return false;
    //if (!m_muonSelectTool->passedHighPtCuts(muon)) return false;
    if (isAccepted.isAvailable(*muon) && !isAccepted(*muon)) return false;
    	  
    if (m_ApplyPtCut)
    {
      if ( muon->muonType() == xAOD::Muon::CaloTagged)
	return (muon->pt()/HG::GeV >  15.); 
      else 
	return (muon->pt()/HG::GeV > m_PtCut);
    }

    if (m_ApplyIPCuts && passIPCut.isAvailable(*muon) && !passIPCut(*muon)) return false;

    return true;
  }
  
  //______________________________________________________________________________
  void MuonHandler::decorateIPCut(xAOD::Muon &muon)
  {
    passIPCut(muon) = false;

    const xAOD::TrackParticle* tp = muon.primaryTrackParticle();
    if (tp==nullptr) return;
    
    const xAOD::EventInfo *eventInfo = 0;
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) fatal("Cannot access EventInfo");

    // d0significance apply for all muon types
    double d0sig = xAOD::TrackingHelpers::d0significance(tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), eventInfo->beamPosSigmaXY() );
    if ( fabs(d0sig) > m_d0BySigd0Cut) return;

     // d0/z0 apply for all muon types except standalone muons
    if (!(muon.muonType() == xAOD::Muon::MuonStandAlone))
    {     
      const xAOD::VertexContainer* vertexCont =0;
      if (m_event->retrieve(vertexCont,"PrimaryVertices").isFailure()) return;

      const xAOD::Vertex* pvx = xAOD::PVHelpers::getHardestVertex(vertexCont);
      if (pvx == nullptr) return;

      double z0 = tp->z0() + tp->vz() - pvx->z();
      z0 = z0 * sin(tp->theta());

      if ( fabs(z0) > m_z0Cut ) return;
    }
    
    // Still here? It passes!
    passIPCut(muon) = true;
  }


  //______________________________________________________________________________
  bool MuonHandler::passIsoCut(const xAOD::Muon *muon, HG::Iso::IsolationType iso)
  {
    /// applies PID cut specified in config file
    if (iso == HG::Iso::Undefined) 
    {
      if (!m_isoDecorators[m_defaultIso]->isAvailable(*muon)) return true;
      return m_isoTools[m_defaultIso]->accept(*muon);
    }
        
    if (m_isoTools.find(iso) != m_isoTools.end()) 
    {
      if (!m_isoDecorators[iso]->isAvailable(*muon))
        return true;
      return (*m_isoDecorators[iso])(*muon);
    }

    fatal("Muon isolation cut requested that wasn't specified in config file. Exiting.");
    return false;
  }


  //______________________________________________________________________________
  void MuonHandler::decorateIso(xAOD::Muon &muon)
  {
    for (auto dec: m_isoDecorators) {
      if (m_isoTools[dec.first]->accept(muon))
        (*dec.second)(muon) = true;
      else
        (*dec.second)(muon) = false;
    }
  }
  
  
  //______________________________________________________________________________
  void MuonHandler::printMuon(const xAOD::Muon *muon, TString comment)
  {
    // prints details about the photon
    printf("Muon %2zu  %s\n", muon->index(), comment.Data());
    
    // print the 4-vector
    printf("   (pT,eta,phi,m) = (%5.1f GeV,%6.3f,%6.3f,%4.4f GeV)\n", muon->pt()/GeV, muon->eta(), muon->phi(), muon->m()/GeV);
    
    // print some more information
    TString str;
    if (muon->isAvailable<float>("Ecalib_ratio"))
      str+=Form("   calibFactor = %.3f", muon->auxdata<float>("Ecalib_ratio"));
    
    if (muon->isAvailable<float>("EfficiencyScaleFactor"))
      str+=Form("   scalefactor = %5.5f",muon->auxdata<float>("EfficiencyScaleFactor"));
    if (str.Sizeof()) printf("%s\n",str.Data());
  }

  //______________________________________________________________________________
  HG::Iso::IsolationType MuonHandler::getIsoType(TString isoName) {
    if      (isoName == "LooseTrackOnly") return HG::Iso::LooseTrackOnly;
    else if (isoName == "Loose") return HG::Iso::Loose;
    else if (isoName == "Gradient") return HG::Iso::Gradient;
    else if (isoName == "GradientLoose") return HG::Iso::GradientLoose;
    else if (isoName == "FixedCutTightTrackOnly") return HG::Iso::FixedCutTightTrackOnly;
    else if (isoName == "FixedCutLoose") return HG::Iso::FixedCutLoose;
    else if (isoName == "UserDefined") return HG::Iso::UserDefined;
    else fatal("Isolation "+isoName+" read from: "+
               m_name+".Selection.IsolationCriteria is not Tight, Gradient, Loose, or UserDefined. Exiting.");
    return HG::Iso::Undefined;
  }
    
}
