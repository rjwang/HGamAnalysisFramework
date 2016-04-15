#include "HGamAnalysisFramework/TruthHandler.h"

#include "HGamAnalysisFramework/PhotonHandler.h"
#include "HGamAnalysisFramework/ElectronHandler.h"
#include "HGamAnalysisFramework/MuonHandler.h"
#include "HGamAnalysisFramework/JetHandler.h"
#include "HGamAnalysisFramework/ETmissHandler.h"
#include "HGamAnalysisFramework/HGamVariables.h"

#include "HGamAnalysisFramework/HgammaIncludes.h"

#include "MCTruthClassifier/MCTruthClassifier.h"

typedef std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> ClassifierResult;

namespace HG {

  //______________________________________________________________________________
  SG::AuxElement::Accessor<char> TruthHandler::isIsolated("isIsolated");
  SG::AuxElement::Accessor<float> TruthHandler::etcone20("etcone20");
  SG::AuxElement::Accessor<float> TruthHandler::etcone40("etcone40");
  SG::AuxElement::Accessor<float> TruthHandler::ptcone20("ptcone20");
  SG::AuxElement::Accessor<float> TruthHandler::ptcone40("ptcone40");
  SG::AuxElement::Accessor<float> TruthHandler::pt("pt");
  SG::AuxElement::Accessor<float> TruthHandler::eta("eta");

  //______________________________________________________________________________
  TruthHandler::TruthHandler(xAOD::TEvent *event, xAOD::TStore *store)
  : m_event(event)
  , m_store(store)
  , m_truthClass(nullptr)
  { }

  //______________________________________________________________________________
  TruthHandler::~TruthHandler()
  {
    SafeDelete(m_truthClass);
  }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::initialize(Config &config)
  {
    m_MxAODName     = "HGam";

    //Production mode
    TString key("mcChannelNumbers");
    for (auto id:config.getNumV("ggF."+key)) m_prodModeMap[int(id)]=GGF;
    for (auto id:config.getNumV("VBF."+key)) m_prodModeMap[int(id)]=VBF;
    for (auto id:config.getNumV("WH."+key))  m_prodModeMap[int(id)]=WH;
    for (auto id:config.getNumV("ZH."+key))  m_prodModeMap[int(id)]=ZH;
    for (auto id:config.getNumV("ttH."+key)) m_prodModeMap[int(id)]=TTH;
    for (auto id:config.getNumV("bbH."+key)) m_prodModeMap[int(id)]=BBH;
    
    // Container names
    m_particleName = config.getStr("TruthHandler.ParticleContainerName"   , "TruthParticles"  );
    m_eventName    = config.getStr("TruthHandler.EventContainerName"      , "TruthEvents"     );
    m_photonName   = config.getStr("TruthHandler.PhotonContainerName"     , "TruthPhotons"    );
    m_electronName = config.getStr("TruthHandler.ElectronContainerName"   , "TruthElectrons"  );
    m_muonName     = config.getStr("TruthHandler.MuonContainerName"       , "TruthMuons"      );
    m_jetName      = config.getStr("TruthHandler.JetContainerName"        , "AntiKt4TruthJets");
    m_metName      = config.getStr("TruthHandler.MissingETContainerName"  , "MET_Truth"       );
    m_higgsName    = config.getStr("TruthHandler.HiggsBosonContainerName" , "TruthHiggsBosons");

    // Photon selections
    m_phMaxEta      = config.getNum ("TruthHandler.Photons.MaxAbsEta"            , 2.37);
    m_phMinPt       = config.getNum ("TruthHandler.Photons.PtPreCutGeV"          , 25.0)*HG::GeV;
    m_phRejectCrack = config.getBool("TruthHandler.Photons.ApplyCrackRejection"  , true);
    m_phMinCrack    = config.getNum ("TruthHandler.Photons.BarrelMaxAbsEta"      , 1.37);
    m_phMaxCrack    = config.getNum ("TruthHandler.Photons.EndcapMinAbsEta"      , 1.52);
    m_phIsoCone     = config.getNum ("TruthHandler.Photons.IsolationCone"        , 0.40);
    m_phIsoCutSlope = config.getNum ("TruthHandler.Photons.IsolationCutSlope"    , 0.10);
    m_phIsoCutConst = config.getNum ("TruthHandler.Photons.IsolationCutConstGeV" , 1.00)*HG::GeV;

    // Electron selections
    m_elMaxEta      = config.getNum ("TruthHandler.Electrons.MaxAbsEta"          , 2.47);
    m_elMinPt       = config.getNum ("TruthHandler.Electrons.PtPreCutGeV"        , 25.0)*HG::GeV;
    m_elRejectCrack = config.getBool("TruthHandler.Electrons.ApplyCrackRejection", true);
    m_elMinCrack    = config.getNum ("TruthHandler.Electrons.BarrelMaxAbsEta"    , 1.37);
    m_elMaxCrack    = config.getNum ("TruthHandler.Electrons.EndcapMinAbsEta"    , 1.52);
    m_elIsoCone     = config.getNum ("TruthHandler.Electrons.IsolationCone"      , 0.40);
    m_elIsoCut      = config.getNum ("TruthHandler.Electrons.IsolationCutGeV"    , 14.0)*HG::GeV;

    // Muon selections
    m_muMaxEta  = config.getNum("TruthHandler.Muons.MaxAbsEta"      , 2.50);
    m_muMinPt   = config.getNum("TruthHandler.Muons.PtPreCutGeV"    , 10.0)*HG::GeV;
    m_muIsoCone = config.getNum("TruthHandler.Muons.IsolationCone"  , 0.40);
    m_muIsoCut  = config.getNum("TruthHandler.Muons.IsolationCutGeV", 14.0)*HG::GeV;

    // Jet selections
    m_jetMaxRapidity = config.getNum("TruthHandler.Jets.MaxAbsRapidity", 4.40);
    m_jetMinPt       = config.getNum("TruthHandler.Jets.PtPreCutGeV"   , 25.0)*HG::GeV;

    // MET selections
    m_metTypes       = config.getStrV("TruthHandler.MissingET.METTypes");

    // MC truth classifier
    m_truthClass = new MCTruthClassifier("MCTruthClassifier");
    if (m_truthClass->initialize().isFailure())
      HG::fatal("Couldn't initialize MCTruthClassifier");

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setP4<xAOD::Muon>(const xAOD::IParticle *p, xAOD::Muon &copy)
  { copy.setP4(p->pt(), p->eta(), p->phi()); }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setP4<xAOD::TruthParticle>(const xAOD::IParticle *p, xAOD::TruthParticle &copy)
  {
    const xAOD::TruthParticle *truth = static_cast<const xAOD::TruthParticle*>(p);
    if (truth == nullptr)
      fatal("Couldn't cast passed IParticle to TruthParticle, which is not supported. Exiting.");

    copy.setPx(truth->px());
    copy.setPy(truth->py());
    copy.setPz(truth->pz());
    copy.setE (truth->e ());
    copy.setM (truth->m ());
  }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setOriginalObjectLink<xAOD::MissingET>(const xAOD::MissingET *orig, xAOD::MissingET *copy)
  {
    // ASG function only works for IParticle, do nothing for MissingET
  }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setOriginalObjectLink<xAOD::MissingET>(const DataVector<xAOD::MissingET> *orig, DataVector<xAOD::MissingET> *copy)
  {
    // ASG function only works for DataVector<IParticle>, do nothing for MissingET
  }

  //______________________________________________________________________________
  void TruthHandler::decorateClassification(xAOD::TruthParticle &part)
  {
    // Try to retrieve pointer to original particle, if possible
    const xAOD::TruthParticle *_part = &part;

    static SG::AuxElement::ConstAccessor<ElementLink<xAOD::IParticleContainer> > link("originalObjectLink");
    if (link.isAvailable(part))
      _part = dynamic_cast<const xAOD::TruthParticle*>(xAOD::getOriginalObject(part));

    // Truth type/origin decorations
    static SG::AuxElement::Accessor<int> truthType("truthType");
    static SG::AuxElement::Accessor<int> truthOrigin("truthOrigin");

    // Use MCTruthClassifier to retrieve type/origin info
    ClassifierResult result = m_truthClass->particleTruthClassifier(_part);
    truthType(part) = result.first;
    truthOrigin(part) = result.second;
  }


  //______________________________________________________________________________
  TruthParticles* TruthHandler::getTruthParticles()
  {
    if (!m_event->contains<xAOD::TruthParticleContainer>(m_particleName))
      return nullptr;

    const xAOD::TruthParticleContainer *truthParticles = nullptr;
    if (m_event->retrieve(truthParticles, m_particleName).isFailure())
      HG::fatal("Can't access TruthParticleContainer");

    return truthParticles;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getPhotons()
  {
    // Check if HGam truth photons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_photonName, true))
      return cont;

    // Get all truth particles, then good photons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_photonName+", exiting!");

    TruthContainer  truthcont      = HG::getGoodTruthPhotons(truthParticles);

    // Make a shallow copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_photonName);
    cont.sort(comparePt);

    // Add decorations
    static std::vector<int> ignorePdgIds = {13, 12, 14, 16, 18}; // mu, nus
    for (auto part: cont) {
      etcone20(*part) = HG::getTruthIsolation(part, truthParticles, 0.2, false, ignorePdgIds);
      etcone40(*part) = HG::getTruthIsolation(part, truthParticles, 0.4, false, ignorePdgIds);
      ptcone20(*part) = HG::getTruthIsolation(part, truthParticles, 0.2, true );
      ptcone40(*part) = HG::getTruthIsolation(part, truthParticles, 0.4, true );

      isIsolated(*part) = true;
      if (m_phIsoCone > 0 &&
          HG::getTruthIsolation(part, truthParticles, m_phIsoCone, false, ignorePdgIds) >= m_phIsoCutSlope*part->pt() + m_phIsoCutConst)
        isIsolated(*part) = false;

      pt (*part) = part->pt();
      eta(*part) = part->eta();

      decorateClassification(*part);
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getElectrons()
  {
    // Check if HGam truth electrons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_electronName, true))
      return cont;

    // Get all truth particles, then good electrons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_electronName+", exiting!");

    TruthContainer  truthcont      = HG::getGoodTruthElectrons(truthParticles);

    // Make a shallow copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_electronName);
    cont.sort(comparePt);

    // Add decorations
    for (auto part: cont) {
      isIsolated(*part) = true;
      if (m_elIsoCone > 0 &&
          HG::getTruthIsolation(part, truthParticles, m_elIsoCone) > m_elIsoCut)
        isIsolated(*part) = false;

      pt (*part) = part->pt();
      eta(*part) = part->eta();
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getMuons()
  {
    // Check if HGam truth muons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_muonName, true))
      return cont;

    // Get all truth particles, then good muons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_muonName+", exiting!");

    TruthContainer  truthcont      = HG::getGoodTruthMuons(truthParticles);

    // Make a shallow copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_muonName);
    cont.sort(comparePt);

    // Add decorations
    for (auto part: cont) {
      isIsolated(*part) = true;
      if (m_muIsoCone > 0 &&
          HG::getTruthIsolation(part, truthParticles, m_muIsoCone) > m_muIsoCut)
        isIsolated(*part) = false;

      pt (*part) = part->pt();
      eta(*part) = part->eta();
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::JetContainer TruthHandler::getJets()
  {
    // Check if HGam truth jets are in TEvent/TStore already
    xAOD::JetContainer jets;
    if (checkEventAndStore<xAOD::Jet>(jets, m_jetName, true))
      return jets;

    // Check if raw truth jets are in TEvent/TStore already
    if (!checkEventAndStore<xAOD::Jet>(jets, m_jetName))
      fatal("Couldn't retrieve raw truth jet container, exiting.");

    jets.sort(comparePt);
    return jets;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer TruthHandler::getMissingET()
  {
    // Check if HGam truth MissingETs are in TEvent/TStore already
    xAOD::MissingETContainer met;
    if (checkEventAndStore<xAOD::MissingET>(met, m_MxAODName + m_metName))
      return met;

    // Check if raw truth MissingETs are in TEvent/TStore already
    if (!checkEventAndStore<xAOD::MissingET>(met, m_metName))
      fatal("Couldn't retrieve raw truth MissingET container, exiting.");

    return met;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getHiggsBosons()
  {
    // Check if HGam Higgs Bosons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_MxAODName + m_higgsName))
      return cont;

    // Get all truth particles, then good muons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_higgsName+", exiting!");

    TruthContainer truthcont = HG::getFinalHiggsBosons(truthParticles);

    // Make a deep copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_higgsName);
    cont.sort(comparePt);

    // Add decorations
    for (auto part: cont) {
      pt (*part) = part->pt();
      eta(*part) = part->eta();
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::applyPhotonSelection(xAOD::TruthParticleContainer &photons)
  {
    xAOD::TruthParticleContainer selected(SG::VIEW_ELEMENTS);
    for (auto ph: photons) {
      // Pt cuts
      if (ph->pt() < m_phMinPt) continue;

      // Eta cuts
      double aeta = fabs(ph->eta());
      if (aeta > m_phMaxEta) continue;
      if (m_phRejectCrack &&
          (aeta > m_phMinCrack && aeta < m_phMaxCrack))
        continue;

      // Isolation cuts
      if (isIsolated.isAvailable(*ph) && !isIsolated(*ph))
        continue;

      // Passed cuts, add to selected container
      selected.push_back(ph);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::applyElectronSelection(xAOD::TruthParticleContainer &electrons)
  {
    xAOD::TruthParticleContainer selected(SG::VIEW_ELEMENTS);
    for (auto el: electrons) {
      // Pt cuts
      if (el->pt() < m_elMinPt) continue;

      // Eta cuts
      double aeta = fabs(el->eta());
      if (aeta > m_elMaxEta) continue;
      if (m_elRejectCrack &&
          (aeta > m_elMinCrack && aeta < m_elMaxCrack))
        continue;

      // Isolation cuts
      if (isIsolated.isAvailable(*el) && !isIsolated(*el))
        continue;

      // Passed cuts, add to selected container
      selected.push_back(el);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::applyMuonSelection(xAOD::TruthParticleContainer &muons)
  {
    xAOD::TruthParticleContainer selected(SG::VIEW_ELEMENTS);
    for (auto mu: muons) {
      // Pt cuts
      if (mu->pt() < m_muMinPt) continue;

      // Eta cuts
      if (fabs(mu->eta()) > m_muMaxEta) continue;

      // Isolation cuts
      if (isIsolated.isAvailable(*mu) && !isIsolated(*mu))
        continue;

      // Passed cuts, add to selected container
      selected.push_back(mu);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::JetContainer TruthHandler::applyJetSelection(xAOD::JetContainer &jets)
  {
    xAOD::JetContainer selected(SG::VIEW_ELEMENTS);
    for (auto jet: jets) {
      // Pt cuts
      if (jet->pt() < m_jetMinPt) continue;

      // Eta cuts
      if (fabs(jet->rapidity()) > m_jetMaxRapidity) continue;

      // Passed cuts, add to selected container
      selected.push_back(jet);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer TruthHandler::applyMissingETSelection(xAOD::MissingETContainer &mets)
  {
    xAOD::MissingETContainer selected(SG::VIEW_ELEMENTS);
    for (auto met: mets) {
      // Limit MET to the types specified in config (TST, CST, ...)
      if (std::find(m_metTypes.begin(), m_metTypes.end(), met->name().c_str()) == m_metTypes.end())
        continue;

      // Passed cuts, add to selected container
      selected.push_back(met);
    }

    return selected;
  }

  //______________________________________________________________________________
  void TruthHandler::removeOverlap(xAOD::TruthParticleContainer &photons  ,
                                   xAOD::JetContainer           &jets     ,
                                   xAOD::TruthParticleContainer &electrons,
                                   xAOD::TruthParticleContainer &muons    )
  {
    // jets overlapping a pT>15 GeV photon or electron are removed
    for ( auto jet=jets.rbegin(); jet!=jets.rend(); ++jet) {
      bool overlap=false;
      for (auto gam:photons)
        if (gam->pt()>15.0*HG::GeV && HG::DRrap(gam,*jet)<0.4)
          overlap=true;
      for (auto e:electrons)
        if (e->pt()>15.0*HG::GeV && HG::DRrap(e,*jet)<0.4)
          overlap=true;
      if (overlap) jets.erase(jet.base()-1);
    }
  }

  //______________________________________________________________________________
  double TruthHandler::truthVertexZ()
  {
    if (var::truthVertexZ.exists())
      return var::truthVertexZ();

    const xAOD::TruthEventContainer *truthEvents = nullptr;
    if (m_event->retrieve(truthEvents, m_eventName).isFailure())
      HG::fatal("Can't access TruthEvents");

    if (truthEvents->size() < 1) {
      Warning("TruthHandler::truthVertexZ()","No TruthEvents?");
      return -999;
    }

    static int nNullVtx = 0;
    if (truthEvents->at(0)->signalProcessVertex() == nullptr) {
      nNullVtx++;
      if (nNullVtx < 5)
        Warning("TruthHandler::truthVertexZ()","No signalProcessVertex for event");
      if (nNullVtx == 5)
        Warning("TruthHandler::truthVertexZ()","Supressing WARNING for: No signalProcessVertex for event");
      return -999;
    }

    var::truthVertexZ.setValue(truthEvents->at(0)->signalProcessVertex()->z());
    return var::truthVertexZ();
  }

  //______________________________________________________________________________
  bool TruthHandler::passFiducial(const xAOD::TruthParticleContainer *allPhotons,
                                  const xAOD::TruthParticleContainer *electrons ,
                                  const xAOD::TruthParticleContainer *muons     ,
                                  const xAOD::JetContainer           *jets      )
  {
    if (var::isFiducial.exists())
      return var::isFiducial();

    var::isFiducial.setTruthValue(false);

    // Kinematic cuts
    if (not passFiducialKinOnly(allPhotons, electrons, muons, jets))
      return false;

    const xAOD::TruthParticle *gam1 = (*allPhotons)[0], *gam2 = (*allPhotons)[1];

    // Isolation cuts
    if (not isIsolated(*gam1) || not isIsolated(*gam2))
      return false;

    var::isFiducial.setTruthValue(true);

    // All cuts passed, this event is in fiducial volume
    return true;
  }
  //______________________________________________________________________________
  bool TruthHandler::passFiducialKinOnly(const xAOD::TruthParticleContainer *allPhotons,
                                         const xAOD::TruthParticleContainer *electrons ,
                                         const xAOD::TruthParticleContainer *muons     ,
                                         const xAOD::JetContainer           *jets      )
  {
    if (var::isFiducialKinOnly.exists())
      return var::isFiducialKinOnly();

    var::isFiducialKinOnly.setTruthValue(false);

    // Safety check
    if (allPhotons == nullptr)
      return false;

    // At least two photons
    if (allPhotons->size() < 2)
      return false;

    const xAOD::TruthParticle *gam1 = (*allPhotons)[0], *gam2 = (*allPhotons)[1];

    // Pt cuts
    if (gam1->pt() < m_phMinPt || gam2->pt() < m_phMinPt)
      return false;

    // Eta cuts
    double aeta1 = fabs(gam1->eta()), aeta2 = fabs(gam2->eta());
    if (aeta1 > m_phMaxEta || aeta2 > m_phMaxEta)
      return false;

    if (m_phRejectCrack                                  &&
        ((aeta1 > m_phMinCrack && aeta1 < m_phMaxCrack)  ||
         (aeta2 > m_phMinCrack && aeta2 < m_phMaxCrack)) )
      return false;

    // Relative pT cuts
    TLorentzVector yy = gam1->p4() + gam2->p4();
    if (gam1->pt()/yy.M() < 0.35 || gam2->pt()/yy.M() < 0.25)
      return false;

    // Mass window cut
    if (yy.M() < 105e3 || 160e3 <= yy.M())
      return false;

    var::isFiducialKinOnly.setTruthValue(true);

    // All cuts passed, this event is in kinematic fiducial volume
    return true;
  }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writePhotons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_photonName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeElectrons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_electronName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeMuons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_muonName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeJets(xAOD::JetContainer &container)
  { return writeContainer<xAOD::Jet>(container, m_MxAODName + m_jetName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeMissingET(xAOD::MissingETContainer &container)
  { return writeContainer<xAOD::MissingET>(container, m_MxAODName + m_metName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeHiggsBosons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_higgsName); }

  //______________________________________________________________________________
  bool TruthHandler::comparePt(const xAOD::IParticle *a, const xAOD::IParticle *b)
  { return a->pt() > b->pt(); }
  
  //______________________________________________________________________________
  ProdMode TruthHandler::getProductionMode(int mcID) {

    // return UNKOWN and print a warning the first time if the mcID isn't known
    if (m_prodModeMap.find(mcID)==m_prodModeMap.end()) {
      static bool first = true;
      if (first) 
	Warning("getProductionMode",
		"Cannot figure out production mode of MC channel number %i, %s",
		mcID,"Did you forget to add it to the config?");
      first=false;
      return UNKNOWN;
    }

    return m_prodModeMap[mcID];
  }
  
  //______________________________________________________________________________
  int TruthHandler::vbfTopology(std::vector<TLorentzVector> v_jets, TLorentzVector v_higgs) {
    //  0 for no VBF topology
    //  >=1 for VBF topology
    //  1 for passing 2-jet cut and failing pT(Hjj) cut
    //  2 for passing 2-jet cut and pT(Hjj) cut
    if (v_jets.size()<2) return 0;
    
    //2-jet VBF cut
    double m_jj = (v_jets.at(0)+v_jets.at(1)).M(); 
    double Dy_jj = std::abs(v_jets.at(0).Rapidity()-v_jets.at(1).Rapidity()); 
    if ( m_jj < 400*GeV || Dy_jj < 4 ) return 0;
    
    //pT(Hjj) cut
    double pT_Hjj = (v_higgs+v_jets.at(0)+v_jets.at(1)).Pt();
    if(pT_Hjj < 30*GeV) return 2;
    return 1;
  }
  
   //______________________________________________________________________________
  int TruthHandler::truthCategory()
  {
    if (var::truthCategory.exists())
      return var::truthCategory();
    
    int truthCategory = -999999;
    
     // get the event info 
    const xAOD::EventInfo *eventInfo = 0; 
    if (m_event->retrieve(eventInfo, "EventInfo").isFailure()) { 
      fatal("Cannot access EventInfo"); 
    }
    
    int mcChannelNumber = eventInfo->mcChannelNumber();

    //production mode flags from a config file 
    bool b_isGGF = isGGF(mcChannelNumber); 
    bool b_isVBF = isVBF(mcChannelNumber); 
    bool b_isWH = isWH(mcChannelNumber); 
    bool b_isZH = isZH(mcChannelNumber); 
    bool b_isTTH = isTTH(mcChannelNumber); 
    bool b_isBBH = isBBH(mcChannelNumber);
    
    bool b_isVH=false;
    if(b_isZH ||b_isWH) b_isVH = true;
    
    //if not Higgs sample, return without setting the variables
    if(!b_isGGF && !b_isVBF && !b_isVH && !b_isTTH) return 0;
    
    bool isProdgg = false;
    bool isProdqq = false;
    
    bool isVqq = false;
    bool isVlnu = false;
    bool isVll = false;
    bool isVnunu = false;
    
    bool isW1qq = false;
    bool isW2qq = false;

    bool found_t1 = false; 
    bool found_W1 = false; 
    int ip_t1 = -999; 
    int ip_W1 = -999;
    
    TLorentzVector v_higgs(0.,0.,0.,0.);
    TLorentzVector v_photon1(0.,0.,0.,0.);
    TLorentzVector v_photon2(0.,0.,0.,0.);
    TLorentzVector v_V(0.,0.,0.,0.);
    TLorentzVector v_electron1(0.,0.,0.,0.);
    TLorentzVector v_electron2(0.,0.,0.,0.);
    TLorentzVector v_V_q1(0.,0.,0.,0.); 
    TLorentzVector v_V_q2(0.,0.,0.,0.);
    TLorentzVector v_b_t1(0.,0.,0.,0.); 
    TLorentzVector v_b_t2(0.,0.,0.,0.);
    TLorentzVector v_W1_q1(0.,0.,0.,0.); 
    TLorentzVector v_W1_q2(0.,0.,0.,0.);
    TLorentzVector v_W2_q1(0.,0.,0.,0.); 
    TLorentzVector v_W2_q2(0.,0.,0.,0.);
    
    //==============Access to TruthEvent
    const xAOD::TruthEventContainer *truthEvents = nullptr; 
    if (m_event->retrieve(truthEvents, m_eventName).isFailure()) 
      HG::fatal("Can't access TruthEvents"); 
 
    if (truthEvents->size() < 1) { 
      Warning("TruthHandler::truthVertexZ()","No TruthEvents?"); 
      return -999; 
    }
    
    //get the scale 
    float qcd_scale = 0.;
    for (const auto* const tePtr : *truthEvents) {
      std::string VarName="Q";
      if (truthEvents->isAvailable<float>(VarName))
	qcd_scale = (*tePtr).auxdataConst< float >( VarName );
    }
    
    //==============Access to TruthParticles
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+", exiting!");
    
    int ip = 0;
    for (xAOD::TruthParticleContainer::const_iterator pItr=truthParticles->begin(); pItr!=truthParticles->end(); ++pItr) { 
      const xAOD::TruthParticle* particle = *pItr; 
       
      ip++;
      //std::cout << ip << " " << particle->pdgId() << " "<< particle->status() << " " << particle->barcode() << " "  << particle->pt() << std::endl;
      
      //===VH
      if(b_isVH){
	//gg->VH not used for now, no need to be careful
	if(ip<5){
	  if(particle->pdgId()==21) isProdgg = true;
	  if(fabs(particle->pdgId())<7) isProdqq = true;
	}
	
	//find V for VH mode (cannot link it easily to the Higgs)
	if((particle->pdgId()==23 || fabs(particle->pdgId())==24) && particle->status()==62){
	  
	  v_V = particle->p4();
	  
	  //search for the daughters
	  if ( particle->hasDecayVtx() ) { 
	    const xAOD::TruthVertex* vertex = particle->decayVtx(); 
	    
	    int pdgid_d1 = 0;
	    int pdgid_d2 = 0;
	    
	    if(vertex->nOutgoingParticles()<2) continue;
	    
	    const xAOD::TruthParticle* daughter1 = vertex->outgoingParticle(0);
	    const xAOD::TruthParticle* daughter2 = vertex->outgoingParticle(1);
	    pdgid_d1 = daughter1->pdgId();
	    pdgid_d2 = daughter2->pdgId();
	    
	    bool isd1_q = false;
	    bool isd1_l = false;
	    bool isd1_nu = false;
	    if(fabs(pdgid_d1)>0 && fabs(pdgid_d1)<7) isd1_q = true;
	    if(fabs(pdgid_d1)==11 || fabs(pdgid_d1)==13 || fabs(pdgid_d1)==15) isd1_l = true;
	    if(fabs(pdgid_d1)==12 || fabs(pdgid_d1)==14 || fabs(pdgid_d1)==16) isd1_nu = true;
	    
	    bool isd2_q = false;
	    bool isd2_l = false;
	    bool isd2_nu = false;
	    if(fabs(pdgid_d2)>0 && fabs(pdgid_d2)<7) isd2_q = true;
	    if(fabs(pdgid_d2)==11 || fabs(pdgid_d2)==13 || fabs(pdgid_d2)==15) isd2_l = true;
	    if(fabs(pdgid_d2)==12 || fabs(pdgid_d2)==14 || fabs(pdgid_d2)==16) isd2_nu = true;
	    
	    //set the V decay flag
	    if(isd1_q && isd2_q){
	      isVqq = true;
	      v_V_q1 = daughter1->p4();
	      v_V_q2 = daughter2->p4();
	    }
	    
	    if((isd1_l && isd2_nu) || (isd2_l && isd1_nu)) isVlnu = true;
	    if(isd1_l && isd2_l) isVll = true;
	    //nunu added in the ll category later
	    if(isd1_nu && isd2_nu) isVnunu = true;
	  }
	}//end V boson
      }//end isVH
      
      //===ttH
      if(b_isTTH){
	
	//t1
	if(fabs(particle->pdgId())==6 && particle->status()==62 &&!found_t1){
	  if ( particle->hasDecayVtx() ) {  
            const xAOD::TruthVertex* t1_vertex = particle->decayVtx();  
	    
            int pdgid_t1_d1 = 0; 
            int pdgid_t1_d2 = 0; 
            
            if(t1_vertex->nOutgoingParticles()<2) continue;
	    found_t1=true;
	    ip_t1=ip;
	    
            const xAOD::TruthParticle* t1_daughter1 = t1_vertex->outgoingParticle(0); 
            const xAOD::TruthParticle* t1_daughter2 = t1_vertex->outgoingParticle(1); 
	    
	    if(fabs(t1_daughter1->pdgId())!=5 && fabs(t1_daughter2->pdgId())!=5) continue;
	    if(fabs(t1_daughter1->pdgId())!=24 && fabs(t1_daughter2->pdgId())!=24) continue;
	    
	    if(fabs(t1_daughter1->pdgId())==5){//assume Wb, switch otherwise
	      t1_daughter1 = t1_vertex->outgoingParticle(1);  
	      t1_daughter2 = t1_vertex->outgoingParticle(0);
	    }
	    
            pdgid_t1_d1 = t1_daughter1->pdgId(); 
            pdgid_t1_d2 = t1_daughter2->pdgId();
	    
	    if(fabs(pdgid_t1_d1)!=24 || fabs(pdgid_t1_d2)!=5){
	      HG::fatal("Not a Wb decay, exiting");
            }

	    //b-jet
	    v_b_t1 = t1_daughter2->p4();
	  }
	}
	
	//look directly for the W
	if(fabs(particle->pdgId())==24 && !found_W1){
	  
	  if ( particle->hasDecayVtx() ) {   
	    const xAOD::TruthVertex* W1_vertex = particle->decayVtx();   
	    
	    int pdgid_W1_d1 = 0;  
	    int pdgid_W1_d2 = 0;  
	      
	    if(W1_vertex->nOutgoingParticles()<2) continue;  
	    found_W1 = true;
	    ip_W1 = ip;
	    
	    const xAOD::TruthParticle* W1_daughter1 = W1_vertex->outgoingParticle(0);  
	    const xAOD::TruthParticle* W1_daughter2 = W1_vertex->outgoingParticle(1);  
	    
	    pdgid_W1_d1 = W1_daughter1->pdgId();
	    pdgid_W1_d2 = W1_daughter2->pdgId();
	    
	    bool isd1_W1_q = false; 
	    bool isd1_W1_lnu = false;
	    if(fabs(pdgid_W1_d1)>0 && fabs(pdgid_W1_d1)<7) isd1_W1_q = true; 
	    if(fabs(pdgid_W1_d1)>10 && fabs(pdgid_W1_d1)<17) isd1_W1_lnu = true; 
	    
	    //qq
	    if(isd1_W1_q){
	      isW1qq = true;
	      v_W1_q1 = W1_daughter1->p4();
	      v_W1_q2 = W1_daughter2->p4();
	    }
	    
	    //lnu
	    // no need to keep the electron
	    
	  }//W1 decay
	}//W1
	
	//t2
	if(fabs(particle->pdgId())==6 && particle->status()==62 && ip!=ip_t1){
	  
	  if ( particle->hasDecayVtx() ) {  
            const xAOD::TruthVertex* t2_vertex = particle->decayVtx();  
	    
            int pdgid_t2_d1 = 0; 
            int pdgid_t2_d2 = 0; 
	    
            if(t2_vertex->nOutgoingParticles()<2) continue; 
	    
            const xAOD::TruthParticle* t2_daughter1 = t2_vertex->outgoingParticle(0); 
            const xAOD::TruthParticle* t2_daughter2 = t2_vertex->outgoingParticle(1); 

	    if(fabs(t2_daughter1->pdgId())!=5 && fabs(t2_daughter2->pdgId())!=5) continue;
	    if(fabs(t2_daughter1->pdgId())!=24 && fabs(t2_daughter2->pdgId())!=24) continue;

	    if(fabs(t2_daughter1->pdgId())==5){//assume Wb, switch otherwise
	      t2_daughter1 = t2_vertex->outgoingParticle(1);  
	      t2_daughter2 = t2_vertex->outgoingParticle(0);
	    }
	    
            pdgid_t2_d1 = t2_daughter1->pdgId(); 
            pdgid_t2_d2 = t2_daughter2->pdgId();
	    
	    if(fabs(pdgid_t2_d1)!=24 || fabs(pdgid_t2_d2)!=5){
	      HG::fatal("Not a Wb decay, exiting");  
            }

	    //b-jet
	    v_b_t2 = t2_daughter2->p4();
	  }
	}
	
	//look directly for the W
	if(fabs(particle->pdgId())==24 && ip!=ip_W1){
	  
	  if ( particle->hasDecayVtx() ) {   
	    const xAOD::TruthVertex* W2_vertex = particle->decayVtx();   
	    
	    int pdgid_W2_d1 = 0;  
	    int pdgid_W2_d2 = 0;  
	      
	    if(W2_vertex->nOutgoingParticles()<2) continue;  
	    	    
	    const xAOD::TruthParticle* W2_daughter1 = W2_vertex->outgoingParticle(0);  
	    const xAOD::TruthParticle* W2_daughter2 = W2_vertex->outgoingParticle(1);  
	    
	    pdgid_W2_d1 = W2_daughter1->pdgId();
	    pdgid_W2_d2 = W2_daughter2->pdgId();
	    
	    bool isd1_W2_q = false; 
	    bool isd1_W2_lnu = false;
	    if(fabs(pdgid_W2_d1)>0 && fabs(pdgid_W2_d1)<7) isd1_W2_q = true; 
	    if(fabs(pdgid_W2_d1)>10 && fabs(pdgid_W2_d1)<17) isd1_W2_lnu = true; 
	    
	    //qq
	    if(isd1_W2_q){
	      isW2qq = true;
	      v_W2_q1 = W2_daughter1->p4(); 
              v_W2_q2 = W2_daughter2->p4();
	    }
	    
	    //lnu
	    //no need to keep the electron
	    
	  }//W2 decay
	}//W2
      }//end ttH
      
       //===find Higgs  
      if(particle->pdgId()==25 && particle->status()==62){
	v_higgs = particle->p4();
      }//end Higgs
      
    }//end particle loop

    //Built in derivation TruthParticles collections (dressed leptons)
    //egammaTruthParticles
    //use this collection to access the H photons and the V electrons
    const xAOD::TruthParticleContainer* egammaTruthParticles; 
    if (m_event->retrieve(egammaTruthParticles, "egammaTruthParticles").isFailure()) {   
      fatal("Cannot access egammaTruthParticles");   
    }  
    
    bool photon_1 = false; 
    bool electron_1 = false; 
    for (xAOD::TruthParticleContainer::const_iterator pItr=egammaTruthParticles->begin(); pItr!=egammaTruthParticles->end(); ++pItr) { 
      const xAOD::TruthParticle* particle = *pItr; 
      int pdgID = particle->pdgId();  
      
      //get mother
      const xAOD::TruthParticle* parent = particle->parent(0);
      bool isFromHiggs = false;
      bool isFromV = false;
      
      if (particle->nParents()){
	if(parent->pdgId()==25) isFromHiggs = true;
	if(parent->pdgId()==23 || fabs(parent->pdgId())==24) isFromV = true;
      }
       
      if(isFromHiggs && pdgID==22){
	//get the Higgs photons from there
	if(!photon_1){
	  v_photon1 = particle->p4();   
	  photon_1=true;
	  continue;
	}
	
	v_photon2 = particle->p4();    
      }
      
      if(isFromV && fabs(pdgID)==11){
	if(!electron_1){
	  v_electron1 = particle->p4();    
	  electron_1=true;
	  continue;
	}
	v_electron2 = particle->p4();     
      }  
      
    }//end egammaTruthParticles

    int n_electrons = 0; 
    if(v_electron1.Pt()>0) n_electrons++; 
    if(v_electron2.Pt()>0) n_electrons++;

    /* 
    //MuonTruthParticles ==> NOT USED FOR NOW 
    const xAOD::TruthParticleContainer* MuonTruthParticles;  
    if (m_event->retrieve(egammaTruthParticles, "MuonTruthParticles").isFailure()) {      
    fatal("Cannot access MuonTruthParticles"); 
     
    for (xAOD::TruthParticleContainer::const_iterator pItr=MuonTruthParticles->begin(); pItr!=MuonTruthParticles->end(); ++pItr) {  
    const xAOD::TruthParticle* particle = *pItr;  
    //int pdgID = particle->pdgId();   
    //std::cout << particle->auxdata<int>("truthOrigin") << std::endl; 
    } 
    */
    
    //TLorentz vectors to ease the jet overlap removal already defined
    //WARNING: removal of jets from taus from V decays not yet implemented
    
    //AntiKt4TruthJets (container already sorted in pt)
    const xAOD::JetContainer* AntiKt4TruthJets; 
    if (m_event->retrieve(AntiKt4TruthJets, "AntiKt4TruthJets").isFailure()) {  
      fatal("Cannot access AntiKt4TruthJets");  
    } 
    
    //select the "extra" jets
    std::vector<TLorentzVector> v_jets;
    int njets = 0;
    for (xAOD::JetContainer::const_iterator pItr=AntiKt4TruthJets->begin(); pItr!=AntiKt4TruthJets->end(); ++pItr) { 
      const xAOD::Jet* truthjet = *pItr; 
      
      //kinematic selection
      if(truthjet->pt()<30000. || fabs(truthjet->eta())>4) continue;

      TLorentzVector v_jet(0., 0., 0., 0.); 
      v_jet = truthjet->p4();

      //remove jets from the higgs photons (all modes)
      float delta_R_photon1 = v_jet.DeltaR(v_photon1);
      float delta_R_photon2 = v_jet.DeltaR(v_photon2);
      if(delta_R_photon1<0.4 || delta_R_photon2<0.4) continue;

      //remove jets from the V electrons (VH and ttH) 
      if(b_isVH || b_isTTH){
	float delta_R_electron1 = 999.;
	float delta_R_electron2 = 999.;
	if(n_electrons>0) delta_R_electron1 = v_jet.DeltaR(v_electron1); 
	if(n_electrons>1) delta_R_electron2 = v_jet.DeltaR(v_electron2);
	if(delta_R_electron1<0.2 || delta_R_electron2<0.2) continue;
      }
      
      //remove jets from VH
      if(b_isVH && isVqq){
      float delta_R_q1 = 999.;
      float delta_R_q2 = 999.;
      delta_R_q1 = v_jet.DeltaR(v_V_q1);
      delta_R_q2 = v_jet.DeltaR(v_V_q2);
      if(delta_R_q1<0.4 || delta_R_q2<0.4) continue;
      }
      
      //remove jets from ttH
      if(b_isTTH){
	//b-jets
	float delta_R_b_t1 = 999.; 
	float delta_R_b_t2 = 999.;  
	delta_R_b_t1 = v_jet.DeltaR(v_b_t1);  
	delta_R_b_t2 = v_jet.DeltaR(v_b_t2);  
	if(delta_R_b_t1<0.4 || delta_R_b_t2<0.4) continue;
	
	//Wqq
	if(isW1qq){
	  float delta_R_W1_q1 = 999.;
	  float delta_R_W1_q2 = 999.; 
	  delta_R_W1_q1 = v_jet.DeltaR(v_W1_q1); 
	  delta_R_W1_q2 = v_jet.DeltaR(v_W1_q2); 
	  if(delta_R_W1_q1<0.4 || delta_R_W1_q2<0.4) continue;
	}
	if(isW2qq){
	  float delta_R_W2_q1 = 999.;
	  float delta_R_W2_q2 = 999.; 
	  delta_R_W2_q1 = v_jet.DeltaR(v_W2_q1); 
	  delta_R_W2_q2 = v_jet.DeltaR(v_W2_q2); 
	  if(delta_R_W2_q1<0.4 || delta_R_W2_q2<0.4) continue;
	}
      }//end ttH
      
      //jet selected
      v_jets.push_back(v_jet);
      
    }//end jets
    
    //==============The categorization happens here
    //general variables
    njets = v_jets.size();
    //jet bins: 0, 1, >=2
    int njetBin = njets > 9 ? 9 : njets; 
    int vbfTopo = vbfTopology(v_jets, v_higgs);
    double pT_H = v_higgs.Pt();
    int pT_HBin = 0;
    int BSMbin = 0;
    int WHZHBin = 0;
    int WHZHDecayBin = 0;
    double pT_V = v_V.Pt();
    int pT_VBin = 0;
    //GGF 0-999
    if(b_isGGF){
      pT_HBin = pT_H < 60*GeV ? 0 : pT_H < 200*GeV ? 1 : 2;
      truthCategory = njetBin*100 + pT_HBin*10 + vbfTopo;
    }
    //VBF 1000-1999
    else if(b_isVBF){
      if(njets>0) BSMbin = v_jets.at(0).Pt() < 100*GeV ? 0 : 1; 
      truthCategory = 1000 + njetBin*100 + BSMbin*10 + vbfTopo;
    }
    //WH 2000-2999 and ZH 3000-3999
    else if(b_isVH){ 
      WHZHBin = b_isWH ? 2 : b_isZH ? 3 : -999999;
      WHZHDecayBin = isVqq ? 0 : isVlnu ? 1 :  (isVll || isVnunu) ? 2 : -999999;
      pT_VBin = pT_V < 120*GeV ? 0 : pT_V < 200*GeV ? 1 : 2;
      truthCategory = 1000*WHZHBin + njetBin*100 + WHZHDecayBin*10 + pT_VBin;
    }
    //ttH 4000-4999
    else if(b_isTTH){ 
      if(njets>0) BSMbin = v_jets.at(0).Pt() < 100*GeV ? 0 : 1;
      truthCategory = 4000 + njetBin*100 + BSMbin*10;
    }
    //bbH 5000-5999
    else if(b_isBBH){ 
      truthCategory = 5000 + njetBin*100;
    }
    
    //N.B: truthCategory=-999999 for unclassified events
    //std::cout << "TRUTH CATEGORY = " << truthCategory << std::endl;
    
    //==============Get subprocess
    int truthProcess = 0;
    //enum: ggF (1), VBF (2), WH_qq (3), WH_lnu (4) , ZH_qq (5), ZH_nunu (6), ZH_ll (7), ttH_qqqq (8), ttH_qqlnu (9), ttH_lnulnu (10)
    if(b_isGGF) truthProcess = 1;
    if(b_isVBF) truthProcess = 2;
    if(b_isWH){
      if(isVqq) truthProcess = 3;
      if(isVlnu) truthProcess = 4;
    }
    if(b_isZH){
      if(isVqq) truthProcess = 5;
      if(isVnunu) truthProcess = 6;
      if(isVll) truthProcess = 7;
    }
    if(b_isTTH){
      if(isW1qq && isW2qq) truthProcess = 8;
      if((isW1qq && !isW2qq) || (!isW1qq && isW2qq)) truthProcess = 9;
      if(!isW1qq && !isW2qq) truthProcess = 10;
    }
    //std::cout << "TRUTH PROCESS = " << truthProcess << std::endl;
    var::truthProcess.setValue(truthProcess);
    //==============End Get subprocess
    
    var::truthCategory.setValue(truthCategory);
    return var::truthCategory();
  }
  
} // namespace HG
