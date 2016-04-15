#ifndef HGamAnalysisFramework_TruthHandler_H
#define HGamAnalysisFramework_TruthHandler_H

#include "HGamAnalysisFramework/HgammaIncludes.h"

class MCTruthClassifier;

namespace HG {

  typedef const xAOD::TruthParticleContainer            TruthParticles;
  typedef ConstDataVector<xAOD::TruthParticleContainer> TruthContainer;

  //  Production modes - move to TruthUtils ?
  enum ProdMode {
    GGF = 0, VBF = 1, WH = 2, ZH = 3, TTH = 4, BBH = 5, 
    UNKNOWN = 9
  };

  class TruthHandler {
  private:
    xAOD::TEvent *m_event;
    xAOD::TStore *m_store;

    MCTruthClassifier *m_truthClass;

    std::string m_MxAODName;
    std::string m_particleName;
    std::string m_eventName;
    std::string m_photonName;
    std::string m_electronName;
    std::string m_muonName;
    std::string m_jetName;
    std::string m_metName;
    std::string m_higgsName;

    double m_phMaxEta;
    double m_phMinPt;
    bool   m_phRejectCrack;
    double m_phMinCrack;
    double m_phMaxCrack;
    double m_phIsoCone;
    double m_phIsoCutSlope;
    double m_phIsoCutConst;

    double m_elMaxEta;
    double m_elMinPt;
    bool   m_elRejectCrack;
    double m_elMinCrack;
    double m_elMaxCrack;
    double m_elIsoCone;
    double m_elIsoCut;

    double m_muMaxEta;
    double m_muMinPt;
    double m_muIsoCone;
    double m_muIsoCut;

    double m_jetMaxRapidity;
    double m_jetMinPt;

    std::vector<TString> m_metTypes;

    static SG::AuxElement::Accessor<char>  isIsolated;
    static SG::AuxElement::Accessor<float> etcone20;
    static SG::AuxElement::Accessor<float> etcone40;
    static SG::AuxElement::Accessor<float> ptcone20;
    static SG::AuxElement::Accessor<float> ptcone40;
    static SG::AuxElement::Accessor<float> pt;
    static SG::AuxElement::Accessor<float> eta;

    // produciont mode map. MCchannelNumber -> ProdMode
    std::map<int,ProdMode> m_prodModeMap;

  protected:
    template <class T>
    void
    setOriginalObjectLink(const T *orig, T *copy);

    template <class T>
    void
    setOriginalObjectLink(const DataVector<T> *orig, DataVector<T> *copy);

    template <class T>
    DataVector<T>
    getDeepCopy(const DataVector<T> *parts, std::string name);

    template <class T>
    DataVector<T>
    getShallowCopy(const DataVector<T> *parts, std::string name);

    template <class T>
    DataVector<T>
    getContainer(const TruthContainer *parts, std::string name);

    template <class T>
    bool
    checkEventAndStore(DataVector<T> &cont, std::string name, bool fullCheck = false);

    template <class T>
    void
    setP4(const xAOD::IParticle *p, T &copy);

    template <class T>
    EL::StatusCode
    writeContainer(DataVector<T> &container, std::string name);

    void decorateClassification(xAOD::TruthParticle &part);


  public:
    TruthHandler(xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~TruthHandler();

    virtual EL::StatusCode initialize(Config &config);

    /// Helper functions
    TruthParticles*              getTruthParticles();
    xAOD::TruthParticleContainer getPhotons();
    xAOD::TruthParticleContainer getElectrons();
    xAOD::TruthParticleContainer getMuons();
    xAOD::JetContainer           getJets();
    xAOD::MissingETContainer     getMissingET();

    xAOD::TruthParticleContainer getHiggsBosons();

    xAOD::TruthParticleContainer applyPhotonSelection   (xAOD::TruthParticleContainer &photons  );
    xAOD::TruthParticleContainer applyElectronSelection (xAOD::TruthParticleContainer &electrons);
    xAOD::TruthParticleContainer applyMuonSelection     (xAOD::TruthParticleContainer &muons    );
    xAOD::JetContainer           applyJetSelection      (xAOD::JetContainer           &jets     );
    xAOD::MissingETContainer     applyMissingETSelection(xAOD::MissingETContainer     &mets     );

    void removeOverlap(xAOD::TruthParticleContainer &photons,
                       xAOD::JetContainer           &jets,
                       xAOD::TruthParticleContainer &electrons,
                       xAOD::TruthParticleContainer &muons);

    double truthVertexZ();

    int vbfTopology(std::vector<TLorentzVector> v_jets, TLorentzVector v_higgs);
    
    int truthCategory();

    // Which production mode the current event is
    // (based on MC channel number + configuration)
    ProdMode getProductionMode(int mcID);

    bool isGGF(int mcID) { return getProductionMode(mcID) == GGF; }
    bool isVBF(int mcID) { return getProductionMode(mcID) == VBF; }
    bool isWH(int mcID) { return getProductionMode(mcID) == WH; }
    bool isZH(int mcID) { return getProductionMode(mcID) == ZH; }
    bool isTTH(int mcID) { return getProductionMode(mcID) == TTH; }
    bool isBBH(int mcID) { return getProductionMode(mcID) == BBH; }

    bool   passFiducial(const xAOD::TruthParticleContainer *photons            ,
                        const xAOD::TruthParticleContainer *electrons = nullptr,
                        const xAOD::TruthParticleContainer *muons     = nullptr,
                        const xAOD::JetContainer           *jets      = nullptr);
    bool   passFiducialKinOnly(const xAOD::TruthParticleContainer *photons            ,
                               const xAOD::TruthParticleContainer *electrons = nullptr,
                               const xAOD::TruthParticleContainer *muons     = nullptr,
                               const xAOD::JetContainer           *jets      = nullptr);
    
    EL::StatusCode writePhotons    (xAOD::TruthParticleContainer &parts);
    EL::StatusCode writeElectrons  (xAOD::TruthParticleContainer &parts);
    EL::StatusCode writeMuons      (xAOD::TruthParticleContainer &parts);
    EL::StatusCode writeJets       (xAOD::JetContainer           &parts);
    EL::StatusCode writeMissingET  (xAOD::MissingETContainer     &met  );
    EL::StatusCode writeHiggsBosons(xAOD::TruthParticleContainer &parts);

    static bool comparePt(const xAOD::IParticle *a, const xAOD::IParticle *b);

  };
}

#include "HGamAnalysisFramework/TruthHandler.hpp"

#endif // HGamAnalysisFramework_TruthHandler_H
