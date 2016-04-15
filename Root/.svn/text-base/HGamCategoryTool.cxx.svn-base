#include<cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include "PATInterfaces/SystematicVariation.h"

#include <HGamAnalysisFramework/HgammaAnalysis.h>
#include <HGamAnalysisFramework/HgammaUtils.h>
#include <HGamAnalysisFramework/HGamVariables.h>
#include <HGamAnalysisFramework/HGamCategoryTool.h>

#include "PhotonVertexSelection/PhotonVertexHelpers.h"
#include "PhotonVertexSelection/PhotonPointingTool.h"

#include "HGamAnalysisFramework/VarHandler.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


namespace HG {

  //______________________________________________________________________________
  HGamCategoryTool::HGamCategoryTool(xAOD::TEvent *event, xAOD::TStore *store)
  : m_event(event)
  , m_store(store)
  , m_readerVBF(nullptr)
  , m_readerVH_had(nullptr)
  { }

  //______________________________________________________________________________
  HGamCategoryTool::~HGamCategoryTool()
  {
    SafeDelete(m_readerVBF);
    SafeDelete(m_readerVH_had);
  }

  //______________________________________________________________________________
  EL::StatusCode HGamCategoryTool::initialize(Config &config)
  {
    // Setup VBF MVA reader
    m_readerVBF = new TMVA::Reader( "!Color:!Silent" );
    m_readerVBF->AddVariable("yy_pTt",             &t_pTt_yy);
    m_readerVBF->AddVariable("jj_m",               &t_m_jj);
    m_readerVBF->AddVariable("jj_DeltaEta",        &t_dEta_jj);
    m_readerVBF->AddVariable("yy_jj_DeltaPhi",     &t_dPhi_yy_jj);
    m_readerVBF->AddVariable("abs(eta_Zeppenfeld)",&t_Zepp);
    m_readerVBF->AddVariable("y_j_DeltaR_min",     &t_Drmin_y_j);
    TString readerPathVBF = PathResolverFindCalibFile("HGamAnalysisFramework/MVA_config_VBF.xml");
    m_readerVBF->BookMVA( "BDTG",readerPathVBF);
    // Setup VH_had MVA reader
    m_readerVH_had = new TMVA::Reader( "!Color:!Silent" );
    m_readerVH_had->AddVariable("m_m_jj",               &t_m_jj);
    m_readerVH_had->AddVariable("m_pTt_yy",             &t_pTt_yy);
    m_readerVH_had->AddVariable("m_Dy_yy_jj",           &t_Dy_yy_jj);
    m_readerVH_had->AddVariable("m_costs_yy_jj",        &t_costs);
    TString readerPathVH  = PathResolverFindCalibFile("HGamAnalysisFramework/MVA_config_VH_had.xml");
    m_readerVH_had->BookMVA( "BDTG", readerPathVH);

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  std::pair<int, float> HGamCategoryTool::getCategoryAndWeight(const xAOD::PhotonContainer    *photons  ,
                                             const xAOD::ElectronContainer  *electrons,
                                             const xAOD::MuonContainer      *muons    ,
                                             const xAOD::JetContainer       *jets     ,
                                             const xAOD::MissingETContainer *met      )
  {
    if (not photons || not electrons || not muons || not jets || not met) return std::make_pair(-1, 0.0);
    // Here the defination of different number is :
    // -1 for events without at least 2 photons
    // ggF : 1   // VBF_low : 2  //VBF_high :3　 // VH_had_low 4  // VH_had_high 5 // VH_MET : 6 // VH_LEP : 7  // VH_dilep : 8
    // tth(Had) : 9   // tth(Lep) : 10
    // the selection order is 10987654321//
    // the Pass function for VH_had and VBF is precuts only now.
    //
    if (photons->size()< 2 )                                    return std::make_pair(-1, 0.0);
    if (Passes_ttH_Leptonic(photons,electrons,muons,jets,met))  return std::make_pair( 10, m_weight);
    if (Passes_ttH_Hadronic(photons,electrons,muons,jets,met))  return std::make_pair( 9, m_weight);
    if (Passes_VH_dileptons(photons,electrons,muons,jets,met))  return std::make_pair( 8, m_weight);
    if (Passes_VH_leptonic(photons,electrons,muons,jets,met))   return std::make_pair( 7, m_weight);
    if (Passes_VH_MET(photons,met))                             return std::make_pair( 6, m_weight);
    if (Passes_VH_hadronic(photons,jets)){
       if(getMVAWeight(m_readerVH_had)> 0.56) return std::make_pair( 5, m_weight);
       else if(getMVAWeight(m_readerVH_had)> 0.2 && getMVAWeight(m_readerVH_had)<0.56 ) return std::make_pair( 4, m_weight);
    }

    if (Passes_VBF(photons,jets))  {
      if(getMVAWeight(m_readerVBF)>0.89) return std::make_pair( 3, m_weight);
      else if(getMVAWeight(m_readerVBF)<0.89 &&getMVAWeight(m_readerVBF)>0.5 ) return std::make_pair(2, m_weight);
    }
    return std::make_pair(1, 1.0);
  }

  bool HGamCategoryTool::Passes_VBF(const xAOD::PhotonContainer    *photons,
                  const xAOD::JetContainer       *jets   )

  {
    m_weight = 1.0;
    if (jets->size()<2) return 0 ;
    TLorentzVector gam1,gam2,jet1,jet2;
    gam1=photons->at(0)->p4();//............
    gam2=photons->at(1)->p4();//............
    jet1=jets->at(0)->p4();//............
    jet2=jets->at(1)->p4();//............
    float m_Dy_j_j = fabs(jet1.Eta()-jet2.Eta());
    float m_eta_Zeppenfeld = fabs( (gam1+gam2).Eta() - ((jet1.Eta()+ jet2.Eta())/2.0 ) ); //JV
    if (m_Dy_j_j<2)     return 0;
    if (fabs(m_eta_Zeppenfeld)>5) return 0;
    return 1;
  }

  bool HGamCategoryTool::Passes_VH_hadronic(const xAOD::PhotonContainer    *photons,
                          const xAOD::JetContainer       *jets   )

  {
    m_weight = 1.0;
    if (jets->size()<2) return 0 ;
    TLorentzVector gam1,gam2,jet1,jet2;
    gam1=photons->at(0)->p4();//............
    gam2=photons->at(1)->p4();//............
    jet1=jets->at(0)->p4();//............
    jet2=jets->at(1)->p4();//............
    float m_m_jj = (jet1+jet2).M();
    if (m_m_jj*0.001<50)  return 0;
    if (m_m_jj*0.001>150) return 0;
    return 1;
  }

  bool HGamCategoryTool::Passes_VH_MET(const xAOD::PhotonContainer    *photons,
                     const xAOD::MissingETContainer *met    )
  {
    m_weight = 1.0;
    TLorentzVector gam1,gam2;
    gam1=photons->at(0)->p4();//............
    gam2=photons->at(1)->p4();//............
    float m_pT_yy = (gam1+gam2).Pt();
    float TST_met = (*met)["TST"]->met();//  <-Here?
    float TST_sumet=  (*met)["TST"]->sumet();
    float      m_MET_signi=-99999;
    if (TST_met!=0) m_MET_signi=TST_met*0.001/sqrt(TST_sumet*0.001);
    if (m_pT_yy*0.001<90)     return 0;
    if (m_MET_signi<7)      return 0;
    return 1;
  }

  bool HGamCategoryTool::Passes_VH_leptonic(const xAOD::PhotonContainer    *photons  ,
                          const xAOD::ElectronContainer  *electrons,
                          const xAOD::MuonContainer      *muons    ,
                          const xAOD::JetContainer       *jets     ,
                          const xAOD::MissingETContainer *met      )
  {
    m_weight = 1.0;
    int m_Nelectrons = electrons->size();
    int m_Nmuons = muons->size();
    m_weight *= getLeptonSFs( electrons, muons ); // lepton SFs

    int m_Njets = jets->size();
    int m_Nbjets = 0;
    for (auto jet: *jets) {
      bool MV2c_20_77=jet->auxdata<char>("MV2c20_FixedCutBEff_77");
      m_weight *= jet->auxdata<float>("SF_MV2c20_FixedCutBEff_77"); // b-tagging SFs
      if (MV2c_20_77) m_Nbjets++;
    }

    TLorentzVector gam1,gam2;
    gam1=photons->at(0)->p4();//............
    gam2=photons->at(1)->p4();//............
    float m_pT_yy = (gam1+gam2).Pt();
    float TST_sumet=(*met)["TST"]->sumet();
    float TST_met=(*met)["TST"]->met();
    float m_MET_signi=-99999;
    if (TST_met!=0) m_MET_signi=TST_met*0.001/sqrt(TST_sumet*0.001);

    if (m_pT_yy*0.001<60)       return 0;
    if (m_Nelectrons+m_Nmuons<1)        return 0;
    if ((m_Nelectrons+m_Nmuons)!=1)        return 0;
    if (m_MET_signi<4.5)   return 0;
    if (m_Nbjets>0)  return 0;
    if (m_Njets>=5)      return 0;
    return 1;
  }

  bool HGamCategoryTool::Passes_VH_dileptons(const xAOD::PhotonContainer    *photons  ,
                           const xAOD::ElectronContainer  *electrons,
                           const xAOD::MuonContainer      *muons    ,
                           const xAOD::JetContainer       *jets     ,
                           const xAOD::MissingETContainer *met      )
  {
    m_weight = 1.0;
    int m_Nelectrons = electrons->size();
    int m_Nmuons = muons->size();
    m_weight *= getLeptonSFs( electrons, muons ); // lepton SFs

    int m_Nbjets = 0;
    for (auto jet: *jets) {
      bool MV2c_20_77=jet->auxdata<char>("MV2c20_FixedCutBEff_77");
      m_weight *= jet->auxdata<float>("SF_MV2c20_FixedCutBEff_77"); // b-tagging SFs
      if (MV2c_20_77) m_Nbjets++;
    }

    TLorentzVector gam1,gam2;
    gam1=photons->at(0)->p4();//............
    gam2=photons->at(1)->p4();//............
    float TST_sumet=(*met)["TST"]->sumet();
    float TST_met=(*met)["TST"]->met();
    float m_MET_signi=-99999;
    if (TST_met!=0) m_MET_signi=TST_met*0.001/sqrt(TST_sumet*0.001);
    TLorentzVector el1, el2, mu1,mu2;
    float m_m_mumu = 0, m_m_ee=0;
    if (m_Nelectrons>1){
      el1 = electrons->at(0)->p4();
      el2= electrons->at(1)->p4();
      m_m_ee =( el1+el2).M();
    }
    if (m_Nmuons>1){
      mu1 = muons->at(0)->p4();
      mu2 = muons->at(1)->p4();
      m_m_mumu =( mu1+mu2).M();
    }
    if (!(( m_Nmuons>=2 && m_m_mumu*0.001>70 && m_m_mumu*0.001<110)||( m_Nelectrons>=2 && m_m_ee*0.001>70 && m_m_ee*0.001<110))) return 0;
    if (m_Nbjets>0)    return 0;
    if (m_MET_signi>3)  return 0;
    return 1;
  }

  // *******************************************
  //   ttH Category (jared.vasquez@yale.edu)
  // *******************************************
  int HGamCategoryTool::Passes_ttH(const xAOD::PhotonContainer    *photons  ,
                  const xAOD::ElectronContainer  *electrons,
                  const xAOD::MuonContainer      *muons    ,
                  const xAOD::JetContainer       *jets     ,
                  const xAOD::MissingETContainer *met      )
  {
    m_weight = 1.0;
    int n_electrons = electrons->size();
    int n_muons = muons->size();
    int n_leptons = (n_electrons + n_muons);

    int n_jets(0), n_jets30(0);
    int n_tags77(0), n_tags77_30(0);
    double tagSF77(1.0), tagSF77_30(1.0);

    // count only central jets s.t. |eta| < 2.5
    for ( auto jet: *jets ) {
      if (fabs(jet->eta()) > 2.5) continue;
      n_jets++;

      // count b-tags (leave room for continuous b-tagging)
      bool MV2c20_77 = jet->auxdata<char>("MV2c20_FixedCutBEff_77");
      tagSF77 *= jet->auxdata<float>("SF_MV2c20_FixedCutBEff_77");
      if (MV2c20_77) n_tags77++;

      // Count 30 GeV Jets and Tags
      if (jet->pt() > 30.e3) {
        n_jets30++;
        tagSF77_30 *= jet->auxdata<float>("SF_MV2c20_FixedCutBEff_77");
        if (MV2c20_77) n_tags77_30++;
      }

    }

    // may use met_sig in near future, study needed.
    float TST_met=(*met)["TST"]->met() * HG::invGeV;

    // Leptonic Channel Selection
    if (n_leptons >= 1) {
      bool passJets = (n_jets >= 2);
      bool passTags = (n_tags77 >= 1);
      bool passMET  = ((TST_met >= 20) || (n_tags77 > 1));

      bool passZey(true);
      for ( auto el: *electrons ) {
        double mey1 = (el->p4() + photons->at(0)->p4()).M() * HG::invGeV;
        double mey2 = (el->p4() + photons->at(1)->p4()).M() * HG::invGeV;
        if ( (fabs(mey1-89) < 5) || (fabs(mey2-89) < 5) ) passZey = false;
      }

      m_weight *= getLeptonSFs( electrons, muons ); // lepton SFs
      m_weight *= tagSF77; // b-tag SFs

      return int( passJets && passTags && passZey && passMET ); // Return 1 for Leptonic Category, 0 if fails

    // Hadronic Channel Selection
    } else {
      bool passJets = (n_jets30 >= 5);
      bool passTags = (n_tags77_30 >= 1);
      m_weight *= tagSF77_30; // b-tag SFs

      return int( passJets && passTags )*2; // Return 2 for Hadronic Category, 0 if fails
    }

    // should never reach this point
    return 0;
  }


  double HGamCategoryTool::getLeptonSFs( const xAOD::ElectronContainer *electrons, const xAOD::MuonContainer *muons )
  {
    double lepSF = 1.0;
    for( auto el : *electrons ) lepSF *= el->auxdata< float >("scaleFactor");
    for( auto mu : *muons )     lepSF *= mu->auxdata< float >("scaleFactor");
    return lepSF;
  }


  void HGamCategoryTool::resetReader()
  {
    t_pTt_yy     = var::pTt_yy()     ;
    t_costs      = var::cosTS_yyjj();
    t_m_jj       = var::m_jj()     ;
    t_dEta_jj    = var::Dy_j_j();
    t_dPhi_yy_jj = var::Dphi_yy_jj();
    t_Zepp       = fabs(var::Zepp());
    t_Drmin_y_j  = var::DRmin_y_j();
    t_Dy_yy_jj   = var::Dy_yy_jj();
  }

//get MVA weight
//
  float HGamCategoryTool::getMVAWeight( TMVA::Reader *xReader)
  {
      resetReader();
      float BDTG_weight=  xReader->EvaluateMVA( "BDTG" ) ;
      return BDTG_weight;
  }

  std::pair<int, float> HGamCategoryTool::getCategoryAndWeightMoriond(const xAOD::PhotonContainer    *photons  ,
                                                    const xAOD::ElectronContainer  *electrons,
                                                    const xAOD::MuonContainer      *muons    ,
                                                    const xAOD::JetContainer       *jets     )
  {
    if (not photons || not electrons || not muons || not jets) return std::make_pair(-1, 0.0);
    if(photons->size()< 2 )return std::make_pair(-1, 0.0);
    if( Passes_VBF( photons  , jets ) )return std::make_pair(2, 1.0);
    return std::make_pair(1, 1.0);

  }
}
