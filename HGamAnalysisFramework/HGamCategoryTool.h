#ifndef HGamAnalysisFramework_HGamCategoryTool_H
#define HGamAnalysisFramework_HGamCategoryTool_H

#include "HGamAnalysisFramework/EventHandler.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include <EventLoop/Worker.h>
#include "HGamAnalysisFramework/HgammaHandler.h"
#include "HGamAnalysisFramework/PhotonHandler.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

namespace TMVA {
  class Reader;
}

namespace HG {

  class HGamCategoryTool {
  protected:
    xAOD::TEvent              *m_event;
    xAOD::TStore              *m_store;

    TMVA::Reader              *m_readerVBF; //!
    TMVA::Reader              *m_readerVH_had; //!

    float                      t_pTt_yy;//!
    float                      t_m_jj; //!
    float                      t_dEta_jj;//!
    float                      t_dPhi_yy_jj;//!
    float                      t_Zepp;//!
    float                      t_Drmin_y_j; //!
    float                      t_Dy_yy_jj; //!
    float                      t_costs; //!
    
    double m_weight; 

  protected:
    double getLeptonSFs( const xAOD::ElectronContainer *electrons, 
                         const xAOD::MuonContainer     *muons );

    bool Passes_VBF(const xAOD::PhotonContainer    *photons  ,
                    const xAOD::JetContainer       *jets     );
 
    bool Passes_VH_MET(const xAOD::PhotonContainer    *photons  ,
                       const xAOD::MissingETContainer *met      );

    bool Passes_VH_hadronic(const xAOD::PhotonContainer    *photons  ,
                            const xAOD::JetContainer       *jets     );

    bool Passes_VH_leptonic(const xAOD::PhotonContainer    *photons  ,
                            const xAOD::ElectronContainer  *electrons,
                            const xAOD::MuonContainer      *muons    ,
                            const xAOD::JetContainer       *jets     ,
                            const xAOD::MissingETContainer *met      );

    bool Passes_VH_dileptons(const xAOD::PhotonContainer    *photons  ,
                             const xAOD::ElectronContainer  *electrons,
                             const xAOD::MuonContainer      *muons    ,
                             const xAOD::JetContainer       *jets     ,
                             const xAOD::MissingETContainer *met      );

    int Passes_ttH(const xAOD::PhotonContainer     *photons  ,
                    const xAOD::ElectronContainer  *electrons,
                    const xAOD::MuonContainer      *muons    ,
                    const xAOD::JetContainer       *jets     ,
                    const xAOD::MissingETContainer *met      );

    inline bool Passes_ttH_Leptonic(const xAOD::PhotonContainer      *photons  ,
                                    const xAOD::ElectronContainer    *electrons,
                                    const xAOD::MuonContainer        *muons    ,
                                    const xAOD::JetContainer         *jets     ,
                                    const xAOD::MissingETContainer   *met      ) { 
      return (Passes_ttH(photons, electrons, muons, jets, met) == 1);
    }
    
    inline bool Passes_ttH_Hadronic(const xAOD::PhotonContainer      *photons  ,
                                    const xAOD::ElectronContainer    *electrons,
                                    const xAOD::MuonContainer        *muons    ,
                                    const xAOD::JetContainer         *jets     ,
                                    const xAOD::MissingETContainer   *met      ) { 
      return (Passes_ttH(photons, electrons, muons, jets, met) == 2);
    }
    float getMVAWeight(TMVA::Reader *XReader);
    void resetReader( );

  public:
    HGamCategoryTool(xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~HGamCategoryTool();

    virtual EL::StatusCode initialize(Config &config);

    // Get category for default coupling analysis
    std::pair<int, float> getCategoryAndWeightMoriond(const xAOD::PhotonContainer    *photons  ,
                                                      const xAOD::ElectronContainer  *electrons,
                                                      const xAOD::MuonContainer      *muons    ,
                                                      const xAOD::JetContainer       *jets     ) ;
    std::pair<int, float> getCategoryAndWeight(const xAOD::PhotonContainer    *photons  ,
                                               const xAOD::ElectronContainer  *electrons,
                                               const xAOD::MuonContainer      *muons    ,
                                               const xAOD::JetContainer       *jets     , 
                                               const xAOD::MissingETContainer *met      );

  };

}
#endif // HGamAnalysisFramework_HGamCategoryTool_H
