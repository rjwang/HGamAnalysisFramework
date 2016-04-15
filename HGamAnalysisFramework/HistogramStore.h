#ifndef HISTOGRAM_STORE_H
#define HISTOGRAM_STORE_H

/********************
 * HistogramStore: 
 *   class to create, store, fill and retrieve TH1 and TH2 histograms
 *
 * Usage:
 *   HistogramStore HistoStore;
 *   HistoStore.createTH1F("Nphotons",40,-0.5,39.5,";#it{N}_{photon-clusters}");
 *   vector<TH1*> AllHistos = HistoStore.getListOfHistograms();
 *
 */

#include <map>
#include <vector>

#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TProfile.h>

//JB
#include <TTree.h>

/*! \brief Class that handles creating, filling, and retrieving histograms
 *        through use of an internal map
 *  \author Nathan Readioff
 */
class HistogramStore
{
 private:
  std::map<TString,TH1F*>     m_histoTH1F;
  std::map<TString,TH2F*>     m_histoTH2F;
  std::map<TString,TProfile*> m_histoTProfile;

  //JB
  std::map<TString,TTree*>    m_histoTTree;

  
 public:
  //! \brief Create and store TH1F histogram
  void createTH1F(TString name, int Nbins, double xmin, double xmax, TString title = "");

  //! \brief Create and store TH1F histogram
  void createTH1F(TString name, const std::vector<double> &bins, TString title = "");

  //! \brief Create and store TH2F histogram
  void createTH2F(TString name, int NbinsX, double xmin, double xmax, int NBinsY, double ymin, double ymax, TString title = ""); 

  //! \brief Create and store TH2F histogram
  void createTH2F(TString name, const std::vector<double> &xbins, const std::vector<double> &ybins, TString title = ""); 

  //! \brief Create and store TProfile histogram
  void createTProfile(TString name, int NbinsX, double xmin, double xmax, TString title = ""); 

  //! \brief Create and store TProfile histogram
  void createTProfile(TString name, const std::vector<double> &xbins, TString title = ""); 

  
  
  //! \brief Fill existing TH1F histogram
  inline void fillTH1F(TString name, double x, double w=1.0) {getTH1F(name)->Fill(x,w);}
  
  //! \brief Fill existing TH2F histogram
  inline void fillTH2F(TString name, double x, double y, double w=1.0) {getTH2F(name)->Fill(x, y, w);}
  
  //! \brief Fill existing TH2F histogram
  inline void fillTProfile(TString name, double x, double y, double w=1.0) {getTProfile(name)->Fill(x, y, w);}
  
  
  //! \brief check wether a given TH1F exist in the store
  inline bool hasTH1F(TString name) {return m_histoTH1F.count(name)>0;}

  //! \brief check wether a given TH2F exist in the store
  inline bool hasTH2F(TString name) {return m_histoTH2F.count(name)>0;}

  //! \brief check wether a given TH2F exist in the store
  inline bool hasTProfile(TString name) {return m_histoTProfile.count(name)>0;}
  


  //! \brief Retrieve TH1F histogram from internal store
  TH1F* getTH1F(TString name);

  //! \brief Retrieve TH2F histogram from internal store
  TH2F* getTH2F(TString name);

  //! \briefRetrieve TProfile histogram from internal store
  TProfile* getTProfile(TString name);

  //! \brief Retrieve List of all histograms in internal store
  std::vector<TH1*> getListOfHistograms();

  //JB
  TTree *getTTree(TString name);
  void createTTree(TString name);
  inline void fillTTree(TString name) { m_histoTTree[name]->Fill(); };
  std::vector<TTree*> getListOfTrees();
};

#endif
