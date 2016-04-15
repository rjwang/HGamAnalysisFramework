#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include "PATInterfaces/SystematicVariation.h"

#include <HGamAnalysisFramework/HgammaAnalysis.h>
#include <HGamAnalysisFramework/HgammaUtils.h>
#include <HGamAnalysisFramework/HGamVariables.h>

#include "PhotonVertexSelection/PhotonVertexHelpers.h"
#include "PhotonVertexSelection/PhotonPointingTool.h"

// #include "AsgTools/SgTEventMeta.h"
// #include "xAODMetaData/FileMetaData.h"

#include "TTree.h"
#include "TBranch.h"

typedef ElementLink<xAOD::TruthParticleContainer> TruthLink_t;
typedef ElementLink<xAOD::PhotonContainer> PhotonLink_t;

// this is needed to distribute the algorithm to the workers
ClassImp(HgammaAnalysis)

HgammaAnalysis :: HgammaAnalysis (const char *name)
: m_event(nullptr)
, m_store(nullptr)
, m_histoStore(nullptr)
, m_name(name)
, m_vertexTool(nullptr)
, m_photonHandler(nullptr)
, m_electronHandler(nullptr)
, m_jetHandler(nullptr)
, m_muonHandler(nullptr)
, m_eventHandler(nullptr)
, m_truthHandler(nullptr)
, m_overlapHandler(nullptr)
, m_etmissHandler(nullptr)
, m_catTool(nullptr)
, m_isInit(false)
, m_isAOD(false)
{
  // Must have no pointer initialization, for CINT
}

EL::StatusCode HgammaAnalysis :: setupJob (EL::Job& job)
{
  job.useXAOD ();

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init(m_name.Data()).ignore(); // call before opening first file

  // tell EventLoop about our output ntuple:
  EL::OutputStream out("MxAOD", "xAODNoMeta");
  job.outputAdd(out);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode HgammaAnalysis :: createOutput() {
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HgammaAnalysis :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  m_histoStore = new HistogramStore();

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (MetaData == nullptr)
    HG::fatal("Couldn't find MetaData TTree in event, is this a proper xAOD file? Exiting.");

  m_isAOD  = MetaData->GetBranch("StreamAOD");
  m_isMAOD = !MetaData->GetBranch("TriggerMenu");
  m_isDAOD = !m_isAOD && !m_isMAOD;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  // file name
  Info("changeInput", "Processing file \"%s\"",wk()->inputFile()->GetName());
  Info("changeInput", "This file has %lli entries", wk()->xaodEvent()->getEntries());

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  m_eventCounter = 0;

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  HG::VarHandler::getInstance()->setEventAndStore(event(), store());

  // asg::SgTEventMeta meta(asg::SgTEventMeta::InputStore);
  // const xAOD::FileMetaData *metad = nullptr;
  // if (meta.retrieve(metad, "MetaData").isFailure())
  //   HG::fatal("Couldn't retrieve MetaData");
  // std::string val;
  // metad->value(xAOD::FileMetaData::conditionsTag, val);

  // TEnv uses value from first file it's specified in.
  // If specified, read in additional configuration
  if (m_config.isDefined("Include"))
    for (TString cfg : m_config.getStrV("Include"))
      m_config.addFile(cfg);

  // Fill unspecified values from default config, specified here.
  if (!m_config.isDefined("BaseConfig")) {
     HG::fatal("You must specify a base configuration file, option: BaseConfig. Exiting.");
  } else {
     m_config.addFile(m_config.getStr("BaseConfig"));
  }

    //Check if isAFII is set in MCSample.cfg
    //Override the isAFII option in config if defined
    if(isMC()){
        const xAOD::EventInfo *eventInfo = 0;
        if (m_event->retrieve(eventInfo, "EventInfo").isFailure()){
          HG::fatal("Cannot access EventInfo");
        }
        int mcID = eventInfo->mcChannelNumber();
        if(m_config.isDefined(Form("IsAFII.%d",mcID))){
            m_config.setValue("IsAFII",m_config.getBool(Form("IsAFII.%d",mcID), false)?"YES":"NO");
        }
    }

  // Currently a hack, for passing whether it's an MxAOD to EventHandler
  if (m_isMAOD)
    m_config.setValue("IsMxAOD", "YES");
  else
    m_config.setValue("IsMxAOD", "NO");

  // What systematics to save
  StrV uncComps = config()->getStrV("HgammaAnalysis.UncertaintyComponents", {"*"});
  for (TString unc: uncComps) {
    TRegexp exp(unc, true);
    m_uncComps.push_back(exp);
  }

  // Print configuration database, if requested
  if (m_config.getBool("HgammaAnalysis.PrintConfig", true)) {
    Info("initialize()", "Printing full configuration:");
    m_config.printDB();
    Info("initialize()", " ");
  }

   // sample name
  TString sampleName = wk()->metaData()->castString("sample_name");
  Info("initialize()", "Sample name = %s", sampleName.Data());

  // Vertex selection tool
  CP::PhotonPointingTool *pointTool = new CP::PhotonPointingTool("PointingTool");
  if (pointTool->initialize().isFailure())
    HG::fatal("Failed vertex init");

  ToolHandle<CP::IPhotonPointingTool> tpoint(pointTool);

  m_vertexTool = new CP::PhotonVertexSelectionTool("PhotonVertexSelectionTool");
  CP_CHECK(m_name, m_vertexTool->setProperty("PhotonPointingTool", tpoint));

  if (m_vertexTool->initialize().isFailure())
    HG::fatal("Failed vertex init");

  m_markPhotonCand = m_config.getBool("PhotonHandler.SortCandidatesFirst", false);

  m_photonHandler = new HG::PhotonHandler("PhotonHandler", m_event, m_store);
  m_photonHandler->initialize(m_config);

  m_electronHandler = new HG::ElectronHandler("ElectronHandler", m_event, m_store);
  m_electronHandler->initialize(m_config);

  m_jetHandler = new HG::JetHandler("JetHandler", m_event, m_store);
  m_jetHandler->setAOD(m_isAOD);
  m_jetHandler->initialize(m_config);

  m_muonHandler = new HG::MuonHandler("MuonHandler", m_event, m_store);
  m_muonHandler->initialize(m_config);

  m_eventHandler = new HG::EventHandler(m_event, m_store);
  m_eventHandler->initialize(m_config);

  m_truthHandler = new HG::TruthHandler(m_event, m_store);
  m_truthHandler->initialize(m_config);

  m_overlapHandler = new HG::OverlapRemovalHandler();
  m_overlapHandler->initialize(m_config);

  m_etmissHandler = new HG::ETmissHandler("ETmissHandler", m_event, m_store);
  m_etmissHandler->initialize(m_config);

  m_catTool = new HG::HGamCategoryTool(m_event, m_store);
  m_catTool->initialize(m_config);

  // Check for HgammaAnalysis specific configs
  m_calcCat = m_config.getBool("HgammaAnalysis.CalculateCouplingCategory", true);

  m_doTwoGoodPhotonsCut = m_config.getBool("HgammaAnalysis.CheckTwoGoodPhotons", true);

  m_doRelPtCut = m_config.getBool("HgammaAnalysis.CheckRelativePtCuts", true);
  m_relPtCut1  = m_config.getNum ("HgammaAnalysis.RelPtFractionFirst" , 0.35);
  m_relPtCut2  = m_config.getNum ("HgammaAnalysis.RelPtFractionSecond", 0.25);

  m_doVertex   = m_config.getBool("HgammaAnalysis.SelectVertex", true);
  m_doHardPV   = m_config.getBool("HgammaAnalysis.UseHardestVertex", false);

  m_doMyyCut   = m_config.getBool("HgammaAnalysis.CheckMyyWindowCut", true);
  m_myyLow     = m_config.getNum("HgammaAnalysis.LowMyyGeV",105.0)*HG::GeV;
  m_myyHigh    = m_config.getNum("HgammaAnalysis.HighMyyGeV",160.0)*HG::GeV;

  m_doJetClean  = m_config.getBool("HgammaAnalysis.CheckJetEventCleaning", false);
  m_jetCleanPt  = m_config.getNum("JetHandler.Selection.EventCleanMinPtGeV", 20.0)*HG::GeV;
  if (m_doJetClean && m_config.getStr("JetHandler.Selection.CutLevel", "") != "LooseBad")
    HG::fatal("Currently you must clean jets with LooseBad to check CheckJetEventCleaning.");

  // Set up trigger matching map
  m_doTrigMatch      = m_config.getBool("EventHandler.CheckTriggerMatching", false);
  m_requiredTriggers = m_config.getStrV("EventHandler.RequiredTriggers");
  for (auto trig: m_requiredTriggers) {
    m_trigMatch[trig] = TrigType::Undefined;
    TString temp = m_config.getStr("EventHandler.TriggerMatchType."+trig, "");
    if (temp == "DiPhoton"      ) m_trigMatch[trig] = TrigType::DiPhoton;
    if (temp == "DiMuon"        ) m_trigMatch[trig] = TrigType::DiMuon;
    if (temp == "DiElectron"    ) m_trigMatch[trig] = TrigType::DiElectron;
    if (temp == "SinglePhoton"  ) m_trigMatch[trig] = TrigType::SinglePhoton;
    if (temp == "SingleMuon"    ) m_trigMatch[trig] = TrigType::SingleMuon;
    if (temp == "SingleElectron") m_trigMatch[trig] = TrigType::SingleElectron;
  }

  // Get list of systematic uncertainties
  // Must be done after all helper tools defined
  const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
  auto recommendedSystematics = registry.recommendedSystematics();

  m_sysList.push_back(CP::SystematicSet());
  for (auto sys: recommendedSystematics) {
    // Check if we want to run over this shift
    TString sysname = sys.name().c_str();
    bool considerSys = false;
    for (TRegexp exp: m_uncComps) {
      if (sysname.Contains(exp)) {
        considerSys = true;
        break;
      }
    }
    if (not considerSys) continue;

    // Add the shift to the register
    if (sys.name().find("continuous") != std::string::npos) {
      TString sysname = sys.name();
      sysname.ReplaceAll("__continuous", "");
      sysname.ReplaceAll(" ", "_");

      m_sysList.push_back(CP::SystematicSet());
      m_sysList.back().insert(CP::SystematicVariation(sysname.Data(), 1));

      m_sysList.push_back(CP::SystematicSet());
      m_sysList.back().insert(CP::SystematicVariation(sysname.Data(), -1));
    } else {
      m_sysList.push_back(CP::SystematicSet());
      m_sysList.back().insert(sys);
    }
  }

  // Need this after tools initialized, so that all systematic histograms are made
  // histInitialize();
  createOutput();

  // move to the user class? This is pretty standard
  TFile *file = wk()->getOutputFile ("MxAOD");
  if(!m_event->writeTo(file).isSuccess()){
    Error("initialize()", "Failed to write event to output file!");
    return EL::StatusCode::FAILURE;
  }

  // register all histograms
  for (auto *histo : m_histoStore->getListOfHistograms()) {
    wk()->addOutput(histo);
  }

  m_isInit = true;

  // Get Photon Fake Rate 2D Histogram file
  TString fileLocation = PathResolverFindCalibFile(config()->getStr("HgammaAnalysis.PhotonSelection.PhotonFakeFile").Data());
  m_PhotonFakeRateFile = TFile::Open(fileLocation, "read");

  if(m_PhotonFakeRateFile == 0)
      Fatal("Initialize","Failed to open photon fake rate file...");

  m_PhotonFakes2DHist = (TH2F*)m_PhotonFakeRateFile->Get("PhotonFakeRates");
  if(m_PhotonFakes2DHist == 0)
      Fatal("Initialize","Failed to fetch photon fake histogram...");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if (!m_isInit)
    HG::fatal("HgammaAnalysis was not initiialized. Did you forget to call HgammaAnalysis::initialize() ?");

  if(m_eventCounter==0) m_startTime = time(nullptr); //time in seconds

  static int progressInterval = config()->getInt("OutputMessage.ProcessedEventsInterval",1000);
  if ( m_eventCounter && m_eventCounter % progressInterval == 0 ) {
    Info("execute()","%i events processed so far  <<<===",
         static_cast< int >(m_eventCounter));
    Info("execute()","Processing rate = %.3f Hz",
         float(m_eventCounter)/(time(nullptr)-m_startTime));
  }
  m_eventCounter++;

  // This function will print the errors, no checking is required
  CP_CHECK(m_name, applySystematicVariation(CP::SystematicSet()));

  // Clear containers which point to objects from previous event
  HG::VarHandler::getInstance()->clearContainers();

  if (m_doVertex) selectVertex();

  setWeightInitial();

  return EL::StatusCode::SUCCESS;
}




EL::StatusCode HgammaAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode HgammaAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  Info("finalize()","Finished processing %i events", m_eventCounter);
  double nSecs = time(nullptr)-m_startTime;
  Info("finalize()","Total time elapsed: %dh %dm %ds",
       int(nSecs)/3600,(int(nSecs)%3600)/60,int(nSecs)%60);
  Info("finalize()","Processing rate = %.3f Hz", float(m_eventCounter)/(time(nullptr)-m_startTime));

  //Delete the Photon Fake Rates File
  if(m_PhotonFakeRateFile != 0)
      m_PhotonFakeRateFile->Close();

  SafeDelete(m_PhotonFakeRateFile);
  SafeDelete(m_vertexTool);
  SafeDelete(m_photonHandler);
  SafeDelete(m_electronHandler);
  SafeDelete(m_jetHandler);
  SafeDelete(m_muonHandler);
  SafeDelete(m_histoStore);
  SafeDelete(m_eventHandler);
  SafeDelete(m_truthHandler);
  SafeDelete(m_overlapHandler);
  SafeDelete(m_etmissHandler);
  SafeDelete(m_catTool);

  TFile *file = wk()->getOutputFile ("MxAOD");

  if(!m_event->finishWritingTo( file ).isSuccess() ) {
    Error("finalize()","Failed to finish writing event to output file!");
    return EL::StatusCode::FAILURE;
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode HgammaAnalysis :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}

/// Configures all handlers for systematic variation according to the specified systematic set
CP::SystematicCode HgammaAnalysis::applySystematicVariation(const CP::SystematicSet &sys)
{
  static const char *METHOD = "HgammaAnalysis::applySystematicVariation";
  CP_CHECK( METHOD, m_eventHandler               ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_photonHandler              ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_muonHandler                ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_electronHandler            ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_jetHandler                 ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, m_etmissHandler              ->applySystematicVariation(sys) );
  CP_CHECK( METHOD, HG::VarHandler::getInstance()->applySystematicVariation(sys) );

  setWeightInitial();

  return CP::SystematicCode::Ok;
}



bool HgammaAnalysis::pass(const xAOD::PhotonContainer   *photons,
                          const xAOD::ElectronContainer *electrons,
                          const xAOD::MuonContainer     *muons,
                          const xAOD::JetContainer      *jets)
{
  if (var::isPassed.exists())
    return var::isPassed();

  if (m_doTwoGoodPhotonsCut) {
    if (photons == nullptr) return false;
    if (!passTwoGoodPhotonsCut(*photons)) return false;
  }

  if (m_doTrigMatch &&
      !passTriggerMatch(photons, electrons, muons, jets))
    return false;

  if (m_doRelPtCut) {
    if (photons == nullptr) return false;
    if (!passRelativePtCuts(*photons)) return false;
  }

  if (m_doMyyCut) {
    if (photons == nullptr) return false;
    if (!passMyyWindowCut(*photons)) return false;
  }

  if (m_doJetClean &&
      !passJetEventCleaning())
    return false;

  return true;
}

// Return the generator Higgs mass, in GeV, from config specified by
//   GeneratorHiggsMass.MCCHANNELNUMBER
// if not defined, the code is aborted
double HgammaAnalysis::getGeneratorHiggsMass(int mcID)
{
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_higgsMass.find(mcID) == m_higgsMass.end())
    m_higgsMass[mcID] = config()->getNum(Form("GeneratorHiggsMass.%d",mcID));
  return m_higgsMass[mcID];
}

// Return the generator efficiency from config specified by
//   GeneratorEfficiency.MCCHANNELNUMBER
// if not defined, 1.0 is returned
double HgammaAnalysis::getGeneratorEfficiency(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_genEffs.find(mcID) == m_genEffs.end())
    m_genEffs[mcID] = config()->getNum(Form("GeneratorEfficiency.%d",mcID),1.0);
  return m_genEffs[mcID];
}

// Return the kFactor from config specified by
//   GeneratorEfficiency.MCCHANNELNUMBER
// if not defined, 1.0 is returned
double HgammaAnalysis::getKFactor(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_kFactors.find(mcID) == m_kFactors.end())
    m_kFactors[mcID] = config()->getNum(Form("kFactor.%d",mcID),1.0);
  return m_kFactors[mcID];
}

// Return the cross section, in pb, from config specified by
//   CrossSection.MCCHANNELNUMBER
// if not defined, the code is aborted
double HgammaAnalysis::getCrossSection(int mcID)
{
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_crossSections.find(mcID) == m_crossSections.end())
    m_crossSections[mcID] = config()->getNum(Form("CrossSection.%d", mcID));
  return m_crossSections[mcID];
}

TString HgammaAnalysis::getMCSampleName(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_mcNames.find(mcID) == m_mcNames.end())
    m_mcNames[mcID] = config()->getStr(Form("SampleName.%d",mcID)).ReplaceAll(" ", "");
  return m_mcNames[mcID];
}

int HgammaAnalysis::getNtotalEvents(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_nTotEvents.find(mcID) == m_nTotEvents.end())
    m_nTotEvents[mcID] = config()->getInt(Form("TotalNEvents.%d",mcID));
  return m_nTotEvents[mcID];
}

TH1F* HgammaAnalysis::getCutFlowHistogram(int mcID, TString suffix) {
  // access the initial number of weighed events
  TString cutFlowName(Form("CutFlow_%s%d%s",isData()?"Run":"MC",mcID,suffix.Data()));
  bool hasMCname = config()->isDefined(Form("SampleName.%d",mcID));
  if (hasMCname) cutFlowName = Form("CutFlow_%s%s",
				    getMCSampleName(mcID).Data(),suffix.Data());
  else Warning("","SampleName.%d not specfied in config!",mcID);
  TH1F *cflowHist = (TH1F*)wk()->inputFile()->Get(cutFlowName);
  if (cflowHist==nullptr)
    HG::fatal("Cannot access cut-flow histogram "+cutFlowName+" in input file");
  return cflowHist;
}

// intial sum of events, including pileup weights
double HgammaAnalysis::getIntialSumOfWeights(int mcID) {
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_NevtsInitial.find(mcID) == m_NevtsInitial.end()) {
    // Hard-coding to bin number 3 = ALLEVTS
    m_NevtsInitial[mcID] = getCutFlowHistogram(mcID,"_noDalitz_weighted")->GetBinContent(3);
  }
  return m_NevtsInitial[mcID];
}

double HgammaAnalysis::lumiXsecWeight(double intLumi, int mcID, bool printFirst) {
  if (intLumi < 0) intLumi = config()->getNum("IntegratedLuminosity_fbInv",1.0);
  if (mcID == -1) mcID = eventInfo()->mcChannelNumber();
  if (m_weightXsec.find(mcID) == m_weightXsec.end()) {
    double sigma      = getCrossSection(mcID);
    double gen_eff    = getGeneratorEfficiency(mcID);
    double kFactor    = getKFactor(mcID);
    double sumInitial = getIntialSumOfWeights(mcID);

    // Hard-coding to bin number 1,2
    double NxAOD      = getCutFlowHistogram(mcID,"_weighted")->GetBinContent(1);
    double NDxAOD     = getCutFlowHistogram(mcID,"_weighted")->GetBinContent(2);

    int Ntot = config()->isDefined(Form("TotalNEvents.%d",mcID)) ? getNtotalEvents(mcID) : -1;
    double skim_eff = NDxAOD / NxAOD;

    m_weightXsec[mcID] = intLumi * 1e3 * sigma * gen_eff * skim_eff * kFactor / sumInitial;
    if (printFirst) {
      printf("\nMC sample %d: %s\n",mcID,getMCSampleName(mcID).Data());
      printf("  Cross section:                %10.4e pb\n",sigma);
      if (gen_eff!=1.0) printf("  Generator efficiency:         %10.4e\n",gen_eff);
      if (kFactor!=1.0) printf("  k-factor:                     %10.2f\n",kFactor);
      if (Ntot) printf("  N events in AMI:              %10d\n",Ntot);
      printf("  sum w in xAOD:                %10.2f\n",NxAOD);
      printf("  sum w in DxAOD:               %10.2f\n",NDxAOD);
      if (skim_eff!=1.0)
	printf("  DxAOD efficiency:             %10.2f%%\n",skim_eff*100);
      printf("  Sum of inital event weights:  %10.2f\n\n",sumInitial);
      // L * sigma * eff * kFactor / Nevts
      printf("  Integrated lumi.:             %10.4f fb-1\n",intLumi);
      printf("  N exp. events for analysis:   %10.2e\n",intLumi*1e3*sigma*gen_eff*skim_eff*kFactor);
      printf("  Cross section event weight:   %10.4e\n\n",m_weightXsec[mcID]);
    }
  }
  return m_weightXsec[mcID];
}

void HgammaAnalysis::selectVertex()
{

  // If the event doesn't contain PVs, can't correct anything
  if (!m_event->contains<xAOD::VertexContainer>("PrimaryVertices")) {
    static bool first=true;
    if (first && not m_isMAOD) Warning("selectVertex","No PrimaryVertices container.%s",
                                       " No PV correction can be applied!!");
    first=false;
    return;
  }

  m_photonHandler->setVertexCorrected(false);
  m_jetHandler->setVertexCorrected(false);
  m_electronHandler->setVertexCorrected(false);
  m_muonHandler->setVertexCorrected(false);

  const xAOD::Vertex *vertex = nullptr;
  if (m_doHardPV) {
    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Couldn't retrieve PrimaryVertices, exiting!");

    vertex = xAOD::PVHelpers::getHardestVertex(vertices);
  } else {
    xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer();
    xAOD::PhotonContainer presel  = photonHandler()->applyPreSelection(photons);

    // If there aren't two photons, just use the hardest vertex
    if (presel.size() < 2) {
      const xAOD::VertexContainer *vertices = nullptr;
      if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
        HG::fatal("Couldn't retrieve PrimaryVertices, exiting!");

      vertex = xAOD::PVHelpers::getHardestVertex(vertices);
    } else {
      // Only use the two leading photons
      presel.resize(2);

      // Mark the nominal preselected photons as the Higgs candidate
      if (m_markPhotonCand                                &&
          (!photonHandler()->markAsCandidate(presel[0]) ||
           !photonHandler()->markAsCandidate(presel[1]) )  )
        Warning("selectVertex()", "Couldn't mark leading photons as Higgs candidates");

      // Get the pointed vertex
      m_vertexTool->getVertex(presel, vertex).ignore();
    }
  }

  if (vertex == nullptr) {
    const xAOD::VertexContainer *vertices = nullptr;
    if (m_event->retrieve(vertices, "PrimaryVertices").isFailure())
      HG::fatal("Couldn't retrieve PrimaryVertices, exiting!");

    vertex = xAOD::PVHelpers::getHardestVertex(vertices);
  }

  // Only set the vertex to use for kinematic corrections if it's not a nullptr
  if (vertex != nullptr) {
    ConstDataVector<xAOD::VertexContainer> *vcont = new ConstDataVector<xAOD::VertexContainer>(SG::VIEW_ELEMENTS);
    vcont->push_back(vertex);

    if (m_store->record(vcont, "HGamVertices").isFailure())
      HG::fatal("Couldn't add HGamVertices to TStore, exiting.");
  }

}



void HgammaAnalysis::setWeightInitial()
{
  // If already set, don't need to do it again
  // Also allows the code to run on MxAODs
  if (var::weightInitial.exists()) return;

  // Determine the initial event weight
  double weight = 1.0;

  if (isMC()) {
    // First MC generator weights
    weight *= eventHandler()->mcWeight();

    // Pileup weight
    weight *= eventHandler()->pileupWeight();

    // z-vertex weight
    weight *= eventHandler()->vertexWeight();
  }

  var::weightInitial.setValue(weight);
}


/// Get initial event weight: MC, pileup, z-vertex
double HgammaAnalysis::weightInitial()
{
  if (!var::weightInitial.exists()) {
    HG::fatal("Initial event weight not found, did you call HgammaAnalysis::execute() ?");
    return -1.0; // should never get called
  }

  return var::weightInitial();
}



/// Get event weight: initial weight * leading two photon SFs
double HgammaAnalysis::weight()
{
  if (!var::weight.exists()) {
    HG::fatal("You must call setSelectedObjects before retrieving the weight!");
    return -1.0; // should never get called
  }

  return var::weight();
}



/// Get category weight: weight * objects used for category selection
double HgammaAnalysis::weightCatCoup_dev()
{
  if (!var::weightCatCoup_dev.exists()) {
    HG::fatal("You must call setSelectedObjects before retrieving the category weight!");
    return -1.0; // should never get called
  }

  return var::weightCatCoup_dev();
}

/// Get category weight: weight * objects used for category selection
double HgammaAnalysis::weightCatCoup_Moriond2016()
{
  if (!var::weightCatCoup_Moriond2016.exists()) {
    HG::fatal("You must call setSelectedObjects before retrieving the category weight!");
    return -1.0; // should never get called
  }

  return var::weightCatCoup_Moriond2016();
}

void HgammaAnalysis::addTruthLinks(std::string recoName, std::string truthName)
{
  // Get the containers
  xAOD::IParticleContainer *recoCont = nullptr;
  if (m_event->retrieve(recoCont, recoName).isFailure())
    HG::fatal("Couldn't retrieve reco container");

  xAOD::TruthParticleContainer *truthCont = nullptr;
  if (m_event->retrieve(truthCont, truthName).isFailure())
    HG::fatal("Couldn't retrieve truth container");

  addTruthLinks(recoCont, truthCont);
}

void HgammaAnalysis::addTruthLinks(xAOD::IParticleContainer *recoCont, xAOD::TruthParticleContainer *truthCont)
{

  // Helpful variables
  std::vector<size_t> linked;
  static SG::AuxElement::Accessor<ElementLink<xAOD::IParticleContainer> > accTruthLink("truthLink");
  static SG::AuxElement::Accessor<ElementLink<xAOD::IParticleContainer> > accRecoLink("recoLink");
  static SG::AuxElement::Accessor<int> accPdgId("pdgId");
  static SG::AuxElement::Accessor<int> accParentPdgId("parentPdgId");

  //Loop over each reco particle
  for (auto recoPtcl : *recoCont) {
    // Grab the truth link from the reconstructed particle
    TruthLink_t link = recoPtcl->auxdata<TruthLink_t>("truthParticleLink");

    // Create a null truthPtcl. This is set to the truth particle in the new container
    xAOD::TruthParticle *truthPtcl = nullptr;

    accPdgId(*recoPtcl) = 0;
    accParentPdgId(*recoPtcl) = 0;

    // Only try search if we have a valid link
    if (link.isValid()) {
      // Get barcode of truth particle in original particle container
      int origPtclBarcode = (*link)->barcode();

      accPdgId(*recoPtcl) = (*link)->pdgId();
      if ((*link)->nParents() > 0 && (*link)->parent(0))
        accParentPdgId(*recoPtcl) = (*link)->parent()->pdgId();

      // Search through new container for the truth particle barcode (From old container link)
      for (auto ptcl : *truthCont) {
        if(ptcl->barcode() == origPtclBarcode) {
          truthPtcl = ptcl;
          break;
        }
      }
    }

    if (truthPtcl != nullptr) {
      ElementLink<xAOD::IParticleContainer> truthLink(*truthCont, truthPtcl->index());
      accTruthLink(*recoPtcl) = truthLink;

      ElementLink<xAOD::IParticleContainer> recoLink(*recoCont, recoPtcl->index());
      accRecoLink(*truthPtcl) = recoLink;

      linked.push_back(truthPtcl->index());
    } else {
      ElementLink<xAOD::IParticleContainer> truthLink;
      accTruthLink(*recoPtcl) = truthLink;
    }
  }

  // For all truth particles not linked already, add a null-link
  for (auto ptcl: *truthCont) {
    if (std::find(linked.begin(), linked.end(), ptcl->index()) == linked.end()) {
      ElementLink<xAOD::IParticleContainer> recoLink;
      accRecoLink(*ptcl) = recoLink;
    }
  }
}

xAOD::PhotonContainer HgammaAnalysis::getFakePhotons(double &photonFakeWeight){
    photonFakeWeight = 1;

    xAOD::PhotonContainer fakeCombination(SG::VIEW_ELEMENTS);

    //Do not carry on if we do not have MC (No truth info in data)
    if(!isMC()) return fakeCombination;

    if(event()->contains<xAOD::PhotonContainer>("HGamPhotonsWithFakes")){
        //fetch muons that overlap with jets - Needed for muon correction and they get removed in overlap tool
        const xAOD::PhotonContainer *HGamPhotonsWithFakes = 0;

        event()->retrieve(HGamPhotonsWithFakes,"HGamPhotonsWithFakes");
        fakeCombination = *HGamPhotonsWithFakes;

        photonFakeWeight = eventHandler()->getVar<float>("weightFakePhotons");

        return fakeCombination;
    }

    // Check if we have an MxAOD with no HGamPhotonsWithFakes container. Print warning and return empty container.
    if(isMAOD() && !event()->contains<xAOD::PhotonContainer>("HGamPhotonsWithFakes")){
        Warning("getFakePhotons()", "MxAOD as input but no HGamPhotonsWithFakes container... No fakes calculated");
        return fakeCombination;
    }

    //Note: We dont actually want to apply fakes if we do not pass the VERTEX cut. (No point)
    //No way to check cutflow from inside HgammaAnalysis.cxx so going to put this guard in HGamCutflowMxAOD.cxx when calling addFakePhotonDecoration.
    //if(m_cutFlow<= VERTEX) return fakeCombination;

    //Here fakePhotons are preSelected photons that do not overlap truth photons
    //Here realPhotons are the photons which do overlap and are matched to the applySelected photons
    xAOD::PhotonContainer fakePhotons(SG::VIEW_ELEMENTS);
    xAOD::PhotonContainer realPhotons(SG::VIEW_ELEMENTS);

    xAOD::TruthParticleContainer truthPhotons = truthHandler()->getPhotons();
    xAOD::PhotonContainer all_photons = photonHandler()->getCorrectedContainer();
    xAOD::PhotonContainer selPhotons = photonHandler()->applySelection(all_photons);

    for(auto photon : all_photons){
        //___________________________________________________________________________________________
        //Must pass preselction (Object Quality and pT/Eta Cuts)
        if (!photonHandler()->passPtEtaCuts(photon) ||
            !photonHandler()->passAmbCut(photon) ||
            !photonHandler()->passAuthorCut(photon)) continue;

        //Do not consider some photons that are not from real hard scattering process
        if(photon->auxdata<int>("truthOrigin") == 38 && photon->auxdata<int>("truthType") == 16) continue; //Bkg photon from UndrPhoton
        if(photon->auxdata<int>("truthOrigin") == 0  && photon->auxdata<int>("truthType") == 13) continue; //Unknown Photon from undefined origin
        if(photon->auxdata<int>("truthOrigin") == 5  && photon->auxdata<int>("truthType") == 4) continue;  //Bkg Electron from Photon Conversion
        //Sort into fake and real photons
        if(HG::minDR(photon, truthPhotons) < 0.4)
            realPhotons.push_back(photon);
        else
            fakePhotons.push_back(photon);
    }

    //Count the number of photons that are not fakes and pass Tight,Iso cuts
    int truePhotonCount = realPhotons.size();

    //Count the number of fake photons
    int fakePhotonCount = fakePhotons.size();

    //Do not continue if we can never make two photons
    if(truePhotonCount + fakePhotonCount < 2)return fakeCombination;

    //We should always go through the fake combination code
//    //If we have two real photons passing myy then we dont need to choose any fakes.
//    if(passRelativePtCuts(realPhotons) && passMyyWindowCut(realPhotons)){
//        fakeCombination = realPhotons;
//        return fakeCombination;
//    }

    //Fetch good combinations and choose one that pass myy cut if there is any
    std::vector<std::pair<xAOD::PhotonContainer,double>> goodCombi = getPhotonCombinations(realPhotons,fakePhotons);
    std::pair<xAOD::PhotonContainer,double> selCombination = chooseRandomComb(goodCombi);

    //Now check if we have at least two photon candidates
    if(selCombination.first.size()>=2){
        fakeCombination = selCombination.first;
        photonFakeWeight = selCombination.second;
        return fakeCombination; //No need to continue from here
    }
    return fakeCombination;
}

double HgammaAnalysis::getFPprob(xAOD::Photon *photon){
    TLorentzVector photon4V = photon->p4()*HG::invGeV;
    double photonPt = photon4V.Pt();
    double photonEta = std::fabs(photon4V.Eta());

    int xbin = m_PhotonFakes2DHist->GetXaxis()->FindBin(photonPt);
    int ybin = m_PhotonFakes2DHist->GetYaxis()->FindBin(photonEta);

    double prob = m_PhotonFakes2DHist->GetBinContent(xbin,ybin);

    return prob;
}

std::vector<std::pair<xAOD::PhotonContainer,double>> HgammaAnalysis::getPhotonCombinations(xAOD::PhotonContainer &realPhotons, xAOD::PhotonContainer &fakePhotons){
    std::vector<std::pair<xAOD::PhotonContainer,double>> fakeCombi;

    int fakeCount = fakePhotons.size();

    //Loop over all possible combinations of fake photons being promoted or not.
    //2 fakes have 4 possible unique combinations: (0 = fake, 1 = promoted)
    // 00 01 10 11 => This is given by 2^2 or 2^fakeCount which is what 1 << fakeCount does using bitshifts
    for (int combi = 0; combi < 1 << fakeCount; combi++){
        xAOD::PhotonContainer promotedPhotons = realPhotons;
        double fakeComboWt = 1;

        for(int i = 0; i < fakeCount; i++){
            double fakeWt = getFPprob(fakePhotons[i]);
            if ( combi&(1<<i) )
                promotedPhotons.push_back(fakePhotons[i]);
            else
                fakeWt = 1 - fakeWt;
            fakeComboWt*=fakeWt;
        }

        //Sort the true and promoted photons by pT
        promotedPhotons.sort(photonHandler()->comparePt);
        //Returning all combinations
        fakeCombi.push_back({promotedPhotons,fakeComboWt});
    }
    return fakeCombi;
}

std::pair<xAOD::PhotonContainer,double> HgammaAnalysis::chooseRandomComb(std::vector<std::pair<xAOD::PhotonContainer,double>> &fakeCombi){
    //Create empty selected photons and weight = 1;
    xAOD::PhotonContainer selPhotons (SG::VIEW_ELEMENTS);
    double combiWeight = 1;
    std::pair<xAOD::PhotonContainer,double> selCombi = {selPhotons,combiWeight};

    //First check size of fake combinations. If zero then we return empty photons and weight =1.
    int combiCount = fakeCombi.size();

    if(combiCount == 0)
        return selCombi;

    //Get ready to filter out the good combinations
    std::vector<std::pair<xAOD::PhotonContainer,double>> goodFakeCombi;

    double goodSumWt = 0;
    //Loop over all fake combinations and sum weights for good fake combinations (Ones that pass cuts)
    for(auto comb : fakeCombi)
        if(passRelativePtCuts(comb.first) && passMyyWindowCut(comb.first) && comb.second != 0 ){
            goodSumWt += comb.second;
            goodFakeCombi.push_back(comb);
        }
    //Now we have all good combinations

    //Starting value for random selection (Ranged based on combination weight). Generate random num.
    //Setting seed as event number
    double rndmNum = goodSumWt*gRandom->Rndm(eventInfo()->eventNumber());
    double start = 0.0;

    // Choose a random good combinations. Based on their weights. (Larger = more likely)
    for (int i = 0; i < (int)goodFakeCombi.size(); i++){
        // check to see where the random number falls between [0, goodSumWt]
        // pick the combination in whose range the random number falls
        if (rndmNum > start && rndmNum < (start+goodFakeCombi[i].second)){
            selCombi = {goodFakeCombi[i].first,goodSumWt};
            break;
        }
        start += goodFakeCombi[i].second;
    }
    return selCombi;
}


/// Set selected collections
void HgammaAnalysis::setSelectedTruthObjects(const xAOD::TruthParticleContainer *photons  ,
                                             const xAOD::TruthParticleContainer *electrons,
                                             const xAOD::TruthParticleContainer *muons    ,
                                             const xAOD::JetContainer           *jets     ,
                                             const xAOD::MissingETContainer     *mets     )
{
  HG::VarHandler::getInstance()->setTruthContainers(photons, electrons, muons, jets, mets);

  setWeightInitial();

  if (!var::weight.exists())
    var::weight.setValue(weightInitial());

  if (m_calcCat && !var::catCoup_dev.exists()) {
    var::catCoup_dev.setValue(0);
    var::weightCatCoup_dev.setValue(weightInitial());
  }

  if (m_calcCat && !var::catCoup_Moriond2016.exists()) {
    var::catCoup_Moriond2016.setValue(0);
    var::weightCatCoup_Moriond2016.setValue(weightInitial());
  }
}



/// Set selected collections
void HgammaAnalysis::setSelectedObjects(const xAOD::PhotonContainer    *photons  ,
                                        const xAOD::ElectronContainer  *electrons,
                                        const xAOD::MuonContainer      *muons    ,
                                        const xAOD::JetContainer       *jets     ,
                                        const xAOD::MissingETContainer *mets     )
{
  HG::VarHandler::getInstance()->setContainers(photons, electrons, muons, jets, mets);

  if (!var::weight.exists()) {
    // Determine total weight (leading photonSFs)
    static SG::AuxElement::Accessor<float> scaleFactor("scaleFactor");
    double myweight = weightInitial();
    if (photons != nullptr) {
      for (size_t i = 0; i < photons->size() && i < 2; ++i)
        myweight *= scaleFactor(*photons->at(i));
    }

    var::weight.setValue(myweight);
  }

  if (m_calcCat && !var::catCoup_dev.exists()) {
    // Determine the category and weight
    std::pair<int, float> catCoup_dev = m_catTool->getCategoryAndWeight(photons, electrons, muons, jets, mets);
    var::catCoup_dev.setValue(catCoup_dev.first);
    var::weightCatCoup_dev.setValue(catCoup_dev.second*var::weight());
  }

  if (m_calcCat && !var::catCoup_Moriond2016.exists()) {
    // Determine the category and weight
    std::pair<int, float> catCoup_Moriond2016 = m_catTool->getCategoryAndWeightMoriond(photons, electrons, muons, jets);
    var::catCoup_Moriond2016.setValue(catCoup_Moriond2016.first);
    var::weightCatCoup_Moriond2016.setValue(catCoup_Moriond2016.second*var::weight());
  }
}

/// Checks if event level jet cleaning cut is passed
bool HgammaAnalysis::passJetEventCleaning()
{
  if (var::isPassedJetEventClean.exists())
    return var::isPassedJetEventClean();

  static SG::AuxElement::ConstAccessor<char>  isClean("isClean");

  bool isEventClean = true;

  xAOD::JetContainer jets = m_jetHandler->getCorrectedContainer();
  for (auto jet: jets) {
    if (m_jetCleanJvt > 0.0) {
      // Cut when JVT is used
      if (jet->pt() > m_jetCleanPt      &&
          m_jetHandler->passJVTCut(jet) &&
          not isClean(*jet)             ) {
        isEventClean = false;
        break;
      }
    } else {
      // Cut when JVT is not used
      if (jet->pt() > m_jetCleanPt      &&
          not isClean(*jet)             ) {
        isEventClean = false;
        break;
      }
    }
  }

  var::isPassedJetEventClean.setValue(isEventClean);

  return isEventClean;
}



//______________________________________________________________________________
bool HgammaAnalysis::passTwoGoodPhotonsCut(const xAOD::PhotonContainer &photons)
{
  if (photons.size() < 2)
    return false;

  xAOD::PhotonContainer leading = photons;
  leading.resize(2);

  xAOD::PhotonContainer sel = photonHandler()->applySelection(leading);
  if (sel.size() < 2)
    return false;

  return true;
}

/// Checks if relative pT cuts for photons are passed
bool HgammaAnalysis::passRelativePtCuts(const xAOD::PhotonContainer &photons)
{
  // If there aren't two photons, the cut fails
  if (photons.size() < 2) return false;

  // Assume Higgs mass from two leading photons
  double myy = (photons[0]->p4() + photons[1]->p4()).M();

  // Check if relative pT cuts are satisfied
  if (photons[0]->pt()/myy < m_relPtCut1) return false;
  if (photons[1]->pt()/myy < m_relPtCut2) return false;

  return true;
}

/// Checks if myy is in the required window
bool HgammaAnalysis::passMyyWindowCut(const xAOD::PhotonContainer &photons)
{
  // If there aren't two photons, the cut fails
  if (photons.size() < 2) return false;
  double myy = (photons[0]->p4() + photons[1]->p4()).M();
  return m_myyLow <= myy && myy < m_myyHigh;
}


bool HgammaAnalysis::passTriggerMatch(const xAOD::PhotonContainer   *photons,
                                      const xAOD::ElectronContainer *electrons,
                                      const xAOD::MuonContainer     *muons,
                                      const xAOD::JetContainer      *jets)
{
  // Check whether at least one passing trigger is matched to selected objects
  for (auto trig: m_requiredTriggers) {
    if (m_eventHandler->passTrigger(trig) &&
        passTriggerMatch(trig.Data(), photons, electrons, muons, jets))
      return true;
  }

  return false;
}



//______________________________________________________________________________
bool HgammaAnalysis::passTriggerMatch(const TString &trig,
                                      const xAOD::PhotonContainer   *photons,
                                      const xAOD::ElectronContainer *electrons,
                                      const xAOD::MuonContainer     *muons,
                                      const xAOD::JetContainer      * /*jets*/)
{
  switch (m_trigMatch[trig]) {
    case TrigType::Undefined:
      return true;
    case TrigType::DiPhoton:
      return photons && photons->size() > 1 &&
             m_eventHandler->passTriggerMatch_DiPhoton(trig,
                                                       *photons->at(0),
                                                       *photons->at(1));
    case TrigType::DiMuon:
      return muons && muons->size() > 1 &&
             m_eventHandler->passTriggerMatch_DiMuon(trig,
                                                     *muons->at(0),
                                                     *muons->at(1));
    case TrigType::DiElectron:
      return electrons && electrons->size() > 1 &&
             m_eventHandler->passTriggerMatch_DiElectron(trig,
                                                         *electrons->at(0),
                                                         *electrons->at(1));
    case TrigType::SinglePhoton:
      return photons &&
             ( (photons->size() > 0 && m_eventHandler->passTriggerMatch_SinglePhoton(trig, *photons->at(0))) ||
               (photons->size() > 1 && m_eventHandler->passTriggerMatch_SinglePhoton(trig, *photons->at(1))) );
    case TrigType::SingleMuon:
      return muons &&
             ( (muons->size() > 0 && m_eventHandler->passTriggerMatch_SingleMuon(trig, *muons->at(0))) ||
               (muons->size() > 1 && m_eventHandler->passTriggerMatch_SingleMuon(trig, *muons->at(1))) );
    case TrigType::SingleElectron:
      return electrons &&
             ( (electrons->size() > 0 && m_eventHandler->passTriggerMatch_SingleElectron(trig, *electrons->at(0))) ||
               (electrons->size() > 1 && m_eventHandler->passTriggerMatch_SingleElectron(trig, *electrons->at(1))) );
    default:
      return true;
  }

  // If option isn't recognized above, default to failing match
  return false;
}



enum HgammaAnalysis::TrigType HgammaAnalysis::getTriggerType(TString Trigger)
{
  if (m_trigMatch.count(Trigger)>0)
    return m_trigMatch[Trigger];
  else
    return TrigType::Undefined;
}
