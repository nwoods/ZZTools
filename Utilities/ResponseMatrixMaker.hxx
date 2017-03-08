#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TChain.h"
#include "TROOT.h"

#include<vector>
#include<unordered_map>
#include<map>
#include<string>
#include<memory>

#include<iostream>

// less typing
template <typename T, typename U>
using Map = typename std::map<T,U>;
template <typename T, typename U>
using UMap = typename std::map<T,U>;
template <typename T>
using Vec = typename std::vector<T>;
typedef std::string Str;
template <typename T>
using UPtr = typename std::unique_ptr<T>;


template <typename T>
class ResponseMatrixMakerBase
{
 public:
  ResponseMatrixMakerBase(const Str& channel, const Str& varName,
                          const Vec<float>& binning);
  virtual ~ResponseMatrixMakerBase(){;}

  // And scale factor histograms for PU reweighting and lepton efficiency
  void registerPUWeights(const TH1D& h, const Str& upOrDown = "") {puWeightHists[upOrDown] = h;}
  // void registerLeptonSFs(const TH2F& h, const Str& lep,
  //                        const Str& upOrDown = "");

  // Add files to run over
  void registerFile(const Str& f) {fileNames.push_back(f);}
  // Same but files are associated to a systematic
  void registerFile(const Str& f, const Str& syst)
  {
    if(systFileNames.find(syst) == systFileNames.end())
      systFileNames[syst] = Vec<Str>();
    systFileNames[syst].push_back(f);
  }

  // Get the response histogram for a particular systematic (or the central
  // value with an empty string)
  // The first time this is called, setup() is called
  const TH2D& operator()(const Str& syst="") {return getResponse(syst);}
  const TH2D& getResponse(const Str& syst="");
  // For PDF errors, get the 100 copies of the response matrix
  const TH3D& getPDFResponses();
  // For the scale errors, get the vector of responses for the various shifts
  const Vec<TH2D>& getScaleResponses();

  const Str& getVar() const {return varName;}
  const Str& getChannel() const {return channel;}

  Vec<Str> knownSystematics() const;

  bool hasSystematic(const Str& s) const {return responses.find(s) != responses.end();}

  void setConstantScale(float c) {scale = c;}
  float getConstantScale() const {return scale;}

  // Can be set to calculate central values only
  void setSkipSystematics(bool skipIfTrue=true) {skipSyst = skipIfTrue;}
  bool willSkipSystematics() {return skipSyst;}

 protected:
  typedef T ValType;

  // Event number -> value(s)
  virtual UPtr<UMap<size_t, T> > getTrueValues(TChain& trueTree,
                                               const Vec<Str>& objects,
                                               const Str& syst = "") const = 0;

  // Point branches to correct addresses (child member objects)
  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects) = 0;

  // Use these branches to get the values to put into a response histo
  // Option can be a systematic to apply or whatever
  virtual T getEventResponse(const Str& option = "") const = 0;

  // Fill a histogram with a (set of) response(s)
  virtual void fillResponse(TH2D& h, const T& val, const T& trueVal,
                            float w) const = 0;
  // For PDF responses
  virtual void fillResponse(TH3D& h, const T& val, const T& trueVal,
                            unsigned iteration, float w) const = 0;

  virtual bool selectEvent(const Str& option = "") const = 0;

  virtual float getLepSF(const Vec<Str>& leptons,
                         float electronSystematic = 0.,
                         float muonSystematic = 0.);

  // requires both true Zs to be on-shell unless overridden
  virtual bool selectTrueEvent(float mZ1, float mZ2) const;

  // If there is already a pointer to obj1_obj2_Mass pointing to the right
  // branch in t, return it. Otherwise, set the branch to use maybeUseThis
  // and return a pointer to it.
  virtual float* getmZPtr(TChain& t, const Str& obj1, const Str& obj2,
                          float& maybeUseThis) const;

  UMap<Str,TH2D> responses;
  TH3D pdfResponses;
  Vec<TH2D> scaleResponses;

  // Vec<float> lPt;
  // Vec<float> lEta;

 private:
  // Set up branches from this class (as opposed to subclasses)
  void setCommonBranches(TChain& t, const Vec<Str>& objects);
  // Make and store response histograms for all systematics
  void setup();

  Vec<Str> fileNames;
  UMap<Str, Vec<Str> > systFileNames;
  const Str varName;
  const Str channel;

  Vec<float> lSF;
  Vec<float> lSFErr;
  float truePU;
  unsigned long long evt;
  float genWeight;
  Vec<float> scaleWeights;
  Vec<float> pdfAndAlphaSWeights;
  UMap<Str, TH1D> puWeightHists;
  UMap<Str, UMap<Str, TH2F> > leptonSFHists;

  const Vec<float> binning;
  const size_t nPDFVariations;
  const size_t iAlphaSUp;
  const size_t iAlphaSDn;
  float scale;
  bool skipSyst;
};


template<class T>
class SimpleValueResponseMatrixMakerBase : public ResponseMatrixMakerBase<T>
{
 public:
  SimpleValueResponseMatrixMakerBase(const Str& channel,
                                     const Str& varName,
                                     const Vec<float>& binning) :
    ResponseMatrixMakerBase<T>(channel, varName, binning)
  {;}
  virtual ~SimpleValueResponseMatrixMakerBase(){;}

 protected:
  typedef typename ResponseMatrixMakerBase<T>::ValType ValType;

  void fillResponse(TH2D& h, const T& val, const T& trueVal, float w) const
  {
    h.Fill(val, trueVal, w);
  }
  // For PDF responses
  void fillResponse(TH3D& h, const T& val, const T& trueVal,
                    unsigned iteration, float w) const
  {
    h.Fill(val, trueVal, float(iteration), w);
  }
};


template<typename T>
class BranchValueResponseMatrixMaker : public SimpleValueResponseMatrixMakerBase<T>
{
 public:
  BranchValueResponseMatrixMaker(const Str& channel, const Str& varName,
                                 const Vec<float>& binning);
  virtual ~BranchValueResponseMatrixMaker(){;}

 protected:
  typedef typename SimpleValueResponseMatrixMakerBase<T>::ValType ValType;

  virtual UPtr<UMap<size_t, T> > getTrueValues(TChain& trueTree,
                                               const Vec<Str>& objects,
                                               const Str& syst = "") const;

  // Point branches to correct addresses (child member objects)
  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

  // Use these branches to get the values to put into a response histo
  // Option can be a systematic to apply or whatever
  virtual T getEventResponse(const Str& option = "") const;

  virtual bool selectEvent(const Str& option = "") const {return true;}


 private:
  T value;
};


// Takes any single-value response matrix maker, and does the same thing but
// with abs(value)
template<class R> // R is another ResponseMatrixMaker class
class AbsValueResponseMatrixMaker : public R
{
 protected:
  typedef typename R::ValType T; // pretend we're like the other templates
 public:
  AbsValueResponseMatrixMaker(const Str& channel, const Str& varName,
                              const Vec<float>& binning);
  virtual ~AbsValueResponseMatrixMaker(){;}

 protected:
  typedef typename R::ValType ValType;

  virtual UPtr<UMap<size_t, T> > getTrueValues(TChain& trueTree,
                                               const Vec<Str>& objects,
                                               const Str& syst = "") const;

  virtual T getEventResponse(const Str& option = "") const;

  // val is modified in place and the same reference is returned
  virtual T& doAbs(T& val) const;// {val = std::abs(val); return val;}
};


// For jet-related variables that do not require an nJets cut
template<typename T>
class JetBranchResponseMatrixMakerBase : public BranchValueResponseMatrixMaker<T>
{
 public:
  JetBranchResponseMatrixMakerBase(const Str& channel, const Str& varName,
                                   const Vec<float>& binning);
  virtual ~JetBranchResponseMatrixMakerBase(){;}

 protected:
  typedef typename BranchValueResponseMatrixMaker<T>::ValType ValType;

  virtual T getEventResponse(const Str& syst = "") const;

  // sets jet systematic branches too
  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

 private:
  T value_jesUp;
  T value_jesDn;
  T value_jerUp;
  T value_jerDn;
};


class DijetBranchResponseMatrixMaker : public JetBranchResponseMatrixMakerBase<float>
{
 public:
  DijetBranchResponseMatrixMaker(const Str& channel, const Str& varName,
                                 const Vec<float>& binning);
  virtual ~DijetBranchResponseMatrixMaker(){;}

 protected:
  typedef JetBranchResponseMatrixMakerBase<float>::ValType ValType;

  virtual UPtr<UMap<size_t, float> > getTrueValues(TChain& trueTree,
                                                   const Vec<Str>& objects,
                                                   const Str& syst = "") const;

  // gets nJets branch (and nJets systematic branches) as well as the
  // value branch
  void setRecoBranches(TChain& t, const Vec<Str>& objects);

  // selects for nJets >= 2
  bool selectEvent(const Str& syst = "") const;

 private:
  unsigned int nJets;
  unsigned int nJets_jesUp;
  unsigned int nJets_jesDn;
  unsigned int nJets_jerUp;
  unsigned int nJets_jerDn;
};


class SelectedZResponseMatrixMakerBase : public SimpleValueResponseMatrixMakerBase<float>
{
 public:
  SelectedZResponseMatrixMakerBase(const Str& channel, const Str& varName,
                                   const Vec<float>& binning,
                                   const Str& comparisonVarName);
  virtual ~SelectedZResponseMatrixMakerBase(){;}

 protected:
  typedef typename SimpleValueResponseMatrixMakerBase<float>::ValType ValType;

  virtual UPtr<UMap<size_t, float> > getTrueValues(TChain& trueTree,
                                                   const Vec<Str>& objects,
                                                   const Str& syst = "") const;

  // Point branches to correct addresses (child member objects)
  void setRecoBranches(TChain& t, const Vec<Str>& objects);

  // Use these branches to get the values to put into a response histo
  // Option can be a systematic to apply or whatever
  virtual float getEventResponse(const Str& option = "") const;

  // functions to indicate whether to use the Z1 value or the Z2 value
  virtual bool z1IsBetter(const float z1Comp, const float z2Comp) const = 0;

  // for derived classes to report the name of the comparison variable
  virtual const Str& compVar() const = 0;

  virtual bool selectEvent(const Str& syst = "") const {return true;}

 private:
  Str constructZVarName(const Str& channel, const Str& var, bool wantZ1) const;

  // branches used for the comparison (mass or pt or whatever)
  float z1CompVar;
  float z2CompVar;
  const Str z1CompVarName;
  const Str z2CompVarName;
  // ptrs will point to compVars or respVarValues as necessary
  float z1RespVarValue;
  float z2RespVarValue;
  float* z1RespVarPtr;
  float* z2RespVarPtr;
  const Str z1RespVarName;
  const Str z2RespVarName;

  const bool compIsResp;
};


class Z1ByMassResponseMatrixMaker : public SelectedZResponseMatrixMakerBase
{
 public:
  Z1ByMassResponseMatrixMaker(const Str& channel, const Str& varName,
                              const Vec<float>& binning);
  virtual ~Z1ByMassResponseMatrixMaker(){;}

 protected:
  typedef SelectedZResponseMatrixMakerBase::ValType ValType;

  // functions to indicate whether to use the Z1 value or the Z2 value
  virtual bool z1IsBetter(const float z1Comp, const float z2Comp) const;

  // for derived classes to report the name of the comparison variable
  virtual const Str& compVar() const {return comparisonVarName;}

 private:
  const Str comparisonVarName;
};


class Z2ByMassResponseMatrixMaker : public SelectedZResponseMatrixMakerBase
{
 public:
  Z2ByMassResponseMatrixMaker(const Str& channel, const Str& varName,
                              const Vec<float>& binning);
  virtual ~Z2ByMassResponseMatrixMaker(){;}

 protected:
  typedef SelectedZResponseMatrixMakerBase::ValType ValType;

  // functions to indicate whether to use the Z1 value or the Z2 value
  virtual bool z1IsBetter(const float z1Comp, const float z2Comp) const;

  // for derived classes to report the name of the comparison variable
  virtual const Str& compVar() const {return comparisonVarName;}

 private:
  const Str comparisonVarName;
};


class Z1ByPtResponseMatrixMaker : public SelectedZResponseMatrixMakerBase
{
 public:
  Z1ByPtResponseMatrixMaker(const Str& channel, const Str& varName,
                            const Vec<float>& binning);
  virtual ~Z1ByPtResponseMatrixMaker(){;}

 protected:
  typedef SelectedZResponseMatrixMakerBase::ValType ValType;

  // functions to indicate whether to use the Z1 value or the Z2 value
  virtual bool z1IsBetter(const float z1Comp, const float z2Comp) const;

  // for derived classes to report the name of the comparison variable
  virtual const Str& compVar() const {return comparisonVarName;}

 private:
  const Str comparisonVarName;
};


class Z2ByPtResponseMatrixMaker : public SelectedZResponseMatrixMakerBase
{
 public:
  Z2ByPtResponseMatrixMaker(const Str& channel, const Str& varName,
                            const Vec<float>& binning);
  virtual ~Z2ByPtResponseMatrixMaker(){;}

 protected:
  typedef SelectedZResponseMatrixMakerBase::ValType ValType;

  // functions to indicate whether to use the Z1 value or the Z2 value
  virtual bool z1IsBetter(const float z1Comp, const float z2Comp) const;

  // for derived classes to report the name of the comparison variable
  virtual const Str& compVar() const {return comparisonVarName;}

 private:
  const Str comparisonVarName;
};


class ZZCompositeResponseMatrixMakerBase : public SimpleValueResponseMatrixMakerBase<float>
{
 public:
  ZZCompositeResponseMatrixMakerBase(const Str& channel, const Str& varName,
                                     const Vec<float>& binning);
  virtual ~ZZCompositeResponseMatrixMakerBase(){;}

 protected:
  typedef typename SimpleValueResponseMatrixMakerBase<float>::ValType ValType;

  virtual UPtr<UMap<size_t, float> > getTrueValues(TChain& trueTree,
                                                   const Vec<Str>& objects,
                                                   const Str& syst = "") const;

  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

  virtual float getEventResponse(const Str& option = "") const;

  virtual bool selectEvent(const Str& option = "") const {return true;}

  virtual float calculateZZVar(float z1Var, float z2Var) const = 0;

  Str constructZVarName(const Str& channel, const Str& var, bool wantZ1) const;

 private:
  const Str z1VarName;
  const Str z2VarName;
  float z1InputVar;
  float z2InputVar;
};


class ZZDeltaPhiResponseMatrixMaker : public ZZCompositeResponseMatrixMakerBase
{
 public:
  ZZDeltaPhiResponseMatrixMaker(const Str& channel, const Str& varName,
                                const Vec<float>& binning);
  virtual ~ZZDeltaPhiResponseMatrixMaker(){;}

 protected:
  typedef ZZCompositeResponseMatrixMakerBase::ValType ValType;

  // bog standard delta phi function
  virtual float calculateZZVar(float z1Phi, float z2Phi) const;
};


class ZZDeltaRResponseMatrixMaker : public SimpleValueResponseMatrixMakerBase<float>
{
 public:
  ZZDeltaRResponseMatrixMaker(const Str& channel, const Str& varName,
                              const Vec<float>& binning);
  virtual ~ZZDeltaRResponseMatrixMaker(){;}

 protected:
  typedef typename SimpleValueResponseMatrixMakerBase<float>::ValType ValType;

  virtual UPtr<UMap<size_t, float> > getTrueValues(TChain& trueTree,
                                                   const Vec<Str>& objects,
                                                   const Str& syst = "") const;

  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

  virtual float getEventResponse(const Str& option = "") const;

  virtual bool selectEvent(const Str& option = "") const {return true;}

 private:
  Str constructZVarName(const Str& channel, const Str& var, bool wantZ1) const;

  float calculateDeltaR(float z1Eta, float z1Phi,
                        float z2Eta, float z2Phi) const;

  const Str z1EtaBranchName;
  const Str z2EtaBranchName;
  const Str z1PhiBranchName;
  const Str z2PhiBranchName;
  float z1Eta;
  float z2Eta;
  float z1Phi;
  float z2Phi;
};


template<typename T>
class VectorValueResponseMatrixMakerBase : public ResponseMatrixMakerBase<Vec<T> >
{
 public:
  VectorValueResponseMatrixMakerBase(const Str& channel, const Str& varName,
                                     const Vec<float>& binning) :
    ResponseMatrixMakerBase<Vec<T> >(channel, varName, binning)
  {;}

  virtual ~VectorValueResponseMatrixMakerBase(){;}

 protected:
  typedef typename ResponseMatrixMakerBase<Vec<T> >::ValType ValType;

  // Fill with all values in vector (only works for simple T right now)
  void fillResponse(TH2D& h, const Vec<T>& vals,
                    const Vec<T>& trueVals, float w) const
  {
    for(size_t i = 0; i < vals.size() && i < trueVals.size(); ++i)
      h.Fill(vals[i], trueVals[i], w);
  }
  // For PDF responses
  void fillResponse(TH3D& h, const Vec<T>& vals, const Vec<T>& trueVals,
                    unsigned iteration, float w) const
  {
    for(size_t i = 0; i < vals.size() && i < trueVals.size(); ++i)
      h.Fill(vals[i], trueVals[i], float(iteration), w);
  }
};


template<typename T>
class MultiBranchResponseMatrixMakerBase : public VectorValueResponseMatrixMakerBase<T>
{
 public:
  MultiBranchResponseMatrixMakerBase(const Str& channel, const Str& varName,
                                     const Vec<float>& binning);
  virtual ~MultiBranchResponseMatrixMakerBase(){;}

 protected:
  typedef typename VectorValueResponseMatrixMakerBase<T>::ValType ValType;

  virtual UPtr<UMap<size_t, Vec<T> > > getTrueValues(TChain& trueTree,
                                                     const Vec<Str>& objects,
                                                     const Str& syst = "") const;

  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

  virtual Vec<T> getEventResponse(const Str& option = "") const;

  virtual bool selectEvent(const Str& option = "") const {return true;}

  virtual Vec<Str> constructVarNames(const Str& channel, const Str& var) const;

  // For derived classes to specify objects
  virtual Vec<Str> constructObjectNames(const Str& channel) const = 0;

 private:
  Vec<T> values;
};


class AllLeptonBranchResponseMatrixMaker : public MultiBranchResponseMatrixMakerBase<float>
{
 public:
  AllLeptonBranchResponseMatrixMaker(const Str& channel, const Str& varName,
                                     const Vec<float>& binning);
  virtual ~AllLeptonBranchResponseMatrixMaker(){;}

 protected:
  typedef MultiBranchResponseMatrixMakerBase<float>::ValType ValType;

  Vec<Str> constructObjectNames(const Str& channel) const;
};


class BothZsBranchResponseMatrixMaker : public MultiBranchResponseMatrixMakerBase<float>
{
 public:
  BothZsBranchResponseMatrixMaker(const Str& channel, const Str& varName,
                                  const Vec<float>& binning);
  virtual ~BothZsBranchResponseMatrixMaker(){;}

 protected:
  typedef MultiBranchResponseMatrixMakerBase<float>::ValType ValType;

  Vec<Str> constructObjectNames(const Str& channel) const;
};


class LeptonMaxBranchResponseMatrixMaker : public SimpleValueResponseMatrixMakerBase<float>
{
 public:
  LeptonMaxBranchResponseMatrixMaker(const Str& channel, const Str& varName,
                                     const Vec<float>& binning);
  virtual ~LeptonMaxBranchResponseMatrixMaker(){;}

 protected:
  typedef SimpleValueResponseMatrixMakerBase<float>::ValType ValType;

  virtual UPtr<UMap<size_t, float> > getTrueValues(TChain& trueTree,
                                                   const Vec<Str>& objects,
                                                   const Str& syst = "") const;

  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

  virtual float getEventResponse(const Str& option = "") const;

  virtual bool selectEvent(const Str& option = "") const {return true;}

  virtual Vec<Str> constructVarNames(const Str& channel, const Str& var) const;

 private:
  const Vec<Str> varNames;
  Vec<float> values;
};


template<typename T, size_t _N>
class NthJetResponseMatrixMaker : public SimpleValueResponseMatrixMakerBase<T>
{
 public:
  NthJetResponseMatrixMaker(const Str& channel, const Str& varName,
                            const Vec<float>& binning);
  virtual ~NthJetResponseMatrixMaker(){;}

 protected:
  typedef typename SimpleValueResponseMatrixMakerBase<T>::ValType ValType;

  virtual UPtr<UMap<size_t, T> > getTrueValues(TChain& trueTree,
                                               const Vec<Str>& objects,
                                               const Str& syst = "") const;

  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

  virtual T getEventResponse(const Str& syst = "") const;

  virtual bool selectEvent(const Str& syst = "") const;

 private:
  // some pointer bullshit to appease the gods of ROOT
  UMap<Str,Vec<T>*> allJetValues; // values keyed to systematic
  Vec<T> allJetValues_object;
  Vec<T> allJetValues_jesUp_object;
  Vec<T> allJetValues_jesDn_object;
  Vec<T> allJetValues_jerUp_object;
  Vec<T> allJetValues_jerDn_object;
};


template<class R> // R is the type of the wrapped response maker
class UseSFHists : public R
{
 public:
  UseSFHists(const Str& channel, const Str& varName,
             const Vec<float>& binning);
  virtual ~UseSFHists(){;}

  void registerElectronSelectionSFHist(const TH2F& h);
  void registerElectronSelectionGapSFHist(const TH2F& h);
  void registerElectronRecoSFHist(const TH2F& h);
  void registerMuonSFHist(const TH2F& h);
  void registerMuonSFErrorHist(const TH2F& h);

 protected:
  typedef typename R::ValType ValType;

  virtual float getLepSF(const Vec<Str>& leptons,
                         float eSyst=0., float mSyst=0.);

  virtual void setRecoBranches(TChain& t, const Vec<Str>& objects);

 private:
  void setupOneLepton(TChain& t, const Str& obj,
                      float*& ptPtr, float& ptVal,
                      float*& etaPtr, float& etaVal,
                      bool*& isGapPtr, bool& isGapVal);

  UPtr<TH2F> hEleSelSF;
  UPtr<TH2F> hEleSelGapSF;
  UPtr<TH2F> hEleRecoSF;
  UPtr<TH2F> hMuSF;
  UPtr<TH2F> hMuSFErr;

  // pointers will point to floats below if this class makes the branches,
  // or to existing branch addresses if they're already made by the base class
  Vec<float*> lPtsSF;
  Vec<float*> lEtasSF; // supercluster etas for electrons
  Vec<bool*> lIsGapSF;
  float l1PtSF_value;
  float l2PtSF_value;
  float l3PtSF_value;
  float l4PtSF_value;
  float l1EtaSF_value;
  float l2EtaSF_value;
  float l3EtaSF_value;
  float l4EtaSF_value;
  bool l1IsGapSF_value;
  bool l2IsGapSF_value;
  bool l3IsGapSF_value;
  bool l4IsGapSF_value;
};


template<class R> // R is the type of the wrapped response maker
class RelaxGenZCuts : public R
{
 public:
  RelaxGenZCuts(const Str& channel, const Str& varName,
                const Vec<float>& binning);
  virtual ~RelaxGenZCuts(){;}

 protected:
  typedef typename R::ValType ValType;

  // no cuts
  bool selectTrueEvent(float mZ1, float mZ2) const {return true;}

 private:
  // no copy constructor
  RelaxGenZCuts(const RelaxGenZCuts<R>&);
};
