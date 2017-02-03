#include "ResponseMatrixMaker.hxx"
#include<cstdlib>
#include<iostream>
#include<cmath> // std::sqrt


namespace
{
  float deltaPhi(float phi1, float phi2)
  {
    const float pi = 3.14159265;
    float out = phi1 - phi2;
    while(out > pi)
      out -= 2. * pi;
    while(out < -1. * pi)
      out += 2. * pi;

    return out;
  }

  float deltaR(float eta1, float phi1, float eta2, float phi2)
  {
    float dPhi = deltaPhi(phi1, phi2);
    float dEta = eta1 - eta2;
    return std::sqrt(dPhi * dPhi + dEta * dEta);
  }
};


template<typename T>
ResponseMatrixMakerBase<T>::ResponseMatrixMakerBase(const Str& channel,
                                                    const Str& varName,
                                                    const Vec<float>& binning) :
  varName(varName),
  channel(channel),
  binning(binning),
  nPDFVariations(100),
  iAlphaSUp(nPDFVariations),
  iAlphaSDn(nPDFVariations+1)
{
  scale = 1.;
  skipSyst = false;

  pdfResponses = TH3D("null", "null", 1, 0., 1., 1, 0., 1., 1, 0., 1.);
}


// template<typename T>
// void ResponseMatrixMakerBase<T>::registerLeptonSFs(const TH2F& h, const Str& lep,
//                                               const Str& upOrDown = "")
// {
//   if(leptonSFHists.find(lep) == UMap<Str, TH2F&>::end)
//     leptonSFHists[lep] = UMap<Str, TH2F&>();
//
//   leptonSFHists[lep][upOrDown] = h;
// }


template<typename T>
const TH2D& ResponseMatrixMakerBase<T>::getResponse(const Str& syst)
{
  if(!responses.size())
    setup();

  return responses.at(syst);
}


template<typename T>
const TH3D& ResponseMatrixMakerBase<T>::getPDFResponses()
{
  if(!responses.size())
    setup();

  return pdfResponses;
}


template<typename T>
const Vec<TH2D>& ResponseMatrixMakerBase<T>::getScaleResponses()
{
  if(!responses.size())
    setup();

  return scaleResponses;
}


template<typename T>
Vec<Str> ResponseMatrixMakerBase<T>::knownSystematics() const
{
  Vec<Str> out;
  for(auto&& resp : responses)
    {
      if(!resp.first.empty())
        out.push_back(resp.first);
    }

  return out;
}


template<typename T>
void ResponseMatrixMakerBase<T>::setup()
{
  if(responses.size())
    return;

  UPtr<TChain> recoTree(new TChain((getChannel()+"/ntuple").c_str(),
                                   ("recoChain_"+getVar() + "_" + getChannel()).c_str()));
  for(const auto& fn : fileNames)
    recoTree->Add(fn.c_str());

  // set up lots of things
  Vec<Str> systs = Vec<Str>({"",
        "pu_up","pu_dn",
        });
  const bool isJetVar = (varName.find("jet") != Str::npos ||
                         varName.find("Jet") != Str::npos ||
                         varName.find("jj") != Str::npos);
  if(isJetVar && !skipSyst)
    {
      systs.push_back("jer_up");
      systs.push_back("jer_dn");
      systs.push_back("jes_up");
      systs.push_back("jes_dn");
    }
  bool hasE, hasMu;
  Vec<Str> objects;
  if(channel.find("eeee") != Str::npos)
    {
      objects = Vec<Str>({"e1","e2","e3","e4"});
      hasE = true;
      hasMu = false;
      if(!skipSyst)
        {
          systs.push_back("eEff_up");
          systs.push_back("eEff_dn");
          systs.push_back("eScale_up");
          systs.push_back("eScale_dn");
          systs.push_back("eRhoRes_up");
          systs.push_back("eRhoRes_dn");
          systs.push_back("ePhiRes_up");
        }
    }
  else if(channel.find("eemm") != Str::npos)
    {
      objects = Vec<Str>({"e1","e2","m1","m2"});
      hasE = true;
      hasMu = true;
      if(!skipSyst)
        {
          systs.push_back("eEff_up");
          systs.push_back("eEff_dn");
          systs.push_back("mEff_up");
          systs.push_back("mEff_dn");
          systs.push_back("eScale_up");
          systs.push_back("eScale_dn");
          systs.push_back("eRhoRes_up");
          systs.push_back("eRhoRes_dn");
          systs.push_back("ePhiRes_up");
          systs.push_back("mClosure_up");
          systs.push_back("mClosure_dn");
        }
    }
  else
    {
      objects = Vec<Str>({"m1","m2","m3","m4"});
      hasE = false;
      hasMu = true;
      if(!skipSyst)
        {
          systs.push_back("mEff_up");
          systs.push_back("mEff_dn");
          systs.push_back("mClosure_up");
          systs.push_back("mClosure_dn");
        }
    }

  // Scale and PDF systematics only done for samples that have LHE info (e.g.
  // not MCFM)
  bool hasLHE = bool(recoTree->FindBranch("pdfWeights"));
  const Vec<size_t> scaleIndicesWeCareAbout = {1,2,3,4,6,8};
  if(hasLHE && !skipSyst)
    {
      for(auto i : scaleIndicesWeCareAbout)
        {
          scaleResponses.push_back(
            TH2D(("scaleVariation"+std::to_string(i)).c_str(), "",
                 binning.size()-1, &binning[0],
                 binning.size()-1, &binning[0]));
        }

      if(!skipSyst)
        {
          Vec<float> iterationBins;
          for(size_t i = 0; i <= nPDFVariations; ++i)
            iterationBins.push_back(float(i));
          pdfResponses = TH3D("pdfResponses", "",
                              binning.size()-1, &binning[0],
                              binning.size()-1, &binning[0],
                              iterationBins.size()-1, &iterationBins[0]);

          auto scaleWtPtr = &scaleWeights;
          auto pdfWtPtr = &pdfAndAlphaSWeights;
          recoTree->SetBranchAddress("scaleWeights", &scaleWtPtr);//&&scaleWeights);
          recoTree->SetBranchAddress("pdfWeights", &pdfWtPtr);//&&pdfAndAlphaSWeights);


          systs.push_back("alphaS_up");
          systs.push_back("alphaS_dn");
        }
    }

  for(auto& s : systs)
    {
      responses[s] = TH2D(s.c_str(), "",
                          binning.size()-1, &binning[0],
                          binning.size()-1, &binning[0]);
    }

  responses[""] = TH2D("nominal", "",
                       binning.size()-1, &binning[0],
                       binning.size()-1, &binning[0]);

  // Get gen info
  UPtr<TChain> trueTree = UPtr<TChain>(new TChain((getChannel()+"Gen/ntuple").c_str(),
                                                  ("trueChain_"+getVar() + "_" + getChannel()).c_str()));
  for(const auto& fn : fileNames)
    trueTree->Add(fn.c_str());

  UPtr<UMap<size_t, T> > trueVals = this->getTrueValues(*trueTree);
  trueTree.reset(); // gone -- don't use any more

  // Set up common branches
  // lPt = Vec<float>(objects.size());
  // lEta = Vec<float>(objects.size());
  lSF = Vec<float>(objects.size());
  lSFErr = Vec<float>(objects.size());
  setCommonBranches(*recoTree, objects);
  this->setRecoBranches(*recoTree, objects);

  bool doPUWt = (puWeightHists.find("") != puWeightHists.end());
  bool doPUWtUp = (puWeightHists.find("up") != puWeightHists.end());
  bool doPUWtDn = (puWeightHists.find("dn") != puWeightHists.end());

  // Loop through base reco tree, fill most things
  for(size_t row = 0; row < size_t(std::abs(recoTree->GetEntries())); ++row)
    {
      recoTree->GetEntry(row);

      // elements needed for event weights
      float puWt = (doPUWt ? getBinFromHist(puWeightHists.at(""), truePU) : 1.);
      float puWtUp = (doPUWtUp ? getBinFromHist(puWeightHists.at("up"), truePU) : 1.);
      float puWtDn = (doPUWtDn ? getBinFromHist(puWeightHists.at("dn"), truePU) : 1.);

      float lepSF = 1.;
      float lepSFEUp = 1.;
      float lepSFEDn = 1.;
      float lepSFMUp = 1.;
      float lepSFMDn = 1.;
      for(size_t i = 0; i < lSF.size(); ++i)
        {
          lepSF *= lSF.at(i);
          if(objects.at(i).compare(0,1,"e") == 0)
            {
              lepSFEUp *= lSF.at(i) + lSFErr.at(i);
              lepSFEDn *= lSF.at(i) - lSFErr.at(i);
            }
          if(objects.at(i).compare(0,1,"m") == 0)
            {
              lepSFMUp *= lSF.at(i) + lSFErr.at(i);
              lepSFMDn *= lSF.at(i) - lSFErr.at(i);
            }
        }

      // for(size_t i = 0; i < objects.size(); ++i)
      //   {
      //     if(objects.at(i)[0] == "e")
      //       {
      //         lepSF *= getBinFromHist(leptonSFHists["e"][""], lPt.at(i), std::abs(lEta.at(i)));
      //         lepSFEUp *= getBinFromHist(leptonSFHists["e"]["up"], lPt.at(i), std::abs(lEta.at(i)));
      //         lepSFEDn *= getBinFromHist(leptonSFHists["e"]["dn"], lPt.at(i), std::abs(lEta.at(i)));
      //       }
      //     if(objects.at(i)[0] == "m")
      //       {
      //         lepSF *= getBinFromHist(leptonSFHists["m"][""], lPt.at(i), lEta.at(i));
      //         lepSFMUp *= getBinFromHist(leptonSFHists["m"]["up"], lPt.at(i), lEta.at(i));
      //         lepSFMDn *= getBinFromHist(leptonSFHists["m"]["dn"], lPt.at(i), lEta.at(i));
      //       }
      //   }

      auto iTrue = trueVals->find(evt);
      if(iTrue == trueVals->end())
        continue;
      const T& trueVal = iTrue->second;

      if(this->selectEvent())
        {
          // Nominal value
          const T val = this->getEventResponse();

          float nominalWeight = scale * puWt * lepSF * genWeight;

          // fill histos that use nominal value but with different weights
          this->fillResponse(responses[""], val, trueVal, nominalWeight);

          if(!skipSyst)
            {
              this->fillResponse(responses["pu_up"], val, trueVal, scale * puWtUp * lepSF * genWeight);
              this->fillResponse(responses["pu_dn"], val, trueVal, scale * puWtDn * lepSF * genWeight);

              if(hasE)
                {
                  this->fillResponse(responses["eEff_up"], val, trueVal, scale * puWt * lepSFEUp * genWeight);
                  this->fillResponse(responses["eEff_dn"], val, trueVal, scale * puWt * lepSFEDn * genWeight);
                }
              if(hasMu)
                {
                  this->fillResponse(responses["mEff_up"], val, trueVal, scale * puWt * lepSFMUp * genWeight);
                  this->fillResponse(responses["mEff_dn"], val, trueVal, scale * puWt * lepSFMDn * genWeight);
                }

              if(hasLHE)
                {
                  // fill once for each scale variation
                  for(size_t ind = 0; ind < scaleIndicesWeCareAbout.size(); ++ind)
                    this->fillResponse(scaleResponses.at(ind), val, trueVal,
                                       nominalWeight * scaleWeights.at(scaleIndicesWeCareAbout.at(ind)));

                  // fill the 3-D histogram with one response for each PDF variation
                  for(size_t ind = 0; ind < nPDFVariations; ++ind)
                    this->fillResponse(pdfResponses, val, trueVal, ind,
                                       nominalWeight * pdfAndAlphaSWeights.at(ind));

                  // the last two items in the PDF weight vector are alpha_S variations
                  this->fillResponse(responses["alphaS_up"], val, trueVal,
                                     nominalWeight * pdfAndAlphaSWeights.at(iAlphaSUp));
                  this->fillResponse(responses["alphaS_dn"], val, trueVal,
                                     nominalWeight * pdfAndAlphaSWeights.at(iAlphaSDn));
                }

              // changes to jet scale/resolution actually change numbers
              if(isJetVar)
                {
                  if(this->selectEvent("jer_up"))
                    this->fillResponse(responses["jer_up"], this->getEventResponse("jer_up"),
                                       trueVal, nominalWeight);
                  if(this->selectEvent("jer_dn"))
                    this->fillResponse(responses["jer_dn"], this->getEventResponse("jer_dn"),
                                       trueVal, nominalWeight);
                  if(this->selectEvent("jes_up"))
                    this->fillResponse(responses["jes_up"], this->getEventResponse("jes_up"),
                                       trueVal, nominalWeight);
                  if(this->selectEvent("jes_dn"))
                    this->fillResponse(responses["jes_dn"], this->getEventResponse("jes_dn"),
                                       trueVal, nominalWeight);
                }
            }
        }
    }

  recoTree.reset(); // gone -- don't use any more

  Map<Str,Str> systTreesNeeded;
  if(hasE && !skipSyst)
    {
      if(systFileNames.find("eScaleUp") != systFileNames.end())
        systTreesNeeded["eScale_up"] = "eScaleUp";
      if(systFileNames.find("eScaleDn") != systFileNames.end())
        systTreesNeeded["eScale_dn"] = "eScaleDn";
      if(systFileNames.find("eRhoResUp") != systFileNames.end())
        systTreesNeeded["eRhoRes_up"] = "eRhoResUp";
      if(systFileNames.find("eRhoResDn") != systFileNames.end())
        systTreesNeeded["eRhoRes_dn"] = "eRhoResDn";
      if(systFileNames.find("ePhiResUp") != systFileNames.end())
        systTreesNeeded["ePhiRes_up"] = "ePhiResUp";
    }
  if(hasMu && !skipSyst)
    {
      if(systFileNames.find("mClosureUp") != systFileNames.end())
        systTreesNeeded["mClosure_up"] = "mClosureUp";
      if(systFileNames.find("mClosureDn") != systFileNames.end())
        systTreesNeeded["mClosure_dn"] = "mClosureDn";
    }

  // systematics requiring other ntuples
  for(auto& treeInfo : systTreesNeeded)
    {
      const Str& systName = treeInfo.first;
      const Str& treeName = treeInfo.second;

      UPtr<TChain> t = UPtr<TChain>(new TChain((getChannel()+"/ntuple").c_str(),
                                               ("chain_"+getVar() + "_" + getChannel()).c_str()));
      for(const auto& fn : systFileNames[treeName])
        t->Add(fn.c_str());

      setCommonBranches(*t, objects);
      this->setRecoBranches(*t, objects);
      for(size_t row = 0; row < size_t(std::abs(t->GetEntries())); ++row)
        {
          t->GetEntry(row);

          float puWt = (doPUWt ? getBinFromHist(puWeightHists.at(""), truePU) : 1.);

          float lepSF = 1.;
          for(size_t i = 0; i < objects.size(); ++i)
            lepSF *= lSF.at(i);

           // for(size_t i = 0; i < objects.size(); ++i)
           //   {
           //     if(objects.at(i)[0] == "e")
           //       {
           //         lepSF *= getBinFromHist(leptonSFHists["e"][""], lPt.at(i), std::abs(lEta.at(i)));
           //       }
           //     if(objects.at(i)[0] == "m")
           //       {
           //         lepSF *= getBinFromHist(leptonSFHists["m"][""], lPt.at(i), lEta.at(i));
           //       }
           //   }

          auto iTrue = trueVals->find(evt);
          if(iTrue == trueVals->end())
            continue;
          const T& trueVal = iTrue->second;

          if(this->selectEvent(systName))
            {
              this->fillResponse(responses[systName], this->getEventResponse(),
                                 trueVal, scale * puWt * lepSF * genWeight);
            }
        }
    } // new tree disappears here
}


template<typename T>
void ResponseMatrixMakerBase<T>::setCommonBranches(TChain& t, const Vec<Str>& objects)
{
  // if(lPt.size() != objects.size())
  //   lPt.resize(objects.size());
  // if(lEta.size() != objects.size())
  //   lEta.resize(objects.size());
  if(lSF.size() != objects.size())
    lSF.resize(objects.size());
  if(lSFErr.size() != objects.size())
    lSFErr.resize(objects.size());

  for(size_t i = 0; i < objects.size(); ++i)
    {
      // t.SetBranchAddress((objects[i]+"Pt").c_str(), &(lPt[i]));
      // t.SetBranchAddress((objects[i]+"Eta").c_str(), &(lEta[i]));
      t.SetBranchAddress((objects[i]+"EffScaleFactor").c_str(), &(lSF[i]));
      t.SetBranchAddress((objects[i]+"EffScaleFactorError").c_str(), &(lSFErr[i]));
    }
  t.SetBranchAddress("nTruePU", &truePU);
  t.SetBranchAddress("evt", &evt);
  t.SetBranchAddress("genWeight", &genWeight);
}


template<typename T>
BranchValueResponseMatrixMaker<T>::BranchValueResponseMatrixMaker(const Str& channel, const Str& varName,
                                                                  const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<T>(channel, varName, binning)
{;}


template<typename T>
UPtr<UMap<size_t, T> >
BranchValueResponseMatrixMaker<T>::getTrueValues(TChain& t, const Str& syst) const
{
  UPtr<UMap<size_t, T> > out(new UMap<size_t, T>());

  unsigned long long trueEvt;
  T val;
  t.SetBranchAddress("evt", &trueEvt);
  t.SetBranchAddress(this->getVar().c_str(), &val);

  for(size_t row = 0; row < size_t(std::abs(t.GetEntries())); ++row)
    {
      t.GetEntry(row);

      (*out)[trueEvt] = val;
    }

  return std::move(out);
}


template<typename T> void
BranchValueResponseMatrixMaker<T>::setRecoBranches(TChain& t, const Vec<Str>& objects)
{
  t.SetBranchAddress(this->getVar().c_str(), &value);
}


template<typename T> T
BranchValueResponseMatrixMaker<T>::getEventResponse(const Str& opt) const
{
  return value;
}


template<class R>
AbsValueResponseMatrixMaker<R>::AbsValueResponseMatrixMaker(const Str& channel,
                                                            const Str& varName,
                                                            const Vec<float>& binning) :
  R(channel, varName, binning)
{;}


template<class R>
UPtr<UMap<size_t, typename R::ValType> >
AbsValueResponseMatrixMaker<R>::getTrueValues(TChain& trueTree,
                                              const Str& syst) const
{
  UPtr<UMap<size_t, typename R::ValType> > out = R::getTrueValues(trueTree, syst);
  for(auto& val : *out)
    this->doAbs(val.second);

  return std::move(out);
}


template<class R>
typename R::ValType
AbsValueResponseMatrixMaker<R>::getEventResponse(const Str& option) const
{
  typename R::ValType out = R::getEventResponse(option);
  this->doAbs(out);
  return out;
}


template<class R>
typename R::ValType&
AbsValueResponseMatrixMaker<R>::doAbs(typename R::ValType& val) const
{
  val = std::abs(val);
  return val;
}


template<typename T>
JetBranchResponseMatrixMakerBase<T>::JetBranchResponseMatrixMakerBase(const Str& channel,
                                                                      const Str& varName,
                                                                      const Vec<float>& binning) :
  BranchValueResponseMatrixMaker<T>(channel, varName, binning)
{;}


template<typename T>
T JetBranchResponseMatrixMakerBase<T>::getEventResponse(const Str& syst) const
{
  if(syst.find("jes_up") != Str::npos)
    return value_jesUp;
  if(syst.find("jes_dn") != Str::npos)
    return value_jesDn;
  if(syst.find("jer_up") != Str::npos)
    return value_jerUp;
  if(syst.find("jer_dn") != Str::npos)
    return value_jerDn;

  return BranchValueResponseMatrixMaker<T>::getEventResponse(syst);
}


template<typename T>
void
JetBranchResponseMatrixMakerBase<T>::setRecoBranches(TChain& t,
                                                     const Vec<Str>& objects)
{
  BranchValueResponseMatrixMaker<T>::setRecoBranches(t, objects);

  t.SetBranchAddress((this->getVar()+"_jesUp").c_str(), &value_jesUp);
  t.SetBranchAddress((this->getVar()+"_jesDown").c_str(), &value_jesDn);
  t.SetBranchAddress((this->getVar()+"_jerUp").c_str(), &value_jerUp);
  t.SetBranchAddress((this->getVar()+"_jerDown").c_str(), &value_jerDn);
}


DijetBranchResponseMatrixMaker::DijetBranchResponseMatrixMaker(const Str& channel,
                                                               const Str& varName,
                                                               const Vec<float>& binning) :
  JetBranchResponseMatrixMakerBase<float>(channel, varName, binning)
{;}


UPtr<UMap<size_t, float> >
DijetBranchResponseMatrixMaker::getTrueValues(TChain& trueTree,
                                              const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  float trueVal;
  unsigned int trueNJets;
  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(this->getVar().c_str(), &trueVal);
  trueTree.SetBranchAddress("nJets", &trueNJets);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      if(trueNJets >= 2)
        (*out)[trueEvt] = trueVal;
    }

  return std::move(out);
}


void
DijetBranchResponseMatrixMaker::setRecoBranches(TChain& t, const Vec<Str>& objects)
{
  JetBranchResponseMatrixMakerBase<float>::setRecoBranches(t, objects);

  t.SetBranchAddress("nJets", &nJets);
  t.SetBranchAddress("nJets_jesUp", &nJets_jesUp);
  t.SetBranchAddress("nJets_jesDown", &nJets_jesDn);
  t.SetBranchAddress("nJets_jerUp", &nJets_jerUp);
  t.SetBranchAddress("nJets_jerDown", &nJets_jerDn);
}


bool
DijetBranchResponseMatrixMaker::selectEvent(const Str& syst) const
{
  if(syst.empty())
    return nJets >= 2;

  if(syst.find("jes_up") != Str::npos)
    return nJets_jesUp >= 2;
  if(syst.find("jes_dn") != Str::npos)
    return nJets_jesDn >= 2;
  if(syst.find("jer_up") != Str::npos)
    return nJets_jerUp >= 2;
  if(syst.find("jer_dn") != Str::npos)
    return nJets_jerDn >= 2;

  return nJets >= 2;
}


SelectedZResponseMatrixMakerBase::SelectedZResponseMatrixMakerBase(const Str& channel,
                                                                   const Str& varName,
                                                                   const Vec<float>& binning,
                                                                   const Str& compVar) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  z1CompVarName(constructZVarName(channel, compVar, true)),
  z2CompVarName(constructZVarName(channel, compVar, false)),
  z1RespVarName(constructZVarName(channel, varName, true)),
  z2RespVarName(constructZVarName(channel, varName, false)),
  compIsResp(z1CompVarName == z1RespVarName || z2CompVarName == z2RespVarName)
{;}


UPtr<UMap<size_t, float> >
SelectedZResponseMatrixMakerBase::getTrueValues(TChain& trueTree,
                                                const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  float z1CompTrue;
  float z2CompTrue;
  float z1RespTrueValue;
  float z2RespTrueValue;
  float* z1RespTruePtr;
  float* z2RespTruePtr;

  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(z1CompVarName.c_str(), &z1CompTrue);
  trueTree.SetBranchAddress(z2CompVarName.c_str(), &z2CompTrue);
  if(compIsResp)
    {
      z1RespTruePtr = &z1CompTrue;
      z2RespTruePtr = &z2CompTrue;
    }
  else
    {
      trueTree.SetBranchAddress(z1RespVarName.c_str(), &z1RespTrueValue);
      trueTree.SetBranchAddress(z2RespVarName.c_str(), &z2RespTrueValue);
      z1RespTruePtr = &z1RespTrueValue;
      z2RespTruePtr = &z2RespTrueValue;
    }

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      (*out)[trueEvt] = (this->z1IsBetter(z1CompTrue, z2CompTrue) ?
                         *z1RespTruePtr :
                         *z2RespTruePtr);
    }

  return std::move(out);
}


void
SelectedZResponseMatrixMakerBase::setRecoBranches(TChain& t,
                                                  const Vec<Str>& objects)
{
  t.SetBranchAddress(z1CompVarName.c_str(), &z1CompVar);
  t.SetBranchAddress(z2CompVarName.c_str(), &z2CompVar);
  if(compIsResp)
    {
      z1RespVarPtr = &z1CompVar;
      z2RespVarPtr = &z2CompVar;
    }
  else
    {
      t.SetBranchAddress(z1RespVarName.c_str(), &z1RespVarValue);
      t.SetBranchAddress(z2RespVarName.c_str(), &z2RespVarValue);
      z1RespVarPtr = &z1RespVarValue;
      z2RespVarPtr = &z2RespVarValue;
    }
}


float
SelectedZResponseMatrixMakerBase::getEventResponse(const Str& option) const
{
  if(this->z1IsBetter(z1CompVar, z2CompVar))
    return *z1RespVarPtr;
  return *z2RespVarPtr;
}


Str
SelectedZResponseMatrixMakerBase::constructZVarName(const Str& channel,
                                                    const Str& var,
                                                    bool wantZ1) const
{
  if(channel == "mmmm")
    {
      if(wantZ1)
        return "m1_m2_" + var;
      return "m3_m4_" + var;
    }

  if(wantZ1)
    return "e1_e2_" + var;

  if(channel == "eemm")
    return "m1_m2_" + var;

  return "e3_e4_" + var;
}


Z1ByMassResponseMatrixMaker::Z1ByMassResponseMatrixMaker(const Str& channel,
                                                         const Str& varName,
                                                         const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Mass"),
  comparisonVarName("Mass")
{;}


bool Z1ByMassResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                             const float z2Comp) const
{
  return std::abs(z1Comp - 91.1876) < std::abs(z2Comp - 91.1876);
}


Z2ByMassResponseMatrixMaker::Z2ByMassResponseMatrixMaker(const Str& channel,
                                                         const Str& varName,
                                                         const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Mass"),
  comparisonVarName("Mass")
{;}


bool Z2ByMassResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                             const float z2Comp) const
{
  return std::abs(z1Comp - 91.1876) > std::abs(z2Comp - 91.1876);
}


Z1ByPtResponseMatrixMaker::Z1ByPtResponseMatrixMaker(const Str& channel,
                                                     const Str& varName,
                                                     const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Pt"),
  comparisonVarName("Pt")
{;}


bool Z1ByPtResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                           const float z2Comp) const
{
  return z1Comp > z2Comp;
}


Z2ByPtResponseMatrixMaker::Z2ByPtResponseMatrixMaker(const Str& channel,
                                                     const Str& varName,
                                                     const Vec<float>& binning) :
  SelectedZResponseMatrixMakerBase(channel, varName, binning, "Pt"),
  comparisonVarName("Pt")
{;}


bool Z2ByPtResponseMatrixMaker::z1IsBetter(const float z1Comp,
                                           const float z2Comp) const
{
  return z1Comp < z2Comp;
}


ZZCompositeResponseMatrixMakerBase::ZZCompositeResponseMatrixMakerBase(const Str& channel,
                                                                       const Str& varName,
                                                                       const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  z1VarName(constructZVarName(channel, varName, true)),
  z2VarName(constructZVarName(channel, varName, false))
{;}


UPtr<UMap<size_t, float> >
ZZCompositeResponseMatrixMakerBase::getTrueValues(TChain& trueTree,
                                                  const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  float z1VarTrue;
  float z2VarTrue;

  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(z1VarName.c_str(), &z1VarTrue);
  trueTree.SetBranchAddress(z2VarName.c_str(), &z2VarTrue);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      (*out)[trueEvt] = this->calculateZZVar(z1VarTrue, z2VarTrue);
    }

  return std::move(out);
}


void
ZZCompositeResponseMatrixMakerBase::setRecoBranches(TChain& t,
                                                    const Vec<Str>& objects)
{
  t.SetBranchAddress(z1VarName.c_str(), &z1InputVar);
  t.SetBranchAddress(z2VarName.c_str(), &z2InputVar);
}


float
ZZCompositeResponseMatrixMakerBase::getEventResponse(const Str& option) const
{
  return this->calculateZZVar(z1InputVar, z2InputVar);
}


Str
ZZCompositeResponseMatrixMakerBase::constructZVarName(const Str& channel,
                                                      const Str& var,
                                                      bool wantZ1) const
{
  if(channel == "mmmm")
    {
      if(wantZ1)
        return "m1_m2_" + var;
      return "m3_m4_" + var;
    }

  if(wantZ1)
    return "e1_e2_" + var;

  if(channel == "eemm")
    return "m1_m2_" + var;

  return "e3_e4_" + var;
}


ZZDeltaPhiResponseMatrixMaker::ZZDeltaPhiResponseMatrixMaker(const Str& channel,
                                                             const Str& varName,
                                                             const Vec<float>& binning) :
  ZZCompositeResponseMatrixMakerBase(channel, "Phi", binning)
{;}


float ZZDeltaPhiResponseMatrixMaker::calculateZZVar(float z1Phi,
                                                    float z2Phi) const
{
  return ::deltaPhi(z1Phi, z2Phi);
}


ZZDeltaRResponseMatrixMaker::ZZDeltaRResponseMatrixMaker(const Str& channel,
                                                         const Str& varName,
                                                         const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  z1EtaBranchName(constructZVarName(channel, "Eta", true)),
  z2EtaBranchName(constructZVarName(channel, "Eta", false)),
  z1PhiBranchName(constructZVarName(channel, "Phi", true)),
  z2PhiBranchName(constructZVarName(channel, "Phi", false))
{;}


UPtr<UMap<size_t, float> >
ZZDeltaRResponseMatrixMaker::getTrueValues(TChain& trueTree,
                                           const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  float z1EtaTrue;
  float z2EtaTrue;
  float z1PhiTrue;
  float z2PhiTrue;

  trueTree.SetBranchAddress("evt", &trueEvt);
  trueTree.SetBranchAddress(z1EtaBranchName.c_str(), &z1EtaTrue);
  trueTree.SetBranchAddress(z2EtaBranchName.c_str(), &z2EtaTrue);
  trueTree.SetBranchAddress(z1PhiBranchName.c_str(), &z1PhiTrue);
  trueTree.SetBranchAddress(z2PhiBranchName.c_str(), &z2PhiTrue);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      (*out)[trueEvt] = ::deltaR(z1EtaTrue, z1PhiTrue, z2EtaTrue, z2PhiTrue);
    }

  return std::move(out);
}

void
ZZDeltaRResponseMatrixMaker::setRecoBranches(TChain& t, const Vec<Str>& objects)
{
  t.SetBranchAddress(z1EtaBranchName.c_str(), &z1Eta);
  t.SetBranchAddress(z2EtaBranchName.c_str(), &z2Eta);
  t.SetBranchAddress(z1PhiBranchName.c_str(), &z1Phi);
  t.SetBranchAddress(z2PhiBranchName.c_str(), &z2Phi);
}


float
ZZDeltaRResponseMatrixMaker::getEventResponse(const Str& option) const
{
  return ::deltaR(z1Eta, z1Phi, z2Eta, z2Phi);
}


Str
ZZDeltaRResponseMatrixMaker::constructZVarName(const Str& channel,
                                               const Str& var,
                                               bool wantZ1) const
{
  if(channel == "mmmm")
    {
      if(wantZ1)
        return "m1_m2_" + var;
      return "m3_m4_" + var;
    }

  if(wantZ1)
    return "e1_e2_" + var;

  if(channel == "eemm")
    return "m1_m2_" + var;

  return "e3_e4_" + var;
}


template<typename T>
MultiBranchResponseMatrixMakerBase<T>::MultiBranchResponseMatrixMakerBase(const Str& channel,
                                                                          const Str& varName,
                                                                          const Vec<float>& binning) :
  VectorValueResponseMatrixMakerBase<T>(channel, varName, binning)
{;}


template<typename T>
UPtr<UMap<size_t, Vec<T> > >
MultiBranchResponseMatrixMakerBase<T>::getTrueValues(TChain& trueTree,
                                                     const Str& syst) const
{
  UPtr<UMap<size_t, Vec<T> > > out(new UMap<size_t, Vec<T> >());

  const Vec<Str> varNames = this->constructVarNames(this->getChannel(),
                                                    this->getVar());

  unsigned long long trueEvt;
  Vec<T> trueValues(varNames.size());

  trueTree.SetBranchAddress("evt", &trueEvt);
  for(size_t i = 0; i < varNames.size(); ++i)
    trueTree.SetBranchAddress(varNames[i].c_str(), &trueValues[i]);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      // copies elements so should be safe
      (*out)[trueEvt] = trueValues;
    }

  return std::move(out);
}


template<typename T>
void
MultiBranchResponseMatrixMakerBase<T>::setRecoBranches(TChain& t,
                                                       const Vec<Str>& objects)
{
  const Vec<Str> varNames = this->constructVarNames(this->getChannel(),
                                                    this->getVar());

  values.resize(varNames.size());

  // this is technically dangerous, but I think it should be ok as long as
  // values never changes size between calls to this function
  for(size_t i = 0; i < varNames.size(); ++i)
    t.SetBranchAddress(varNames[i].c_str(), &values[i]);
}


template<typename T>
Vec<T>
MultiBranchResponseMatrixMakerBase<T>::getEventResponse(const Str& option) const
{
  // copies elements so should be safe
  return values;
}


template<typename T>
Vec<Str>
MultiBranchResponseMatrixMakerBase<T>::constructVarNames(const Str& channel,
                                                         const Str& var) const
{
  Vec<Str> objects = this->constructObjectNames(channel);

  Vec<Str> out;
  for(const auto& obj : objects)
    out.push_back(obj + var);

  return out;
}


AllLeptonBranchResponseMatrixMaker::AllLeptonBranchResponseMatrixMaker(const Str& channel,
                                                                       const Str& varName,
                                                                       const Vec<float>& binning) :
  MultiBranchResponseMatrixMakerBase<float>(channel, varName, binning)
{;}


Vec<Str>
AllLeptonBranchResponseMatrixMaker::constructObjectNames(const Str& channel) const
{
  if(channel == "mmmm")
    return Vec<Str>({"m1","m2","m3","m4"});
  else if(channel == "eemm")
    return Vec<Str>({"m1","m2","e1","e2"});
  else
    return Vec<Str>({"e1","e2","e3","e4"});
}


BothZsBranchResponseMatrixMaker::BothZsBranchResponseMatrixMaker(const Str& channel,
                                                                 const Str& varName,
                                                                 const Vec<float>& binning) :
  MultiBranchResponseMatrixMakerBase<float>(channel, varName, binning)
{;}


Vec<Str>
BothZsBranchResponseMatrixMaker::constructObjectNames(const Str& channel) const
{
  if(channel == "mmmm")
    return Vec<Str>({"m1_m2_","m3_m4_"});
  else if(channel == "eemm")
    return Vec<Str>({"m1_m2_","e1_e2_"});
  else
    return Vec<Str>({"e1_e2_","e3_e4_"});
}


LeptonMaxBranchResponseMatrixMaker::LeptonMaxBranchResponseMatrixMaker(const Str& channel,
                                                                       const Str& varName,
                                                                       const Vec<float>& binning) :
  SimpleValueResponseMatrixMakerBase<float>(channel, varName, binning),
  varNames(constructVarNames(channel, varName))
{;}


UPtr<UMap<size_t, float> >
LeptonMaxBranchResponseMatrixMaker::getTrueValues(TChain& trueTree,
                                                  const Str& syst) const
{
  UPtr<UMap<size_t, float> > out(new UMap<size_t, float>());

  unsigned long long trueEvt;
  Vec<float> trueValues(varNames.size());

  trueTree.SetBranchAddress("evt", &trueEvt);
  for(size_t i = 0; i < varNames.size(); ++i)
    trueTree.SetBranchAddress(varNames[i].c_str(), &trueValues[i]);

  for(size_t row = 0; row < size_t(std::abs(trueTree.GetEntries())); ++row)
    {
      trueTree.GetEntry(row);

      float max = -999999999.;
      for(const auto& val : trueValues)
        {
          if(val > max)
            max = val;
        }

      (*out)[trueEvt] = max;
    }

  return std::move(out);
}


void
LeptonMaxBranchResponseMatrixMaker::setRecoBranches(TChain& t,
                                                    const Vec<Str>& objects)
{
  values.resize(varNames.size());

  // this is technically dangerous, but I think it should be ok as long as
  // values never changes size between calls to this function
  for(size_t i = 0; i < varNames.size(); ++i)
    t.SetBranchAddress(varNames[i].c_str(), &values[i]);
}


float
LeptonMaxBranchResponseMatrixMaker::getEventResponse(const Str& option) const
{
  float max = -999999999.;
  for(const auto& val : values)
    {
      if(val > max)
        max = val;
    }

  return max;
}

Vec<Str>
LeptonMaxBranchResponseMatrixMaker::constructVarNames(const Str& channel,
                                                      const Str& var) const
{
  Vec<Str> objects;

  if(channel == "mmmm")
    objects = {"m1","m2","m3","m4"};
  else if(channel == "eemm")
    objects = {"m1","m2","e1","e2"};
  else
    objects = {"e1","e2","e3","e4"};

  Vec<Str> out;
  for(const auto& obj : objects)
    out.push_back(obj + var);

  return out;
}





typedef SimpleValueResponseMatrixMakerBase<float> FloatResponseMatrixMakerBase;
typedef BranchValueResponseMatrixMaker<float> FloatBranchResponseMatrixMaker;
typedef AbsValueResponseMatrixMaker<FloatBranchResponseMatrixMaker> AbsFloatBranchResponseMatrixMaker;
typedef BranchValueResponseMatrixMaker<unsigned int> UIntBranchResponseMatrixMaker;
typedef AbsValueResponseMatrixMaker<DijetBranchResponseMatrixMaker> AbsDijetBranchResponseMatrixMaker;
typedef JetBranchResponseMatrixMakerBase<unsigned int> JetUIntBranchResponseMatrixMaker;
typedef JetBranchResponseMatrixMakerBase<float> JetFloatBranchResponseMatrixMaker;
typedef AbsValueResponseMatrixMaker<ZZDeltaPhiResponseMatrixMaker> ZZAbsDeltaPhiResponseMatrixMaker;


#if defined(__ROOTCLING__)
#pragma link C++ class FloatBranchResponseMatrixMaker;
#pragma link C++ class FloatResponseMatrixMaker;
#pragma link C++ class AbsFloatResponseMatrixMaker;
#pragma link C++ class UIntResponseMatrixMaker;
#pragma link C++ class DijetBranchResponseMatrixMaker;
#pragma link C++ class AbsDijetBranchResponseMatrixMaker;
#pragma link C++ class JetUIntBranchResponseMatrixMaker;
#pragma link C++ class JetFloatBranchResponseMatrixMaker;
#pragma link C++ class Z1ByMassResponseMatrixMaker;
#pragma link C++ class Z2ByMassResponseMatrixMaker;
#pragma link C++ class Z1ByPtResponseMatrixMaker;
#pragma link C++ class Z2ByPtResponseMatrixMaker;
#pragma link C++ class ZZDeltaRResponseMatrixMaker;
#pragma link C++ class ZZAbsDeltaPhiResponseMatrixMaker;
#pragma link C++ class AllLeptonBranchResponseMatrixMaker;
#pragma link C++ class LeptonMaxBranchResponseMatrixMaker;
#pragma link C++ class BothZsBranchResponseMatrixMaker;
#endif