// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONKFPARTICLE_H
#define QAG4SIMULATIONKFPARTICLE_H

#include <CLHEP/Vector/LorentzVector.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/IteratorRange.h>
#include <TH1.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <kfparticle_sphenix/KFParticle_particleList.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/getClass.h>
#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>
#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "QAHistManagerDef.h"

class PHCompositeNode;
class PHG4TruthInfoContainer;
class SvtxTrackMap;

class QAG4SimulationKFParticle : public SubsysReco
{
 public:
  QAG4SimulationKFParticle(const std::string &name = "QAG4SimulationKFParticle");

  virtual ~QAG4SimulationKFParticle() = default;

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  std::string get_histo_prefix();

  void setTrackMapName(const std::string &name) { m_trackMapName = name; }

 protected:
  SvtxClusterEval *clustereval = nullptr;

 private:
  int load_nodes(PHCompositeNode *);

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  PHG4Particle *getTruthTrack(SvtxTrack *thisTrack);
  CLHEP::HepLorentzVector *makeHepLV(PHCompositeNode *topNode, int track_number);

  PHG4TruthInfoContainer *m_truthContainer = nullptr;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;

  SvtxTrackMap *m_trackMap = nullptr;
  PHG4TruthInfoContainer *m_truthInfo = nullptr;

  std::string m_trackMapName = "SvtxTrackMap";
};

#endif  // QAG4SIMULATIONKFPARTICLE_H
