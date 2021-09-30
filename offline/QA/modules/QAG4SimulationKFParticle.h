// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONKFPARTICLE_H
#define QAG4SIMULATIONKFPARTICLE_H

#include <fun4all/SubsysReco.h>

#include <g4eval/SvtxEvalStack.h>

#include <CLHEP/Vector/LorentzVector.h> 

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>
#include <HepMC/IteratorRange.h> 

#include <memory>
#include <set>
#include <string>


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

 private:
  int load_nodes(PHCompositeNode *);

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  CLHEP::HepLorentzVector* makeHepLV(PHCompositeNode *topNode, int track_number);

  PHG4TruthInfoContainer *m_truthContainer = nullptr;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;

  SvtxTrackMap *m_trackMap = nullptr;
  PHG4TruthInfoContainer *m_truthInfo = nullptr;

  std::string m_trackMapName = "SvtxTrackMap";
};

#endif  // QAG4SIMULATIONKFPARTICLE_H
