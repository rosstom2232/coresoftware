// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIOND0_H
#define QAG4SIMULATIOND0_H

#include <fun4all/SubsysReco.h>

#include <g4eval/SvtxEvalStack.h>

#include <memory>
#include <set>
#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHHepMCGenEventMap;

namespace HepMC
{
  class GenEvent;
}

class QAG4SimulationD0 : public SubsysReco
{
 public:
  QAG4SimulationD0(const std::string &name = "QAG4SimulationD0");

  virtual ~QAG4SimulationD0() = default;

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  std::string get_histo_prefix();

  void addEmbeddingID(int embeddingID);

 private:
  int load_nodes(PHCompositeNode *);

  PHG4TruthInfoContainer *m_truthContainer = nullptr;

  unsigned int _nlayers_maps = 3;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;

  SvtxTrackMap *m_trackMap = nullptr;

  std::set<int> m_embeddingIDs;

  PHHepMCGenEventMap *m_Geneventmap = nullptr;
  int _embedding_id = 1;
};

#endif  // QAG4SIMULATIOND0_H
