#include "QAG4SimulationKFParticle.h"

using namespace std;

enum MOTHER_ID
{
  D0 = 421,
};

// Create necessary objects
typedef std::pair<int, float> particle_pair;
KFParticle_particleList kfp_particleList_evtReco;

//Particle masses are in GeV
std::map<std::string, particle_pair> particleMasses = kfp_particleList_evtReco.getParticleList();

QAG4SimulationKFParticle::QAG4SimulationKFParticle(const std::string &name)
  : SubsysReco(name)
{
}

int QAG4SimulationKFParticle::InitRun(PHCompositeNode *topNode)
{
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
    m_svtxEvalStack->set_strict(false);
    m_svtxEvalStack->set_verbosity(Verbosity());
  }
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    cout << __PRETTY_FUNCTION__ << " Fatal Error : missing " << m_trackMapName << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_trackMap)
  {
    cout << __PRETTY_FUNCTION__ << " Fatal Error : missing G4TruthInfo" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationKFParticle::Init(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  h = new TH1F(TString(get_histo_prefix()) + "D0_InvMass",  //
               ";m_{K^{-}#pi^{+}} [GeV/c^{2}];Entries", 200, 1, 2);
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationKFParticle::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    cout << "QAG4SimulationKFParticle::process_event() entered" << endl;

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1F *h_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "D0_InvMass"));
  assert(h_mass);

  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);

  std::vector<CLHEP::HepLorentzVector> daughters;
  for (auto &[key, track] : *m_trackMap)
  {
    SvtxTrack *thisTrack = getTrack(key, m_trackMap);
    CLHEP::HepLorentzVector *theVector = makeHepLV(topNode, thisTrack->get_id());
    if (theVector) daughters.push_back(*theVector);
  }

  CLHEP::HepLorentzVector mother;
  cout << "Number of dauughters: " << daughters.size() << endl;
  if (daughters.size() >= 2)
  {
    for (CLHEP::HepLorentzVector daughter : daughters)
    {
      mother += daughter;
    }
  }

  h_mass->Fill(mother.m());

  cout << "Mother mass = " << mother.m() << endl;

  daughters.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

SvtxTrack *QAG4SimulationKFParticle::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = NULL;

  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    if (iter->first == track_id) matched_track = iter->second;
  }

  return matched_track;
}

PHG4Particle *QAG4SimulationKFParticle::getTruthTrack(SvtxTrack *thisTrack)
{
  if (!clustereval)
  {
    clustereval = m_svtxEvalStack->get_cluster_eval();
  }

  TrkrDefs::cluskey clusKey = *thisTrack->begin_cluster_keys();
  PHG4Particle *particle = clustereval->max_truth_particle_by_cluster_energy(clusKey);

  return particle;
}

CLHEP::HepLorentzVector *QAG4SimulationKFParticle::makeHepLV(PHCompositeNode *topNode, int track_number)
{
  SvtxTrack *track = getTrack(track_number, m_trackMap);
  PHG4Particle *g4particle = getTruthTrack(track);
  CLHEP::HepLorentzVector *lvParticle = NULL;

  PHHepMCGenEventMap *m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "Missing node PHHepMCGenEventMap" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    return NULL;
  }

  PHHepMCGenEvent *m_genevt = m_geneventmap->get(1);
  if (!m_genevt)
  {
    std::cout << "Missing node PHHepMCGenEvent" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    return NULL;
  }

  HepMC::GenEvent *theEvent = m_genevt->getEvent();

  bool breakOut = false;
  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
  {
    assert((*p));
    if ((*p)->barcode() == g4particle->get_barcode())
    {
      for (HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()->particles_begin(HepMC::parents);
           mother != (*p)->production_vertex()->particles_end(HepMC::parents); ++mother)
      {
        if (abs((*mother)->pdg_id()) == MOTHER_ID::D0)
        {
          double mass = 0;
          for (auto it = particleMasses.begin(); it != particleMasses.end(); ++it)
          {
            if (it->second.first == abs((*p)->pdg_id()))
            {
              mass = it->second.second;
            }
          }
          lvParticle = new CLHEP::HepLorentzVector();
          //lvParticle->setVectM(CLHEP::Hep3Vector(track->get_px(), track->get_py(), track->get_pz()), mass);
          lvParticle->setVectM(CLHEP::Hep3Vector(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz()), mass);
        }
        else
          continue;
        break;
      }
      breakOut = true;
    }
    if (breakOut) break;
  }

  return lvParticle;
}

int QAG4SimulationKFParticle::load_nodes(PHCompositeNode *topNode)
{
  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthContainer)
  {
    cout << "QAG4SimulationTracking::load_nodes - Fatal Error - "
         << "unable to find DST node "
         << "G4TruthInfo" << endl;
    assert(m_truthContainer);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

string QAG4SimulationKFParticle::get_histo_prefix()
{
  return string("h_") + Name() + string("_") + m_trackMapName + string("_");
}
