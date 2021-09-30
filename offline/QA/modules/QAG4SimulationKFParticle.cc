#include "QAG4SimulationKFParticle.h"

#include "QAHistManagerDef.h"

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxVertexEval.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>  // for pair
#include <vector>

using namespace std;

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

int QAG4SimulationKFParticle::Init(PHCompositeNode */*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  h = new TH1F(TString(get_histo_prefix()) + "D0_InvMass",  //
               "D0_InvMass;InvMass;Count", 50, 1.75, 1.95);
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

  TH1D *h_mass = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "D0_InvMass"));
  assert(h_mass);

  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);

  std::vector<CLHEP::HepLorentzVector> daughters;
  for(auto& [key, track] : *m_trackMap) 
  {
    SvtxTrack* thisTrack = getTrack(key, m_trackMap);
    daughters.push_back(*makeHepLV(topNode, thisTrack->get_id()));
  } 
  for(auto& output : daughters)
  { 
    cout << output.px() << ", " << output.py() << ", " << output.pz() << ", " << output.e() << endl;
  }
 
  daughters.clear();
  return Fun4AllReturnCodes::EVENT_OK;
}

 enum MOTHER_ID
  {  
    D0 = 421,
  };

SvtxTrack* QAG4SimulationKFParticle::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
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

CLHEP::HepLorentzVector* QAG4SimulationKFParticle::makeHepLV(PHCompositeNode *topNode, int track_number)
{
  /*
  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode*>(nodeIter.findFirst("SvtxTrackMap"));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  }
  */

  SvtxTrack *track = getTrack(track_number, m_trackMap);
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
    std::cout <<  "Missing node PHHepMCGenEvent" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    return NULL;
  }

  HepMC::GenEvent* theEvent = m_genevt->getEvent();

  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
  {
    if ((unsigned int)((*p)->barcode()) == track->get_id())
    {
      for (HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()->particles_begin(HepMC::parents);
           mother != (*p)->production_vertex()->particles_end(HepMC::parents); ++mother)
      {
        if (abs((*mother)->barcode()) == MOTHER_ID::D0)
        {
          std::cout << "This track came from the correct mother!!!" << std::endl;
	  double mass = abs((*p)->pdg_id()) == 321 ? 0.4937 : 0.139;
          lvParticle = new CLHEP::HepLorentzVector(track->get_px()
                                           ,track->get_py()
                                           ,track->get_pz()
                                           ,mass);
        }
        else
        {
          std::cout << "This track is background" << std::endl;
        }
      }
    }
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

string
QAG4SimulationKFParticle::get_histo_prefix()
{
  return string("h_") + Name() + string("_") + m_trackMapName + string("_");
}
