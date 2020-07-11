///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "heatherh/ana/scripts/modules.hh"
#include "Stntuple/scripts/global_vars.h"

def_name trigger_001("trigger_ana001");
//def_name track_ce ("track_ana001");
//def_name track_dio("track_ana001");
///////////////////////////////////////////////////////////////////////////////



//-----------------------------------------------------------------------------
void  trigger_ana001(int IsBkg=0, int GenCode=28, double NEvents=-1) {
//-----------------------------------------------------------------------------
// configure analysis module
//-----------------------------------------------------------------------------
  m_trig = (TTriggerAna001Module*) g.x->AddModule("TTriggerAna001Module",0);  
  //  m_trig->SetDebugBit(54, 1);
  m_trig->SetNormalization(NEvents);
  m_trig->SetBackgroundFlag(IsBkg);
  m_trig->SetGeneratorCode(GenCode);

}
