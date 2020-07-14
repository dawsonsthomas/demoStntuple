//
#include "TInterpreter.h"
#include "demoStntuple/ana/scripts/modules.hh"
//-----------------------------------------------------------------------------
int load_stnana_scripts_demoStntuple() {
  char        macro[200];
  const char* script[] = { 
    //    "global_vars.cc",
//    "init_geometry.C",
//    "cosmics.C",
//    "dio_calib.C",
//    "ghist.C",
    "track.C",
    0 
  };

  const char* work_dir = gSystem->Getenv("MU2E_BASE_RELEASE");

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i++) {
    sprintf(macro,"%s/demoStntuple/ana/scripts/%s",work_dir,script[i]);
    if (! cint->IsLoaded(macro)) {
      cint->LoadMacro(macro);
    }
  }
  
  return 0;
}
