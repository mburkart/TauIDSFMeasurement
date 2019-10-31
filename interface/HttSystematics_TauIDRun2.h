#ifndef TauID2017_HttSystematics_SMRun2_h
#define TauID2017_HttSystematics_SMRun2_h
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
// Run2 SM analysis systematics
// Implemented in src/HttSystematics_SMRun2.cc
void AddTauIDRun2Systematics(CombineHarvester& cb, bool jetfakes, bool embedding, bool regional_jec, int era);
}

#endif
