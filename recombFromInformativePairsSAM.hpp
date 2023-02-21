//
//  recombFromInformativePairsSAM.hpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 02.09.22.
//

#ifndef recombFromInformativePairsSAM_hpp
#define recombFromInformativePairsSAM_hpp

#include <stdio.h>
#include "recombUtils.hpp"

void parseRecombFromSAMOptions(int argc, char** argv);
int RecombFromSAMMain(int argc, char** argv);

void collectRateStatsBasedOnInsertLength(const std::vector<PhaseSwitch*>& phaseSwitches, const std::vector<std::vector<int>>& phaseConcordanceCoords);

#endif /* recombFromInformativePairsSAM_hpp */
