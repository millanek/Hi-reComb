//
//  TrioPhase.hpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 23.08.23.
//

#ifndef TrioPhase_hpp
#define TrioPhase_hpp

#include "generalUtils.hpp"

int trioPhaseMain(int argc, char** argv);
void parseTrioPhaseOptions(int argc, char** argv);

inline void printPhasedLine(std::ostream& outFilePhasedHets, const int numPhasedHets, const VariantInfo& v, const char offspringAllele1, const char offspringAllele2) {
    outFilePhasedHets << numPhasedHets << "\t0\t1\t" << v.chr << "\t" << v.posInt << "\t" << offspringAllele1 << "\t" << offspringAllele2 << "\t0/1\t0\t.\t100.00\t100" << std::endl;
}

#endif /* TrioPhase_hpp */
