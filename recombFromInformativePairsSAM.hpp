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

class RecombReadPairs {
    public:
    // Parse the samtools file to find reads that match records from the pairstools file
    // and therefore  can be informative about the phasing and recombination
    RecombReadPairs(string readFileName) {
        std::ifstream* samFile = new std::ifstream(readFileName.c_str()); assertFileOpen(*samFile, readFileName);
        string line; int readN = 0;
        std::vector<RecombRead*> informativeReads;
        while (getline(*samFile,line)) {
            // std::cerr << line << std::endl;
            readN++;
            std::vector<string> samRecVec = split(line, '\t'); //assert(pairVec.size() == 8);
            RecombRead* thisRead = new RecombRead(samRecVec);
            informativeReads.push_back(thisRead);
            if (readN % 2 == 0) {
                RecombReadPair* thisReadPair = new RecombReadPair(informativeReads[0], informativeReads[1]);
                readPairs.push_back(thisReadPair);
                informativeReads.clear();
            }
        }
    };
    
    std::vector<RecombReadPair*> readPairs;
    std::vector<RecombReadPair*> informativeReadPairs;
    
    // linking stats
    int num0het = 0; int num1het = 0; int num2plusHets = 0;
    int totalUsedReadLengthBp = 0;
    
    // Base quality stats
    int numMatch = 0; int numMismatch = 0;
    std::vector<double> matchBaseScores; std::vector<double> mismatchBaseScores;
    
    // Recombination stats
    int numConcordant = 0; int numDiscordant = 0;
    long long int totalEffectiveLength = 0;
    
    
    void linkWithHets(std::map<int,PhaseInfo*>& posToPhase, int minBQ) {
        for (int i = 0; i < readPairs.size(); i++) {
            RecombReadPair* thisReadPair = readPairs[i];
            thisReadPair->findAndCombinePairHets(posToPhase);
            thisReadPair->filterHetsByQuality(minBQ);
            
            // Collecting stats
            if (thisReadPair->hetSites.size() == 0) num0het++;
            else if (thisReadPair->hetSites.size() == 1) num1het++;
            else num2plusHets++;
            
            // Check if pairs of hets are in the same phase-block
            if (thisReadPair->hetSites.size() > 1) {
                thisReadPair->read1->linkHetsWithPhaseBlock();
                thisReadPair->read2->linkHetsWithPhaseBlock();
                for (std::map<int, std::vector<int>>::iterator it = thisReadPair->read1->BlockIDsToHetPos.begin();
                    it != thisReadPair->read1->BlockIDsToHetPos.end(); it++) {
                    if (thisReadPair->read2->BlockIDsToHetPos.count(it->first) == 1) {
                        informativeReadPairs.push_back(thisReadPair);
                        totalUsedReadLengthBp += thisReadPair->read1->usedLength;
                        totalUsedReadLengthBp += thisReadPair->read2->usedLength;
                    }
                }
            }
            
        }
    }
    
    void printReadPairStats() {
        std::cout << "Initial Read Pairs.size(): " << readPairs.size() << std::endl;
        std::cout << "num0het: " << num0het << std::endl;
        std::cout << "num1het: " << num1het << std::endl;
        std::cout << "num2plusHets: " << num2plusHets << std::endl;
        std::cout << "informativeReadPairs.size(): " << informativeReadPairs.size() << std::endl;
        std::cout << std::endl;
    }
    
    void printBaseQualityStats() {
        std::cout << "numMatch: " << numMatch << std::endl;
        std::cout << "numMismatch: " << numMismatch << std::endl;
        std::cout << "Mean mismatchBaseScores: " << vector_average(mismatchBaseScores) << std::endl;
        std::cout << "Mean matchBaseScores: " << vector_average(matchBaseScores) << std::endl;
    }
    
    void printRecombinationBaseStats() {
        std::cout << "Effective coverage (bp): " << totalEffectiveLength << std::endl;
        std::cout << "numConcordant: " << numConcordant << std::endl;
        std::cout << "numDiscordant: " << numDiscordant << std::endl;
    }
    
};


#endif /* recombFromInformativePairsSAM_hpp */
