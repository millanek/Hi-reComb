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

void collectRateStatsBasedOnInsertLength(const std::vector<DefiningRecombInfo*>& phaseSwitches, const std::vector<DefiningRecombInfo*>& phaseConcordanceCoords);

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
    
    // Read Pairs
    std::vector<RecombReadPair*> readPairs;
    std::vector<RecombReadPair*> informativeReadPairs;
    
    // Base quality stats
    int numMatch = 0; int numMismatch = 0;
    std::vector<double> matchBaseScores; std::vector<double> mismatchBaseScores;
//    std::vector<double> concordantBaseScores; std::vector<double> discordantBaseScores;
    
    // Recombination stats
    int numConcordant = 0; int numDiscordant = 0;
    long long int totalEffectiveLength = 0;
    // Details about discordant read pairs to output by the printSwitchInfoIntoFile() method
    std::vector<DefiningRecombInfo*> phaseSwitches;
    std::vector<DefiningRecombInfo*> concordantPairs;
    std::vector<DefiningRecombInfo*> allInformativePairs;
    
    // Records about recombination-informative het sites
    std::vector<int> coveredHetPos;
    
    void linkWithHets(const std::map<int,PhaseInfo*>& posToPhase, const std::map <int, bool>& subsetLoci, const int minBQ) {
        for (int i = 0; i < readPairs.size(); i++) {
            RecombReadPair* thisReadPair = readPairs[i];
            thisReadPair->findAndCombinePairHets(posToPhase);
            thisReadPair->filterHetsByQuality(minBQ);
            
            // Check if het-bases in the reads match the genotype calls from the VCF/hapcut2 file
            categoriseBaseMatchMismatch(thisReadPair);
            
            // Find the number of informative heterozygous sites on this read-pair
            collectHetNumStats(thisReadPair);
            
            bool bContainsSubsetLoci = false; // Ignore read pairs with excluded SNPs
            
            for (int j = 0; j < thisReadPair->hetSites.size(); j++) {
                //std::cout << "thisReadPair->hetSites[j]->pos: " << thisReadPair->hetSites[j]->pos << std::endl;
                if (subsetLoci.count(thisReadPair->hetSites[j]->pos) == 1) {
                    bContainsSubsetLoci = true;
                  //  std::cout << "thisReadPair->hetSites[j]->pos: " << thisReadPair->hetSites[j]->pos << std::endl;
                }
            }
            
            // Check if there is a pair of hets on the two reads belonging to the same phase-block
            // If yes, the read-pair is then categorised as phase-informative
            if (thisReadPair->hetSites.size() > 1 && !bContainsSubsetLoci) {
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
    
    void printSwitchInfoIntoFile(string fileName) {
        std::ofstream* phaseSwitchFile = new std::ofstream(fileName);
        for (int i = 0; i != phaseSwitches.size(); i++) {
            *phaseSwitchFile << phaseSwitches[i]->posLeft << "\t" << phaseSwitches[i]->posRight << "\t" << phaseSwitches[i]->dist << "\t" << phaseSwitches[i]->phaseQualLeft << "\t" << phaseSwitches[i]->phaseQualRight << std::endl;
        }
    }
    
    
    void printBaseQualityStats() {
        std::cout << "Initial Read Pairs.size(): " << readPairs.size() << std::endl;
        std::cout << "numMatch: " << numMatch << std::endl;
        std::cout << "numMismatch: " << numMismatch << std::endl;
        std::cout << "Mean mismatchBaseScores: " << vector_average(mismatchBaseScores) << std::endl;
        std::cout << "Mean matchBaseScores: " << vector_average(matchBaseScores) << std::endl;
    }
    
    void printReadPairStats() {
        std::cout << "num0het: " << num0het << std::endl;
        std::cout << "num1het: " << num1het << std::endl;
        std::cout << "num2plusHets: " << num2plusHets << std::endl;
        std::cout << "informativeReadPairs.size(): " << informativeReadPairs.size() << std::endl;
        std::cout << std::endl;
    }
    
    void printConcordDiscordStats() {
        std::cout << "Effective coverage (bp): " << totalEffectiveLength << std::endl;
        std::cout << "numConcordant: " << numConcordant << std::endl;
        std::cout << "numDiscordant: " << numDiscordant << std::endl;
    }
    
    void findUniqueHetsCoveredByReadsAndSortThem() {
        std::sort(coveredHetPos.begin(), coveredHetPos.end());
        std::vector<int>::iterator it = std::unique(coveredHetPos.begin(), coveredHetPos.end());
        coveredHetPos.resize(distance(coveredHetPos.begin(),it));
        std::cout << "coveredHetPos.size() " << coveredHetPos.size() << std::endl;
        for (int i = 0; i != coveredHetPos.size(); i++) {
            hetPosToIndex[coveredHetPos[i]] = i;
        }
        
    }
    
    void findBoundingHetIndicesForEachReadPair() {
        for (int i = 0; i != allInformativePairs.size(); i++) {
            allInformativePairs[i]->indexLeft = hetPosToIndex.at(allInformativePairs[i]->posLeft);
            allInformativePairs[i]->indexRight = hetPosToIndex.at(allInformativePairs[i]->posRight);
        }
    }
    
    
    std::vector<DefiningRecombInfo*> getBootstrapSample() {
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<> distr(0, (int)allInformativePairs.size() - 1); // define the range
        
        std::vector<DefiningRecombInfo*> bootstrapSample(allInformativePairs.size());
        for (int i = 0; i < allInformativePairs.size(); i++) {
            int s = distr(gen);
            bootstrapSample[i] = allInformativePairs[s];
        }
        return bootstrapSample;
    }
    
    
private:
    // linking stats
    int num0het = 0; int num1het = 0; int num2plusHets = 0;
    int totalUsedReadLengthBp = 0;
    std::map<int,int> hetPosToIndex;
    
    void categoriseBaseMatchMismatch(RecombReadPair* thisReadPair) {
        std::vector<HetInfo*>::iterator it = thisReadPair->hetSites.begin();
        while(it != thisReadPair->hetSites.end()) {
            HetInfo* thisHet = *it;
            if (thisHet->readPhaseBaseMismatch) {
                numMismatch++;
                mismatchBaseScores.push_back(thisHet->thisBaseQuality);
                it = thisReadPair->hetSites.erase(it);
            } else {
                // std::cout << "Fine: " << std::endl;
                matchBaseScores.push_back(thisHet->thisBaseQuality);
                numMatch++;
                coveredHetPos.push_back(thisHet->pos);
                ++it;
            }
        }
    }
    
    void collectHetNumStats(RecombReadPair* thisReadPair) {
        if (thisReadPair->hetSites.size() == 0) num0het++;
        else if (thisReadPair->hetSites.size() == 1) num1het++;
        else num2plusHets++;
    }
    
};

class RecombInterval {
public:
    RecombInterval() {};
    
    RecombInterval(int leftIndex, double meanRate, int lc, int rc): sumRecombFractionPerBP(0), sumConcordantFraction(0) {
        j = leftIndex;
        recombFractionPerBp = meanRate;
        
        leftCoord = lc; rightCoord = rc;
        dj = rightCoord - leftCoord + 1;
        rj = recombFractionPerBp * dj;
    };
    
    int j;
    int leftCoord;
    int rightCoord;
    double recombFractionPerBp;
    int dj; double rj;
    
    std::vector<DefiningRecombInfo*> coveringReadPairs;
    double sum_P_ij;
    double sumConcordantFraction;
    double sumRecombFractionPerBP;

    
    void initialiseInterval(RecombReadPairs* rp) {
        for (int i = 0; i != rp->allInformativePairs.size(); i++) {
            if(rp->allInformativePairs[i]->posLeft <= leftCoord && rp->allInformativePairs[i]->posRight >= rightCoord){
                coveringReadPairs.push_back(rp->allInformativePairs[i]);
                if (rp->allInformativePairs[i]->isRecombined) {
                    double recombFractionPerBpThisPair = (double)1.0/(double)rp->allInformativePairs[i]->dist;
                    sumRecombFractionPerBP += recombFractionPerBpThisPair;
                } else {
                    sumConcordantFraction++;
                }
            }
        }
        
        if (coveringReadPairs.size() > 10) recombFractionPerBp = sumRecombFractionPerBP/sumConcordantFraction;
        rj = recombFractionPerBp * dj;
    }
    
};


class RecombMap {
public:
    RecombMap(RecombReadPairs* rp) {
        // Get the mean recombination rate per bp
        meanRate = (double)rp->numDiscordant/(double)rp->totalEffectiveLength;
        std::cout << "meanRecombinationRate " << meanRate << std::endl;
        
        
        recombIntervals.resize(rp->coveredHetPos.size() - 1);
        for (int j = 0; j < recombIntervals.size(); j++) {
            // Initialise the recombination fractions assuming a uniform recombination rate along the genome (this uses the per-bp rate)
            RecombInterval ri(j, meanRate, rp->coveredHetPos[j], rp->coveredHetPos[j + 1]);
            ri.initialiseInterval(rp); // Then calculate the first iteration
            recombIntervals[j] = ri;
            if (j % 5000 == 0) printUpdateOnPrompt(j);
        }
        std::cout << std::endl;
    }; 
    
    double meanRate;
    std::vector<RecombInterval> recombIntervals;
    
    
    void printPcUpdateOnPrompt(int pc) {
        if (pc == 100) std::cout << pc << "% (DONE)" << std::endl;
        else { std::cout << pc << "%..."; std::cout.flush(); }
    }
    
    double EMiteration(int EMiterationNum) {
        std::cout << "Iteration " << EMiterationNum << ": " << std::endl;
        updatePijs();
        double delta = updateRecombFractions();
        std::cout << std::endl;
        return delta;
    }
    
    void printUpdateOnPrompt(int j) {
        std::cout << "numProcessedHets: " << j << " ("<< (double)j/recombIntervals.size() << "%)"<< std::endl;
        std::cout << "pos: " << recombIntervals[j].leftCoord << "bp"<< std::endl;
    }
    
    void outputMapToFile(string fileName) {
        std::ofstream* f = new std::ofstream(fileName);
        // That's just approximating the chromosome beginning; do we want it?
        *f << "0\t" << recombIntervals[0].leftCoord << "\t" << recombIntervals[0].recombFractionPerBp << std::endl;
        // Now the actual map
        for (int i = 0; i != recombIntervals.size(); i++) {
            *f << recombIntervals[i].leftCoord << "\t" << recombIntervals[i].rightCoord << "\t" << recombIntervals[i].recombFractionPerBp << std::endl;
        }
    }
    
    void calculateAndPrintPerHetCoverageStats(string fn, RecombReadPairs* rp) {
        std::ofstream* depthFile = new std::ofstream(fn + ".txt");
        for (int i = 0; i != recombIntervals.size() - 1; i++) {
            int coverageDirectlyOnTheHet = 0;
            int thisHetPos = recombIntervals[i].leftCoord;
            for (int j = 0; j != rp->allInformativePairs.size(); j++) {
                if(rp->allInformativePairs[j]->posLeft == thisHetPos || rp->allInformativePairs[j]->posRight == thisHetPos){
                    coverageDirectlyOnTheHet++;
                }
            }
            int coveredHetEffectiveDepth = (int)recombIntervals[i].coveringReadPairs.size();
            *depthFile << thisHetPos << "\t" << recombIntervals[i].recombFractionPerBp << "\t" << coveredHetEffectiveDepth << "\t" << coverageDirectlyOnTheHet << std::endl;
        }
    }
    
private:
    void getSumPijForInterval(int intervalIndex) {
        RecombInterval* i = &recombIntervals[intervalIndex];
        
       /* if (intervalIndex == 40) {
            std::cout << "i.coveringReadPairs.size(): " << i.coveringReadPairs.size() << std::endl;
        } */
        i->sum_P_ij = 0;
        for (int r = 0; r < i->coveringReadPairs.size(); r++) {
            if (!i->coveringReadPairs[r]->isRecombined) continue;
            double sum_r_k = 0;
            /* if (intervalIndex == 40) {
                std::cout << "i.coveringReadPairs[r]->indexLeft: " << i.coveringReadPairs[r]->indexLeft << std::endl;
                std::cout << "i.coveringReadPairs[r]->indexRight: " << i.coveringReadPairs[r]->indexRight << std::endl;
            } */
            for (int k = i->coveringReadPairs[r]->indexLeft; k < i->coveringReadPairs[r]->indexRight; k++) {
                sum_r_k += recombIntervals[k].rj;
            //    if (r == 5) {
            //    std::cout << "k: " << k << std::endl;
            //    std::cout << "recombIntervals[k].rj: " << recombIntervals[k].rj << std::endl;
            //    }
            }
            double p_ij = i->rj/sum_r_k;  // eq. (2) from proposal
             /* if (intervalIndex == 40) {
                std::cout << "r: " << r << std::endl;
                std::cout << "i.rj: " << i.rj << std::endl;
                std::cout << "sum_r_k: " << sum_r_k << std::endl;
                std::cout << "p_ij: " << p_ij << std::endl;
                std::cout << std::endl;
            } */
            i->sum_P_ij += p_ij;
        }
        
            /*if (intervalIndex == 40) {
                std::cout << "intervalIndex: " << intervalIndex << std::endl;
            std::cout << "i.sum_P_ij: " << i->sum_P_ij << std::endl;
        } */
    }
    
    void updatePijs() {
       // std::cout << "Updating p_ij values: " << std::endl;
        int pc = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            getSumPijForInterval(j);
            
            int update5pcInterval = (int)recombIntervals.size() / 20;
            if (j > 0 && j % update5pcInterval == 0) { pc += 5; printPcUpdateOnPrompt(pc); }
        }
    }
    
    double updateRecombFractions() {
        // std::cout << "Updating recombination fractions: " << std::endl;
        double delta = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            double newRj = recombIntervals[j].sum_P_ij / recombIntervals[j].sumConcordantFraction;
            if (!isnan(newRj)) delta += abs(newRj - recombIntervals[j].rj);
            recombIntervals[j].rj = newRj;
            recombIntervals[j].recombFractionPerBp = recombIntervals[j].rj / recombIntervals[j].dj;
        }
        std::cout << "delta: " << delta << std::endl;
        return delta;
    }
    
};

#endif /* recombFromInformativePairsSAM_hpp */
