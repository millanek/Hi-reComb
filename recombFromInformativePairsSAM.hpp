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

class RecombReadPairsStats {
public:
    
    RecombReadPairsStats(): numConcordant(0), numDiscordant(0), totalEffectiveLength(0), meanRate(0.0) {
        windowSizeMins = {0,1000,2000,5000,10000,100000,1000000};
        windowSizeMax = {1000,2000,5000,10000,100000,1000000,1000000000};
        numRecombsInSizeWindows.resize(windowSizeMins.size(),0.0);
        numNonRecombsInSizeWindows.resize(windowSizeMins.size(),0.0);
        lengthOfInformativeSequenceWindows.resize(windowSizeMins.size(),0);
        windowRates.resize(windowSizeMins.size(),0.0);
    }
    
    void collectStats(const std::vector<DefiningRecombInfo*>& allInformativePairs, const bool bUptate = false) {
        
        if (bUptate) { resetVals();}
        
        for (int j = 0; j < allInformativePairs.size(); j++) {
            int l = allInformativePairs[j]->dist;
          //  std::cout << "l = " << l << std::endl;
          //  std::cout << "windowSizeMins.size() = " << windowSizeMins.size() << std::endl;
            for (int k = 0; k < windowSizeMins.size(); k++) {
          //      std::cout << "k = " << k << std::endl;
                if (l > windowSizeMins[k] && l <= windowSizeMax[k]) {
                    numRecombsInSizeWindows[k] += allInformativePairs[j]->probabilityRecombined;
                    numDiscordant += allInformativePairs[j]->probabilityRecombined;
                    numNonRecombsInSizeWindows[k] += (1 - allInformativePairs[j]->probabilityRecombined);
                    numConcordant += (1 - allInformativePairs[j]->probabilityRecombined);
                    lengthOfInformativeSequenceWindows[k] += l; totalEffectiveLength += l;
                }
            }
        }
        
        for (int j = 0; j < windowSizeMins.size(); j++) {
            double thisWindowRate = (double)numRecombsInSizeWindows[j]/lengthOfInformativeSequenceWindows[j];
            windowRates[j] = thisWindowRate;
        }
        
        meanLength = totalEffectiveLength / (double)allInformativePairs.size();
        meanRate = numDiscordant/totalEffectiveLength;
    };
    
    double numConcordant; double numDiscordant;
    long long int totalEffectiveLength; double meanLength;
    double meanRate;
    std::vector<double> windowRates;
    std::vector<int> windowSizeMins;
    
    void printRecombReadPairStats() {
        for (int j = 0; j != windowSizeMins.size(); j++) {
            std::cout << "window: " << windowSizeMins[j] << " - " << windowSizeMax[j] <<
                "; rate = " << windowRates[j] << "; n recomb = " << numRecombsInSizeWindows[j] << "; n non-recomb = " << numNonRecombsInSizeWindows[j] <<
                "; seqLength = " << lengthOfInformativeSequenceWindows[j] << std::endl;
        }
        std::cout << "total rate = " << meanRate << "; n recomb = " << numDiscordant << "; n non-recomb = " << numConcordant << "; seqLength = " << totalEffectiveLength << std::endl;
    }
    
    void printBasicStats() { // CURRENTLY NOT USED
        std::cout << "Effective coverage (bp): " << totalEffectiveLength << std::endl;
        std::cout << "numConcordant: " << numConcordant << std::endl;
        std::cout << "numDiscordant: " << numDiscordant << std::endl;
    }
    
    
private:
    // Size windows are 0 - 1000bp, 1001 - 2000 bp, 2000 - 5000bp, 5000 - 10000, 10000 - 100000, 100000 - 1000000, 1M+
    std::vector<int> windowSizeMax;
    std::vector<double> numRecombsInSizeWindows;
    std::vector<double> numNonRecombsInSizeWindows;
    std::vector<long long int> lengthOfInformativeSequenceWindows;
    
    void resetVals() {
        std::fill(numRecombsInSizeWindows.begin(), numRecombsInSizeWindows.end(), 0.0);
        std::fill(numNonRecombsInSizeWindows.begin(), numNonRecombsInSizeWindows.end(), 0.0);
        std::fill(windowRates.begin(), windowRates.end(), 0.0);
        std::fill(lengthOfInformativeSequenceWindows.begin(), lengthOfInformativeSequenceWindows.end(), 0);
        numDiscordant = 0; numConcordant = 0; totalEffectiveLength = 0;
    }
    
};

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
    std::vector<DefiningRecombInfo*> allInformativePairs;
    std::vector<DefiningRecombInfo*> bootstrappedInformativePairs;
    RecombReadPairsStats* stats = new RecombReadPairsStats();
    
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
    
    void printReadPairFileInfoIntoFile(const string& fileName, const bool switches) {
        std::ofstream* readPairFile = new std::ofstream(fileName);
        for (int i = 0; i != allInformativePairs.size(); i++) {
            if (allInformativePairs[i]->isRecombined == switches) {
                *readPairFile << allInformativePairs[i]->posLeft << "\t" << allInformativePairs[i]->posRight << "\t" << allInformativePairs[i]->dist << "\t" << allInformativePairs[i]->phaseErrorP_left << "\t" << allInformativePairs[i]->phaseErrorP_right << std::endl;
            }
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
    
    
    std::vector<DefiningRecombInfo*> getBootstrapSample(const std::vector<DefiningRecombInfo*> orginalInformativePairs) {
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); // seed the generator
        std::uniform_int_distribution<> distr(0, (int)orginalInformativePairs.size() - 1); // define the range
        
        std::vector<DefiningRecombInfo*> bootstrapSample(orginalInformativePairs.size());
        for (int i = 0; i < orginalInformativePairs.size(); i++) {
            int s = distr(gen);
            bootstrapSample[i] = orginalInformativePairs[s];
        }
        return bootstrapSample;
    }
    
    void considerDoubleCrossovers() {
        for (int j = 0; j < allInformativePairs.size(); j++) {
            double lambda = stats->meanRate * allInformativePairs[j]->dist;
            if(lambda > 3) {
                std::cerr << "WARNING: This looks very suspicious..." << std::endl;
                std::cerr << "stats->meanRate: " << stats->meanRate << std::endl;
                std::cerr << "allInformativePairs[j]->dist: " << allInformativePairs[j]->dist << std::endl;
                std::cerr << "lambda: " << lambda << std::endl;
                std::cerr << std::endl;
            }
                        
            double doubleCrossoverProbability;
            if (allInformativePairs[j]->dist > 1000000) {
                // Poisson Probability mass function using the mean rate
                // not worried for now about triple and more crossover
                doubleCrossoverProbability = (pow(lambda, 2.0) * exp(-lambda)) / 2.0;
            } else {
                // If distance is less than 1Mb, we set double crossover to zero because of crossover interference
                doubleCrossoverProbability = 0;
            }
            
            if (!allInformativePairs[j]->isRecombined) {
                // p(ph1=T) * p(ph2=T) * p(b1=A) * p(b2=C)
                // and a double cross-over probability:
                allInformativePairs[j]->probabilityRecombined += (1 - allInformativePairs[j]->phaseErrorP_left) * (1 - allInformativePairs[j]->phaseErrorP_right) * (1 - allInformativePairs[j]->baseErrorP_left) * (1 - allInformativePairs[j]->baseErrorP_right) * (doubleCrossoverProbability);
            }
        }
    }
    
    
    void adjustRecombinationProbabilities() {
       // std::cout << "meanL: " << meanL << std::endl;
        
        int maxLengthWindow = stats->windowSizeMins.back();
     //   double meanL_At1MbpPlus = getMeanLengthAboveThreshold(maxLengthWindow);
        
        //std::cout << "meanL: " << meanL_At1MbpPlus << std::endl;
        // std::cout << "maxLengthWindow: " << maxLengthWindow << std::endl;
        
        for (int j = 0; j < allInformativePairs.size(); j++) {
            if (allInformativePairs[j]->isRecombined) {
                if (allInformativePairs[j]->dist < maxLengthWindow) {
                    double scaleFactor = getScalingFactorFromInterpolation(allInformativePairs[j]->dist, (double)maxLengthWindow);
                    allInformativePairs[j]->probabilityRecombined = allInformativePairs[j]->probabilityRecombined * scaleFactor;
                    // allInformativePairs[j]->probabilityRecombined = allInformativePairs[j]->dist / stats->meanLength;
                }
            }
        }
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
    
    double getScalingFactorFromInterpolation(const int thisLength, const double meanL_At1MbpPlus) {
        double rateAt1000bp = stats->windowRates[1];
        double rateAt1MbpPlus = stats->windowRates.back();
        double rateRatio = rateAt1MbpPlus/rateAt1000bp;
        int minLength = stats->windowSizeMins[1];
        
        double scaleFactor = rateRatio + ((1-(rateRatio))/(meanL_At1MbpPlus-minLength)) * (thisLength - minLength); // Linear interpolation between 1Kb and 1Mb
        return scaleFactor;
    }
    
    double getMeanLengthAboveThreshold(const int threshold) {
        unsigned long long int sumLengths = 0; int nAboveThreshold = 0;
        for (int j = 0; j < allInformativePairs.size(); j++) {
            if (allInformativePairs[j]->dist >= threshold) {
                sumLengths += allInformativePairs[j]->dist; nAboveThreshold++;
            }
        }
        double meanLengthAboveThreshold = (double)sumLengths/(double)nAboveThreshold;
        return meanLengthAboveThreshold;
    }
    
};

class RecombInterval {
public:
    RecombInterval() {};
    
    RecombInterval(int leftIndex, double meanRate, RecombReadPairs* rp): directReadCoverageConcord(0), directReadCoverageDiscord(0) {
        j = leftIndex;
        recombFractionPerBp = meanRate;
        leftCoord = rp->coveredHetPos[j]; rightCoord = rp->coveredHetPos[j + 1];
        dj = rightCoord - leftCoord + 1;
        rj = recombFractionPerBp * dj;
        initialiseInterval(rp);
        effectiveCoverage = (int)coveringReadPairs.size();
        getDirectCoverageOnLeftCoord();
    };
    
    int j;
    int leftCoord;
    int rightCoord;
    double recombFractionPerBp;
    std::vector<double> bootstrapRecombFractions;
    int dj; double rj;
    
    std::vector<DefiningRecombInfo*> coveringReadPairs;
    int effectiveCoverage;
    int directReadCoverageConcord;
    int directReadCoverageDiscord;
    double sum_P_ij;
    double sumConcordantFraction;
    double sumRecombFractionPerBP;
    
    void updateVals(const double minCoverage) {
        if (coveringReadPairs.size() > minCoverage) {
            recombFractionPerBp = sumRecombFractionPerBP/sumConcordantFraction;
            rj = recombFractionPerBp * dj;
        }
    }
    
    void initialiseInterval(RecombReadPairs* rp) {
        sumRecombFractionPerBP = 0; sumConcordantFraction = 0; coveringReadPairs.clear();
        for (int i = 0; i != rp->allInformativePairs.size(); i++) {
            if(rp->allInformativePairs[i]->posLeft <= leftCoord && rp->allInformativePairs[i]->posRight >= rightCoord){
                coveringReadPairs.push_back(rp->allInformativePairs[i]);
                //  double recombFractionPerBpThisPair = (double)1.0/(double)rp->allInformativePairs[i]->dist;
                sumRecombFractionPerBP += (double)rp->allInformativePairs[i]->probabilityRecombined/(double)rp->allInformativePairs[i]->dist;
                sumConcordantFraction += 1 - rp->allInformativePairs[i]->probabilityRecombined;
            }
        }
    }

private:
    
    void getDirectCoverageOnLeftCoord() {
        for (int j = 0; j != coveringReadPairs.size(); j++) {
            if(coveringReadPairs[j]->posLeft == leftCoord || coveringReadPairs[j]->posRight == leftCoord){
                if (coveringReadPairs[j]->isRecombined) {
                    directReadCoverageDiscord++;
                } else {
                    directReadCoverageConcord++;
                }
            }
        }
    }
};


class RecombMap {
public:
    RecombMap(RecombReadPairs* rp, double minCoverageFraction) {
        // Get the mean recombination rate per bp
        meanRate = rp->stats->meanRate;
        std::cout << "meanRecombinationRate " << meanRate << std::endl;
        
        cummulativeRates.resize(rp->coveredHetPos.size() - 1);
        recombIntervals.resize(rp->coveredHetPos.size() - 1);
        int pc = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            // Initialise the recombination fractions assuming a uniform recombination rate along the genome (this uses the per-bp rate)
            RecombInterval ri(j, meanRate, rp); // Calculate the first iteration
            recombIntervals[j] = ri;
            int update5pcInterval = (int)recombIntervals.size() / 20;
            if (j > 0 && j % update5pcInterval == 0) { pc += 5; printPcUpdateOnPrompt(pc); }
            // if (j % 5000 == 0) printUpdateOnPrompt(j);
        }
        
        meanEffectiveCoverage = calculateMeanCoverage();
        
        std::cout << "recombIntervals.size() " << recombIntervals.size() << std::endl;
        
        for (int j = 0; j < recombIntervals.size(); j++) {
            recombIntervals[j].updateVals(minCoverageFraction * meanEffectiveCoverage);
        }
        mapLength = calculateCummulativeRates();
        std::cout << std::endl;
    }; 
    
    std::vector<RecombInterval> recombIntervals;
    double meanRate;
    std::vector<double> cummulativeRates; // Cummulative rates in cM
    double mapLength; // Map length in cM
    double meanEffectiveCoverage;
    
    void printPcUpdateOnPrompt(int pc) {
        if (pc == 100) std::cout << pc << "% (DONE)" << std::endl;
        else { std::cout << pc << "%..."; std::cout.flush(); }
    }
    
    double EMiteration(const int EMiterationNum, const double minCoverageFraction, bool bLoud = true) {
        if (bLoud) std::cout << "Iteration " << EMiterationNum << ": " << std::endl;
        updatePijs(bLoud);
        double delta = updateRecombFractions(minCoverageFraction * meanEffectiveCoverage, bLoud);
        mapLength = calculateCummulativeRates();
        if (bLoud) std::cout << std::endl;
        return delta;
    }
    
    void printUpdateOnPrompt(int j) {
        std::cout << "numProcessedHets: " << j << " ("<< (double)j/recombIntervals.size() << "%)"<< std::endl;
        std::cout << "pos: " << recombIntervals[j].leftCoord << "bp"<< std::endl;
    }
    
    void outputMapToFile(string fileName) {
        std::cout << "Writing recombination map into: " << fileName;
        std::ofstream* f = new std::ofstream(fileName);
        // That's just approximating the chromosome beginning; do we want it?
        *f << "1\t" << recombIntervals[0].leftCoord << "\t" << recombIntervals[0].recombFractionPerBp << std::endl;
        // Now the actual map
        for (int i = 0; i != recombIntervals.size(); i++) {
            *f << recombIntervals[i].leftCoord << "\t" << recombIntervals[i].rightCoord << "\t" << recombIntervals[i].recombFractionPerBp << std::endl;
        }
        std::cout << "...DONE!" << std::endl;
    }
    
    void outputBootstrapToFile(string fileName) {
        std::cout << "Writing bootstrapped recombination maps into:" << fileName;
        std::ofstream* f = new std::ofstream(fileName);
        // That's just approximating the chromosome beginning; do we want it?
        *f << "1\t" << recombIntervals[0].leftCoord << "\t"; print_vector(recombIntervals[0].bootstrapRecombFractions, *f);
        // Now the actual map
        for (int i = 0; i != recombIntervals.size(); i++) {
            *f << recombIntervals[i].leftCoord << "\t" << recombIntervals[i].rightCoord << "\t";
            print_vector(recombIntervals[i].bootstrapRecombFractions, *f);
        }
        std::cout << "...DONE!" << std::endl;
    }
    
    void outputMapToFileFixedWindowSizes(string fileName, int windowSize) {
        std::cout << "Writing recombination map with a fixed window size of " << windowSize << "bp into:" << fileName;
        std::ofstream* f = new std::ofstream(fileName);
        
        vector< vector<int> > allCoordsVectors; allCoordsVectors.resize(2);
        allCoordsVectors[0].resize(recombIntervals.size()); allCoordsVectors[1].resize(recombIntervals.size());
        vector<double> perBPrecombinationValueVector(recombIntervals.size(),0.0);
        for (int i = 0; i != recombIntervals.size(); i++) {
            allCoordsVectors[0][i] = recombIntervals[i].leftCoord;
            allCoordsVectors[1][i] = recombIntervals[i].rightCoord;
            perBPrecombinationValueVector[i] = recombIntervals[i].recombFractionPerBp;
        }
        int lastCoord = allCoordsVectors[1].back();
        int maxCoordToOutput = roundToNearestValue(lastCoord, windowSize);
        for (int start = 1; start < maxCoordToOutput; start = start + windowSize) {
            int physicalWindowEnd = start + windowSize - 1;
            *f << start << "\t" << physicalWindowEnd << "\t";
            *f << getAverageForPhysicalWindow(allCoordsVectors, perBPrecombinationValueVector, start, physicalWindowEnd);
            *f << std::endl;
        }
        std::cout << "...DONE!" << std::endl;
    }
    
    void printPerHetCoverageStats(string fn, RecombReadPairs* rp) {
        std::ofstream* depthFile = new std::ofstream(fn);
        *depthFile << "pos" << "\t" << "ratePerBp" << "\t" << "EffectiveDepth" << "\t" << "directReadCoverageConcord" << "\t" <<
                      "directReadCoverageDiscord" << std::endl;
        for (int i = 0; i != recombIntervals.size() - 1; i++) {
            int thisHetPos = recombIntervals[i].leftCoord;
            int coveredHetEffectiveDepth = (int)recombIntervals[i].coveringReadPairs.size();
            *depthFile << thisHetPos << "\t" << recombIntervals[i].recombFractionPerBp << "\t" << coveredHetEffectiveDepth << "\t" << recombIntervals[i].directReadCoverageConcord << "\t" << recombIntervals[i].directReadCoverageDiscord << std::endl;
        }
    }
    
private:
    void getSumPijForInterval(int intervalIndex) {
        RecombInterval* i = &recombIntervals[intervalIndex];
        
        i->sum_P_ij = 0;
        for (int r = 0; r < i->coveringReadPairs.size(); r++) {
            if (!i->coveringReadPairs[r]->isRecombined) continue;
            if (std::isnan(i->coveringReadPairs[r]->sum_r_k)) {
                i->coveringReadPairs[r]->sum_r_k = 0;
                for (int k = i->coveringReadPairs[r]->indexLeft; k < i->coveringReadPairs[r]->indexRight; k++) {
                    i->coveringReadPairs[r]->sum_r_k += recombIntervals[k].rj;
                }
            }
            double p_ij = (i->rj *  i->coveringReadPairs[r]->probabilityRecombined) / i->coveringReadPairs[r]->sum_r_k;  // eq. (2) from proposal
            i->sum_P_ij += p_ij;
        }
    }
    
    void updatePijs(bool bLoud = true) {
       // std::cout << "Updating p_ij values: " << std::endl;
        int pc = 0; int update5pcInterval = (int)recombIntervals.size() / 20;
        if (bLoud) std::cout << "Updating p_ij values: ";
        for (int j = 0; j < recombIntervals.size(); j++) {
            getSumPijForInterval(j);
            
            if (bLoud) if (j > 0 && j % update5pcInterval == 0) { pc += 5; printPcUpdateOnPrompt(pc); }
        }
        
        // Reset the pre-calculated sums before the next round
        for (int j = 0; j < recombIntervals.size(); j++) {
            RecombInterval* i = &recombIntervals[j];
            for (int r = 0; r < i->coveringReadPairs.size(); r++) {
                i->coveringReadPairs[r]->sum_r_k = NAN;
            }
        }
    }
    
    double updateRecombFractions(const double minCoverage, bool bLoud = true) {
        // std::cout << "Updating recombination fractions: " << std::endl;
        double delta = 0;
        int pc = 0; int update5pcInterval = (int)recombIntervals.size() / 20;
        if (bLoud) std::cout << "Updating recombination fractions: ";
        for (int j = 0; j < recombIntervals.size(); j++) {
            double newRj;
            if (bLoud) if (j > 0 && j % update5pcInterval == 0) { pc += 5; printPcUpdateOnPrompt(pc); }
            
            if (recombIntervals[j].effectiveCoverage > minCoverage) {
                newRj = recombIntervals[j].sum_P_ij / recombIntervals[j].sumConcordantFraction;
            } else if (recombIntervals[j].effectiveCoverage > 5) {
                newRj = weighedAverageRateUnderCoveringReads(recombIntervals[j].coveringReadPairs) * recombIntervals[j].dj;
               // std::cout << "newRj: " << newRj << std::endl;
               // std::cout << "newRj/ recombIntervals[j].dj: " << newRj/ recombIntervals[j].dj << std::endl;
            } else {
                newRj = meanRate * recombIntervals[j].dj;
            }
            delta += std::abs(newRj - recombIntervals[j].rj);
            recombIntervals[j].rj = newRj;
            recombIntervals[j].recombFractionPerBp = recombIntervals[j].rj / recombIntervals[j].dj;
        }
        if (bLoud) std::cout << "delta: " << delta << std::endl;
        return delta;
    }
    
    double weighedAverageRateUnderCoveringReads(std::vector<DefiningRecombInfo*>& coveringReadPairs) {
        std::vector<int> coverageCounts; coverageCounts.resize(recombIntervals.size(),0);
        for (int i = 0; i < coveringReadPairs.size(); i++) {
            for (int j = coveringReadPairs[i]->indexLeft; j < coveringReadPairs[i]->indexRight; j++) {
                coverageCounts[j]++;
            }
        }
        
        double sumRecomb = 0; double sumLengths = 0;
        for (int j = 0; j < recombIntervals.size();j++) {
            sumRecomb += recombIntervals[j].rj * coverageCounts[j];
            sumLengths += recombIntervals[j].dj * coverageCounts[j];
        }
        double weighedAvgRate = sumRecomb/sumLengths;
       // std::cout << "vector_sum(coverageCounts): " << vector_sum(coverageCounts) << std::endl;
       // std::cout << "weighedAvgRate: " << weighedAvgRate << std::endl;
        return weighedAvgRate;
    }
    
    // Returns total map length in cM
    double calculateCummulativeRates() {
        double cumSum = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            cummulativeRates[j] = cumSum * 100;
            cumSum += recombIntervals[j].rj;
        }
        return (cumSum * 100);
    }
    
    double calculateMeanCoverage() {
        int total = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            total += recombIntervals[j].effectiveCoverage;
        }
        double meanCoverage = (double)total/recombIntervals.size();
        return meanCoverage;
    }
    
    double getAverageForPhysicalWindow(const vector<vector<int>>& intervalCoordinates, const vector<double>& values, const int start, const int end) {
        // Binary search to find the first interval whose end coordinate is greater
        // or equal to the start of the region in question
        vector<int>::const_iterator itStart = lower_bound(intervalCoordinates[1].begin(),intervalCoordinates[1].end(),start);
        int numBPtotal = 0; int numBPthisInterval = 0;
        double sumPerBPvalue = 0; double meanValue = NAN;
        
        if (itStart != intervalCoordinates[1].end()) {  // if (start < f[1])    ---  excludind case 1)
            vector<int>::size_type index = std::distance(intervalCoordinates[1].begin(), itStart);
            // Sum the lengths
            while (intervalCoordinates[0][index] <= end && index < intervalCoordinates[0].size()) { // if (f[0] >= end)   ---    excluding case 2)
                double valueThisInterval = values[index];

              //  std::cerr << "valuesThisFeature[0]\t" << valuesThisFeature[0] << std::endl;
                if (intervalCoordinates[0][index] < start && intervalCoordinates[1][index] <= end)
                    numBPthisInterval = (intervalCoordinates[1][index] - start) + 1;
                else if (intervalCoordinates[0][index] >= start && intervalCoordinates[1][index] <= end)
                    numBPthisInterval = (intervalCoordinates[1][index] - intervalCoordinates[0][index]);
                else if (intervalCoordinates[0][index] >= start && intervalCoordinates[1][index] > end)
                    numBPthisInterval = (end - intervalCoordinates[0][index]);
                else if (intervalCoordinates[0][index] < start && intervalCoordinates[1][index] > end)
                    numBPthisInterval = (end - start) + 1;
                
                numBPtotal += numBPthisInterval;
                sumPerBPvalue += valueThisInterval * numBPthisInterval;
                index++;
            }
            
                meanValue = (double)sumPerBPvalue/numBPtotal;
            //std::cerr << "meanValue: " << meanValue << std::endl;
        }
        return meanValue;
        
    }
    
};

#endif /* recombFromInformativePairsSAM_hpp */
