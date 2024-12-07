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
    
    RecombReadPairsStats(const int minDistanceToConsider): numConcordant(0), numDiscordant(0), totalEffectiveLength(0), meanRate(0.0) {
        windowSizeMins = {minDistanceToConsider,1000,5000,10000,100000,1000000};
        windowSizeMax = {1000,5000,10000,100000,1000000,1000000000};
        numRecombsInSizeWindows.resize(windowSizeMins.size(),0.0);
        numNonRecombsInSizeWindows.resize(windowSizeMins.size(),0.0);
        lengthOfInformativeSequenceWindows.resize(windowSizeMins.size(),0);
        windowRates.resize(windowSizeMins.size(),0.0);
    }
    
    void collectAndPrintStats(const string& msg, const std::vector<DefiningRecombInfo*>& informativeReadPairs, const bool bReset, bool bLoud = true) {
        collectStats(informativeReadPairs, bReset);
        if (bLoud) std::cout << msg << std::endl;
        if (bLoud) printRecombReadPairStats();
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
                    numNonRecombsInSizeWindows[k] += (1.0 - allInformativePairs[j]->probabilityRecombined);
                    numConcordant += (1.0 - allInformativePairs[j]->probabilityRecombined);
                    lengthOfInformativeSequenceWindows[k] += l; totalEffectiveLength += l;
                 //   if (k == 3 || k == 4) {
                  //      std::cout << "read pair j = " << j << "; probabilityRecombined = " << allInformativePairs[j]->probabilityRecombined << std::endl;
                 //   }
                    
                    
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
            std::cout << "Size window: " << windowSizeMins[j] << " - " << windowSizeMax[j] <<
                "bp; rate = " << windowRates[j] << "; n recomb = " << numRecombsInSizeWindows[j] << "; n non-recomb = " << numNonRecombsInSizeWindows[j] <<
                "; seqLength = " << lengthOfInformativeSequenceWindows[j] << std::endl;
        }
        std::cout << "total rate = " << meanRate << "; n recomb = " << numDiscordant << "; n non-recomb = " << numConcordant << "; seqLength = " << totalEffectiveLength << std::endl;
        std::cout << std::endl;
    }
    
    void printBasicStats() { // CURRENTLY NOT USED
        std::cout << "Effective coverage (bp): " << totalEffectiveLength << std::endl;
        std::cout << "numConcordant: " << numConcordant << std::endl;
        std::cout << "numDiscordant: " << numDiscordant << std::endl;
    }
    
    
    std::vector<int> windowSizeMax;
    std::vector<double> numRecombsInSizeWindows;
    std::vector<double> numNonRecombsInSizeWindows;
    std::vector<long long int> lengthOfInformativeSequenceWindows;

private:
    void resetVals() {
        std::fill(numRecombsInSizeWindows.begin(), numRecombsInSizeWindows.end(), 0.0);
        std::fill(numNonRecombsInSizeWindows.begin(), numNonRecombsInSizeWindows.end(), 0.0);
        std::fill(windowRates.begin(), windowRates.end(), 0.0);
        std::fill(lengthOfInformativeSequenceWindows.begin(), lengthOfInformativeSequenceWindows.end(), 0);
        numDiscordant = 0; numConcordant = 0; totalEffectiveLength = 0;
    }
    
};

class CoveredHetInfo {
public:
    CoveredHetInfo(): coverageDiscord(0), coverageConcord(0) {}
    
    int coverageDiscord;
    int coverageConcord;
};


class RecombReadPairs {
public:
    // Parse the samtools file to find reads that match records from the pairstools file
    // and therefore  can be informative about the phasing and recombination
    RecombReadPairs(string readFileName, const int minDistanceToConsider) {
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
        stats = new RecombReadPairsStats(minDistanceToConsider);
    };
    
    // Empty base constructor
    RecombReadPairs() {
        stats = new RecombReadPairsStats(200); // 200 = minDistanceToConsider
    };
    
    // Read Pairs
    std::vector<RecombReadPair*> readPairs;
    int numDiscordant = 0; int numConcordant = 0; int numAmbiguous = 0; int numTooShort = 0;
    
    std::vector<RecombReadPair*> informativeReadPairs;
    
    // Base quality stats
    int numMatch = 0; int numMismatch = 0;
    std::vector<double> matchBaseScores; std::vector<double> mismatchBaseScores;
//    std::vector<double> concordantBaseScores; std::vector<double> discordantBaseScores;
    
    // Recombination stats
    std::vector<DefiningRecombInfo*> allInformativePairs;
    std::vector<DefiningRecombInfo*> bootstrappedInformativePairs;
    RecombReadPairsStats* stats;
    
    // Records about recombination-informative het sites
    std::map<int, CoveredHetInfo*> coveredHets;
    std::vector<int> coveredHetPos;
    std::vector<int> coveragePerHetDiscord;
    std::vector<int> coveragePerHetConcord;
    
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
    
    void categorisePairs(const int minDistanceToDefinePairs) {
        DefiningRecombInfo* thisPairInformation;
        for (std::vector<RecombReadPair*>::iterator it = informativeReadPairs.begin(); it != informativeReadPairs.end(); it++) {
            RecombReadPair* thisReadPair = *it;
            
            thisReadPair->findIndicesOfConcordantAndDiscordantPairsOfHets(minDistanceToDefinePairs);
            thisReadPair->determineIfReadPairConcordantOrDiscordant();
            
            if (thisReadPair->pairRecombinationStatus == PAIR_DISCORDANT) {
                thisPairInformation = thisReadPair->getDefiningHetPair(thisReadPair->switchPairI, thisReadPair->switchPairJ);
                numDiscordant++;
            } else if (thisReadPair->pairRecombinationStatus == PAIR_CONCORDANT) {
                thisPairInformation = thisReadPair->getDefiningHetPair(thisReadPair->concordPairI, thisReadPair->concordPairJ);
                numConcordant++;
            } else if (thisReadPair->pairRecombinationStatus == PAIR_TOO_SHORT) {
                numTooShort++;
            } else {
                numAmbiguous++;
                
                /*std::cout << "Ambiguous read pair: " << std::endl;
                for (int h = 0; h < thisReadPair->hetSites.size(); h++) {
                    std::cout << "thisReadPair->hetSites[h]->thisHetPhase01: " << thisReadPair->hetSites[h]->thisHetPhase01 << std::endl;
                    std::cout << "thisReadPair->hetSites[h]->pos: " << thisReadPair->hetSites[h]->pos << std::endl;
                }
                std::cout << "thisReadPair->concordPairI.size(): " << thisReadPair->concordPairI.size() << std::endl;
                std::cout << "thisReadPair->switchPairI.size(): " << thisReadPair->switchPairI.size() << std::endl;
                exit(1);
                 */
            }
            if (thisReadPair->pairRecombinationStatus != PAIR_AMBIGUOUS && thisReadPair->pairRecombinationStatus != PAIR_TOO_SHORT) {
                allInformativePairs.push_back(thisPairInformation);
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
        std::cout << "Initial number of read pairs = " << readPairs.size() << std::endl;
        std::cout << std::endl;
        std::cout << "BASE QUALITY METRICS:" << std::endl;
        std::cout << "numMatch: " << numMatch << "; mean base score = " << vector_average(matchBaseScores) << std::endl;
        std::cout << "numMismatch: " << numMismatch << "; mean base score = " << vector_average(mismatchBaseScores) << std::endl;
        std::cout << std::endl;
    }
    
    void printReadPairStatsByHetNumbers() {
        std::cout << "HET SITE LINKING METRICS:" << std::endl;
        std::cout << "Number of read pairs with 2+ valid hets = " << num2plusHets;
        std::cout << "; (0 het = " << num0het << ", 1 het = " << num1het << ")" << std::endl;
        std::cout << "Number of read pairs with at least one informative het (in the same phase block) on each fragment = " << informativeReadPairs.size() << std::endl;
        std::cout << std::endl;
    }
    
    void printReadPairStatsByCategory(const int minD) {
        std::cout << "READ PAIR CATEGORY METRICS:" << std::endl;
        std::cout << "Read pairs that are too short to use (below " << minD << "bp) = " << numTooShort << std::endl;
        std::cout << "Read pairs that have multiple hets with ambiguous phasing  = " << numAmbiguous << std::endl;
        std::cout << "Read pairs to be directly used in map estimation: " << std::endl;
        std::cout << "Discordant (apparent recombined) read pairs = " << numDiscordant << std::endl;
        std::cout << "Concordant read pairs = " << numConcordant << std::endl;
        std::cout << std::endl;
    }
    
    
    void findUniqueHetsCoveredByReadsAndSortThem() {
        std::sort(coveredHetPos.begin(), coveredHetPos.end());
        std::vector<int>::iterator it = std::unique(coveredHetPos.begin(), coveredHetPos.end());
        coveredHetPos.resize(distance(coveredHetPos.begin(),it));
        for (int i = 0; i != coveredHetPos.size(); i++) {
            coveredHets[coveredHetPos[i]] = new CoveredHetInfo();
        }
        std::cout << "Number of heterozygous sites covered by all informative read pairs = " << coveredHetPos.size() << std::endl;
    }
    
    void findWhichHetsAreUsedByFinalReadPairSet() {
        for (int i = 0; i != coveredHetPos.size(); i++) {
            hetPosToIndex[coveredHetPos[i]] = i;
        }
        std::vector<int> usedIndices;
        for (int i = 0; i != allInformativePairs.size(); i++) {
            usedIndices.push_back(hetPosToIndex.at(allInformativePairs[i]->posLeft));
            usedIndices.push_back(hetPosToIndex.at(allInformativePairs[i]->posRight));
        }
        std::sort(usedIndices.begin(), usedIndices.end());
        std::vector<int>::iterator it = std::unique(usedIndices.begin(), usedIndices.end());
        usedIndices.resize(distance(usedIndices.begin(),it));
        std::cout << "Number of heterozygous sites covered by used informative read pairs = " << usedIndices.size() << std::endl;
        
        std::vector<int> usedCoveredHetPos;
        for (const int& index : usedIndices) {
            if (index >= 0 && index < coveredHetPos.size()) {
                usedCoveredHetPos.push_back(coveredHetPos[index]);
            }
        }
        usedCoveredHetPos.erase( std::unique( usedCoveredHetPos.begin(), usedCoveredHetPos.end() ), usedCoveredHetPos.end() );
        
        coveredHetPos = usedCoveredHetPos;
    }
    
    
    void updateBoundingHetIndicesForEachReadPair() {
        hetPosToIndex.clear();
        for (int i = 0; i != coveredHetPos.size(); i++) {
            hetPosToIndex[coveredHetPos[i]] = i;
        }
        for (int i = 0; i != allInformativePairs.size(); i++) {
            allInformativePairs[i]->indexLeft = hetPosToIndex.at(allInformativePairs[i]->posLeft);
            allInformativePairs[i]->indexRight = hetPosToIndex.at(allInformativePairs[i]->posRight);
        }
    }
    
    void calculateDirectCoverageOnEachHet() {
        coveragePerHetDiscord.resize(coveredHetPos.size());
        coveragePerHetConcord.resize(coveredHetPos.size());
        for (int i = 0; i != coveredHetPos.size(); i++) {
            for (int j = 0; j != allInformativePairs.size(); j++) {
                if(allInformativePairs[j]->posLeft == coveredHetPos[i] || allInformativePairs[j]->posRight == coveredHetPos[i]){
                    if (allInformativePairs[j]->isRecombined) {
                        coveragePerHetDiscord[i]++;
                    } else {
                        coveragePerHetConcord[i]++;
                    }
                }
            }
        }
    }
    
    void calculateDirectCoverageOnEachHetMap() {
        for (int j = 0; j != allInformativePairs.size(); j++) {
            int posLeft = allInformativePairs[j]->posLeft;
            int posRight = allInformativePairs[j]->posRight;
            if (allInformativePairs[j]->isRecombined) {
                coveredHets.at(posLeft)->coverageDiscord++;
                coveredHets.at(posRight)->coverageDiscord++;
            } else {
                coveredHets.at(posLeft)->coverageConcord++;
                coveredHets.at(posRight)->coverageConcord++;
            }
        }
    }
    
    void findAndRemoveReadPairsCoveringMultiHets(const string& rn) {
        
        std::vector<int> problematicSNPs;
    /*    for (int i = 0; i != coveredHetPos.size(); i++) {
            if (coveragePerHetDiscord[i] >= 2 && coveragePerHetConcord[i] == 0) problematicSNPs.push_back(coveredHetPos[i]);
        } */
        
        for (int i = 0; i != coveredHetPos.size(); i++) {
            if (coveredHets.at(coveredHetPos[i])->coverageDiscord >= 2 && coveredHets.at(coveredHetPos[i])->coverageConcord == 0) problematicSNPs.push_back(coveredHetPos[i]);
        }
        
        
       /* std::ofstream* depthFile1 = new std::ofstream("recombMap_wDepth_rp" + rn + ".txt");
        *depthFile1 << "pos" << "\t" << "directReadCoverageConcord" << "\t" << "directReadCoverageDiscord" << std::endl;
        for (int i = 0; i != coveredHetPos.size(); i++) {
            *depthFile1 << coveredHetPos[i] << "\t" << coveragePerHetConcord[i] << "\t" << coveragePerHetDiscord[i] << std::endl;
        }
        */
        
        std::vector<int> readPairsToRemove;
        for (int i = 0; i != problematicSNPs.size(); i++) {
            for (int j = 0; j != allInformativePairs.size(); j++) {
                if((allInformativePairs[j]->posLeft == problematicSNPs[i] || allInformativePairs[j]->posRight == problematicSNPs[i])){
                    readPairsToRemove.push_back(i);
                }
            }
        }
        
        std::cout << "Removing " << readPairsToRemove.size() << " read pairs directly covering " << problematicSNPs.size() << " problematic SNPs " << std::endl;
        
        // Sort the indices in descending order so that removing elements doesn't affect the subsequent indices
        std::sort(readPairsToRemove.begin(), readPairsToRemove.end(), std::greater<int>());

        // Remove elements at specified indices
        for (const int& index : readPairsToRemove) {
            if (index >= 0 && index < allInformativePairs.size()) {
                allInformativePairs.erase(allInformativePairs.begin() + index);
            }
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
            if (!allInformativePairs[j]->isRecombined) { // If the read isn't recombined
                double lambda = stats->meanRate * allInformativePairs[j]->dist;
                if(lambda > 3) {
                    std::cerr << "WARNING: This looks very suspicious..." << std::endl;
                    std::cerr << "stats->meanRate: " << stats->meanRate << std::endl;
                    std::cerr << "allInformativePairs[j]->dist: " << allInformativePairs[j]->dist << std::endl;
                    std::cerr << "lambda: " << lambda << std::endl;
                    std::cerr << std::endl;
                }
                
                
                // If distance is less than 1Mb, we set double crossover to zero because of crossover interference
                if (allInformativePairs[j]->dist > 1000000) {
                    // Poisson Probability mass function using the mean rate
                    // not worried for now about triple and more crossover
                    allInformativePairs[j]->doubleCrossoverProbability = (pow(lambda, 2.0) * exp(-lambda)) / 2.0;
                }
                
                // p(ph1=T) * p(ph2=T) * p(b1=A) * p(b2=C)
                // and a double cross-over probability:
                allInformativePairs[j]->probabilityRecombined += (1 - allInformativePairs[j]->phaseErrorP_left) * (1 - allInformativePairs[j]->phaseErrorP_right) * (1 - allInformativePairs[j]->baseErrorP_left) * (1 - allInformativePairs[j]->baseErrorP_right) * allInformativePairs[j]->doubleCrossoverProbability;
            }
        
        }
    }
    
    
    void adjustRecombinationProbabilities() {
        int maxLengthWindow = stats->windowSizeMins.back();
        for (int j = 0; j < allInformativePairs.size(); j++) {
           // if (allInformativePairs[j]->isRecombined) {
                if (allInformativePairs[j]->dist < maxLengthWindow) {
                    double scaleFactor = getScalingFactorFromInterpolation(allInformativePairs[j]->dist, (double)maxLengthWindow);
                    if (scaleFactor > 1) continue;
                    allInformativePairs[j]->probabilityRecombined = allInformativePairs[j]->probabilityRecombined * scaleFactor;
                    // allInformativePairs[j]->probabilityRecombined = allInformativePairs[j]->dist / stats->meanLength;
                }
                checkProbabilityValid(j, allInformativePairs[j]);
          // }
        }
    }
    
    void adjustRecombinationProbabilitiesBayes() {
       // std::cout << ": " << meanL << std::endl;
        
        // CALCULATE FALSE POSITIVE RATE FROM RECOMBINATION RATE DIFFERENCES IN DIFFERENT READ LENGTH BRACKETS
        double longreadRate = stats->windowRates.back();
        double shortreadLength = stats->lengthOfInformativeSequenceWindows.front();
        double shortreadRecombs = stats->numRecombsInSizeWindows.front();
        double shortreadNonRecombs = stats->numNonRecombsInSizeWindows.front();
        
        double expectedShortreadRecombs = longreadRate * shortreadLength;
        
        double falsePositiveRate = (shortreadRecombs - expectedShortreadRecombs) / (shortreadRecombs + shortreadNonRecombs);
        if (falsePositiveRate < 0) falsePositiveRate = 0;
        std::cout << "falsePositiveRate estimate " << falsePositiveRate << std::endl;
        std::cout << "(shortreadRecombs + shortreadNonRecombs) = " << shortreadRecombs + shortreadNonRecombs << std::endl;
        
        if (falsePositiveRate > 0) {
            //  std::cout << "allInformativePairs.size() " << allInformativePairs.size() << std::endl;
            
            // ADJUST RECOMBINATION PROBABILITIES USING THE MEAN RATE FROM LONG FRAGMENTS (1Mb+) AS A PRIOR FOR ALL THE READS
            // ALSO INCORPORATING THE OVERALL FALSE POSITIVE RATE ESTIMATE
            for (int j = 0; j < allInformativePairs.size(); j++) {
                if (allInformativePairs[j]->isRecombined) {
                    //if (allInformativePairs[j]->isRecombined) {
                    double PriorPrR = longreadRate * allInformativePairs[j]->dist; if (PriorPrR > 1.0) PriorPrR = 1.0;
                    double PrORgivenR = allInformativePairs[j]->probabilityRecombined - falsePositiveRate;
                    double PrORgivenNR = allInformativePairs[j]->PrORgivenNR;
                    double BayesNumerator = PrORgivenR * PriorPrR;
                    double BayesDenominator = BayesNumerator + ( (falsePositiveRate+PrORgivenNR) * (1-PriorPrR) );
                    double BayesPosterior = BayesNumerator / BayesDenominator;
                    
                    //std::cout << "PrORgivenR: " << PrORgivenR << std::endl;
                    //std::cout << "BayesPosterior: " << BayesPosterior << std::endl;
                    //std::cout << std::endl;
                    
                    allInformativePairs[j]->probabilityRecombined = BayesPosterior;
                    
                    checkProbabilityValid(j, allInformativePairs[j]);
                }
            }
        }
    }
    
    void removeReadPairsAboveAndBelowGivenLength(const int minLength, const double maxLengthFractionOfChromosome = 0.5) {
        int firstHetPos = coveredHetPos.front(); int lastHetPos = coveredHetPos.back();
        int chrLength = lastHetPos - firstHetPos; assert(chrLength > 0);
        std::vector<int> indicesToRemove;
        for (int i = (int)allInformativePairs.size() - 1; i >= 0; i--) {
            int readPairLength = allInformativePairs[i]->posRight - allInformativePairs[i]->posLeft;
            if (readPairLength > (maxLengthFractionOfChromosome * chrLength)) indicesToRemove.push_back(i);
        }
        std::cout << "Removing " << indicesToRemove.size() << " read pairs because of length above " << maxLengthFractionOfChromosome;
        std::cout << " of the chromosome (" << maxLengthFractionOfChromosome * chrLength << "bp)" << std::endl;
        // Remove elements at specified indices
        for (const int& index : indicesToRemove) {
            if (index >= 0 && index < allInformativePairs.size()) {
                allInformativePairs.erase(allInformativePairs.begin() + index);
            }
        }
        
        indicesToRemove.clear();
        for (int i = (int)allInformativePairs.size() - 1; i >= 0; i--) {
            int readPairLength = allInformativePairs[i]->posRight - allInformativePairs[i]->posLeft;
            if (readPairLength < minLength) indicesToRemove.push_back(i);
        }
        std::cout << "Removing " << indicesToRemove.size() << " read pairs because of length below " << minLength << "bp" << std::endl;
        // Remove elements at specified indices
        for (const int& index : indicesToRemove) {
            if (index >= 0 && index < allInformativePairs.size()) {
                allInformativePairs.erase(allInformativePairs.begin() + index);
            }
        }
        
        // Update stats:
        numConcordant = 0; numDiscordant = 0;
        for (int i = (int)allInformativePairs.size() - 1; i >= 0; i--) {
            if (allInformativePairs[i]->isRecombined) numDiscordant++;
            else numConcordant++;
        }
        
    }
    
    void adjustReadPairsToMapEdges() {
        for (int i = 0; i != coveredHetPos.size(); i++) {
            hetPosToIndex[coveredHetPos[i]] = i;
        }
        
        //std::cerr << "coveredHetPos.size(): " << coveredHetPos.size() << std::endl;
        
        for (int i = 0; i != allInformativePairs.size(); i++) {
            if (allInformativePairs[i]->posLeft < coveredHetPos.front()) {
                allInformativePairs[i]->posLeft = coveredHetPos.front();
                //std::cerr << "coveredHetPos.front(): " << coveredHetPos.front() << std::endl;
            }
            if (allInformativePairs[i]->posRight > coveredHetPos.back()) {
                allInformativePairs[i]->posRight = coveredHetPos.back();
            }
            allInformativePairs[i]->indexLeft = hetPosToIndex.at(allInformativePairs[i]->posLeft);
            allInformativePairs[i]->indexRight = hetPosToIndex.at(allInformativePairs[i]->posRight);
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
        double rateAt1000bp = stats->windowRates.front();
        double rateAt1MbpPlus = stats->windowRates.back();
        double rateRatio = rateAt1MbpPlus/rateAt1000bp;
        int minLength = stats->windowSizeMins.front();
        
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
    
    void checkProbabilityValid(const int j, const DefiningRecombInfo* dri) {
        if (dri->probabilityRecombined < 0 || dri->probabilityRecombined > 1) {
            std::cout << "ERROR: Recombination probability out of bounds for this read pair:" << std::endl;
            std::cout << "read pair j = " << j << "; probabilityRecombined = " << dri->probabilityRecombined << std::endl;
            std::cout << "allInformativePairs[j]->phaseErrorP_left = " << dri->phaseErrorP_left << std::endl;
            std::cout << "allInformativePairs[j]->phaseErrorP_right = " << dri->phaseErrorP_right << std::endl;
            std::cout << "allInformativePairs[j]->baseErrorP_left = " << dri->baseErrorP_left << std::endl;
            std::cout << "allInformativePairs[j]->baseErrorP_right = " << dri->baseErrorP_right << std::endl;
            std::cout << "allInformativePairs[j]->doubleCrossoverProbability = " << dri->doubleCrossoverProbability << std::endl;
            std::cout << "allInformativePairs[j]->dist = " << dri->dist << std::endl;
            std::cout << "allInformativePairs[j]->posLeft = " << dri->posLeft << std::endl;
            std::cout << "allInformativePairs[j]->posRight = " << dri->posRight << std::endl;
            std::cout << "allInformativePairs[j]->isRecombined = " << dri->isRecombined << std::endl;
            //std::cout << "doubleCrossoverProbability = " << doubleCrossoverProbability << std::endl;
            exit(1);
        }
    }
    
};

class RecombInterval: public RecombIntervalBase {
public:
    RecombInterval() {};
    
    RecombInterval(int leftIndex, double meanRate, RecombReadPairs* rp): directReadCoverageConcord(0), directReadCoverageDiscord(0), bSufficientCoverage(false) {
        j = leftIndex;
        recombFractionPerBp = meanRate;
        leftCoord = rp->coveredHetPos[j]; rightCoord = rp->coveredHetPos[j + 1];
        dj = rightCoord - leftCoord + 1;
        rj = recombFractionPerBp * dj;
        initialiseInterval(rp);
        effectiveCoverage = (int)coveringReadPairs.size();
    };
    
    std::vector<double> bootstrapRecombFractions;
    
    std::vector<DefiningRecombInfo*> coveringReadPairs;
    int effectiveCoverage;
    int directReadCoverageConcord;
    int directReadCoverageDiscord;
    double sum_P_ij;
    double sumConcordantFraction;
    double sumRecombFractionPerBP;
    bool bSufficientCoverage;
    
    void updateVals() {
        recombFractionPerBp = sumRecombFractionPerBP/sumConcordantFraction;
        rj = recombFractionPerBp * dj;
        if(rj > 1) {
            std::cerr << "ERROR: rj needs to be a fraction < 1" << std::endl;
            std::cerr << "rj = " << rj << std::endl;
            std::cerr << "sumRecombFractionPerBP = " << sumRecombFractionPerBP << std::endl;
            std::cerr << "sumConcordantFraction = " << sumConcordantFraction << std::endl;
            std::cerr << "dj = " << dj << std::endl;
            exit(1);
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
            
            getDirectCoverageOnLeftCoord(rp->allInformativePairs[i]);
        }
    }

private:
    
    void getDirectCoverageOnLeftCoord(DefiningRecombInfo* rInf) {
        
            if(rInf->posLeft == leftCoord || rInf->posRight == leftCoord){
                if (rInf->isRecombined) {
                    directReadCoverageDiscord++;
                } else {
                    directReadCoverageConcord++;
                }
            }
        }
};



class RecombMap : public RecombMapBase {
public:
    RecombMap(RecombReadPairs* rp, double minCoverageFraction): delta(std::numeric_limits<double>::max()), EMiterationNum(0)  {
        // Get the mean recombination rate per bp
        readPairs = rp;
        meanRate = readPairs->stats->meanRate;
        std::cout << "Mean Recombination Rate = " << meanRate << std::endl;
        
        cummulativeRates.resize(readPairs->coveredHetPos.size() - 1);
        recombIntervals.resize(readPairs->coveredHetPos.size() - 1);
        int pc = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            // Initialise the recombination fractions assuming a uniform recombination rate along the genome (this uses the per-bp rate)
            RecombInterval ri(j, meanRate, readPairs); // Calculate the first iteration
            recombIntervals[j] = ri;
            int update5pcInterval = (int)recombIntervals.size() / 20;
            if (j > 0 && j % update5pcInterval == 0) { pc += 5; printPcUpdateOnPrompt(pc); }
        }
        
        meanEffectiveCoverage = calculateMeanCoverage();
        std::cout << "Mean Effective Coverage = " << meanEffectiveCoverage << std::endl;
        
        mapPhysicalLength = recombIntervals.back().rightCoord - recombIntervals.front().leftCoord + 1;
        // edgeMinimumCoverage = (mapPhysicalLength / (double) readPairs->stats->meanLength) * 5;
        edgeMinimumCoverage = minCoverageFraction * meanEffectiveCoverage;
        // std::cout << "edgeMinimumCoverage Coverage v2 = " << edgeMinimumCoverage << std::endl;
        std::cout << "edgeMinimumCoverage = " << minCoverageFraction * meanEffectiveCoverage << std::endl;
        
        removeIntervalsWithInsufficientEffectiveCoverage();
        
        // Filling these helper variables for the getAverageRateForPhysicalWindow function
        intervalCoordsVectors.resize(2);
        intervalCoordsVectors[0].resize(recombIntervals.size()); intervalCoordsVectors[1].resize(recombIntervals.size());
        intervalPerBPrVector.resize(recombIntervals.size(),0.0);
        for (int i = 0; i != recombIntervals.size(); i++) {
            intervalCoordsVectors[0][i] = recombIntervals[i].leftCoord;
            intervalCoordsVectors[1][i] = recombIntervals[i].rightCoord;
            intervalPerBPrVector[i] = recombIntervals[i].recombFractionPerBp;
        }
    };
    
    std::vector<RecombInterval> recombIntervals;
    double edgeMinimumCoverage;
    RecombReadPairs* readPairs;
    
    // Variables for the EM algorithm
    int EMiterationNum = 0; double delta;
    
    bool debug = false;
    
    void firstUpdate() {
        std::vector<int> cannotGetAverage;
        for (int j = 0; j < recombIntervals.size(); j++) {
            assert(recombIntervals[j].coveringReadPairs.size() > edgeMinimumCoverage);
            //     if (recombIntervals[j].coveringReadPairs.size() > edgeMinimumCoverage) {
            recombIntervals[j].updateVals();
        }
         /*      if (j > 6000 && recombIntervals[j].recombFractionPerBp < 5e-08) {
                   std::cout << "Interval: " << j << std::endl;
                   std::cout << "weighedAvgRate: " << recombIntervals[j].recombFractionPerBp << std::endl;
                   std::cout << "recombIntervals[j].coveringReadPairs.size(): " << recombIntervals[j].coveringReadPairs.size() << std::endl;
                   //std::cout << "rightMostIndex: " << rightMostIndex << std::endl;
                   //print_vector(coverageCounts, std::cout, '\n');
                   exit(1);
                }
          
                
                
            } else {
                std::cout << "Removing interval: " << j << std::endl;
                std::cout << "Coverage: " << recombIntervals[j].coveringReadPairs.size() << std::endl;
                std::cout << "Effective coverage: " << recombIntervals[j].effectiveCoverage << std::endl;
                std::cout << "leftCoord: " << recombIntervals[j].leftCoord << "; rightCoord: " << recombIntervals[j].rightCoord << std::endl;
                cannotGetAverage.push_back(j);
            }
        }
        
       std::vector<int> cannotGetAverage;
        for (int j = 0; j < recombIntervals.size(); j++) {
            if (recombIntervals[j].coveringReadPairs.size() <= edgeMinimumCoverage) {
                recombIntervals[j].recombFractionPerBp = weighedAverageRateUnderCoveringReads(recombIntervals[j].coveringReadPairs, j);
                if (recombIntervals[j].recombFractionPerBp == -1) {
                    cannotGetAverage.push_back(j);
                } else {
                    recombIntervals[j].rj = recombIntervals[j].recombFractionPerBp * recombIntervals[j].dj;
                }
            }
        }
        
        std::sort(cannotGetAverage.begin(), cannotGetAverage.end(), std::greater<int>());
        // Remove elements at specified indices
        for (const int& index : cannotGetAverage) {
            if (index >= 0 && index < recombIntervals.size()) {
                recombIntervals.erase(recombIntervals.begin() + index);
                intervalCoordsVectors[0].erase(intervalCoordsVectors[0].begin() + index);
                intervalCoordsVectors[1].erase(intervalCoordsVectors[1].begin() + index);
                intervalPerBPrVector.erase(intervalPerBPrVector.begin() + index);
            }
        }
    */
        
        mapLength = calculateCummulativeRates();
        std::cout << std::endl;
    }
    
    void printPcUpdateOnPrompt(int pc) {
        if (pc == 100) std::cout << pc << "% (DONE)" << std::endl;
        else { std::cout << pc << "%..."; std::cout.flush(); }
    }
    
    void EMiteration(bool bLoud = true) {
        EMiterationNum++;
        if (bLoud) std::cout << "Iteration " << EMiterationNum << ": " << std::endl;
        updatePijs(bLoud);
        delta = updateRecombFractions(bLoud);
        double previousMapLength = mapLength;
        mapLength = calculateCummulativeRates();
        if (previousMapLength - mapLength > 10) {
            debug = true;
        }
        if (bLoud) std::cout << std::endl;
    }
    
    void printUpdateOnPrompt(int j) {
        std::cout << "numProcessedHets: " << j << " ("<< (double)j/recombIntervals.size() << "%)"<< std::endl;
        std::cout << "pos: " << recombIntervals[j].leftCoord << "bp"<< std::endl;
    }
    
    void outputMapToFile(string fileName) {
        std::cout << "Writing recombination map into: " << fileName << std::endl;
        std::ofstream* f = new std::ofstream(fileName);
        // That's just approximating the chromosome beginning; do we want it?
        *f << "1\t" << recombIntervals[0].leftCoord << "\t" << recombIntervals[0].recombFractionPerBp << std::endl;
        // Now the actual map
        for (int i = 0; i != recombIntervals.size(); i++) {
            *f << recombIntervals[i].leftCoord << "\t" << recombIntervals[i].rightCoord << "\t" << recombIntervals[i].recombFractionPerBp << std::endl;
        }
        // std::cout << "...DONE!" << std::endl;
    }
    
    void outputBootstrapToFile(string fileName) {
        std::cout << "Writing bootstrapped recombination maps into:" << fileName << std::endl;
        std::ofstream* f = new std::ofstream(fileName);
        // That's just approximating the chromosome beginning; do we want it?
        *f << "1\t" << recombIntervals[0].leftCoord << "\t"; print_vector(recombIntervals[0].bootstrapRecombFractions, *f);
        // Now the actual map
        for (int i = 0; i != recombIntervals.size(); i++) {
            *f << recombIntervals[i].leftCoord << "\t" << recombIntervals[i].rightCoord << "\t";
            print_vector(recombIntervals[i].bootstrapRecombFractions, *f);
        }
        // std::cout << "...DONE!" << std::endl;
    }
    
    void calculateMapForFixedWindowSizes(int windowSize) {
        
        for (int i = 0; i != recombIntervals.size(); i++) {
            intervalCoordsVectors[0][i] = recombIntervals[i].leftCoord;
            intervalCoordsVectors[1][i] = recombIntervals[i].rightCoord;
            intervalPerBPrVector[i] = recombIntervals[i].recombFractionPerBp;
           // std::cout << "intervalPerBPrVector[i]: " << intervalPerBPrVector[i] << std::endl;
        }
        
        int lastCoord = intervalCoordsVectors[1].back();
        int maxCoordToOutput = roundToNearestValue(lastCoord, windowSize);
        for (int start = 1; start < maxCoordToOutput; start = start + windowSize) {
            int physicalWindowEnd = start + windowSize - 1;
            physicalWindowStartEnd[0].push_back(start); physicalWindowStartEnd[1].push_back(physicalWindowEnd);
            double avg = getAverageRateForPhysicalWindow(start, physicalWindowEnd);
            physicalWindowR.push_back(avg);
        }
    }
    
    void outputMapToFileFixedWindowSizes(string fileName, int windowSize) {
        std::cout << "Writing recombination map with a fixed window size of " << windowSize << "bp into:" << fileName;
        std::ofstream* f = new std::ofstream(fileName);
        
       // std::cout << "physicalWindowStartEnd[0].size(): " << physicalWindowStartEnd[0].size() << fileName;
       // std::cout << "physicalWindowR.size(): " << physicalWindowR.size() << fileName;
        
        assert(physicalWindowStartEnd[0].size() == physicalWindowR.size());
        assert(physicalWindowStartEnd[0].size() == physicalWindowStartEnd[1].size());
        
        for (int i = 0; i != physicalWindowStartEnd[0].size(); i++) {
            int start = physicalWindowStartEnd[0][i];
            int physicalWindowEnd = physicalWindowStartEnd[1][i];
            *f << start << "\t" << physicalWindowEnd << "\t" << physicalWindowR[i] << std::endl;
            //std::cout << "avg: " << avg << std::endl;
        }
        std::cout << "...DONE!" << std::endl;
    }
    
    void printPerHetCoverageStats(string fn) {
        std::vector<int> problematicSNPs;
        std::ofstream* depthFile = new std::ofstream(fn);
        *depthFile << "pos" << "\t" << "ratePerBp" << "\t" << "EffectiveDepth" << "\t" << "directReadCoverageConcord" << "\t" <<
                      "directReadCoverageDiscord" << std::endl;
        for (int i = 0; i != recombIntervals.size(); i++) {
            int thisHetPos = recombIntervals[i].leftCoord;
            int coveredHetEffectiveDepth = (int)recombIntervals[i].coveringReadPairs.size();
            *depthFile << thisHetPos << "\t" << recombIntervals[i].recombFractionPerBp << "\t" << coveredHetEffectiveDepth << "\t" << recombIntervals[i].directReadCoverageConcord << "\t" << recombIntervals[i].directReadCoverageDiscord << std::endl;
            if (recombIntervals[i].directReadCoverageDiscord >= 2) problematicSNPs.push_back(thisHetPos);
        }
        std::cout << "There are = " << problematicSNPs.size() << " problematic intervals with directly covering read pairs" << std::endl;
    }
    
private:
    void getSumPijForInterval(int intervalIndex) {
        RecombInterval* i = &recombIntervals[intervalIndex];
        
        i->sum_P_ij = 0;
        for (int r = 0; r < i->coveringReadPairs.size(); r++) {
            if (i->coveringReadPairs[r]->probabilityRecombined == 0) continue;
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
            if(debug) {
                std::cerr << "DEBUGGING:" << std::endl;
                std::cerr << "rj = " << recombIntervals[j].rj << std::endl;
                std::cerr << "leftCoord = " << recombIntervals[j].leftCoord << std::endl;
                std::cerr << "recombIntervals[j].sum_P_ij = " << recombIntervals[j].sum_P_ij << std::endl;
                std::cerr << "sumRecombFractionPerBP = " << recombIntervals[j].sumRecombFractionPerBP << std::endl;
                std::cerr << "sumConcordantFraction = " << recombIntervals[j].sumConcordantFraction << std::endl;
                std::cerr << "dj = " << recombIntervals[j].dj << std::endl;
                std::cerr << "coveringReadPairs.size() = " << recombIntervals[j].coveringReadPairs.size() << std::endl;
                //std::cerr << "coveringReadPairs.size() = " << recombIntervals[j].s << std::endl;
                for (int r = 0; r < recombIntervals[j].coveringReadPairs.size(); r++) {
                  //  if (!i->coveringReadPairs[r]->isRecombined) continue;
                    if (std::isnan(recombIntervals[j].coveringReadPairs[r]->sum_r_k)) {
                        recombIntervals[j].coveringReadPairs[r]->sum_r_k = 0;
                        for (int k = recombIntervals[j].coveringReadPairs[r]->indexLeft; k < recombIntervals[j].coveringReadPairs[r]->indexRight; k++) {
                            recombIntervals[j].coveringReadPairs[r]->sum_r_k += recombIntervals[k].rj;
                        }
                    }
                    double p_ij = (recombIntervals[j].rj *  recombIntervals[j].coveringReadPairs[r]->probabilityRecombined) / recombIntervals[j].coveringReadPairs[r]->sum_r_k;  // eq. (2) from proposal
                    if(p_ij > 1) {
                        std::cerr << "dist = " << recombIntervals[j].coveringReadPairs[r]->dist << "; p_ij = " << p_ij << std::endl;
                        
                        std::cerr << "coveringReadPairs[r]->probabilityRecombined = " << recombIntervals[j].coveringReadPairs[r]->probabilityRecombined << std::endl;
                        std::cerr << "coveringReadPairs[r]->sum_r_k = " << recombIntervals[j].coveringReadPairs[r]->sum_r_k << std::endl;
                        std::cerr << "coveringReadPairs[r]->posLeft = " << recombIntervals[j].coveringReadPairs[r]->posLeft << std::endl;
                        std::cerr << "coveringReadPairs[r]->posRight = " << recombIntervals[j].coveringReadPairs[r]->posRight << std::endl;
                        std::cerr << "Intervals covered by this funky read:" << std::endl;
                        
                        for (int k = recombIntervals[j].coveringReadPairs[r]->indexLeft; k < recombIntervals[j].coveringReadPairs[r]->indexRight; k++) {
                            std::cerr << "leftCoord:" << recombIntervals[k].leftCoord << std::endl;
                            std::cerr << "rightCoord:" << recombIntervals[k].rightCoord << std::endl;
                            std::cerr << "rj:" << recombIntervals[k].rj << std::endl;
                            std::cerr << std::endl;
                        }
                        std::cerr << std::endl;
                    }
                    
                   // i->sum_P_ij += p_ij;
                }
                
                
                
                exit(1);
            }
            
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
    
    double updateRecombFractions(bool bLoud = true) {
        // std::cout << "Updating recombination fractions: " << std::endl;
        double delta = 0;
        int pc = 0; int update5pcInterval = (int)recombIntervals.size() / 20;
        if (bLoud) std::cout << "Updating recombination fractions: ";
        
        std::vector<double> updatedRecombFractionsPerBp; updatedRecombFractionsPerBp.resize(recombIntervals.size());
        for (int j = 0; j < recombIntervals.size(); j++) {
            if (bLoud) if (j > 0 && j % update5pcInterval == 0) { pc += 5; printPcUpdateOnPrompt(pc); }
            
            if (recombIntervals[j].effectiveCoverage > edgeMinimumCoverage) {
                updatedRecombFractionsPerBp[j] = (recombIntervals[j].sum_P_ij / recombIntervals[j].sumConcordantFraction)/recombIntervals[j].dj;
            }
        }
        
        for (int j = 0; j < recombIntervals.size(); j++) {
            if (recombIntervals[j].effectiveCoverage <= edgeMinimumCoverage) {
              //  std::cout << "j: " << j << std::endl;
                updatedRecombFractionsPerBp[j] = weighedAverageRateUnderCoveringReads(recombIntervals[j].coveringReadPairs, j);
               // std::cout << "newRj: " << newRj << std::endl;
               // std::cout << "newRj/ recombIntervals[j].dj: " << newRj/ recombIntervals[j].dj << std::endl;
            } //else {
            //    newRj = meanRate * recombIntervals[j].dj;
            //}
        }
        
        for (int j = 0; j < recombIntervals.size(); j++) {
            recombIntervals[j].recombFractionPerBp = updatedRecombFractionsPerBp[j];
            double newRj = updatedRecombFractionsPerBp[j] * recombIntervals[j].dj;
            if (!std::isnan(newRj) && !std::isnan(recombIntervals[j].rj)) {
                delta += std::abs(newRj - recombIntervals[j].rj);
            }
            recombIntervals[j].rj = newRj;
        }
        
        if (bLoud) std::cout << "delta: " << delta << std::endl;
        return delta;
    }
    
    
    // Weighed average recombination rate for reads that cover this interval
    double weighedAverageRateUnderCoveringReads(std::vector<DefiningRecombInfo*>& coveringReadPairs, int rpJ) {
        std::vector<int> coverageCounts; coverageCounts.resize(recombIntervals.size(),0);
        int leftMostIndex = std::numeric_limits<int>::max(); int rightMostIndex = 0;
        for (int i = 0; i < coveringReadPairs.size(); i++) {
            if (coveringReadPairs[i]->indexLeft < leftMostIndex) leftMostIndex = coveringReadPairs[i]->indexLeft;
            if (coveringReadPairs[i]->indexRight > rightMostIndex) rightMostIndex = coveringReadPairs[i]->indexRight;
            for (int j = coveringReadPairs[i]->indexLeft; j < coveringReadPairs[i]->indexRight; j++) {
                coverageCounts[j]++;
            }
        }
        
        
        double sumRecomb = 0; double sumLengths = 0;
        for (int j = leftMostIndex; j <= rightMostIndex;j++) {
            if(recombIntervals[j].bSufficientCoverage) {
                sumRecomb += recombIntervals[j].rj * coverageCounts[j];
                sumLengths += recombIntervals[j].dj * coverageCounts[j];
            }
        }
        //if (sumLengths == 0) { return -1;}
        double weighedAvgRate = sumRecomb/sumLengths;
        
       /* std::cout << "Interval: " << rpJ << std::endl;
        std::cout << "weighedAvgRate: " << weighedAvgRate << std::endl;
       if (weighedAvgRate < 5e-08) {
            std::cout << "leftMostIndex: " << leftMostIndex << std::endl;
            std::cout << "rightMostIndex: " << rightMostIndex << std::endl;
            print_vector(coverageCounts, std::cout, '\n');
            exit(1);
        }
        
        */
     
       // std::cout << "vector_sum(coverageCounts): " << vector_sum(coverageCounts) << std::endl;
       // std::cout << "weighedAvgRate: " << weighedAvgRate << std::endl;
        return weighedAvgRate;
    }
    
    // Returns total map length in cM
    double calculateCummulativeRates() {
        double cumSum = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            if (std::isnan(recombIntervals[j].rj)) {
                if (j > 0) cummulativeRates[j] = cummulativeRates[j-1];
                else cummulativeRates[j] = 0;
            } else {
                cummulativeRates[j] = cumSum * 100;
                cumSum += recombIntervals[j].rj;
            }
        }
        return (cumSum * 100);
    }
    
    double calculateMeanCoverage() {
        long long int total = 0;
        for (int j = 0; j < recombIntervals.size(); j++) {
            total += recombIntervals[j].effectiveCoverage;
            assert(total > 0); // check for overflow
        }
        double meanCoverage = total/(double)recombIntervals.size();
        return meanCoverage;
    }
    
    void removeIntervalsWithInsufficientEffectiveCoverage() {
        std::vector<int> indicesToRemove;
      //  std::cerr << "readPairs->coveredHetPos.size(): " << readPairs->coveredHetPos.size() << std::endl;
      //  std::cerr << "recombIntervals.size(): " << recombIntervals.size() << std::endl;
      //  assert(readPairs->coveredHetPos.size() == recombIntervals.size());
        for (int j = 0; j != recombIntervals.size(); j++) {
            if (recombIntervals[j].coveringReadPairs.size() <= edgeMinimumCoverage) {
                indicesToRemove.push_back(j);
            } else {
                recombIntervals[j].bSufficientCoverage = true;
            }
        }
        
        std::sort(indicesToRemove.begin(), indicesToRemove.end(), std::greater<int>());
        // Remove elements at specified indices
        for (const int& index : indicesToRemove) {
            if (index >= 0 && index < recombIntervals.size()) {
                recombIntervals.erase(recombIntervals.begin() + index);
                readPairs->coveredHetPos.erase(readPairs->coveredHetPos.begin() + index);
            }
        }
        readPairs->adjustReadPairsToMapEdges();
    }
};



#endif /* recombFromInformativePairsSAM_hpp */
