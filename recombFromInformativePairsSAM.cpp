//
//  recombFromInformativePairsSAM.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 02.09.22.
//


// Sam flags:
// 65 -     Paired, first in pair, + strand
// 81 -     Paired, first in pair, - strand
// 97 -     Paired, first in pair, + strand
// 113 -    Paired, first in pair, - strand
// 129 - Paired, second in pair, + strand
// 145 - Paired, second in pair, - strand
// 161 - Paired, second in pair, + strand
// 177 - Paired, second in pair, - strand
// >2000 - secondary alignment


#include "recombFromInformativePairsSAM.hpp"

#define SUBPROGRAM "RecombMap"

#define DEBUG 1

static const char *DISCORDPAIRS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] hapcutBlockFile.txt INFORMATVE_PAIRS.sam\n"
"Select reads/fragments from the PAIRTOOLS_FILE which could be informative about recombination:\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -m, --min-MQ                            (default: 20) the minimum mapping quality for a read to be considered\n"
"       -b, --min-BQ                            (default: 30) the minimum base quality for assesssing discordant phase\n"
"       -d, --min-Dist                          (default: 500) the minimum distance (bp) to consider discordant phase a recombination\n"
"                                               as opposed to gene conversion\n"
"       -p, --min-PQ                            (default: 30) the minimum phase quality for assesssing discordant phase (relevant with the --hapCut option)\n"
"       --hapCut                                the het positions come from HapCut output\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:b:m:p:d:";

//enum { OPT_ANNOT, OPT_AF  };
enum { OPT_HAPCUT  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "min-MQ",   required_argument, NULL, 'm' },
    { "min-BQ",   required_argument, NULL, 'b' },
    { "min-PQ",   required_argument, NULL, 'p' },
    { "min-Dist",   required_argument, NULL, 'd' },
    { "hapCut",   no_argument, NULL, OPT_HAPCUT },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static bool hapcutFormat = false;
    static string hetsFile;
    static string samFile;
    static string runName = "";
    static int minMQ = 20;
    static int minBQ = 30;
    static int minPQ = 30;
    static int minDist = 500;
}

int RecombFromSAMMain(int argc, char** argv) {
    parseRecombFromSAMOptions(argc, argv);
    string line; // for reading the input files
        
    std::ofstream* recombFile = new std::ofstream("recombMap" + opt::runName + ".txt");
    std::ofstream* depthFile = new std::ofstream("recombMap_wDepth" + opt::runName + ".txt");
    
    std::cout << "1) Processing hets..." << std::endl;
    AllPhaseInfo* p = new AllPhaseInfo(opt::hetsFile);
    
    std::cout << std::endl;
    std::cout << "2) Loading read-pairs... " << std::endl;
    RecombReadPairs* rp = new RecombReadPairs(opt::samFile);
    
    std::cout << std::endl;
    std::cout << "3) Linking read-pairs and phased hets... " << std::endl;
    rp->linkWithHets(p->posToPhase, opt::minBQ);
    rp->printBaseQualityStats();
    rp->printReadPairStats();
    
    std::cout << "4) Categorising concordant-discordant read-pairs... " << std::endl;
    // Now consider all HiC pairs i that map to the chromosome and cover a het at each end - lets call the physical locations of the hets x_i, y_i ( x < y) and let Z_i = 1 if they are in phase, 0 if out of phase (recombined)
    // Prod_i ((1-Z_i)(G(y_i)-G(x_i)) + Z_i*(1 -(G(y_i)-G(x_i))
    
    // std::vector<PhaseSwitch*> thisPairSwitches;
    int readPairsProcessed = 0;
    for (std::vector<RecombReadPair*>::iterator it = rp->informativeReadPairs.begin(); it != rp->informativeReadPairs.end(); it++) {
        readPairsProcessed++;
        RecombReadPair* thisReadPair = *it;
        
        thisReadPair->findIndicesOfConcordantAndDiscordantPairsOfHets(opt::minDist);
        thisReadPair->determineIfReadPairConcordantOrDiscordant();
        
        if (thisReadPair->pairRecombinationStatus == PAIR_DISCORDANT) {
            rp->numDiscordant++;
            DefiningRecombInfo* thisSwitch = thisReadPair->getDefiningHetPair(thisReadPair->switchPairI, thisReadPair->switchPairJ);
            rp->phaseSwitches.push_back(thisSwitch);
            rp->totalEffectiveLength = rp->totalEffectiveLength + thisSwitch->dist;
        } else if (thisReadPair->pairRecombinationStatus == PAIR_CONCORDANT) {
            rp->numConcordant++;
            DefiningRecombInfo* thisConcordantInfo = thisReadPair->getDefiningHetPair(thisReadPair->concordPairI, thisReadPair->concordPairJ);
            rp->concordantPairs.push_back(thisConcordantInfo);
            rp->totalEffectiveLength = rp->totalEffectiveLength + thisConcordantInfo->dist;
        }
        
    }
    
    rp->printConcordDiscordStats();
    rp->printSwitchInfoIntoFile("switches" + opt::runName + ".txt");
    
    std::cout << std::endl;
    std::cout << "5) Making a genetic map... " << std::endl;
    std::sort(rp->coveredHetPos.begin(), rp->coveredHetPos.end());
    std::vector<int>::iterator it = std::unique(rp->coveredHetPos.begin(), rp->coveredHetPos.end());
    rp->coveredHetPos.resize(distance(rp->coveredHetPos.begin(),it));
    std::cout << "coveredHetPos.size() " << rp->coveredHetPos.size() << std::endl;
    
    double meanRecombinationRate = (double)rp->numDiscordant/(double)rp->totalEffectiveLength;
    std::cout << "meanRecombinationRate " << meanRecombinationRate << std::endl;
    std::vector<double> recombFractions(rp->coveredHetPos.size()+1, meanRecombinationRate);
    
    int numProcessedHets = 0;
    for (int i = 0; i != rp->coveredHetPos.size() - 1; i++) {
        int left = rp->coveredHetPos[i];
        int right = rp->coveredHetPos[i + 1];
        //int distSNPs = (right - left) + 1;
        
        int coveringReadPairs = 0;
        
        double totalRecombFractionPerBP = 0;
        for (int j = 0; j != rp->phaseSwitches.size(); j++) {
            if(rp->phaseSwitches[j]->posLeft <= left && rp->phaseSwitches[j]->posRight >= right){
                coveringReadPairs++;
                double recombFractionPerBP = (double)1.0/(double)rp->phaseSwitches[j]->dist;
                totalRecombFractionPerBP += recombFractionPerBP;
            }
        }
      //  std::cout << "totalRecombFraction: " << totalRecombFraction << std::endl;
        
        double totalConcordantFraction = 0;
        for (int j = 0; j != rp->concordantPairs.size(); j++) {
            if(rp->concordantPairs[j]->posLeft <= left && rp->concordantPairs[j]->posRight >= right){
                coveringReadPairs++;
                //double concordPairFraction = (double)distSNPs/(double)(phaseConcordanceCoords[j][2]);
                totalConcordantFraction++;
            }
        }
       // std::cout << "totalConcordantFraction: " << totalConcordantFraction << std::endl;
        
        if (coveringReadPairs > 10) {
            recombFractions[i+1] = totalRecombFractionPerBP/totalConcordantFraction;
        }
        
        numProcessedHets++;
        if (numProcessedHets % 10000 == 0) {
            std::cout << "numProcessedHets: " << numProcessedHets << " ("<< (double)numProcessedHets/rp->coveredHetPos.size() << "%)"<< std::endl;
            std::cout << "pos: " << left << "bp"<< std::endl;
            
        }
        rp->coveredHetEffectiveDepth.push_back(coveringReadPairs);
    }
    
    *recombFile << "0\t" << rp->coveredHetPos[0] << "\t" << recombFractions[0] << std::endl;
    for (int i = 1; i != rp->coveredHetPos.size(); i++) {
        *recombFile << rp->coveredHetPos[i-1] << "\t" << rp->coveredHetPos[i] << "\t" << recombFractions[i] << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "6) Collecting stats... " << std::endl;
    
    // a) Prints out the distribition of recombination rates based on read pairs with variable-length inserts
    collectRateStatsBasedOnInsertLength(rp->phaseSwitches, rp->concordantPairs);
    
    int coverageDirectlyOnTheHet = 0;
    for (int i = 0; i != rp->coveredHetPos.size() - 1; i++) {
        int thisHetPos = rp->coveredHetPos[i];
        for (int j = 0; j != rp->phaseSwitches.size(); j++) {
            if(rp->phaseSwitches[j]->posLeft == thisHetPos || rp->phaseSwitches[j]->posRight == thisHetPos){
                coverageDirectlyOnTheHet++;
            }
        }
        for (int j = 0; j != rp->concordantPairs.size(); j++) {
            if(rp->concordantPairs[j]->posLeft == thisHetPos || rp->concordantPairs[j]->posRight == thisHetPos){
                coverageDirectlyOnTheHet++;
            }
        }
        rp->coveredHetDirectDepth.push_back(coverageDirectlyOnTheHet);
        *depthFile << thisHetPos << "\t" << recombFractions[i] << "\t" << rp->coveredHetEffectiveDepth[i] << "\t" << coverageDirectlyOnTheHet << std::endl;
    }
    
    return 0;
}

void parseRecombFromSAMOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'b': arg >> opt::minBQ; break;
            case 'm': arg >> opt::minMQ; break;
            case 'd': arg >> opt::minDist; break;
            case OPT_HAPCUT: opt::hapcutFormat = true; break;
            case 'h':
                std::cout << DISCORDPAIRS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DISCORDPAIRS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hetsFile = argv[optind++];
    opt::samFile = argv[optind++];
}

void collectRateStatsBasedOnInsertLength(const std::vector<DefiningRecombInfo*>& phaseSwitches, const std::vector<DefiningRecombInfo*>& phaseConcordanceCoords) {
    // Size windows are 0 - 1000bp, 1001 - 2000 bp, 2000 - 5000bp, 5000 - 10000, 10000 - 100000, 100000 - 1000000, 1M+
    std::vector<int> numRecombsInSizeWindows(7,0); std::vector<int> numNonRecombsInSizeWindows(7,0);
    std::vector<long long int> lengthOfInformativeSequenceWindows(7,0);
    std::vector<int> windowSizeMins = {0,1000,2000,5000,10000,100000,1000000};
    std::vector<int> windowSizeMax = {1000,2000,5000,10000,100000,1000000,1000000000};
    int totalRecombs = 0; int totalNonRecombs = 0; long long int totalL = 0;
    for (int j = 0; j != phaseSwitches.size(); j++) {
        int l = phaseSwitches[j]->posRight - phaseSwitches[j]->posLeft + 1;
        //std::cout << "l = " << l << std::endl;
        for (int k = 0; k != lengthOfInformativeSequenceWindows.size(); k++) {
            if (l > windowSizeMins[k] && l <= windowSizeMax[k]) {
                numRecombsInSizeWindows[k]++; totalRecombs++; lengthOfInformativeSequenceWindows[k] += l; totalL += l;
            }
        }
    }
    for (int j = 0; j != phaseConcordanceCoords.size(); j++) {
        int l = phaseConcordanceCoords[j]->dist;
        
        for (int k = 0; k != lengthOfInformativeSequenceWindows.size(); k++) {
            if (l > windowSizeMins[k] && l <= windowSizeMax[k]) {
                numNonRecombsInSizeWindows[k]++; totalNonRecombs++; lengthOfInformativeSequenceWindows[k] += l; totalL += l;
            }
        }
    }

    for (int j = 0; j != lengthOfInformativeSequenceWindows.size(); j++) {
        double thisWindowRate = (double)numRecombsInSizeWindows[j]/lengthOfInformativeSequenceWindows[j];
        std::cout << "window: " << windowSizeMins[j] << " - " << windowSizeMax[j] <<
            "; rate = " << thisWindowRate << "; n recomb = " << numRecombsInSizeWindows[j] << "; n non-recomb = " << numNonRecombsInSizeWindows[j] <<
            "; seqLength = " << lengthOfInformativeSequenceWindows[j] << std::endl;
    }
    std::cout << "total rate = " << (double)totalRecombs/totalL << "; n recomb = " << totalRecombs << "; n non-recomb = " << totalNonRecombs << "; seqLength = " << totalL << std::endl;

}
