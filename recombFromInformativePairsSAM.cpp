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
"Generate a recombination map from a phased hapcut2 file of heterozygous sites and a sam file with read pairs covering the hets\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -m, --min-MQ                            (default: 20) the minimum mapping quality for a read to be considered\n"
"       -b, --min-BQ                            (default: 30) the minimum base quality for assesssing discordant phase\n"
"       -d, --min-Dist                          (default: 1000) the minimum distance (bp) to consider discordant phase a recombination\n"
"                                               as opposed to gene conversion\n"
"       -p, --min-PQ                            (default: 0) the minimum phase quality for assesssing discordant phase\n"
"       -s, --subsetHets=FILE.txt               (optional) Exclude the sites specified in this file\n"
"\n"
"OUTPUT OPTIONS:\n"
"       -n, --run-name=RN                       (optional) Will be included in the output file name: recombMap_RN.txt\n"
"       -c, --coverageStats=OUTFILE             (optional) Output coverage over each het site into OUTFILE.txt\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:b:m:p:d:s:c:";

//enum { OPT_ANNOT, OPT_AF  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "min-MQ",   required_argument, NULL, 'm' },
    { "min-BQ",   required_argument, NULL, 'b' },
    { "min-PQ",   required_argument, NULL, 'p' },
    { "min-Dist",   required_argument, NULL, 'd' },
    { "subsetHets",   required_argument, NULL, 's' },
    { "coverageStats",   required_argument, NULL, 'c' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string hetsFile;
    static string samFile;
    static string runName = "";
    static int minMQ = 20;
    static int minBQ = 30;
    static int minPQ = 0;
    static int minDist = 1000;
    static string hetsSubset;
    static string coverageStatsFile;
}

int RecombFromSAMMain(int argc, char** argv) {
    parseRecombFromSAMOptions(argc, argv);
    
    std::cout << "1) Processing hets..." << std::endl;
    AllPhaseInfo* p = new AllPhaseInfo(opt::hetsFile, opt::minPQ, opt::hetsSubset);
    std::cout << std::endl;
    
    std::cout << "2) Loading read-pairs... " << std::endl;
    RecombReadPairs* rp = new RecombReadPairs(opt::samFile);
    std::cout << std::endl;
    
    std::cout << "3) Linking read-pairs and phased hets... " << std::endl;
    rp->linkWithHets(p->posToPhase, p->subsetLoci, opt::minBQ);
    rp->printBaseQualityStats();
    rp->printReadPairStats();
    
    std::cout << "4) Categorising concordant-discordant read-pairs... " << std::endl;
    // Now consider all HiC pairs i that map to the chromosome and cover a het at each end - lets call the physical locations of the hets x_i, y_i ( x < y) and let Z_i = 1 if they are in phase, 0 if out of phase (recombined)
    // Prod_i ((1-Z_i)(G(y_i)-G(x_i)) + Z_i*(1 -(G(y_i)-G(x_i))
    
    // std::vector<PhaseSwitch*> thisPairSwitches;
    int readPairsProcessed = 0; DefiningRecombInfo* thisPairInformation;
    for (std::vector<RecombReadPair*>::iterator it = rp->informativeReadPairs.begin(); it != rp->informativeReadPairs.end(); it++) {
        readPairsProcessed++;
        RecombReadPair* thisReadPair = *it;
        
        thisReadPair->findIndicesOfConcordantAndDiscordantPairsOfHets(opt::minDist);
        thisReadPair->determineIfReadPairConcordantOrDiscordant();
        
        if (thisReadPair->pairRecombinationStatus == PAIR_DISCORDANT) {
            rp->numDiscordant++;
            thisPairInformation = thisReadPair->getDefiningHetPair(thisReadPair->switchPairI, thisReadPair->switchPairJ);
            rp->phaseSwitches.push_back(thisPairInformation);
        } else if (thisReadPair->pairRecombinationStatus == PAIR_CONCORDANT) {
            rp->numConcordant++;
            thisPairInformation = thisReadPair->getDefiningHetPair(thisReadPair->concordPairI, thisReadPair->concordPairJ);
            rp->concordantPairs.push_back(thisPairInformation);
        }
        if (thisReadPair->pairRecombinationStatus != PAIR_AMBIGUOUS) {
            rp->allInformativePairs.push_back(thisPairInformation);
            rp->totalEffectiveLength = rp->totalEffectiveLength + thisPairInformation->dist;
        }
    }
    rp->printConcordDiscordStats();
    rp->printSwitchInfoIntoFile("switches" + opt::runName + ".txt");
    std::cout << std::endl;
    
    
    
    std::cout << "5) Making a genetic map... " << std::endl;
    rp->findUniqueHetsCoveredByReadsAndSortThem(); // Find and sort informative SNPs
    rp->findBoundingHetIndicesForEachReadPair(); // For each informative read-pair, find the indices of the bounding SNPs in the sorted het vector
    
    // Now do the initial calculation of recombination fraction for each interval
    RecombMap* rm = new RecombMap(rp);
    
    double delta = std::numeric_limits<double>::max(); int EMiterationNum = 0;
    std::cout << "Starting EM iterations..." << std::endl;
    while (delta > 0.1) {
        EMiterationNum++; delta = rm->EMiteration(EMiterationNum);
    }
    rm->outputMapToFile("recombMap" + opt::runName + ".txt");
    
  //  std::cout << "5a) Bootstrap... " << std::endl;
    
    
    
    std::cout << "6) Collecting stats... " << std::endl;
    // a) Prints out the distribition of recombination rates based on read pairs with variable-length inserts
    collectRateStatsBasedOnInsertLength(rp->phaseSwitches, rp->concordantPairs);
    // b) Calculate and print coverage stats (effective coverage and direct coverage) per het site
    if(!opt::coverageStatsFile.empty()) rm->calculateAndPrintPerHetCoverageStats(opt::coverageStatsFile, rp);
    
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
            case 'p': arg >> opt::minPQ; break;
            case 's': arg >> opt::hetsSubset; break;
            case 'c': arg >> opt::coverageStatsFile; break;
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
        int l = phaseSwitches[j]->dist;
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
