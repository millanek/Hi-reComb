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
#define epsilon 0.0001

static const char *DISCORDPAIRS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] hapcutBlockFile.txt INFORMATVE_PAIRS.sam\n"
"Generate a recombination map from a phased hapcut2 file of heterozygous sites and a sam file with read pairs covering the hets\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -m, --min-MQ                            (default: 20) the minimum mapping quality for a read to be considered\n"
"       -q, --min-BQ                            (default: 30) the minimum base quality for assesssing discordant phase\n"
"       -d, --min-Dist                          (default: 1000) the minimum distance (bp) to consider discordant phase a recombination\n"
"                                               as opposed to gene conversion\n"
"       -p, --min-PQ                            (default: 0) the minimum phase quality for assesssing discordant phase\n"
"       -s, --subsetHets=FILE.txt               (optional) Exclude the sites specified in this file\n"
"\n"
"OUTPUT OPTIONS:\n"
"       -n, --run-name=RN                       (optional) Will be included in the output file name: recombMap_RN.txt\n"
"       -f, --fixed-window=SIZE                 (optional) Output additional file with recombination map in windows of given SIZE (in bp)\n"
"                                               it will be output in file recombMap_RN_fW_SIZE.txt; the SIZE should be at least 1000bp\n"
"       -c, --coverageStats                     (optional) Output coverage over each het site into recombMap_wDepth_RN.txt\n"
"       -r, --readPairInfo                      (optional) Output coverage over each het site into discordantPairs_RN.txt and concordantPairs_RN.txt\n"
"       -b, --bootstrap=N                       (optional) Output N bootstrap replicates into bootstrap_RN.txt\n"

"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:q:m:p:d:s:f:crb:";

//enum { OPT_ANNOT, OPT_AF  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "fixed-window",   required_argument, NULL, 'f' },
    { "min-MQ",   required_argument, NULL, 'm' },
    { "min-BQ",   required_argument, NULL, 'q' },
    { "min-PQ",   required_argument, NULL, 'p' },
    { "min-Dist",   required_argument, NULL, 'd' },
    { "subsetHets",   required_argument, NULL, 's' },
    { "coverageStats",   no_argument, NULL, 'c' },
    { "readPairInfo",   no_argument, NULL, 'r' },
    { "bootstrap",   required_argument, NULL, 'b' },
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
    static bool outputCoverageStats = false;
    static bool outputReadPairInfo = false;
    static int physicalWindowSize = NAN;
    static int nBootstrap = 0;
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
            thisPairInformation = thisReadPair->getDefiningHetPair(thisReadPair->switchPairI, thisReadPair->switchPairJ);
        } else if (thisReadPair->pairRecombinationStatus == PAIR_CONCORDANT) {
            thisPairInformation = thisReadPair->getDefiningHetPair(thisReadPair->concordPairI, thisReadPair->concordPairJ);
        }
        if (thisReadPair->pairRecombinationStatus != PAIR_AMBIGUOUS) {
            rp->allInformativePairs.push_back(thisPairInformation);
        }
    }
    rp->stats->collectStats(rp->allInformativePairs);
    rp->stats->printRecombReadPairStats();
    rp->adjustRecombinationProbabilities(); // Adjust probabilities based on read-length distrubutions
    if (opt::outputReadPairInfo) {
        rp->printReadPairFileInfoIntoFile("switches" + opt::runName + ".txt", true);
        rp->printReadPairFileInfoIntoFile("concords" + opt::runName + ".txt", false);
    }
    std::cout << std::endl;
        
    std::cout << "5) Making a genetic map... " << std::endl;
    rp->findUniqueHetsCoveredByReadsAndSortThem(); // Find and sort informative SNPs
    rp->findBoundingHetIndicesForEachReadPair(); // For each informative read-pair, find the indices of the bounding SNPs in the sorted het vector
    
    // Now do the initial calculation of recombination fraction for each interval
    RecombMap* rm = new RecombMap(rp);
    
    double delta = std::numeric_limits<double>::max(); int EMiterationNum = 0;
    std::cout << "Starting EM iterations..." << std::endl;
    while (delta > (epsilon * rp->stats->numDiscordant)) {
        EMiterationNum++; delta = rm->EMiteration(EMiterationNum);
        std::cout << "Map length = " << rm->mapLength << std::endl;
    }
    std::cout << "DONE.... Map length = " << rm->mapLength << std::endl;
    
    rm->outputMapToFile("recombMap" + opt::runName + ".txt");
    if (!isnan(opt::physicalWindowSize)) rm->outputMapToFileFixedWindowSizes("recombMap" + opt::runName + "_FW_" + numToString(opt::physicalWindowSize) + ".txt", opt::physicalWindowSize);
    
     
    if (opt::nBootstrap > 0) {
        std::cout << std::endl;
        std::cout << "5a) Bootstrap... " << std::endl;
        std::cout << "Sample: 1" << std::endl;
        
        for (int j = 0; j < rm->recombIntervals.size(); j++) {
            rm->recombIntervals[j].bootstrapRecombFractions.push_back(rm->recombIntervals[j].recombFractionPerBp);
        }
        
        std::vector<DefiningRecombInfo*> orginalInformativePairs = rp->allInformativePairs;
        for (int i = 0; i < opt::nBootstrap; i++) {

            rp->allInformativePairs = rp->getBootstrapSample(orginalInformativePairs);
            for (int j = 0; j < rm->recombIntervals.size(); j++) {
                rm->recombIntervals[j].initialiseInterval(rp);
                rm->recombIntervals[j].updateVals(0.1 * rm->meanEffectiveCoverage);
            }
            double delta = std::numeric_limits<double>::max(); int EMiterationNum = 0;
            while (delta > (epsilon * rp->stats->numDiscordant)) {
                EMiterationNum++; delta = rm->EMiteration(EMiterationNum, false);
                //std::cout << "Map length = " << rm->mapLength << std::endl;
            }
            std::cout << "DONE.... Map length = " << rm->mapLength << std::endl;
            for (int j = 0; j < rm->recombIntervals.size(); j++) {
                rm->recombIntervals[j].bootstrapRecombFractions.push_back(rm->recombIntervals[j].recombFractionPerBp);
            }
            std::cout << "Sample: " << i + 1 << std::endl;
        }
        rm->outputBootstrapToFile("bootstrap_" + opt::runName + ".txt");
    }
    // Optional: calculate and print coverage stats (effective coverage and direct coverage) per het site
    if(opt::outputCoverageStats) rm->printPerHetCoverageStats("recombMap_wDepth" + opt::runName + ".txt", rp);
     
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
            case 'q': arg >> opt::minBQ; break;
            case 'm': arg >> opt::minMQ; break;
            case 'd': arg >> opt::minDist; break;
            case 'p': arg >> opt::minPQ; break;
            case 's': arg >> opt::hetsSubset; break;
            case 'f': arg >> opt::physicalWindowSize; break;
            case 'c': opt::outputCoverageStats = true; break;
            case 'r': opt::outputReadPairInfo = true; break;
            case 'b': arg >> opt::nBootstrap; break;
            case 'h':
                std::cout << DISCORDPAIRS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (!isnan(opt::physicalWindowSize) && opt::physicalWindowSize < 1000) {
        std::cerr << "Error: the -f parameter should be at least 1000 bp\n";
        die = true;
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
