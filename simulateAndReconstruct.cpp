//
//  simulateAndReconstruct.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 12.09.23.
//

#include "simulateAndReconstruct.hpp"
#include "recombFromInformativePairsSAM.hpp"

#define SUBPROGRAM "Simulation"

#define DEBUG 1

static const char *SIMULATION_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] RecombMap.txt simulationInput.txt\n"
"Simulates input data corresponding to the parameters and the map and then uses this to attempt to reconstruct the map\n"
"\n"
HelpOption RunNameOption
"       -q, --min-BQ                            (default: 30) the minimum base quality for assesssing discordant phase\n"
"       -d, --min-Dist                          (default: 1000) the minimum distance (bp) to consider discordant phase a recombination\n"
"                                               as opposed to gene conversion\n"
"       -e, --epsilon=NUM                       (default: 0.0001) sets when the EM algorithm is deemed to have converged\n"
"                                               the smaller the epsilon the more EM iterations will be run\n"
"       -p, --min-PQ                            (default: 10) the minimum phase quality for assesssing discordant phase\n"
"       -s, --subsetHets=FILE.txt               (optional) Exclude the sites specified in this file\n"
"\n"
"OUTPUT OPTIONS:\n"
"       -f, --fixed-window=SIZE                 (optional) Output additional file with recombination map in windows of given SIZE (in bp)\n"
"                                               it will be output in file recombMap_RN_fW_SIZE.txt; the SIZE should be at least 1000bp\n"
"       -c, --coverageStats                     (optional) Output coverage over each het site into recombMap_wDepth_RN.txt\n"
"       -r, --readPairInfo                      (optional) Output coverage over each het site into discordantPairs_RN.txt and concordantPairs_RN.txt\n"
"       -b, --bootstrap=N                       (optional) Output N bootstrap replicates into bootstrap_RN.txt\n"
"       --simulationInput                       (optional) Output information needed for simulations to simulationInput_RN.txt\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:q:m:p:d:s:f:crb:e:";

enum { OPT_SIM_INPUT  };

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
    { "epsilon",   required_argument, NULL, 'e' },
    { "simulationInput",   required_argument, NULL, OPT_SIM_INPUT },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string hetsFile;
    static string samFile;
    static string runName = "";
    static int minMQ = 20;
    static int minBQ = 30;
    static int minPQ = 10;
    static int minDist = 1000;
    static string hetsSubset;
    static bool outputCoverageStats = false;
    static bool outputReadPairInfo = false;
    static bool bOutputForSimulations = false;
    static int physicalWindowSize = -1;
    static int nBootstrap = 0;
    static double epsilon = 0.00001;

    // These are fixed for now
    static int maxEMiterations = 10;
    static double minCoverage = 0.33;
    static int minDistanceToDefinePairs = 200;
}

int SimulationMain(int argc, char** argv) {
    parseSimulationOptions(argc, argv);
    
    std::cout << "1) Processing hets..." << std::endl;
    AllPhaseInfo* p = new AllPhaseInfo(opt::hetsFile, opt::minPQ, opt::hetsSubset);
    std::cout << std::endl;
    
    std::cout << "2) Loading read-pairs... " << std::endl;
    RecombReadPairs* rp = new RecombReadPairs(opt::samFile, opt::minDistanceToDefinePairs);
    std::cout << std::endl;
    
    std::cout << "3) Linking read-pairs and phased hets... " << std::endl;
    rp->linkWithHets(p->posToPhase, p->subsetLoci, opt::minBQ);
    rp->printBaseQualityStats();
    rp->printReadPairStatsByHetNumbers();
    
    std::cout << "4) Categorising concordant-discordant read-pairs... " << std::endl;
    rp->categorisePairs(opt::minDistanceToDefinePairs);
    rp->printReadPairStatsByCategory(opt::minDistanceToDefinePairs); rp->stats->collectAndPrintStats("Initial rates/stats by read pair distances: ", rp->allInformativePairs, true);
    rp->considerDoubleCrossovers(); rp->stats->collectAndPrintStats("Adjusted for double-crossovers: ", rp->allInformativePairs, true);
    //rp->adjustRecombinationProbabilities(); // Adjust probabilities based on read-length distributions
    rp->adjustRecombinationProbabilitiesBayes(); rp->adjustRecombinationProbabilities(); rp->stats->collectAndPrintStats("Adjusted for false positive rates: ", rp->allInformativePairs, true);
    rp->findUniqueHetsCoveredByReadsAndSortThem(); // Find and sort informative SNPs
    rp->removeReadPairsAboveAndBelowGivenLength(opt::minDist,0.5);
    
    rp->calculateDirectCoverageOnEachHet();
    rp->findAndRemoveReadPairsCoveringMultiHets(opt::runName);
    
    rp->findWhichHetsAreUsedByFinalReadPairSet(); // For each informative read-pair, find the indices of the bounding SNPs in the sorted het vector
    rp->findBoundingHetIndicesForEachReadPair();
    
    std::cout << "Left with " << rp->allInformativePairs.size() << " read pairs to build the map" << std::endl;
    
    rp->stats->collectAndPrintStats("Final read pair set detailed stats: ", rp->allInformativePairs, true);
    
    if (opt::outputReadPairInfo) {
        rp->printReadPairFileInfoIntoFile("discordantPairs" + opt::runName + ".txt", true);
        rp->printReadPairFileInfoIntoFile("concordantPairs" + opt::runName + ".txt", false);
    }
    std::cout << std::endl;
    
    std::cout << "Probability of two discordant pairs landing on the same SNP = " << binomialPMF(2, rp->numDiscordant*2, 1.0/p->posToPhase.size())  << std::endl;
        
    std::cout << "5) Making a genetic map... " << std::endl;
    
    // Now do the initial calculation of recombination fraction for each interval
    RecombMap* rm = new RecombMap(rp, opt::minCoverage);
    rm->firstUpdate();
    
    rm->outputMapToFile("recombMap_initial" + opt::runName + ".txt");
    
    std::cout << "Starting EM iterations..." << std::endl; int EMi = 0;
    while (rm->delta > (opt::epsilon * rp->stats->numDiscordant) && rm->EMiterationNum < opt::maxEMiterations) {
        EMi++;
        rm->EMiteration();
        std::cout << "Map length = " << rm->mapLength << std::endl;
        rm->outputMapToFile("recombMap_iteration_" + numToString(EMi) + opt::runName + ".txt");
    }
    std::cout << "DONE.... Map length = " << rm->mapLength << std::endl;
    
    rm->outputMapToFile("recombMap" + opt::runName + ".txt");
    if (opt::physicalWindowSize != -1) rm->outputMapToFileFixedWindowSizes("recombMap" + opt::runName + "_FW_" + numToString(opt::physicalWindowSize) + ".txt", opt::physicalWindowSize);
    
    // Optional: calculate and print coverage stats (effective coverage and direct coverage) per het site
    if(opt::outputCoverageStats) rm->printPerHetCoverageStats("recombMap_wDepth" + opt::runName + ".txt", rp);
    
    if (opt::bOutputForSimulations) { rp->outputInfoForSimulation("simulationInput" + opt::runName + ".txt"); }
     
    if (opt::nBootstrap > 0) {
        std::cout << std::endl;
        std::cout << "5a) Bootstrap... " << std::endl;
     //   std::cout << "Sample: 1" << std::endl;
       
    // Add the original ma as the first sample
       for (int j = 0; j < rm->recombIntervals.size(); j++) {
            rm->recombIntervals[j].bootstrapRecombFractions.push_back(rm->recombIntervals[j].recombFractionPerBp);
       }
        
        std::vector<DefiningRecombInfo*> orginalInformativePairs = rp->allInformativePairs;
        for (int i = 0; i < opt::nBootstrap; i++) {

            rp->allInformativePairs = rp->getBootstrapSample(orginalInformativePairs);
            for (int j = 0; j < rm->recombIntervals.size(); j++) {
                rm->recombIntervals[j].initialiseInterval(rp);
                // rm->recombIntervals[j].updateVals(opt::minCoverage * rm->meanEffectiveCoverage);
            }
            rm->firstUpdate();
            rm->delta = std::numeric_limits<double>::max(); rm->EMiterationNum = 0;
            while (rm->delta > (opt::epsilon * rp->stats->numDiscordant) && rm->EMiterationNum < opt::maxEMiterations) {
                rm->EMiteration(false);
                //std::cout << "Map length = " << rm->mapLength << std::endl;
            }
            std::cout << "DONE.... Map length = " << rm->mapLength << std::endl;
            for (int j = 0; j < rm->recombIntervals.size(); j++) {
                rm->recombIntervals[j].bootstrapRecombFractions.push_back(rm->recombIntervals[j].recombFractionPerBp);
            }
            std::cout << "Sample: " << i + 1 << std::endl;
        }
        rm->outputBootstrapToFile("bootstrap" + opt::runName + ".txt");
    }
     
    return 0;
}

void parseSimulationOptions(int argc, char** argv) {
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
            case 'e': arg >> opt::epsilon; break;
            case OPT_SIM_INPUT: opt::bOutputForSimulations = true; break;
            case 'h':
                std::cout << SIMULATION_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::physicalWindowSize != -1 && opt::physicalWindowSize < 1000) {
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
        std::cout << "\n" << SIMULATION_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hetsFile = argv[optind++];
    opt::samFile = argv[optind++];
}
