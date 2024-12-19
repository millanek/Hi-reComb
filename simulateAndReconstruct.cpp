//
//  simulateAndReconstruct.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 12.09.23.
//

#include "simulateAndReconstruct.hpp"

#define SUBPROGRAM "Simulate"

#define DEBUG 1

static const char *SIMULATION_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] ReferenceMap.txt INFORMATVE_PAIRS.sam\n"
"Simulates Hi-C read pairs corresponding to the parameters and the reference map, with insert size distribution sampled from 'INFORMATVE_PAIRS.sam'. Then uses the simulated data to reconstruct the map using the same algorithm as in the 'RecombMap' command.\n"
"\n"
HelpOption RunNameOption
"       -v, --verbose                           Verbose info to std_out\n"
"       -r, --replicates=NUM                    (default: 10) Number of times to simulate and reconstruct\n"
"\n"
"SIMULATION OPTIONS:\n"
"       -t, --targetCoverage                    (default: 1000) effective coverage of the simulated data\n"
"       -x, --errorRate                         (default: 0.01) error rate per read pair (false positive or false negative)\n"
"\n"
"RECONSTRUCTION OPTIONS:\n"
"       -q, --min-BQ                            (default: 30) the minimum base quality for assesssing discordant phase\n"
"       -d, --min-Dist                          (default: 1000) the minimum distance (bp) to consider discordant phase a recombination\n"
"                                               as opposed to gene conversion\n"
"       -e, --epsilon=NUM                       (default: 0.0001) sets when the EM algorithm is deemed to have converged\n"
"                                               the smaller the epsilon the more EM iterations will be run\n"
"       -m, --maxEM=NUM                         (default: 10) maximum number of EM iterations to run\n"
//"       -p, --min-PQ                            (default: 10) the minimum phase quality for assesssing discordant phase\n"
"\n"
"OUTPUT OPTIONS:\n"
"       -f, --fixed-window=SIZE                 (default: 2000) Fixed window size (in bp) for the output; the SIZE should be at least 1000bp\n"
"       --readPairInfo                          (optional) Output location of reads pairs for the first run to discordantPairs_RN.txt and concordantPairs_RN.txt\n"
"       -i, --intermediateMaps                  (optional) Output maps at each EM iteration for the first run\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:q:m:p:d:f:r:e:vt:x:i";

enum { OPT_PAIR_INFO  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "verbose",   no_argument, NULL, 'v' },
    
    { "replicates",   required_argument, NULL, 'r' },
    
    { "targetCoverage",   no_argument, NULL, 't' },
    { "errorRate",   no_argument, NULL, 'x' },
    
    { "fixed-window",   required_argument, NULL, 'f' },
    { "min-BQ",   required_argument, NULL, 'q' },
    { "min-PQ",   required_argument, NULL, 'p' },
    { "min-Dist",   required_argument, NULL, 'd' },
    { "coverageStats",   no_argument, NULL, 'c' },
    { "epsilon",   required_argument, NULL, 'e' },
    { "maxEM",   required_argument, NULL, 'm' },
    { "readPairInfo",   no_argument, NULL, OPT_PAIR_INFO },
    { "intermediateMaps",   no_argument, NULL, 'i' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string mapFile;
    static string samFile;
    static string runName = "";
    static bool v = false;

    static int numReplicates = 10;

    static long long int targetCoverage = 1000;
    static double errorRate = 0.01;

    static bool outputReadPairInfo = false;

    static int minBQ = 30;
    static int minPQ = 10;
    static int minDist = 1000;
    static int physicalWindowSize = 2000;
    static double epsilon = 0.00001;
    static int maxEMiterations = 10;
    static bool bIntermediateMaps = false;

    // These are fixed for now
    static double minCoverage = 0.33;
    static int minDistanceToDefinePairs = 200;
}

int SimulationMain(int argc, char** argv) {
    parseSimulationOptions(argc, argv);
    
    std::cout << "1) Reading the recombination map from which to simulate data ... " << std::endl;
    RecombMapForSimulation* rs = new RecombMapForSimulation(opt::mapFile);
    std::cout << std::endl;
    
    std::cout << "2) Loading read-pairs to determine read-pair distance distribution... " << std::endl;
    RecombReadPairs* rpRef = new RecombReadPairs(opt::samFile, opt::minDistanceToDefinePairs);
    std::cout << std::endl;
    
    std::vector<int> readDistances;
    for (std::vector<RecombReadPair*>::iterator it = rpRef->readPairs.begin(); it != rpRef->readPairs.end(); it++) {
        RecombReadPair* r = *it; int p1 = r->read1->readPos; int p2 = r->read2->readPos; int l = abs(p2 - p1);
        if (l > 200) readDistances.push_back(l);
    }
    
    std::cout << "3) Simulating read pairs and reconstructing maps... " << std::endl;
    std::cout << "Target Effective Coverage: " << opt::targetCoverage << std::endl;
    std::cout << "Error Rate (false positive and false negative crossovers): " << opt::errorRate << std::endl;
    std::cout << std::endl;
    
    std::vector<std::vector<double>> simulationResults2kb; simulationResults2kb.resize(opt::numReplicates + 1);
    rs->calculateMapForFixedWindowSizes(opt::physicalWindowSize);
    simulationResults2kb[0] = rs->physicalWindowR;
    
    for (int r = 1; r != opt::numReplicates + 1; r++) {
        std::cout << "Simulation replicate " << r << std::endl;
    
        std::cout << "3a) Simulating read-pairs: " << std::endl;
        std::cout << std::endl;
        RecombReadPairs* rp = new RecombReadPairs();
        rs->simuateReadPairs(rp, readDistances, opt::targetCoverage, opt::errorRate);
    
        std::cout << "3b) Categorising concordant-discordant read-pairs... " << std::endl;
        rp->categorisePairs(opt::minDistanceToDefinePairs);
        rp->printReadPairStatsByCategory(opt::minDistanceToDefinePairs);
        rp->stats->collectAndPrintStats("Initial rates/stats by read pair distances: ", rp->allInformativePairs, true, opt::v);
        
        rp->adjustRecombinationProbabilitiesBayes();
        rp->adjustRecombinationProbabilities();
        rp->stats->collectAndPrintStats("Adjusted for false positive rates: ", rp->allInformativePairs, true, opt::v);
        
        rp->findUniqueHetsCoveredByReadsAndSortThem(); // Find and sort informative SNPs
        rp->removeReadPairsAboveAndBelowGivenLength(opt::minDist,0.5);
        
        // std::cout << "Removal done... " << std::endl;
        rp->calculateDirectCoverageOnEachHetMap();
        // std::cout << "Direct coverage calculated... " << std::endl;
        rp->findAndRemoveReadPairsCoveringMultiHets(opt::runName);
        
        rp->findWhichHetsAreUsedByFinalReadPairSet(); // For each informative read-pair, find the indices of the bounding SNPs in the sorted het vector
        rp->updateBoundingHetIndicesForEachReadPair();
        
        std::cout << "Left with " << rp->allInformativePairs.size() << " read pairs to build the map" << std::endl;
        
        rp->stats->collectAndPrintStats("Final read pair set detailed stats: ", rp->allInformativePairs, true);
        
        if (opt::outputReadPairInfo && r == 1) {
            rp->printReadPairFileInfoIntoFile("discordantPairs" + opt::runName + ".txt", true);
            rp->printReadPairFileInfoIntoFile("concordantPairs" + opt::runName + ".txt", false);
        }
    
        std::cout << "3c) Map reconstruction... " << std::endl;
   
        // Now do the initial calculation of recombination fraction for each interval
        RecombMap* rm = new RecombMap(rp, opt::minCoverage);
        rm->firstUpdate();
        
        if (opt::bIntermediateMaps && r == 1) rm->outputMapToFile("recombMap_initial" + opt::runName + ".txt");
        
        std::cout << "Starting EM iterations..." << std::endl; int EMi = 0;
        while (rm->delta > (opt::epsilon * rp->stats->numDiscordant) && rm->EMiterationNum < opt::maxEMiterations) {
            EMi++;
            rm->EMiteration(opt::v);
            if (opt::v) std::cout << "Map length = " << rm->mapLength << std::endl;
            if (opt::bIntermediateMaps && r == 1 && (EMi == 3 || EMi == 5 || EMi == 7)) rm->outputMapToFile("recombMap_iteration_" + numToString(EMi) + opt::runName + ".txt");
        }
        std::cout << "DONE.... Map length = " << rm->mapLength << std::endl;
        std::cout << std::endl;
       // rm->outputMapToFile("recombMap" + opt::runName + ".txt");
        rm->calculateMapForFixedWindowSizes(opt::physicalWindowSize);
   //   rm->outputMapToFileFixedWindowSizes("recombMap" + opt::runName + "_FW_" + numToString(opt::physicalWindowSize) + ".txt", opt::physicalWindowSize);
        simulationResults2kb[r] = rm->physicalWindowR;
        simulationResults2kb[r].resize(simulationResults2kb[0].size(),NAN); // Sometimes the reconstructed map can be a bit shorter if there were no read pairs at the end
    }
    
    string fileName = "simulatedMaps_" + opt::runName + "_FW_" + numToString(opt::physicalWindowSize) + ".txt";
    std::cout << "Writing simulated recombination maps with a fixed window size of " << opt::physicalWindowSize << "bp into:" << fileName << std::endl;
    std::ofstream* f = new std::ofstream(fileName);
    
    for (int i = 0; i != simulationResults2kb[0].size(); i++) {
        *f << rs->physicalWindowStartEnd[0][i] << "\t"; *f << rs->physicalWindowStartEnd[1][i] << "\t";
        *f << simulationResults2kb[0][i];
       // std::cout << "and here" << std::endl;
        for (int r = 1; r != opt::numReplicates + 1; r++) {
            *f << "\t" << simulationResults2kb[r][i];
            //std::cout << "avg: " << avg << std::endl;
        }
        *f << std::endl;
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
            case 'v': opt::v = true; break;
                
            case 'r': arg >> opt::numReplicates; break;
            
            case 't': arg >> opt::targetCoverage; break;
            case 'x': arg >> opt::errorRate; break;
                
            case 'q': arg >> opt::minBQ; break;
            case 'm': arg >> opt::maxEMiterations; break;
            case 'd': arg >> opt::minDist; break;
            case 'p': arg >> opt::minPQ; break;
            case 'f': arg >> opt::physicalWindowSize; break;
            case 'e': arg >> opt::epsilon; break;
            
            case 'i': opt::bIntermediateMaps = true; break;
            case OPT_PAIR_INFO: opt::outputReadPairInfo = true; break;
            case 'h':
                std::cout << SIMULATION_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::errorRate < 0 || opt::errorRate > 1) {
        std::cerr << "Error: The error rate (-x parameter) should be between 0 and 1\n";
        die = true;
    }
    
    if (opt::physicalWindowSize != -1 && opt::physicalWindowSize < 1000) {
        std::cerr << "Error: The -f parameter should be at least 1000 bp\n";
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
    opt::mapFile = argv[optind++];
    opt::samFile = argv[optind++];
}
