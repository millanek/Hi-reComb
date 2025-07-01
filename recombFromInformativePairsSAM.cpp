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
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] HAPCUT2_PHASE.txt INFORMATVE_PAIRS.sam\n"
"Generate a recombination map from a phased hapcut2 file of heterozygous sites and a sam file with read pairs covering the hets\n"
"\n"
"COMMON OPTIONS:\n"
HelpOption RunNameOption
"       -v, --verbose                           Verbose info to std_out\n"
"       -q, --min-BQ                            (default: 30) minimum base quality at informative hets\n"
"       -p, --min-PQ                            (default: 10) minimum phase quality for informative hets\n"
"       -d, --minDist                           (default: 1000) minimum distance (bp) between fragments for detecting recombination\n"
"       -m, --maxEM=NUM                         (default: 10) maximum number of EM iterations to run\n"
"       -x, --minCoverageF=NUM                  (default: 0.2) Minimum coverage fraction around chromosome edges\n"
"\n"
"OUTPUT OPTIONS:\n"
"       -b, --bootstrap=N                       Output N bootstrap replicates into bootstrap_RN.txt\n"
"       -f, --fixedWindow=SIZE                  Output additional file with recombination map in windows of given SIZE (in bp)\n"
"                                               it will be output in file recombMap_RN_fW_SIZE.txt; the SIZE should be at least 1000bp\n"
"       -c, --coverageStats                     Output coverage over each het site into recombMap_wDepth_RN.txt\n"
"       -r, --readPairInfo                      Output coverage over each het site into discordantPairs_RN.txt and concordantPairs_RN.txt\n"
"       -i, --intermediateMaps                  Output maps at each EM iteration\n"
"\n"
"ADVANCED FINE TUNING OPTIONS:\n"
"       --minDforDoubleCrossovers=NUM           (default: 1000000) minimum distance (bp) to consider for double crossovers\n"
"                                               should be be based on the extent of crossover interference in the species\n"
"                                               if set to > chromosome length, double crossovers will not be considered\n"
"       --imbalancedSNPs=N_DISCORD,N_CONCORD    (default: 2,0) minimum number of discordant and maximum concordant SNPs to consider a site as imbalanced\n"
"       --FPrateThresholds=SHORT,LONG           (default: 1000,1000000) insert size thresholds for estimating false positive rate (i.e. the correction factor f)\n"
"                                               the comparison is for <=SHORT vs >=LONG insert sizes; the values should be in bp\n"
"       -s, --subsetHets=FILE.txt               (optional) Exclude the het sites whose coordinates are listed in this file (one het per line)\n"
"       -e, --epsilon=NUM                       (default: 0.0001) sets when the EM algorithm is deemed to have converged\n"
"                                               the smaller the epsilon the more EM iterations will be run, up to --maxEM\n"
"                                               i.e. the --maxEM option takes precedence over this one\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_DOUBLE_CROSS, OPT_IMBALANCED_SNPs, OPT_FPRATE_THRESHOLDS };

static const char* shortopts = "hn:q:m:p:d:s:f:crb:e:ivx:";

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "run-name",   required_argument, NULL, 'n' },
    { "fixed-window",   required_argument, NULL, 'f' },
    { "maxEM",   required_argument, NULL, 'm' },
    { "min-BQ",   required_argument, NULL, 'q' },
    { "min-PQ",   required_argument, NULL, 'p' },
    { "minDist",   required_argument, NULL, 'd' },
    { "subsetHets",   required_argument, NULL, 's' },
    { "coverageStats",   no_argument, NULL, 'c' },
    { "readPairInfo",   no_argument, NULL, 'r' },
    { "bootstrap",   required_argument, NULL, 'b' },
    { "epsilon",   required_argument, NULL, 'e' },
    { "intermediateMaps",   no_argument, NULL, 'i' },
    { "minCoverageF",   required_argument, NULL, 'x' },
    { "verbose",   no_argument, NULL, 'v' },
    { "minDforDoubleCrossovers", required_argument, NULL, OPT_DOUBLE_CROSS },
    { "imbalancedSNPs", required_argument, NULL, OPT_IMBALANCED_SNPs },
    { "FPrateThresholds", required_argument, NULL, OPT_FPRATE_THRESHOLDS },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string hetsFile;
    static string samFile;
    static string runName = "";
    static int minBQ = 30;
    static int minPQ = 10;
    static int minDist = 1000;
    static string hetsSubset;
    static bool outputCoverageStats = false;
    static bool outputReadPairInfo = false;
    static bool bIntermediateMaps = false;
    static bool v = false;
    static int physicalWindowSize = -1;
    static int nBootstrap = 0;
    static double epsilon = 0.01;
    static int maxEMiterations = 10;
    static double minDistForDoubleCrossovers = 1000000; // Default is 1Mb, can be set to a larger value to disable double crossovers
    static int minDiscordantSNPs = 2; // Minimum number of discordant SNPs to consider a site as imbalanced
    static int maxConcordantSNPs = 0; // Maximum number of concordant SNPs to consider a site as imbalanced
    static int FPrateShortThreshold = 1000; // Threshold for short insert size to estimate false positive rate
    static int FPrateLongThreshold = 1000000; // Threshold for long insert size to estimate false positive rate
    // These are fixed for now
    static double minCoverage = 0.2;
    static int minDistanceToDefinePairs = 200;
}

int RecombFromSAMMain(int argc, char** argv) {
    parseRecombFromSAMOptions(argc, argv);
    
    std::cout << "1) Processing hets..." << std::endl;
    AllPhaseInfo* p = new AllPhaseInfo(opt::hetsFile, opt::minPQ, opt::hetsSubset);
    std::cout << std::endl;
    
    std::cout << "2) Loading read-pairs... " << std::endl;
    RecombReadPairs* rp = new RecombReadPairs(opt::samFile, opt::minDistanceToDefinePairs, opt::FPrateShortThreshold, opt::FPrateLongThreshold);
    std::cout << std::endl;
    
    std::cout << "3) Linking read-pairs and phased hets... " << std::endl;
    rp->linkWithHets(p->posToPhase, p->subsetLoci, opt::minBQ);
    rp->printBaseQualityStats();
    rp->printReadPairStatsByHetNumbers();
    
    std::cout << "4) Categorising concordant-discordant read-pairs... " << std::endl;
    rp->categorisePairs(opt::minDistanceToDefinePairs);
    rp->printReadPairStatsByCategory(opt::minDistanceToDefinePairs); 
    rp->stats->collectAndPrintStats("Initial rates/stats by read pair distances: ", rp->allInformativePairs, true, opt::v);
    
    rp->considerDoubleCrossovers(opt::minDistForDoubleCrossovers); rp->stats->collectAndPrintStats("Adjusted for double-crossovers: ", rp->allInformativePairs, true, opt::v);
    //rp->adjustRecombinationProbabilities(); // Adjust probabilities based on read-length distributions
    rp->adjustRecombinationProbabilitiesBayes();
    rp->adjustRecombinationProbabilities(); 
    rp->stats->collectAndPrintStats("Adjusted for false positive rates: ", rp->allInformativePairs, true, opt::v);
    rp->findUniqueHetsCoveredByReadsAndSortThem(); // Find and sort informative SNPs
    rp->removeReadPairsAboveAndBelowGivenLength(opt::minDist,0.5);
    std::cout << "Removal done... " << std::endl;
    rp->calculateDirectCoverageOnEachHetMap();
    std::cout << "Direct coverage calculated... " << std::endl;
    rp->findAndRemoveReadPairsCoveringMultiHets(opt::minDiscordantSNPs, opt::maxConcordantSNPs); 
    
    rp->findWhichHetsAreUsedByFinalReadPairSet(); // For each informative read-pair, find the indices of the bounding SNPs in the sorted het vector
    rp->updateBoundingHetIndicesForEachReadPair();
    
    std::cout << "Left with " << rp->allInformativePairs.size() << " read pairs to build the map" << std::endl;
    
    rp->stats->collectAndPrintStats("Final read pair set detailed stats: ", rp->allInformativePairs, true);
    
    if (opt::outputReadPairInfo) {
        rp->printReadPairFileInfoIntoFile("discordantPairs" + opt::runName + ".txt", true);
        rp->printReadPairFileInfoIntoFile("concordantPairs" + opt::runName + ".txt", false);
    }
    std::cout << std::endl;
        
  //  double lambda = rp->numDiscordant*2 / (double)p->posToPhase.size();
  //  std::cout << "Probability of two discordant pairs landing on the same SNP Poisson = " << (pow(lambda, 2.0) * exp(-lambda)) / 2.0  << std::endl;

         
    
        
    std::cout << "5) Making a genetic map... " << std::endl;
    
    // Now do the initial calculation of recombination fraction for each interval
    RecombMap* rm = new RecombMap(rp, opt::minCoverage);
    rm->firstUpdate();
    
    if (opt::bIntermediateMaps) rm->outputMapToFile("recombMap_initial" + opt::runName + ".txt");
    
    std::cout << "Starting EM iterations..." << std::endl; int EMi = 0;
  //  while (rm->delta > (opt::epsilon * rp->stats->numDiscordant) && rm->EMiterationNum < opt::maxEMiterations) {
    while (rm->delta > opt::epsilon && rm->EMiterationNum < opt::maxEMiterations) {
        EMi++;
        rm->EMiteration(opt::v);
        if (opt::v) std::cout << "Map length = " << rm->mapLength << std::endl;
        if (opt::bIntermediateMaps) rm->outputMapToFile("recombMap_iteration_" + numToString(EMi) + opt::runName + ".txt");
    }
    std::cout << "DONE.... Map length = " << rm->mapLength << std::endl;
    
    rm->outputMapToFile("recombMap" + opt::runName + ".txt");
    if (opt::physicalWindowSize != -1) {
        rm->calculateMapForFixedWindowSizes(opt::physicalWindowSize);
        rm->outputMapToFileFixedWindowSizes("recombMap" + opt::runName + "_FW_" + numToString(opt::physicalWindowSize) + ".txt", opt::physicalWindowSize);
    }
    
    // Optional: calculate and print coverage stats (effective coverage and direct coverage) per het site
    if(opt::outputCoverageStats) rm->printPerHetCoverageStats("recombMap_wDepth" + opt::runName + ".txt");
    
    
    if (opt::nBootstrap > 0) {
        std::cout << std::endl;
        std::cout << "5a) Bootstrap... " << std::endl;
     //   std::cout << "Sample: 1" << std::endl;
       
    // Add the original map as the first sample
        rm->physicalWindowBootstraps.push_back(rm->physicalWindowR);
        
        vector<DefiningRecombInfo*> orginalInformativePairs = rp->allInformativePairs;
        for (int i = 0; i < opt::nBootstrap; i++) {
            RecombReadPairs* thisRp = new RecombReadPairs(*rp);
            thisRp->allInformativePairs = thisRp->getBootstrapSample(orginalInformativePairs);
            RecombMap* thisRm = new RecombMap(thisRp, opt::minCoverage, true);
            thisRm->firstUpdate();
    
            while (thisRm->delta > (opt::epsilon) && thisRm->EMiterationNum < opt::maxEMiterations) {
                thisRm->EMiteration(opt::v);
                //std::cout << "Map length = " << rm->mapLength << std::endl;
            }
            std::cout << "DONE.... Map length = " << thisRm->mapLength << std::endl;
            thisRm->calculateMapForFixedWindowSizes(opt::physicalWindowSize);
            thisRm->physicalWindowR.resize(rm->physicalWindowR.size(),NAN);
            rm->physicalWindowBootstraps.push_back(thisRm->physicalWindowR);
            std::cout << "Sample: " << i + 1 << std::endl;
            delete thisRm; delete thisRp;
        }
        rm->outputBootstrapToFilePhysical("bootstrap" + opt::runName + ".txt");
    }
     
    return 0;
}

void parseRecombFromSAMOptions(int argc, char** argv) {
    bool die = false;
    vector<string> nDiscordConcord;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 'q': arg >> opt::minBQ; break;
            case 'm': arg >> opt::maxEMiterations; break;
            case 'd': arg >> opt::minDist; break;
            case 'p': arg >> opt::minPQ; break;
            case 's': arg >> opt::hetsSubset; break;
            case 'f': arg >> opt::physicalWindowSize; break;
            case 'c': opt::outputCoverageStats = true; break;
            case 'r': opt::outputReadPairInfo = true; break;
            case 'i': opt::bIntermediateMaps = true; break;
            case 'v': opt::v = true; break;
            case 'b': arg >> opt::nBootstrap; break;
            case 'e': arg >> opt::epsilon; break;
            case 'x': arg >> opt::minCoverage; break;
            case OPT_DOUBLE_CROSS:
                arg >> opt::minDistForDoubleCrossovers;
                if (opt::minDistForDoubleCrossovers < 0) {
                    std::cerr << "Error: the --minDforDoubleCrossovers parameter should be non-negative\n";
                    die = true;
                }
                break;
            case OPT_IMBALANCED_SNPs:
                nDiscordConcord = split(arg.str(), ',');
                if (nDiscordConcord.size() != 2) {
                    std::cerr << "Error: the --imbalancedSNPs parameter requires two comma-separated values: <nDiscordant>,<nConcordant>\n";
                    die = true;
                }
                int nDiscordant, nConcordant;
                nDiscordant = safeStringToInt(nDiscordConcord[0]);
                if (nDiscordant < 0) {
                    std::cerr << "Error: the --imbalancedSNPs parameter requires a non-negative number of discordant SNPs\n";
                    die = true;
                }
                nConcordant = safeStringToInt(nDiscordConcord[1]);
                if (nConcordant < 0) {
                    std::cerr << "Error: the --imbalancedSNPs parameter requires a non-negative number of concordant SNPs\n";
                    die = true;
                }
                if (nDiscordant == 0 && nConcordant == 0) {
                    std::cerr << "Error: the --imbalancedSNPs parameter requires at least one of the two numbers to be non-zero\n";
                    die = true;
                }
                opt::minDiscordantSNPs = nDiscordant;
                opt::maxConcordantSNPs = nConcordant;
                break;
            case OPT_FPRATE_THRESHOLDS:
                {
                    vector<string> thresholds = split(arg.str(), ',');
                    if (thresholds.size() != 2) {
                        std::cerr << "Error: the --FPrateThresholds parameter requires two comma-separated values: <shortThreshold>,<longThreshold>\n";
                        die = true;
                    }
                    opt::FPrateShortThreshold = safeStringToInt(thresholds[0]);
                    if (opt::FPrateShortThreshold < 0) {
                        std::cerr << "Error: the --FPrateThresholds parameter requires a non-negative short threshold\n";
                        die = true;
                    }
                    opt::FPrateLongThreshold = safeStringToInt(thresholds[1]);
                    if (opt::FPrateLongThreshold < 0) {
                        std::cerr << "Error: the --FPrateThresholds parameter requires a non-negative long threshold\n";
                        die = true;
                    }
                    if (opt::FPrateShortThreshold >= opt::FPrateLongThreshold) {
                        std::cerr << "Error: the --FPrateThresholds parameter requires the short threshold to be less than the long threshold\n";
                        die = true;
                    }
                    if (opt::FPrateShortThreshold >= 4000) {
                        std::cerr << "Error: the --FPrateThresholds parameter requires the short threshold to be less than 5000 bp\n";
                        die = true;
                    }
                    if (opt::FPrateLongThreshold <= 200000) {
                        std::cerr << "Error: the --FPrateThresholds parameter requires the long threshold to be more than 200000 bp\n";
                        die = true;
                    }
                }
                break;
            case 'h':
                std::cout << DISCORDPAIRS_USAGE_MESSAGE;
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
        std::cout << "\n" << DISCORDPAIRS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hetsFile = argv[optind++];
    opt::samFile = argv[optind++];
}
