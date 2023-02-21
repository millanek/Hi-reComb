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
        
    std::ofstream* phaseSwitchFile = new std::ofstream("switches" + opt::runName + ".txt");
    std::ofstream* recombFile = new std::ofstream("recombMap" + opt::runName + ".txt");
    std::ofstream* depthFile = new std::ofstream("recombMap_wDepth" + opt::runName + ".txt");
    
    std::cout << "1) Processing hets..." << std::endl;
    AllPhaseInfo* p = new AllPhaseInfo(opt::hetsFile);
    
    std::cout << std::endl;
    std::cout << "2) Loading read-pairs... " << std::endl;
    RecombReadPairs* rp = new RecombReadPairs(opt::samFile);
    
    std::cout << std::endl;
    std::cout << "3) Linking read-pairs and hets... " << std::endl;
    rp->linkWithHets(p->posToPhase, opt::minBQ);
    rp->printReadPairStats();
    
    std::cout << "4) Categorising concordant-discordant read-pairs... " << std::endl;

    std::vector<double> concordantBaseScores; std::vector<double> discordantBaseScores;
    std::vector<PhaseSwitch*> phaseSwitches;
    std::vector<std::vector<int>> phaseConcordanceCoords;
        
    std::vector<int> coveredHetPos; std::vector<int> coveredHetEffectiveDepth;
    std::vector<int> coveredHetDirectDepth;

    int readPairsProcessed = 0;
    for (int r = 0; r < rp->informativeReadPairs.size(); r++) {
        readPairsProcessed++;
        
       for (int i = 0; i < rp->informativeReadPairs[r]->hetSites.size(); i++) {
            HetInfo* thisHet = rp->informativeReadPairs[r]->hetSites[i];
            if (thisHet->readPhaseBaseMismatch) {
                rp->numMismatch++;
                rp->mismatchBaseScores.push_back(thisHet->thisBaseQuality);
            } else {
                // std::cout << "Fine: " << std::endl;
                rp->matchBaseScores.push_back(thisHet->thisBaseQuality);
                rp->numMatch++;
                coveredHetPos.push_back(thisHet->pos);
            }
        }
        
        // Now consider all HiC pairs i that map to the chromosome and cover a het at each end - lets call the physical locations of the hets x_i, y_i ( x < y) and let Z_i = 1 if they are in phase, 0 if out of phase (recombined)
        
        // Prod_i ((1-Z_i)(G(y_i)-G(x_i)) + Z_i*(1 -(G(y_i)-G(x_i))
        
        // std::vector<PhaseSwitch*> thisPairSwitches;
        std::vector<int> switchPairI; std::vector<int> switchPairJ;
        std::vector<int> concordPairI; std::vector<int> concordPairJ;
        for (int i = 0; i < rp->informativeReadPairs[r]->hetSites.size() - 1; i++) {
            for (int j = 1; j < rp->informativeReadPairs[r]->hetSites.size(); j++) {
                if (rp->informativeReadPairs[r]->hetSites[i]->phaseBlock == rp->informativeReadPairs[r]->hetSites[j]->phaseBlock) {
                    int phaseI = rp->informativeReadPairs[r]->hetSites[i]->thisHetPhase01;
                    int phaseJ = rp->informativeReadPairs[r]->hetSites[j]->thisHetPhase01;
                    int iPos = rp->informativeReadPairs[r]->hetSites[i]->pos;
                    int jPos = rp->informativeReadPairs[r]->hetSites[j]->pos;
                    if (phaseI != phaseJ && abs(jPos - iPos) > opt::minDist) {
                        switchPairI.push_back(i); switchPairJ.push_back(j);
                        //int iPos = informativeReadPairs[r]->hetSites[i]->pos;
                        //int jPos = informativeReadPairs[r]->hetSites[j]->pos;
                        //int iQual = informativeReadPairs[r]->hetSites[i]->thisPhaseQuality;
                        //int jQual = informativeReadPairs[r]->hetSites[j]->thisPhaseQuality;
                        // PhaseSwitch* thisSwitch = new PhaseSwitch(iPos, jPos, iQual, jQual);
                        // thisPairSwitches.push_back(thisSwitch);
                    } else {
                        concordPairI.push_back(i); concordPairJ.push_back(j);
                    }
                }
            }
        }
        
        // TO DO:
        // Select the 'right' switch pair if there are multiple options:
        // The shortest one? Needs more thought....
        if (switchPairI.size() > 0) {
            rp->numDiscordant++;
            int iPos = rp->informativeReadPairs[r]->hetSites[switchPairI[0]]->pos;
            int jPos = rp->informativeReadPairs[r]->hetSites[switchPairJ[0]]->pos;
            int iQual = rp->informativeReadPairs[r]->hetSites[switchPairI[0]]->thisPhaseQuality;
            int jQual = rp->informativeReadPairs[r]->hetSites[switchPairJ[0]]->thisPhaseQuality;
            if (jPos - iPos < 0) {
                int tmp = iPos; iPos = jPos; jPos = tmp;
                tmp = iQual; iQual = jQual; jQual = tmp;
            }
            PhaseSwitch* thisSwitch = new PhaseSwitch(iPos, jPos, iQual, jQual);
            phaseSwitches.push_back(thisSwitch);
            switchPairI.empty(); switchPairJ.empty();
            rp->totalEffectiveLength = rp->totalEffectiveLength + (jPos - iPos);
            
        } else {
            rp->numConcordant++;
            std::vector<int> thisConcordantCoords;
            int maxD = 0; int maxDindex = 0;
            for (int i = 0; i != concordPairI.size(); i++) {
                int iPos = rp->informativeReadPairs[r]->hetSites[concordPairI[i]]->pos;
                int jPos = rp->informativeReadPairs[r]->hetSites[concordPairJ[i]]->pos;
                if (abs(jPos - iPos) > maxD) {
                    maxDindex = i;
                }
            }
            int iPosDindex = rp->informativeReadPairs[r]->hetSites[concordPairI[maxDindex]]->pos;
            int jPosDindex = rp->informativeReadPairs[r]->hetSites[concordPairJ[maxDindex]]->pos;
            if (jPosDindex - iPosDindex < 0) {
                int tmp = iPosDindex; iPosDindex = jPosDindex; jPosDindex = tmp;
            }
            
            int concordantPairDist = jPosDindex - iPosDindex;
            rp->totalEffectiveLength = rp->totalEffectiveLength + concordantPairDist;
            thisConcordantCoords.push_back(iPosDindex);
            thisConcordantCoords.push_back(jPosDindex);
            thisConcordantCoords.push_back(concordantPairDist);
            phaseConcordanceCoords.push_back(thisConcordantCoords);
        }
      /*  if (readPairsProcessed % 10000 == 0) {
            std::cout << "readPairsProcessed: " << readPairsProcessed << std::endl;
           // std::cout << "informativeReadPairs[r]->hetSites.size(): " << informativeReadPairs[r]->hetSites.size() << std::endl;
            std::cout << "phaseSwitches.size(): " << phaseSwitches.size() << std::endl;
            std::cout << "Effective coverage (bp): " << totalEffectiveLength << std::endl;
            std::cout << std::endl;
        } */
    }
    
    rp->printRecombinationBaseStats();
    std::cout << "phaseConcordanceCoords.size(): " << phaseConcordanceCoords.size() << std::endl;
    
    rp->printBaseQualityStats();
    
    for (int i = 0; i != phaseSwitches.size(); i++) {
        *phaseSwitchFile << phaseSwitches[i]->posLeft << "\t" << phaseSwitches[i]->posRight << "\t" << phaseSwitches[i]->dist << "\t" << phaseSwitches[i]->phaseQualLeft << "\t" << phaseSwitches[i]->phaseQualRight << std::endl;
    }
    
    
    std::cout << std::endl;
    std::cout << "5) Making a genetic map... " << std::endl;
    std::sort(coveredHetPos.begin(), coveredHetPos.end());
    std::vector<int>::iterator it = std::unique(coveredHetPos.begin(), coveredHetPos.end());
    coveredHetPos.resize(distance(coveredHetPos.begin(),it));
    std::cout << "coveredHetPos.size() " << coveredHetPos.size() << std::endl;
    
    double meanRecombinationRate = (double)rp->numDiscordant/(double)rp->totalEffectiveLength;
    std::cout << "meanRecombinationRate " << meanRecombinationRate << std::endl;
    std::vector<double> recombFractions(coveredHetPos.size()+1, meanRecombinationRate);
    
    int numProcessedHets = 0;
    for (int i = 0; i != coveredHetPos.size() - 1; i++) {
        int left = coveredHetPos[i];
        int right = coveredHetPos[i + 1];
        //int distSNPs = (right - left) + 1;
        
        int coveringReadPairs = 0;
        
        double totalRecombFractionPerBP = 0;
        for (int j = 0; j != phaseSwitches.size(); j++) {
            if(phaseSwitches[j]->posLeft <= left && phaseSwitches[j]->posRight >= right){
                coveringReadPairs++;
                double recombFractionPerBP = (double)1.0/(double)phaseSwitches[j]->dist;
                totalRecombFractionPerBP += recombFractionPerBP;
            }
        }
      //  std::cout << "totalRecombFraction: " << totalRecombFraction << std::endl;
        
        double totalConcordantFraction = 0;
        for (int j = 0; j != phaseConcordanceCoords.size(); j++) {
            if(phaseConcordanceCoords[j][0] <= left && phaseConcordanceCoords[j][1] >= right){
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
            std::cout << "numProcessedHets: " << numProcessedHets << " ("<< (double)numProcessedHets/coveredHetPos.size() << "%)"<< std::endl;
            std::cout << "pos: " << left << "bp"<< std::endl;
            
        }
        coveredHetEffectiveDepth.push_back(coveringReadPairs);
    }
    
    *recombFile << "0\t" << coveredHetPos[0] << "\t" << recombFractions[0] << std::endl;
    for (int i = 1; i != coveredHetPos.size(); i++) {
        *recombFile << coveredHetPos[i-1] << "\t" << coveredHetPos[i] << "\t" << recombFractions[i] << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "6) Collecting stats... " << std::endl;
    
    // a) Prints out the distribition of recombination rates based on read pairs with variable-length inserts
    collectRateStatsBasedOnInsertLength(phaseSwitches, phaseConcordanceCoords);
    
    int coverageDirectlyOnTheHet = 0;
    for (int i = 0; i != coveredHetPos.size() - 1; i++) {
        int thisHetPos = coveredHetPos[i];
        for (int j = 0; j != phaseSwitches.size(); j++) {
            if(phaseSwitches[j]->posLeft == thisHetPos || phaseSwitches[j]->posRight == thisHetPos){
                coverageDirectlyOnTheHet++;
            }
        }
        for (int j = 0; j != phaseConcordanceCoords.size(); j++) {
            if(phaseConcordanceCoords[j][0] == thisHetPos || phaseConcordanceCoords[j][1] == thisHetPos){
                coverageDirectlyOnTheHet++;
            }
        }
        coveredHetDirectDepth.push_back(coverageDirectlyOnTheHet);
        *depthFile << thisHetPos << "\t" << recombFractions[i] << "\t" << coveredHetEffectiveDepth[i] << "\t" << coverageDirectlyOnTheHet << std::endl;
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

void collectRateStatsBasedOnInsertLength(const std::vector<PhaseSwitch*>& phaseSwitches, const std::vector<std::vector<int>>& phaseConcordanceCoords) {
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
        int l = phaseConcordanceCoords[j][1] - phaseConcordanceCoords[j][0] + 1;
        
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
