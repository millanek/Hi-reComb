//
//  findInformativePairs.cpp
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
// 73 -     Paired, mate unmapped, first in pair, + strand
// 137 -    Paired, mate unmapped, second in pair, + strand
// 185 -    Paired, mate unmapped, mate - strand, second in pair, - strand
// 121 -    Paired, mate unmapped, mate - strand, first in pair, - strand

#include "findInformativePairs.hpp"
#include <unordered_set>

#define SUBPROGRAM "FindInfoPairs"

#define DEBUG 1

static const char *INFOREADS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] HAPCUT2_PHASE.txt\n"
"Select reads from SAMTOOLS_FILE which could be informative about recombination\n"
"Expects SAM input on STDIN, e.g.,:\n"
"samtools view ALIGNEMENT.bam | Hi-reComb FindInfoPairs HAPCUT2_PHASE.txt > INFORMATIVE_READS.sam\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -m, --min-MQ                            (default: 20) the minimum mapping quality for a read to be considered\n"
"       --hapCut                                the het positions come from HapCut output"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* shortopts = "hn:m:";

enum { OPT_HAPCUT  };

static const struct option longopts[] = {
    { "help",   no_argument, NULL, 'h' },
    { "min-MQ",   required_argument, NULL, 'm' },
    { "run-name",   required_argument, NULL, 'n' },
    { "hapCut",   no_argument, NULL, OPT_HAPCUT },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string hetsFile;
    static bool hapcutFormat = false;
    static string runName = "";
    static int minMQ = 20;
}




int InfoReadsMain(int argc, char** argv) {
    parseInfoReadsOptions(argc, argv);
    string line; // for reading the input files
    
    std::istream* samtoolsFile = &std::cin;
    
    std::map<string,std::vector<string> > readNameToSamRecords;
    
    AllPhaseInfo* p = new AllPhaseInfo(opt::hetsFile, 0);
    
    std::cerr << "Finished reading het sites. There are " << p->posToPhase.size() << " hets." << std::endl;
    std::cerr << "Processing reads:" << std::endl;
    
    int readsProcessed = 0;
    // Now parse the pairtools file to find read pairs that can be informative about the phasing and recombination
    while (getline(*samtoolsFile,line)) {
        readsProcessed++;
        if (readsProcessed % 100000 == 0) {
            std::cerr << "Processed " << readsProcessed << " reads" << std::endl;
        }
        // std::cerr << line << std::endl;
        std::vector<string> samRecVec = split(line, '\t'); //assert(pairVec.size() == 8);
        int flag = atoi(samRecVec[1].c_str());
        if (flag > 2000) continue;
        
        int MQ = atoi(samRecVec[4].c_str());
        if (MQ < opt::minMQ) continue;
        
        
        RecombRead* thisRead = new RecombRead(samRecVec);
        
        thisRead->findHetsInRead(p->posToPhase);
        
        if (thisRead->hetSites.size() > 0) {
            readNameToSamRecords[thisRead->readName].push_back(line);
        }
        delete thisRead;
    }
    
    for (std::map<string,std::vector<string>>::iterator it = readNameToSamRecords.begin(); it != readNameToSamRecords.end(); it++) {
        if (it->second.size() > 1) {
            std::cout << it->second[0] << std::endl;
            std::cout << it->second[1] << std::endl;
        }
    }
    
    return 0;
    
}



void parseInfoReadsOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case OPT_HAPCUT: opt::hapcutFormat = true; break;
            case 'm': arg >> opt::minMQ; break;
            case 'h':
                std::cout << INFOREADS_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 1)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << INFOREADS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::hetsFile = argv[optind++];
}
