//
//  recombUtils.hpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 02.09.22.
//

#ifndef recombUtils_hpp
#define recombUtils_hpp

#include <stdio.h>
#include "generalUtils.hpp"

#define SOFT_CLIP_CIGAR 'S'
#define HARD_CLIP_CIGAR 'H'
#define MATCH_CIGAR 'M'
#define INSERTION_CIGAR 'I'
#define DELETION_CIGAR 'D'

#define PAIR_CONCORDANT 0
#define PAIR_DISCORDANT 1
#define PAIR_AMBIGUOUS 2
#define PAIR_TOO_SHORT 3

double transformFromPhred(const double phredScore);

class HetInfo {
    public:
    HetInfo(int position, char readBase, int readBaseQuality, char phase0, char phase1, double phaseQuality, int phaseBlockIn) {
        pos = position;
        thisBase = readBase; thisBaseQuality = readBaseQuality;
        thisHetPhase0 = phase0; thisHetPhase1 = phase1;
        thisPhaseQuality = phaseQuality;
        phaseBlock = phaseBlockIn;
        
        if (thisBase != thisHetPhase0 && thisBase != thisHetPhase1) {
            readPhaseBaseMismatch = true;
        } else {
            readPhaseBaseMismatch = false;
        }
        
        if (thisBase == thisHetPhase0) {
            thisHetPhase01 = 0;
        } else if (thisBase == thisHetPhase1) {
            thisHetPhase01 = 1;
        }
        
    };
    
    int pos;
    int phaseBlock;
    
    char thisBase;
    int thisBaseQuality;
    
    char thisHetPhase0;
    char thisHetPhase1;
    double thisPhaseQuality;
    
    int thisHetPhase01;
    
    bool readPhaseBaseMismatch;

};


class PhaseInfo {
    public:
    PhaseInfo() {};
    
    PhaseInfo(vector<string>& phasedSNPdetails, int blockNumIn): valid(false) {
        int H1phase = atoi(phasedSNPdetails[1].c_str()); int H2phase = atoi(phasedSNPdetails[2].c_str());
        assert(phasedSNPdetails[5].length() == 1); assert(phasedSNPdetails[6].length() == 1);
        char refBase = phasedSNPdetails[5][0]; char altBase = phasedSNPdetails[6][0];
        if (H1phase == 0 && H2phase == 1) {
            phasedVars.push_back(refBase); phasedVars.push_back(altBase); valid = true;
        } else if (H1phase == 1 && H2phase == 0) {
            phasedVars.push_back(altBase); phasedVars.push_back(refBase); valid = true;
        }
        pos = atoi(phasedSNPdetails[4].c_str()); quality = stringToDouble(phasedSNPdetails[10].c_str());
        coverage = atoi(phasedSNPdetails[11].c_str());
        blockNum = blockNumIn;
    };
    
    int pos;
    vector<char> phasedVars;
    double quality;
    int coverage;
    int blockNum;
    bool valid;
    
};

class PhasePair {
    public:
    PhasePair(): phase1(NULL), phase2(NULL) {};
    
    PhaseInfo* phase1;
    PhaseInfo* phase2;
    
};

class ReadLinkSNPpair {
    public:
    ReadLinkSNPpair(int p1, int p2, char b1, char b2) {
        pos1 = p1; pos2 = p2;
        base1.push_back(b1); base2.push_back(b2);
    };
    
    int pos1; int pos2;
    vector<char> base1; vector<char> base2;
    
};




class RecombRead {
    public:
    RecombRead() {};
    
    RecombRead(const vector<string> samRecVec) : usedLength(0) {
        
        flag = atoi(samRecVec[1].c_str());
        readStrand = assignStrandFromFlag();
        readName = samRecVec[0]; readPos = atoi(samRecVec[3].c_str());
        MQ = atoi(samRecVec[4].c_str()); CIGAR = samRecVec[5];
        readSeq = samRecVec[9]; readQual = samRecVec[10];
        
        generateCIGARvectors();
    };
    
    bool operator< (const RecombRead &other) const {
            return readPos < other.readPos;
    }
    
    
    int flag;
    string readStrand; string readName; int readPos; int adjustedReadPos;
    string readSeq; string readQual;
    int MQ; string CIGAR;
    vector<int> GIGARnums; vector<char> GIGARtypes;
    vector<int> GIGARnumsNoSI;
    
    int usedLength;
    vector<HetInfo*> hetSites;
    map<int,vector<int>> BlockIDsToHetPos;
    
    void findHetsInRead(const map<int,PhaseInfo*>& positionToPhase);
    void linkHetsWithPhaseBlock();
    
    private:
        string assignStrandFromFlag();
        void generateCIGARvectors();
        void findHetsInMatchingString(const string& matchSeq, const int startPos, const map<int,PhaseInfo*>& positionToPhase);

};



class DefiningRecombInfo {
    public:
    
    DefiningRecombInfo(HetInfo* iHet, HetInfo* jHet, int pairRecombinationStatus): indexLeft(-1), indexRight(-1), sum_r_k(NAN), doubleCrossoverProbability(0) {
        posLeft = iHet->pos; posRight = jHet->pos; dist = abs(posRight - posLeft) + 1;
        phaseErrorP_left = transformFromPhred(iHet->thisPhaseQuality); phaseErrorP_right = transformFromPhred(jHet->thisPhaseQuality);
        baseErrorP_left = transformFromPhred(iHet->thisBaseQuality); baseErrorP_right = transformFromPhred(jHet->thisBaseQuality);
        isRecombined = pairRecombinationStatus;
        // Given haplotypes:
        // -------A------------------C--------------
        // -------T------------------G--------------
        // And a read pair called as:
        //        A------------------G
        // With BQ=20 and PQ=20 on both sides, we have the following probabilities:
        // p(ph1=T) = 0.99; p(ph1=F) = 0.01
        // p(ph2=T) = 0.99; p(ph2=F) = 0.01
        // p(b1=A) = 0.99; p(b1=T) = 0.01 / 3    Assuming equal probability of any base
        // p(b2=G) = 0.99; p(b2=C) = 0.01 / 3    Assuming equal probability of any base
        if (isRecombined) {
            // 1) p(ph1=T) * p(ph2=T) * p(b1=A) * p(b2=G)
            probabilityRecombined = (1 - phaseErrorP_left) * (1 - phaseErrorP_right) * (1 - baseErrorP_left) * (1 - baseErrorP_right);
            // 2) p(ph1=T) * p(ph2=T) * p(b1=T) * p(b2=C)
            probabilityRecombined += (1 - phaseErrorP_left) * (1 - phaseErrorP_right) * (baseErrorP_left / 3) * (baseErrorP_right / 3);
            // 3a) p(ph1=T) * p(ph2=F) * p(b1=A) * p(b2=C)
            probabilityRecombined += (1 - phaseErrorP_left) * phaseErrorP_right * (1 - baseErrorP_left) * (baseErrorP_right / 3);
            // 3b) p(ph1=F) * p(ph2=T) * p(b1=T) * p(b2=G)
            probabilityRecombined += phaseErrorP_left * (1 - phaseErrorP_right) * (baseErrorP_left / 3) * (1 - baseErrorP_right);
            // 4) p(ph1=F) * p(ph2=F) * p(b1=A) * p(b2=G)
            probabilityRecombined += phaseErrorP_left * phaseErrorP_right * (1 - baseErrorP_left) * (1 - baseErrorP_right);
            
            
            PrORgivenNR = 1 - probabilityRecombined;
            
            /*
            // Calculating the probability of observing recombination if there is no recombination
            // p(ph1=T) * p(ph2=T) * p(b1=A) * p(b2=C) -- truth is the read pair is A------------------C or A------------------T or A------------------A
            PrORgivenNR = (1 - phaseErrorP_left) * (1 - phaseErrorP_right) * (1 - baseErrorP_left) * (baseErrorP_right);
            // p(ph1=T) * p(ph2=T) * p(b1=T) * p(b2=G) -- truth is the read pair is T------------------G or C------------------G or G------------------G
            PrORgivenNR += (1 - phaseErrorP_left) * (1 - phaseErrorP_right) * (baseErrorP_left) * (1 - baseErrorP_right);
            // p(ph1=F) * p(ph2=T) * p(b1=A) * p(b2=G) -- truth is the read pair is A------------------G, but the left phase is wrong
            PrORgivenNR += phaseErrorP_left * (1 - phaseErrorP_right) * (1 - baseErrorP_left) * (1 - baseErrorP_right);
            // p(ph1=F) * p(ph2=T) * p(b1=A) * p(b2=G) -- truth is the read pair is A------------------G, but the right phase is wrong
            PrORgivenNR += (1 - phaseErrorP_left) * phaseErrorP_right * (1 - baseErrorP_left) * (1 - baseErrorP_right);
            if(probabilityRecombined + PrORgivenNR != 1) {
                std::cout << "probabilityRecombined = " << probabilityRecombined << "; PrORgivenNR = " << PrORgivenNR << std::endl;
                std::cout << "probabilityRecombined + PrORgivenNR = " << probabilityRecombined + PrORgivenNR << std::endl;
            }
             */
        } else {
            /*
            // With a read pair called as:
            //        A------------------C
            // We have the following probabilities:
            // p(ph1=T) * p(ph2=T) * p(b1=A) * p(b2=G)
            probabilityRecombined = (1 - phaseErrorP_left) * (1 - phaseErrorP_right) * (1 - baseErrorP_left) * (baseErrorP_right / 3);
            // p(ph1=T) * p(ph2=T) * p(b1=T) * p(b2=C)
            probabilityRecombined += (1 - phaseErrorP_left) * (1 - phaseErrorP_right) * (baseErrorP_left / 3) * (1 - baseErrorP_right);
            // p(ph1=F) * p(ph2=T) * p(b1=A) * p(b2=C)
            probabilityRecombined += phaseErrorP_left * (1 - phaseErrorP_right) * (1 - baseErrorP_left) * (1 - baseErrorP_right);
            // p(ph1=T) * p(ph2=F) * p(b1=A) * p(b2=C)
            probabilityRecombined += (1 - phaseErrorP_left) * phaseErrorP_right * (1 - baseErrorP_left) * (1 - baseErrorP_right);
            
            //probabilityRecombined = 0;
             */
        }
    }
    
    int posLeft; int posRight;
    double phaseErrorP_left; double phaseErrorP_right;
    double baseErrorP_left; double baseErrorP_right;
    int dist;
    int isRecombined;
    double probabilityRecombined;
    double PrORgivenNR;
    double sum_r_k;
    double doubleCrossoverProbability;
    
    // These are the indices in a sorted vector of all informative hets long the chromosome
    // it cannot be assigned in contruction but is filled in later; -1 is just a placeholder for missing data
    int indexLeft; int indexRight;
};




class RecombReadPair {
    public:
    RecombReadPair(RecombRead* r1, RecombRead* r2) {
        if (r1 < r2) {
            read1 = r1;
            read2 = r2;
        } else {
            read1 = r2;
            read2 = r1;
        }
    };
    
    bool operator< (const RecombReadPair &other) const {
            return read1->readPos < other.read1->readPos;
    }
    
    
    
    RecombRead* read1;
    RecombRead* read2;
    
    int pairSpan;
    int pairRecombinationStatus;
    vector<HetInfo*> hetSites;
    
    vector<int> switchPairI; vector<int> switchPairJ;
    vector<int> concordPairI; vector<int> concordPairJ;
    
    
    //map<int,vector<int>> BlockIDsToHetPos;
    //map<int> HetPosToBlockIDs;
    
    void findAndCombinePairHets(const map<int,PhaseInfo*> & positionToPhase);
    void filterHetsByQuality(int minQuality);
    void filterHetsByBlock(int blockNum);
    void findIndicesOfConcordantAndDiscordantPairsOfHets(const int minDistance);
    void determineIfReadPairConcordantOrDiscordant();
    DefiningRecombInfo* getDefiningHetPair(const vector<int>& indicesI, const vector<int>& indicesJ);
    
private:
    DefiningRecombInfo* initialiseRecombInfo(const vector<int>& indicesI, const vector<int>& indicesJ, const int index);
    
};

class AllPhaseInfo {
    public:
    AllPhaseInfo(const string& hapcutFileName, int minPhaseQual, string subsetFileName = "") {
        string line;
        
        std::ifstream* hapcutFile = new std::ifstream(hapcutFileName.c_str()); assertFileOpen(*hapcutFile, hapcutFileName);
        if (!subsetFileName.empty()) {
            std::ifstream* subsetFile = new std::ifstream(subsetFileName.c_str()); assertFileOpen(*subsetFile, subsetFileName);
            while (getline(*subsetFile, line)) {
                subsetLoci[atoi(line.c_str())] = true;
            }
        }
        int blockNum = 0;
        // Parse the Hapcut blocks file
        while (getline(*hapcutFile, line)) {
            if (line[0] == '*') {
            
            } else if (line[0] == 'B' && line[1] == 'L') { // New block - separating the hets by blocks
                blockNum++; phaseBlockSNPnums.push_back(0);
            } else {
                vector<string> phasedSNPdetails = split(line, '\t');
                PhaseInfo* thisPhase = new PhaseInfo(phasedSNPdetails,blockNum);
                if (thisPhase->valid && thisPhase->coverage > 2
                    //&& (subsetLoci.size() == 0 || subsetLoci.count(thisPhase->pos) == 0)
                    && thisPhase->quality >= minPhaseQual) {
                    posToPhase[thisPhase->pos] = thisPhase;
                    phaseBlockSNPnums[blockNum-1]++;
                }
            }
        }
        hapcutFile->close();
    }
    
    map<int,PhaseInfo*> posToPhase;
    map<int,bool> subsetLoci;
    vector<int> phaseBlockSNPnums;
};

template <typename T> int roundToNearestValue(T num, int roundingValue)
{
    int d_i = num;
    int halfRoundingValue = roundingValue/2;
    return ((d_i % roundingValue) < halfRoundingValue) ? d_i - (d_i % roundingValue) : d_i + (roundingValue - (d_i % roundingValue));
}


#endif /* recombUtils_hpp */
