//
//  recombUtils.cpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 02.09.22.
//

#include "recombUtils.hpp"

void RecombRead::findHetsInMatchingString(const string& matchSeq, const int startPos, const map<int,PhaseInfo*>& positionToPhase) {
    for (int i = 0; i < matchSeq.length(); i++) {
        if (positionToPhase.count(startPos + i) == 1) {
            PhaseInfo* thisHetPhase = positionToPhase.at(startPos + i);
            vector<char> phasedSNPbases = thisHetPhase->phasedVars;
            char readBase = matchSeq[i];
            int snpPos = startPos + i;
            HetInfo* het = new HetInfo(snpPos, readBase, int(readQual[i])-33, phasedSNPbases[0], phasedSNPbases[1], thisHetPhase->quality, thisHetPhase->blockNum);
            hetSites.push_back(het);
        }
    }
}

void RecombRead::findHetsInRead(const map<int,PhaseInfo*>& positionToPhase) {
    
    int startPos = readPos; string processedReadSeq = readSeq;
    while (GIGARtypes.size() > 0) {
        switch (GIGARtypes[0])
        {
            case SOFT_CLIP_CIGAR:
                processedReadSeq = processedReadSeq.substr(GIGARnums[0]);
                break;
                
            case HARD_CLIP_CIGAR:
                break;
                
            case DELETION_CIGAR:
                startPos = startPos + GIGARnums[0];
                break;
                
            case INSERTION_CIGAR:
                processedReadSeq = processedReadSeq.substr(GIGARnums[0]);
                break;
                
            case MATCH_CIGAR:
                string matchSeq = processedReadSeq.substr(0, GIGARnums[0]);
                findHetsInMatchingString(matchSeq, startPos, positionToPhase);
                startPos = startPos + GIGARnums[0];
                usedLength = usedLength + GIGARnums[0];
                processedReadSeq = processedReadSeq.substr(GIGARnums[0]);
                
                break;
        }
        
        GIGARtypes.erase(GIGARtypes.begin());
        GIGARnums.erase(GIGARnums.begin());
        
    }
}

void RecombRead::linkHetsWithPhaseBlock() {
    for (int i = 0; i < hetSites.size(); i++) {
        int thisHetPhaseBlock = hetSites[i]->phaseBlock;
        if (BlockIDsToHetPos.count(thisHetPhaseBlock) == 1) {
            BlockIDsToHetPos.at(thisHetPhaseBlock).push_back(hetSites[i]->pos);
        } else {
            BlockIDsToHetPos[thisHetPhaseBlock].push_back(hetSites[i]->pos);
        }
    }
}



string RecombRead::assignStrandFromFlag() {
    string strand = "?";
    if (flag == 81 || flag == 113 || flag == 145 || flag == 177 || flag == 185 || flag == 121) {
        strand = "-";
    } else if (flag == 65 || flag == 73 || flag == 97 || flag == 129 || flag == 161 || flag == 137) {
        strand = "+";
    } else {
        std::cerr << "Unexpected read flag: " << flag << std::endl;
        exit(1);
    }
    return strand;
}

void RecombRead::generateCIGARvectors() {
    string CIGARnum = "";
    for (int i = 0; i < CIGAR.length(); i++) {
       // std::cout << "CIGAR[i]: " << CIGAR[i] << std::endl;
       // std::cout << "isdigit(CIGAR[i]): " << isdigit(CIGAR[i]) << std::endl;
        if (isdigit(CIGAR[i])) {
            CIGARnum += CIGAR[i];
         //   std::cout << "CIGARnum: " << CIGARnum << std::endl;
        } else {
           // std::cout << "CIGARnum: " << CIGARnum << std::endl;
            int CIGARnumInt = atoi(CIGARnum.c_str());
            GIGARnums.push_back(CIGARnumInt);
            if (CIGAR[i] != SOFT_CLIP_CIGAR && CIGAR[i] != INSERTION_CIGAR) GIGARnumsNoSI.push_back(CIGARnumInt);
            GIGARtypes.push_back(CIGAR[i]);
            CIGARnum = "";
        }
    }
}

void RecombReadPair::findAndCombinePairHets(const map<int,PhaseInfo*> & positionToPhase) {
    read1->findHetsInRead(positionToPhase);
    read2->findHetsInRead(positionToPhase);
    
    if (read1->hetSites.size() > 0) {
        hetSites = read1->hetSites;
    }
    
    if (read2->hetSites.size() > 0) {
        if (hetSites.size() == 0) hetSites = read2->hetSites;
        else hetSites.insert(hetSites.end(), read2->hetSites.begin(), read2->hetSites.end());
    }
}

void RecombReadPair::filterHetsByQuality(int minQuality) {
    
    vector<HetInfo*> goodHets;
    for (vector<HetInfo*>::iterator it = hetSites.begin(); it != hetSites.end(); it++) {
        if ((*it)->thisBaseQuality >= minQuality) {
            goodHets.push_back((*it));
        }
    }
    hetSites = goodHets;
}

void RecombReadPair::filterHetsByBlock(int blockNum) {
    
    vector<HetInfo*> goodHets;
    for (vector<HetInfo*>::iterator it = hetSites.begin(); it != hetSites.end(); it++) {
        if ((*it)->phaseBlock == blockNum) {
            goodHets.push_back((*it));
        }
    }
    hetSites = goodHets;
}

void RecombReadPair::findIndicesOfConcordantAndDiscordantPairsOfHets(const int minDistance) {
    for (int i = 0; i < hetSites.size(); ++i) {
        for (int j = i + 1; j < hetSites.size(); ++j) {
            if (hetSites[i]->phaseBlock == hetSites[j]->phaseBlock) {
                int phaseI = hetSites[i]->thisHetPhase01;
                int phaseJ = hetSites[j]->thisHetPhase01;
                int iPos = hetSites[i]->pos;
                int jPos = hetSites[j]->pos;
                if(abs(jPos - iPos) > minDistance) {
                    if (phaseI != phaseJ) {
                        switchPairI.push_back(i); switchPairJ.push_back(j);
                    } else {
                        concordPairI.push_back(i); concordPairJ.push_back(j);
                    }
                }
            }
        }
    }
}

void RecombReadPair::determineIfReadPairConcordantOrDiscordant() {
    if (switchPairI.size() > 0 && concordPairI.size() == 0)
        pairRecombinationStatus = PAIR_DISCORDANT;
    else if (concordPairI.size() > 0 && switchPairI.size() == 0)
        pairRecombinationStatus = PAIR_CONCORDANT;
    else if (concordPairI.size() == 0 && switchPairI.size() == 0)
        pairRecombinationStatus = PAIR_TOO_SHORT;
    else
        pairRecombinationStatus = PAIR_AMBIGUOUS;
}

DefiningRecombInfo* RecombReadPair::getDefiningHetPair(const vector<int>& indicesI, const vector<int>& indicesJ) {
    int maxD = 0; int maxDindex = 0;
    int minD = std::numeric_limits<int>::max(); int minDindex = 0;
    for (int i = 0; i != indicesI.size(); i++) {
        int iPos = hetSites[indicesI[0]]->pos;
        int jPos = hetSites[indicesJ[0]]->pos;
        
        if (abs(jPos - iPos) > maxD) maxDindex = i;
        if (abs(jPos - iPos) < minD) minDindex = i;
    }
    
    DefiningRecombInfo* thisRecombInfo;
    if (pairRecombinationStatus == PAIR_CONCORDANT) {
        thisRecombInfo = initialiseRecombInfo(indicesI, indicesJ, maxDindex);
    } else if (pairRecombinationStatus == PAIR_DISCORDANT) {
        thisRecombInfo = initialiseRecombInfo(indicesI, indicesJ, minDindex);
    } else {
        assert(false);
    }
    return thisRecombInfo;
}

DefiningRecombInfo* RecombReadPair::initialiseRecombInfo(const vector<int>& indicesI, const vector<int>& indicesJ, const int index) {
    
    int iPos = hetSites[indicesI[index]]->pos;
    int jPos = hetSites[indicesJ[index]]->pos;
    
    DefiningRecombInfo* recombInfo;
    if (jPos - iPos > 0) {
        recombInfo = new DefiningRecombInfo(hetSites[indicesI[index]], hetSites[indicesJ[index]], pairRecombinationStatus);
    } else {
        recombInfo = new DefiningRecombInfo(hetSites[indicesJ[index]], hetSites[indicesI[index]], pairRecombinationStatus);
    }
    return recombInfo;
}

double transformFromPhred(const double phredScore) {
    double p = pow(10,-(phredScore/10.0));
    return p;
}
