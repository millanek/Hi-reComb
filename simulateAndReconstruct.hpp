//
//  simulateAndReconstruct.hpp
//  Hi-reComb
//
//  Created by Milan Malinsky on 12.09.23.
//

#ifndef simulateAndReconstruct_hpp
#define simulateAndReconstruct_hpp

#include <stdio.h>
#include "recombUtils.hpp"
#include "recombFromInformativePairsSAM.hpp"

void parseSimulationOptions(int argc, char** argv);
int SimulationMain(int argc, char** argv);


class RecombMapForSimulation : public RecombMapBase {
public:
    RecombMapForSimulation(string& mapFileName) {
        
        std::ifstream* mapFile = new std::ifstream(mapFileName.c_str()); assertFileOpen(*mapFile, mapFileName);
        string line; int mapRecordN = 0;
        while (getline(*mapFile,line)) {
            // std::cerr << line << std::endl;
            mapRecordN++;
            vector<string> mapRecVec = split(line, '\t');
            if(mapRecVec.size() != 3) {
                std::cerr << "ERROR: The recombination map input file " << mapFileName << " must have three tab-separated columns." << std::endl;
                std::cerr << "This condition failed on line: " << mapRecordN << ". Exiting...." << std::endl;
                exit(1);
            }
            RecombIntervalBase* thisInterval = new RecombIntervalBase(mapRecordN, atoi(mapRecVec[0].c_str()), atoi(mapRecVec[1].c_str()), stringToDouble(mapRecVec[2]));
            recombIntervals.push_back(*thisInterval);
        }
        
        // Filling these helper variables for the getAverageRateForPhysicalWindow function
        intervalCoordsVectors.resize(2);
        intervalCoordsVectors[0].resize(recombIntervals.size()); intervalCoordsVectors[1].resize(recombIntervals.size());
        intervalPerBPrVector.resize(recombIntervals.size(),0.0);
        for (int i = 0; i != recombIntervals.size(); i++) {
            intervalCoordsVectors[0][i] = recombIntervals[i].leftCoord;
            intervalCoordsVectors[1][i] = recombIntervals[i].rightCoord;
            intervalPerBPrVector[i] = recombIntervals[i].recombFractionPerBp;
        }
        
        
        // Get the mean recombination rate per bp
        //meanRate = rp->stats->meanRate;
        //std::cout << "Mean Recombination Rate = " << meanRate << std::endl;
        
        cummulativeRates.resize(mapRecordN - 1);
        
        mapPhysicalStart = recombIntervals.front().leftCoord;
        mapPhysicalEnd = recombIntervals.back().rightCoord;
        mapPhysicalLength = (mapPhysicalEnd - mapPhysicalStart) + 1;

    };
    
    vector<RecombIntervalBase> recombIntervals;
    
    void simuateReadPairs(RecombReadPairs* rp, vector<int>& readDistances, long long int targetCoverage, double errorRate) {
        std::random_device rd; // obtain a random number from hardware
        std::mt19937 gen(rd()); std::mt19937 gen2(rd()); std::mt19937 gen3(rd()); std::mt19937 gen4(rd()); // seed the generator
        std::uniform_int_distribution<> midPosDist(mapPhysicalStart, mapPhysicalEnd); // define the range
        std::uniform_int_distribution<> lengthIndexDist(0, (int)readDistances.size() - 1); // define the range
        std::uniform_real_distribution<> rUnif(0, 1); // define the range

        long long int totalSimulatedEffectiveLength = 0; long long int totalTargetEffectiveLength = (targetCoverage * mapPhysicalLength);
        //  std::cout << "totalSimulatedEffectiveLength: " << totalSimulatedEffectiveLength << std::endl;
        //  std::cout << "rsMapLength: " << rsMapLength << std::endl;
        //  std::cout << "totalTargetEffectiveLength: " << totalTargetEffectiveLength << std::endl;
        while (totalSimulatedEffectiveLength <= totalTargetEffectiveLength ) {
            int midPos = midPosDist(gen);
            int thisLength = readDistances[lengthIndexDist(gen2)];
            int startPos = midPos - (thisLength/2); int endPos = midPos + (thisLength/2);
            if (startPos < mapPhysicalStart || endPos > mapPhysicalEnd) continue;
            
            totalSimulatedEffectiveLength += thisLength;
            
            // std::cout << "thisLength: " << thisLength << std::endl;
            
            
            // Determine crossover probability based on the supplied (reference) genetic map
            double meanRate = getAverageRateForPhysicalWindow(startPos, startPos + thisLength);
            
            
            double crossoverProb = meanRate * thisLength;
            //  std::cout << std::endl;
            
            
            // First read with one het
            RecombRead* r1 = new RecombRead(); HetInfo* h1 = new HetInfo(startPos, 'A', 1000, 'A', 'C', 1000, 1); r1->hetSites.push_back(h1);
            
            // Second read
            RecombRead* r2 = new RecombRead(); HetInfo* h2;
            // Need to decide if concordant or discordant
            double cpRandom = rUnif(gen3); double errorRandom = rUnif(gen4);
            /*
            std::cout << "startPos: " << startPos << std::endl;
            std::cout << "startPos + thisLength: " << startPos + thisLength << std::endl;
            std::cout << "thisLength: " << thisLength << std::endl;
            std::cout << "meanRate: " << meanRate << std::endl;
            std::cout << "crossoverProb: " << crossoverProb << std::endl;
            std::cout << "cpRandom: " << cpRandom << std::endl;
            std::cout << "errorRandom: " << errorRandom << std::endl;
            std::cout << "errorRate: " << errorRate << std::endl;
            exit(0);
            */
            if (crossoverProb >= cpRandom && errorRate < errorRandom)
                h2 = new HetInfo(startPos + thisLength, 'C', 1000, 'A', 'C', 1000, 1); // This would be discordant
            else if (crossoverProb >= cpRandom && errorRate >= errorRandom)
                h2 = new HetInfo(startPos + thisLength, 'A', 1000, 'A', 'C', 1000, 1); // This would be in phase - false negative
            else if (crossoverProb < cpRandom && errorRate < errorRandom)
                h2 = new HetInfo(startPos + thisLength, 'A', 1000, 'A', 'C', 1000, 1); // This would be in phase
            else if (crossoverProb < cpRandom && errorRate >= errorRandom)
                h2 = new HetInfo(startPos + thisLength, 'C', 1000, 'A', 'C', 1000, 1); // This would be discordant - false positive
            
            r2->hetSites.push_back(h2);
            
            RecombReadPair* thisRp = new RecombReadPair(r1, r2);
            thisRp->hetSites.push_back(r1->hetSites[0]); thisRp->hetSites.push_back(r2->hetSites[0]);
            thisRp->pairSpan = thisLength;
            
            rp->informativeReadPairs.push_back(thisRp);
            rp->coveredHetPos.push_back(startPos); rp->coveredHetPos.push_back(startPos + thisLength);
        }
    }
    
    void calculateMapForFixedWindowSizes(int windowSize) {
        
        for (int i = 0; i != recombIntervals.size(); i++) {
            intervalCoordsVectors[0][i] = recombIntervals[i].leftCoord;
            intervalCoordsVectors[1][i] = recombIntervals[i].rightCoord;
            intervalPerBPrVector[i] = recombIntervals[i].recombFractionPerBp;
           // std::cout << "intervalPerBPrVector[i]: " << intervalPerBPrVector[i] << std::endl;
        }
        
        int lastCoord = intervalCoordsVectors[1].back();
        int maxCoordToOutput = roundToNearestValue(lastCoord, windowSize);
        for (int start = 1; start < maxCoordToOutput; start = start + windowSize) {
            int physicalWindowEnd = start + windowSize - 1;
            physicalWindowStartEnd[0].push_back(start); physicalWindowStartEnd[1].push_back(physicalWindowEnd);
            double avg = getAverageRateForPhysicalWindow(start, physicalWindowEnd);
            physicalWindowR.push_back(avg);
        }
    }
    
};

class RecombMapReconstructedFromSimulations : public RecombMap {
public:
    
    vector<vector<double>> simulationResults2kb;
};

#endif /* simulateAndReconstruct_hpp */
