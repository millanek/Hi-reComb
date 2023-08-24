
CXXFLAGS=-std=c++11
CXX=g++
BIN := Build
LDFLAGS=-lz

all: $(BIN)/Hi-reComb

$(BIN)/Hi-reComb: $(BIN) $(BIN)/generalUtils.o $(BIN)/recombFromInformativePairsSAM.o $(BIN)/findInformativePairs.o $(BIN)/recombUtils.o $(BIN)/Hi-reComb.o $(BIN)/gzstream.o $(BIN)/countMendelianViolations.o $(BIN)/TrioPhase.o
	$(CXX) $(CXXFLAGS) $(BIN)/generalUtils.o $(BIN)/recombFromInformativePairsSAM.o $(BIN)/findInformativePairs.o $(BIN)/recombUtils.o $(BIN)/countMendelianViolations.o $(BIN)/Hi-reComb.o $(BIN)/gzstream.o $(BIN)/TrioPhase.o -o $@ $(LDFLAGS)

$(BIN)/%.o: %.cpp
		$(CXX) -c $(CXXFLAGS) $< -o $@

$(BIN):
	mkdir -p $@

clean:
	rm $(BIN)/*.o $(BIN)/Hi-reComb
