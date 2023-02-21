
CXXFLAGS=-std=c++11
CXX=g++
BIN := Build
LDFLAGS=-lz

all: $(BIN)/Hi-reComb

$(BIN)/Hi-reComb: $(BIN) $(BIN)/generalUtils.o $(BIN)/recombFromInformativePairsSAM.o $(BIN)/findInformativePairs.o $(BIN)/recombUtils.o $(BIN)/Hi-reComb.o
	$(CXX) $(CXXFLAGS) $(BIN)/generalUtils.o $(BIN)/recombFromInformativePairsSAM.o $(BIN)/findInformativePairs.o $(BIN)/recombUtils.o $(BIN)/Hi-reComb.o -o $@ $(LDFLAGS)

$(BIN)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(BIN):
	mkdir -p $@

clean:
	rm $(BIN)/*.o $(BIN)/Hi-reComb
