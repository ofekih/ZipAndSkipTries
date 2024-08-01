CXX=c++
#CXXFLAGS=-std=c++17 -march=native -O0 -ggdb -g
CXXFLAGS=-std=c++17 -march=native -O3 -w
# CXXFLAGS=-std=c++17 -march=native -g -Wall
CPPFLAGS=-DNDEBUG
LDFLAGS=
#LIBS=-lcityhash

SOURCE=CTriePP.hpp AlphabetAwareZFastTrie.hpp LongString.hpp Fast.hpp any.hpp Factory.hpp CuckooHash.hpp Makefile
MAIN=main.cpp
TEST_MAIN=test.cpp

BIN=c_trei_pp
TEST_BIN=test

# TEXT_PATH=~/data/hoge.txt
# TEXT_PATH=~/data/cere.10KB.txt
# TEXT_PATH=~/data/einstein.en.50MB.txt

# TEXT_PATH=~/data/proteins.50MB
# TEXT_PATH=~/data/proteins.100MB
# TEXT_PATH=~/data/proteins.200MB
# TEXT_PATH=~/data/dna.50MB
# TEXT_PATH=~/data/dna.100MB
# TEXT_PATH=~/data/dna.200MB
# TEXT_PATH=~/data/english.50MB
# TEXT_PATH=~/data/english.100MB
# TEXT_PATH=~/data/english.200MB
# TEXT_PATH=~/data/dblp.xml.50MB
# TEXT_PATH=~/data/dblp.xml.100MB
# TEXT_PATH=~/data/dblp.xml.200MB
TEXT_PATH=~/data/geographic_names
# TEXT_PATH=~/data/proteins.distinct
# TEXT_PATH=~/data/pitches
# TEXT_PATH=~/data/enwiki-20181220-all-titles-in-ns0

all: run

run: ${BIN} test_ok
	@echo "  [RUN] start"
	@./${BIN} ${TEXT_PATH}
	@echo "  [RUN] end"

test_ok : ${TEST_BIN}
	@echo "  [TEST] start"
	@./${TEST_BIN}
	@echo "OK" > test_ok
	@echo "  [TEST] OK"

${BIN}:	${MAIN} ${SOURCE}
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) ${LIBS} -o $(BIN) ${MAIN}

${TEST_BIN}: ${TEST_MAIN} ${SOURCE}
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) ${LIBS} -o $(TEST_BIN) ${TEST_MAIN}

clean:
	@echo " [CLN] Removing object files"
	@rm -f $(BIN) ${TEST_BIN} *~

query : make_query
	./make_query ~/data/proteins.200MB ~/data/proteins.200MB_query
	./make_query ~/data/dna.200MB ~/data/dna.200MB_query
	./make_query ~/data/english.200MB ~/data/english.200MB_query
	./make_query ~/data/dblp.xml.200MB ~/data/dblp.xml.200MB_query

make_query : make_query.cpp
	$(CPP) $(CPPFLAGS) $(LDFLAGS) $(FLAGS) ${DEBUGFLAGS} -o make_query make_query.cpp
