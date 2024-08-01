CXX=c++
CXXFLAGS=-std=c++17 -march=native -g -Wall
# CXXFLAGS=-std=c++17 -march=native -O3
CPPFLAGS=-DNDEBUG
LDFLAGS=
LIBS=-lcityhash

SOURCE = ZFastTrie.hpp BitString.hpp Fast.hpp any.hpp rigtorp/HashMap.h Makefile
MAIN=main.cpp
TEST_MAIN=test.cpp

BIN=ZFastTrie
TEST_BIN=Test

TEXTS=("proteins.distinct" "uk-2002.urls.distinct" "dblp.xml.distinct" "geographic_names" "commoncrawl.ascii.distinct" "all_vital.txt.distinct")
TEXT_DIR=../data
TEXT_PATH=${TEXT_DIR}/proteins.distinct
# TEXT_PATH=${TEXT_DIR}/uk-2002.urls.distinct
# TEXT_PATH=${TEXT_DIR}/dblp.xml.distinct
# TEXT_PATH=${TEXT_DIR}/geographic_names
# TEXT_PATH=${TEXT_DIR}/commoncrawl.ascii.distinct
# TEXT_PATH=${TEXT_DIR}/all_vital.txt.distinct


all: run

run: test_ok ${BIN}
	@echo "  [RUN] start"
	@./${BIN} ${TEXT_PATH} ${TEXT_PATH}.query
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
