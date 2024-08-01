CXX=c++

# GCC 8.2.0
CXXFLAGS=-std=c++1z -std=gnu++1z -march=native -O3
# clang-1000.11.45.2
#CXXFLAGS=-std=c++17 -march=native
# CXXFLAGS=-std=c++17 -march=native -O0 -ggdb
# CPPFLAGS=-g3 -O0 -ggdb
CPPFLAGS=-DNDEBUG
INCLUDES=-I./
# LDFLAGS=-L${HOME}/prefix_search/env/lib -L${HOME}/prefix_search/env/lib64
# LDFLAGS=-I ./CTriePP/external/separate_chaining/external/bit_span/include -I ./CTriePP/external/separate_chaining/external/tudostats/include
LIBS=-lcityhash

BIN=ctrie cedarpp_prefix tsl_htrie_map tsl_htrie_map_32 packed_ctrie_xor packed_ctrie_hash z_fast_trie


all: run

c_trei_pp: bench.cpp Makefile
	$(CXX) $(CXXFLAGS) ${CPPFLAGS} ${INCLUDES} ${LDFLAGS} bench.cpp -o c_trei_pp -DC_TRIE_PP

${BIN}: bench.cpp Makefile
	#$(CXX) $(CXXFLAGS) ${CPPFLAGS} bench.cpp -o ctrie -DC_TRIE
	#$(CXX) $(CXXFLAGS) ${CPPFLAGS} bench.cpp -o packed_ctrie_xor -DPACKED_C_TRIE_XOR
	#$(CXX) $(CXXFLAGS) ${CPPFLAGS} bench.cpp -o packed_ctrie_hash -DPACKED_C_TRIE_HASH
	#$(CXX) $(CXXFLAGS) ${CPPFLAGS} ${INCLUDES} ${LIBS} bench.cpp -o tsl_htrie_map -DBENCH_TSL_HTRIE_MAP
	#$(CXX) $(CXXFLAGS) ${CPPFLAGS} ${INCLUDES} ${LIBS} bench.cpp -o tsl_htrie_map_32 -DBENCH_TSL_HTRIE_MAP -DBENCH_KEY_SIZE_32
	#$(CXX) $(CXXFLAGS) ${CPPFLAGS} ${INCLUDES} ${LIBS} bench.cpp -o cedarpp_prefix                    -DBENCH_CEDARPP
	$(CXX) $(CXXFLAGS) ${CPPFLAGS} ${LIBS} -o z_fast_trie -DZ_FAST_TRIE bench.cpp

run: ${BIN} c_trei_pp
	./bench.sh
	# ./bench.sh > bench.log

clean:
	@echo "  [CLEAN] Removing object files"
	@rm -f $(BIN) c_trei_pp *~
