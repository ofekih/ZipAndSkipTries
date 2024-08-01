CXX=c++
CXXFLAGS=-std=c++17 -march=native -O3

all:

shuf:
	# shuf ~/data/commoncrawl.ascii.distinct -o ~/data/commoncrawl.ascii.distinct.shuf
	# shuf ~/data/all_vital.txt.distinct -o ~/data/all_vital.txt.distinct.shuf
	shuf ~/data/proteins.distinct -o ~/data/proteins.distinct.shuf.query
	shuf ~/data/uk-2002.urls.distinct -o ~/data/uk-2002.urls.distinct.shuf.query
	shuf ~/data/geographic_names -o ~/data/geographic_names.shuf.query
	shuf ~/data/dblp.xml.distinct -o ~/data/dblp.xml.distinct.shuf.query
	shuf ~/data/commoncrawl.ascii.distinct -o ~/data/commoncrawl.ascii.distinct.shuf.query
	shuf ~/data/all_vital.txt.distinct -o ~/data/all_vital.txt.distinct.shuf.query


query: make_query
	# ./make_query ~/data/proteins.distinct ~/data/proteins.distinct.query
	# ./make_query ~/data/uk-2002.urls.distinct ~/data/uk-2002.urls.distinct.query
	# ./make_query ~/data/geographic_names ~/data/geographic_names.query
	# ./make_query ~/data/dblp.xml.distinct ~/data/dblp.xml.distinct.query
	# ./make_query ~/data/commoncrawl.ascii.distinct ~/data/commoncrawl.ascii.distinct.query
	# ./make_query ~/data/all_vital.txt.distinct ~/data/all_vital.txt.distinct.query
	# ./make_query ~/data/proteins.distinct.shuf.query ~/data/proteins.distinct.shuf.query
	# ./make_query ~/data/uk-2002.urls.distinct.shuf.query ~/data/uk-2002.urls.distinct.shuf.query
	# ./make_query ~/data/dblp.xml.distinct.shuf.query ~/data/dblp.xml.distinct.shuf.query
	# ./make_query ~/data/geographic_names.shuf.query ~/data/geographic_names.shuf.query
	# ./make_query ~/data/commoncrawl.ascii.distinct.shuf.query ~/data/commoncrawl.ascii.distinct.shuf.query
	# ./make_query ~/data/all_vital.txt.distinct.shuf.query ~/data/all_vital.txt.distinct.shuf.query
	./make_query ~/data/proteins.distinct.shuf ~/data/proteins.distinct.shuf_s.query
	./make_query ~/data/uk-2002.urls.distinct.shuf ~/data/uk-2002.urls.distinct.shuf_s.query
	./make_query ~/data/dblp.xml.distinct.shuf ~/data/dblp.xml.distinct.shuf_s.query
	./make_query ~/data/geographic_names.shuf ~/data/geographic_names.shuf_s.query
	./make_query ~/data/commoncrawl.ascii.distinct.shuf ~/data/commoncrawl.ascii.distinct.shuf_s.query
	./make_query ~/data/all_vital.txt.distinct.shuf ~/data/all_vital.txt.distinct.shuf_s.query

info: text_info
	./text_info ~/data/proteins.distinct
	./text_info ~/data/uk-2002.urls.distinct
	./text_info ~/data/dblp.xml.distinct
	./text_info ~/data/geographic_names
	./text_info ~/data/commoncrawl.ascii.distinct
	./text_info ~/data/all_vital.txt.distinct

	# ./text_info ~/data/proteins.200MB
	# ./text_info ~/data/dna.200MB
	# ./text_info ~/data/english.200MB
	# ./text_info ~/data/dblp.xml.200MB

distinct: make_distinct
	# ./make_distinct ~/data/proteins ~/data/proteins.distinct
	# ./make_distinct ~/data/uk-2002.urls ~/data/uk-2002.urls.distinct
	# ./make_distinct ~/data/dblp.xml ~/data/dblp.xml.distinct
	# ./make_distinct ~/data/commoncrawl.ascii ~/data/commoncrawl.ascii.distinct
	# ./make_distinct ~/data/all_vital.txt ~/data/all_vital.txt.distinct

wseq :make_wseq
	./make_wseq ~/data/europarl-v7.es-en.en ~/data/word_sequences.en
	./make_wseq ~/data/europarl-v7.es-en.es ~/data/word_sequences.es

geo: make_geo
	./make_geo ~/data/allCountries.txt ~/data/geographic_names

make_wseq: make_wseq.cpp
	$(CXX) $(CXXFLAGS) make_wseq.cpp -o make_wseq

make_query: make_query.cpp
	$(CXX) $(CXXFLAGS) make_query.cpp -o make_query

text_info: text_info.cpp
	$(CXX) $(CXXFLAGS) text_info.cpp -o text_info

make_distinct: make_distinct.cpp
	$(CXX) $(CXXFLAGS) make_distinct.cpp -o make_distinct

make_geo: make_geo.cpp
	$(CXX) $(CXXFLAGS) make_geo.cpp -o make_geo

# run: compile
# 	./bench.sh

clean:
	@echo "  [CLEAN] Removing object files"
	@rm -f $(BIN) *~
