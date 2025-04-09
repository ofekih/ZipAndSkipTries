#!/usr/bin/env bash

DEFAULT_VALUE=$((2**10))
MAX_NUM_WORDS=$((2**19))
MAX_WORD_LENGTH=$((2**22))

NUM_REPETITIONS=${1:-10}

run_variable_lcp_benchmarks() {
	local num_words=$1
	local word_length=$2
	local num_repetitions=${3:-1000}

	local mean_lcp=4

	echo "Running variable LCP benchmarks..."

	while (( mean_lcp <= word_length )); do
		./bin/synthetic_benchmark $num_words $word_length $mean_lcp $num_repetitions

		mean_lcp=$(( mean_lcp * 2 ))
	done
}

run_variable_word_length_benchmarks() {
	local num_words=$1
	local max_word_length=$2
	local mean_lcp=$3
	local num_repetitions=${4:-1000}
	local word_length=4

	echo "Running variable word length benchmarks..."

	while (( word_length <= max_word_length )); do
		./bin/synthetic_benchmark $num_words $word_length $mean_lcp $num_repetitions

		word_length=$(( word_length * 2 ))
	done
}

run_variable_num_words_benchmarks() {
	local max_num_words=$1
	local word_length=$2
	local mean_lcp=$3
	local num_repetitions=${4:-1000}
	local num_words=4

	echo "Running variable num words benchmarks..."

	while (( num_words <= max_num_words )); do
		./bin/synthetic_benchmark $num_words $word_length $mean_lcp $num_repetitions

		num_words=$(( num_words * 2 ))
	done
}

for i in {1..1000}
do
	echo "Iteration $i"

	run_variable_lcp_benchmarks $DEFAULT_VALUE $MAX_WORD_LENGTH $NUM_REPETITIONS
	run_variable_word_length_benchmarks $DEFAULT_VALUE $MAX_WORD_LENGTH $DEFAULT_VALUE $NUM_REPETITIONS
	run_variable_num_words_benchmarks $MAX_NUM_WORDS $DEFAULT_VALUE $DEFAULT_VALUE $NUM_REPETITIONS

	echo "Sleeping for 60 seconds"
	sleep 60
done


