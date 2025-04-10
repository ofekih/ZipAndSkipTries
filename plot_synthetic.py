#!/usr/bin/env python3
"""
Plotting script for synthetic benchmark data.

This script generates plots from benchmark data collected using synthetic_benchmark.cu.
It creates various plots comparing the performance of different trie implementations
on synthetic data with controlled properties, focusing on construction and search times.

The script reads CSV files from the data-synthetic directory and generates plots in the
figures directory. It supports plotting construction time and search time against
different parameters (number of keys, key length, LCP length).

The plots include trend lines with log-log regression to analyze asymptotic behavior.
"""

from enum import Enum, auto, unique
import matplotlib.pyplot as plt

from pathlib import Path
from collections import defaultdict
from functools import cache

from typing import NamedTuple

import numpy as np

import csv
import socket
import warnings

DATA_DIRECTORY = Path('data-synthetic')
HOSTNAME = socket.gethostname()
MAX_WORD_LENGTH = 2**22
CONSTRUCTION_FILENAME = 'construction-data'
REMOVAL_FILENAME = 'removal-data'
SEARCH_FILENAME = 'search-data'
START_PARALLEL = 1

FIGURE_DIRECTORY = Path('figures')
FIGURE_DIRECTORY.mkdir(parents=True, exist_ok=True)

OKABE_COLORS = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=OKABE_COLORS) # type: ignore
plt.rcParams['text.usetex'] = True

def get_dpi(savefig: bool):
	"""
	Set the DPI (dots per inch) for figures based on whether they will be saved or displayed.

	Args:
		savefig: If True, use higher DPI (600) for saved figures, otherwise use lower DPI (300) for display
	"""
	plt.rcParams['figure.dpi'] = 600 if savefig else 300

# Plot configuration
POINT_SIZE = 1  # Default size for data points in plots

@unique
class DataType(Enum):
	"""
	Enumeration of benchmark data types.

	Used to specify which type of benchmark data to load and plot.
	"""
	CONSTRUCTION = auto()  # Construction time benchmarks
	REMOVAL = auto()       # Removal time benchmarks
	SEARCH = auto()        # Search time benchmarks

def get_file(method_name: str, data_type: DataType) -> Path:
	"""
	Get the file path for benchmark data based on method name and data type.

	Args:
		method_name: Name of the method/implementation being benchmarked
		data_type: Type of benchmark data (construction, removal, or search)

	Returns:
		Path object pointing to the CSV file containing the benchmark data
	"""

	match data_type:
		case DataType.CONSTRUCTION:
			return DATA_DIRECTORY / f'{HOSTNAME}-{CONSTRUCTION_FILENAME}-{method_name}.csv'
		case DataType.REMOVAL:
			return DATA_DIRECTORY / f'{HOSTNAME}-{REMOVAL_FILENAME}-{method_name}.csv'
		case DataType.SEARCH:
			return DATA_DIRECTORY / f'{HOSTNAME}-{SEARCH_FILENAME}-{method_name}.csv'

class FitLine(NamedTuple):
	"""
	Represents a fitted trend line for log-log plots.

	Stores the parameters of a linear regression on log-transformed data,
	as well as the points to plot and the R² value indicating goodness of fit.
	"""
	m: float              # Slope of the fitted line in log-log space
	b: float              # Y-intercept of the fitted line in log-log space
	x: tuple[int, ...]    # X coordinates for plotting the trend line
	y: tuple[float, ...]  # Y coordinates for plotting the trend line
	r2: float             # R² value indicating goodness of fit

def binary_search(x: list[int], key: int) -> int:
	"""
	Perform binary search to find the insertion point for key in a sorted list.

	Used to find the index where data points should be skipped when fitting trend lines.

	Args:
		x: Sorted list of integers to search in
		key: Value to find the insertion point for

	Returns:
		Index where key should be inserted to maintain sorted order
	"""

	low = 0
	high = len(x)

	while low < high:
		mid = (high + low) // 2

		if x[mid] < key:
			low = mid + 1
		else:
			high = mid

	return low

def fit_log_line(x: list[int], y: list[float], skip_until: int = 0, all_points: bool = False, add_next: list[int] = []) -> FitLine | None:
	"""
	Fit a line to log-transformed data points for trend line analysis.

	Performs linear regression on log2-transformed data to analyze asymptotic behavior.
	Can skip initial data points to focus on asymptotic behavior at scale.

	Args:
		x: List of x-coordinates (typically problem sizes)
		y: List of y-coordinates (typically running times)
		skip_until: Skip data points with x values less than this threshold
		all_points: If True, include all points in the result; if False, only include endpoints
		add_next: Additional x values to extrapolate to

	Returns:
		A FitLine object containing the fitted line parameters, or None if fitting failed
	"""

	if skip_until > 0:
		skip = binary_search(x, skip_until)
	else:
		skip = 0

	if skip == len(x):
		warnings.warn(f'The trend line starts after the max data, so it is not shown')
		return None

	logx, logy = np.log2(x[skip:]), np.log2(y[skip:])
	m, b = np.polyfit(logx, logy, 1)

	fit = np.poly1d((m, b))

	for next_x in add_next:
		next_logx = np.log2(next_x)
		x.append(next_x)
		logx = np.append(logx, next_logx)


	expected_y = fit(logx)
	average_y = np.sum(logy) / len(logy)
	r2 = np.sum((expected_y[:len(expected_y) - len(add_next)] - average_y) ** 2) / np.sum((logy - average_y) ** 2)

	return FitLine(m = m, b = b,
		x = tuple(x[skip::] if all_points else x[skip::len(x) - skip - 1]),
		y = tuple(2 ** (expected_y if all_points else expected_y[::len(expected_y) - 1])),
		r2 = r2)

def plot(figure_name: str, save: bool, ylabel: str, xlabel: str = 'Num nodes (n)', legend_loc: str = 'best'):
	"""
	Finalize and display or save a plot.

	Adds labels, legend, and either displays the plot or saves it to a file.

	Args:
		figure_name: Base name for the saved figure file (without extension)
		save: If True, save the figure to a file; if False, display it
		ylabel: Label for the y-axis
		xlabel: Label for the x-axis (defaults to 'Num nodes (n)')
		legend_loc: Position for the legend (defaults to 'best')
	"""

	plt.xlabel(xlabel) # type: ignore
	plt.ylabel(ylabel) # type: ignore

	plt.legend(loc = legend_loc) # type: ignore
	if save:
		path = (FIGURE_DIRECTORY / figure_name).with_suffix('.png')
		plt.savefig(path) # type: ignore
	else:
		plt.show() # type: ignore

@cache
def load_agg_data(data_structure: str, data_type: DataType) -> dict[tuple[int, int, int], tuple[int, int]]:
	"""
	Load and aggregate benchmark data from a CSV file.

	Reads benchmark data for a specific method and data type, aggregating duplicate entries.

	Args:
		data_structure: Name of the data structure/implementation being benchmarked
		data_type: Type of benchmark data (construction, removal, or search)

	Returns:
		Dictionary mapping (n, m, l) tuples to (total_time, repetition_count) tuples,
		where n is number of keys, m is key length, and l is LCP length
	"""

	data = defaultdict[tuple[int, int, int], tuple[int, int]](lambda: (0, 0))

	with open(get_file(data_structure, data_type), newline='') as file:
		reader = csv.reader(file)

		for row in reader:
			if data_type == DataType.REMOVAL:
				n, m, time_ns, num_repetitions, par = map(int, row)
				l = m
			else:
				n, m, l, time_ns, num_repetitions, par = map(int, row)

			if not par in (0, START_PARALLEL):
				continue

			entry = data[(n, m, l)]
			data[(n, m, l)] = (entry[0] + time_ns, entry[1] + num_repetitions)

	return data

@cache
def load_data(data_structure: str, data_type: DataType, match: tuple[int, int, int]) -> tuple[tuple[int], tuple[float]]:
	"""
	Extract a specific data series from the benchmark data.

	Loads aggregated data and extracts a specific parameter series based on the match tuple.

	Args:
		data_structure: Name of the data structure/implementation being benchmarked
		data_type: Type of benchmark data (construction, removal, or search)
		match: Tuple of (n, m, l) where -1 indicates the parameter to vary
		       For example, (-1, 1024, 1024) means vary n with m=1024 and l=1024

	Returns:
		Tuple of (x values, y values) sorted by x, where x is the selected parameter and y is elapsed time
	"""

	x: list[int] = []
	y: list[float] = []

	match_n, match_m, match_l = match

	for (n, m, l), (time_ns, num_repetitions) in load_agg_data(data_structure, data_type).items():
		if match_n == -1 and match_m == m and match_l == l:
			x.append(n)
			y.append(time_ns / num_repetitions)
		elif match_n == n and match_m == -1 and match_l == l:
			x.append(m)
			y.append(time_ns / num_repetitions)
		elif match_n == n and match_m == m and match_l == -1:
			x.append(l)
			y.append(time_ns / num_repetitions)

	A, B = zip(*sorted(zip(x, y)))
	return A, B

@cache
def moving_average(x: list[int], y: list[float], window: float = 1.9) -> tuple[list[int], list[float]]:
	"""
	Apply a moving average smoothing to data points.

	First consolidates duplicate x values, then applies a geometric window smoothing
	to reduce noise in the data while preserving trends.

	Args:
		x: List of x-coordinates
		y: List of y-coordinates
		window: Window size factor (points within x/window to x*window are averaged)

	Returns:
		Tuple of (smoothed x values, smoothed y values)
	"""

	new_x: list[int] = []
	new_y: list[float] = []

	i = 0
	while i < len(x):
		total_sum = y[i]
		total_count = 1
		while i < len(x) - 1 and x[i] == x[i + 1]:
			total_sum += y[i + 1]
			total_count += 1
			i += 1

		new_x.append(x[i])
		new_y.append(total_sum / total_count)
		i += 1

	x, y = new_x, new_y

	# Convert x and y to numpy arrays for efficient computation
	x_array = np.array(x)
	y_array = np.array(y)

	# Prepare an array to store the moving averages
	moving_averages = np.zeros_like(y_array)

	# Calculate the moving average for each element in x
	for i, x_val in enumerate(x_array):
		# Define the range for the moving average
		lower_bound = x_val / window
		upper_bound = x_val * window

		# Get indices of x values that fall within the defined window
		window_indices = np.where((x_array >= lower_bound) & (x_array <= upper_bound))[0]

		# Calculate the average of y values within this window
		moving_averages[i] = np.mean(y_array[window_indices])

	# Return the original x values and the calculated moving averages
	return x, moving_averages.tolist()

def add_data_point_to_plot(method_name: str, data_type: DataType, match: tuple[int, int, int], skip_until: int = 0, fitlabel: str = '', max_value: int = -1, markersize: float = POINT_SIZE, draw_best_fit: bool = True, label: str = ''):
	"""
	Add a data series to the current plot with optional trend line.

	Loads data for a specific method and parameter combination, applies smoothing, and adds it to the plot.
	Optionally fits and adds a trend line with equation and R² value in the legend.

	Args:
		method_name: Name of the method/implementation being benchmarked
		data_type: Type of benchmark data (construction, removal, or search)
		match: Tuple of (n, m, l) where -1 indicates the parameter to vary
		skip_until: Skip data points with x values less than this threshold for trend line fitting
		fitlabel: Label to use in the trend line equation (e.g., 'n', 'm', '\\ell')
		max_value: Maximum x value to include (-1 for no limit)
		markersize: Size of the markers for data points
		draw_best_fit: If True, fit and draw a trend line; if False, only plot data points
		label: Custom label for the data series (defaults to method_name if empty)
	"""

	x, y = load_data(method_name, data_type, match)

	if max_value > 0:
		x, y = zip(*[(x[i], y[i]) for i in range(len(x)) if x[i] <= max_value])

	x, y = moving_average(x, y) # type: ignore

	label = label if label else method_name

	fit = fit_log_line(x, y, skip_until=skip_until)

	if draw_best_fit and fit:
		p = plt.loglog(x, y, base=2, marker = 'o', markersize = markersize, linestyle = 'None') # type: ignore
		pcolor = p[0].get_color() # type: ignore

		sign = '+' if fit.b >= 0 else '-'

		plt.loglog(fit.x, fit.y, base=2, color=pcolor, linestyle='--', # type: ignore
			label=f'{label}: $T({fitlabel}) = {fit.m:.1f} {fitlabel} {sign} {abs(fit.b):.1f}, R^2 = {fit.r2:.2f}$')
	else:
		plt.loglog(x, y, base=2, marker = 'o', markersize = markersize, linestyle = 'None', label=label) # type: ignore

def compare_construction_data(savefig: bool = False):
	"""
	Create plots comparing construction time for different implementations.

	Generates two plots:
	1. Construction time vs. LCP length (ℓ)
	2. Construction time vs. key length (m)

	Some plots are commented out in the code but can be uncommented if needed.

	Args:
		savefig: If True, save the figures to files; if False, display them
	"""

	# plt.figure(num = 0, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	# add_data_point_to_plot('ctrie++', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='\\texttt{c-trie++}', skip_until=2**8)
	# # add_data_point_to_plot('skip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='ST', skip_until=2**8)
	# add_data_point_to_plot('memory-intensive-zip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='MI-ZT', skip_until=2**8)
	# add_data_point_to_plot('zip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='ZT', skip_until=2**8)
	# # add_data_point_to_plot('parallel-skip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='PST', skip_until=2**8)
	# add_data_point_to_plot('parallel-memory-intensive-zip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='MI-PZT', skip_until=2**8)
	# add_data_point_to_plot('parallel-zip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='PZT', skip_until=2**8)

	# plot(figure_name='construction-time-num-keys-comparison', save=savefig, ylabel='Time (ns)', xlabel='Number of Keys ($n$)')

	plt.figure(num = 1, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	skip_until = 2**14

	add_data_point_to_plot('ctrie++', DataType.CONSTRUCTION, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='\\texttt{c-trie++}', skip_until=skip_until, draw_best_fit=False)
	# add_data_point_to_plot('skip-trie', DataType.CONSTRUCTION, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='ST', skip_until=2**14)
	add_data_point_to_plot('memory-intensive-zip-trie', DataType.CONSTRUCTION, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='MI-ZT', skip_until=skip_until, draw_best_fit=False)
	add_data_point_to_plot('zip-trie', DataType.CONSTRUCTION, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='ZT', skip_until=skip_until, draw_best_fit=False)
	# add_data_point_to_plot('parallel-skip-trie', DataType.CONSTRUCTION, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='PST', skip_until=2**14)
	add_data_point_to_plot('parallel-memory-intensive-zip-trie', DataType.CONSTRUCTION, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='MI-PZT', skip_until=skip_until, draw_best_fit=False)
	add_data_point_to_plot('parallel-zip-trie', DataType.CONSTRUCTION, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='PZT', skip_until=skip_until, draw_best_fit=False)

	plot(figure_name='construction-time-lcp-comparison', save=savefig, ylabel='Time (ns)', xlabel='Mean LCP Length ($\\ell$)')

	plt.figure(num = 2, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	add_data_point_to_plot('ctrie++', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='\\texttt{c-trie++}', skip_until=2**11)
	# add_data_point_to_plot('skip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='ST', skip_until=2**11)
	add_data_point_to_plot('memory-intensive-zip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='MI-ZT', skip_until=2**11)
	add_data_point_to_plot('zip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='ZT', skip_until=2**11)
	# add_data_point_to_plot('parallel-skip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='PST', skip_until=2**11)
	add_data_point_to_plot('parallel-memory-intensive-zip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='MI-PZT', skip_until=2**11)
	add_data_point_to_plot('parallel-zip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='PZT', skip_until=2**11)

	plot(figure_name='construction-time-key-length-comparison', save=savefig, ylabel='Time (ns)', xlabel='Key Length ($m$)')

def compare_search_data(savefig: bool = False):
	"""
	Create plots comparing search time for different implementations.

	Generates two plots:
	1. Search time vs. key length (m)
	2. Search time vs. LCP length (ℓ)

	Some plots are commented out in the code but can be uncommented if needed.

	Args:
		savefig: If True, save the figures to files; if False, display them
	"""

	# plt.figure(num = 100, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	# add_data_point_to_plot('ctrie++', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='\\texttt{c-trie++}', skip_until=2**13)
	# # add_data_point_to_plot('skip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='ST', skip_until=2**13)
	# add_data_point_to_plot('memory-intensive-zip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='MI-ZT', skip_until=2**13)
	# add_data_point_to_plot('zip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='ZT', skip_until=2**13)
	# # add_data_point_to_plot('parallel-skip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='PST', skip_until=2**13)
	# add_data_point_to_plot('parallel-memory-intensive-zip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='MI-PZT', skip_until=2**13)
	# add_data_point_to_plot('parallel-zip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='PZT', skip_until=2**13)

	# plot(figure_name='search-time-num-keys-comparison', save=savefig, ylabel='Time (ns)', xlabel='Number of Keys ($n$)')

	plt.figure(num = 101, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	add_data_point_to_plot('ctrie++', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='\\texttt{c-trie++}', draw_best_fit=False)
	# add_data_point_to_plot('skip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='ST', draw_best_fit=False)
	add_data_point_to_plot('memory-intensive-zip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='MI-ZT', draw_best_fit=False)
	add_data_point_to_plot('zip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='ZT', draw_best_fit=False)
	# add_data_point_to_plot('parallel-skip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='PST', draw_best_fit=False)
	add_data_point_to_plot('parallel-memory-intensive-zip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='MI-PZT', draw_best_fit=False)
	add_data_point_to_plot('parallel-zip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='PZT', draw_best_fit=False)

	plot(figure_name='search-time-key-length-comparison', save=savefig, ylabel='Time (ns)', xlabel='Key Length ($m$)')

	plt.figure(num = 102, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	add_data_point_to_plot('ctrie++', DataType.SEARCH, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='\\texttt{c-trie++}', skip_until=2**14, draw_best_fit=False)
	# add_data_point_to_plot('skip-trie', DataType.SEARCH, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='ST', skip_until=2**14)
	add_data_point_to_plot('memory-intensive-zip-trie', DataType.SEARCH, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='MI-ZT', skip_until=2**14, draw_best_fit=False)
	add_data_point_to_plot('zip-trie', DataType.SEARCH, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='ZT', skip_until=2**14, draw_best_fit=False)
	# add_data_point_to_plot('parallel-skip-trie', DataType.SEARCH, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='PST', skip_until=2**14)
	add_data_point_to_plot('parallel-memory-intensive-zip-trie', DataType.SEARCH, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='MI-PZT', skip_until=2**14, draw_best_fit=False)
	add_data_point_to_plot('parallel-zip-trie', DataType.SEARCH, (2**10, MAX_WORD_LENGTH, -1), fitlabel='\\ell', label='PZT', skip_until=2**14, draw_best_fit=False)

	plot(figure_name='search-time-lcp-comparison', save=savefig, ylabel='Time (ns)', xlabel='Mean LCP Length ($\\ell$)')

if __name__ == '__main__':
	savefig = False
	# savefig = True

	compare_construction_data(savefig)
	compare_search_data(savefig)
