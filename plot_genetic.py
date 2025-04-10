#!/usr/bin/env python3
"""
Plotting script for genetic benchmark data.

This script generates plots from benchmark data collected using genetic_benchmark.cu.
It creates various plots comparing the performance of different trie implementations
on genetic data from the ABC-HuMi dataset, focusing on construction and search times.

The script reads CSV files from the data-genetic directory and generates plots in the
figures directory. It supports plotting construction time and search time against
different parameters (number of keys, total key length, LCP length).

The plots include trend lines with log-log regression to analyze asymptotic behavior.
"""

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from pathlib import Path

import numpy as np

from typing import NamedTuple
from enum import Enum, auto, unique

from functools import cache

import csv
import socket
import warnings

# Constants and configuration
HOSTNAME = socket.gethostname()

# Directory and file paths
DATA_DIRECTORY = Path('data-genetic')  # Directory containing genetic benchmark data
CONSTRUCTION_FILENAME = 'construction-data'  # Prefix for construction benchmark files
REMOVAL_FILENAME = 'removal-data'  # Prefix for removal benchmark files
SEARCH_FILENAME = 'search-data'  # Prefix for search benchmark files

# Create figures directory if it doesn't exist
FIGURE_DIRECTORY = Path('figures')
FIGURE_DIRECTORY.mkdir(exist_ok = True)

# Color-blind friendly palette (Okabe-Ito)
OKABE_COLORS = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
plt.rcParams['axes.prop_cycle'] = plt.cycler(color = OKABE_COLORS) # type: ignore

# Plot configuration
POINT_SIZE = 0.1  # Default size for data points in plots

def get_dpi(savefig: bool):
	"""Set the DPI (dots per inch) for figures based on whether they will be saved or displayed.

	Args:
		savefig: If True, use higher DPI (600) for saved figures, otherwise use lower DPI (300) for display
	"""
	plt.rcParams['figure.dpi'] = 600 if savefig else 300

@unique
class DataType(Enum):
	"""Enumeration of benchmark data types.

	Used to specify which type of benchmark data to load and plot.
	"""
	CONSTRUCTION = auto()  # Construction time benchmarks
	REMOVAL = auto()       # Removal time benchmarks
	SEARCH = auto()        # Search time benchmarks

class FitLine(NamedTuple):
	"""Represents a fitted trend line for log-log plots.

	Stores the parameters of a linear regression on log-transformed data,
	as well as the points to plot and the R² value indicating goodness of fit.
	"""
	m: float              # Slope of the fitted line in log-log space
	b: float              # Y-intercept of the fitted line in log-log space
	x: tuple[int, ...]    # X coordinates for plotting the trend line
	y: tuple[float, ...]  # Y coordinates for plotting the trend line
	r2: float             # R² value indicating goodness of fit

def binary_search(x: list[int], key: int) -> int:
	"""Perform binary search to find the insertion point for key in a sorted list.

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
	"""Fit a line to log-transformed data points for trend line analysis.

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

def get_file(method_name: str, data_type: DataType) -> Path:
	"""Get the file path for benchmark data based on method name and data type.

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

def plot(figure_name: str, save: bool, ylabel: str, xlabel: str, legend_loc: str | None = 'best'):
	"""Finalize and display or save a plot.

	Adds labels, legend, and either displays the plot or saves it to a file.

	Args:
		figure_name: Base name for the saved figure file (without extension)
		save: If True, save the figure to a file; if False, display it
		ylabel: Label for the y-axis
		xlabel: Label for the x-axis
		legend_loc: Position for the legend (e.g., 'best', 'upper left', etc.)
	"""
	plt.xlabel(xlabel) # type: ignore
	plt.ylabel(ylabel) # type: ignore

	if legend_loc:
		plt.legend(loc = legend_loc) # type: ignore

	if save:
		path = (FIGURE_DIRECTORY / figure_name).with_suffix('.png')
		plt.savefig(path) # type: ignore
	else:
		plt.show() # type: ignore

@cache
def load_data(method_name: str, data_type: DataType, cutoff: int = 0) -> tuple[list[int], list[int], list[int], list[float]]:
	"""Load benchmark data from a CSV file.

	Reads benchmark data for a specific method and data type, filtering by cutoff value.

	Args:
		method_name: Name of the method/implementation being benchmarked
		data_type: Type of benchmark data (construction, removal, or search)
		cutoff: Filter value for the parallel cutoff parameter

	Returns:
		Tuple containing lists of (n values, m values, l values, elapsed times)
		where n is number of keys, m is key length, l is LCP length, and elapsed times are in nanoseconds
	"""
	ns: list[int] = []
	ms: list[int] = []
	ls: list[int] = []
	es: list[float] = []

	with open(get_file(method_name, data_type), 'r') as f:
		reader = csv.reader(f)

		for row in reader:
			if int(row[5]) != cutoff:
				continue

			ns.append(int(row[0]))
			ms.append(int(row[1]))
			ls.append(int(row[2]))
			es.append(int(row[3]) / int(row[4]))

	return ns, ms, ls, es

@cache
def get_data_point(method_name: str, data_type: DataType, index: int, cutoff: int = 0) -> tuple[tuple[int], tuple[float]]:
	"""Extract a specific data series from the benchmark data.

	Loads data and extracts a specific parameter (based on index) and the corresponding times.

	Args:
		method_name: Name of the method/implementation being benchmarked
		data_type: Type of benchmark data (construction, removal, or search)
		index: Which parameter to extract (0 for n, 1 for m, 2 for l)
		cutoff: Filter value for the parallel cutoff parameter

	Returns:
		Tuple of (x values, y values) sorted by x, where x is the selected parameter and y is elapsed time
	"""
	data = load_data(method_name, data_type, cutoff)

	x, y = zip(*sorted(zip(data[index], data[3])))
	return x, y

def moving_average(x: list[int], y: list[float], window: float = 1.9) -> tuple[list[int], list[float]]:
	"""Apply a moving average smoothing to data points.

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

def add_data_point_to_plot(method_name: str, data_type: DataType, index: int, cutoff: int = 0, skip_until: int=0, fitlabel: str='', max_value: int=-1, markersize: float = POINT_SIZE, draw_best_fit: bool = True, label: str = ''):
	"""Add a data series to the current plot with optional trend line.

	Loads data for a specific method and parameter, applies smoothing, and adds it to the plot.
	Optionally fits and adds a trend line with equation and R² value in the legend.

	Args:
		method_name: Name of the method/implementation being benchmarked
		data_type: Type of benchmark data (construction, removal, or search)
		index: Which parameter to plot against (0 for n, 1 for m, 2 for l)
		cutoff: Filter value for the parallel cutoff parameter
		skip_until: Skip data points with x values less than this threshold for trend line fitting
		fitlabel: Label to use in the trend line equation (e.g., 'n', 'm', 'L')
		max_value: Maximum x value to include (-1 for no limit)
		markersize: Size of the markers for data points
		draw_best_fit: If True, fit and draw a trend line; if False, only plot data points
		label: Custom label for the data series (defaults to method_name if empty)
	"""
	x, y = get_data_point(method_name, data_type, index, cutoff)

	if max_value > 0:
		x, y = zip(*[(x[i], y[i]) for i in range(len(x)) if x[i] <= max_value])

	x, y = moving_average(list(x), list(y))

	label = label if label else method_name

	fit = fit_log_line(x, y, skip_until=skip_until)

	if draw_best_fit and fit:
		p = plt.loglog(x, y, base=2, marker = 'o', markersize = markersize, linestyle = 'None') # type: ignore
		pcolor = p[0].get_color() # type: ignore

		sign = '+' if fit.b >= 0 else '-'

		plt.loglog(fit.x, fit.y, base=2, color=pcolor, linestyle='--', alpha=0.7, linewidth=0.8, # type: ignore
			label=f'{label}: $T({fitlabel}) = {fit.m:.1f} {fitlabel} {sign} {abs(fit.b):.1f}, R^2 = {fit.r2:.2f}$')
	else:
		plt.loglog(x, y, base=2, marker = 'o', markersize = markersize, linestyle = 'None', label=label) # type: ignore



def compare_construction_data(savefig: bool):
	"""Create plots comparing construction time for different implementations.

	Generates three plots:
	1. Construction time vs. number of keys (n)
	2. Construction time vs. total key length (N)
	3. Construction time vs. total LCP length (L)

	Args:
		savefig: If True, save the figures to files; if False, display them
	"""
	plt.figure(num = 0, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	skip_until = 2 ** 8
	cutoff = 1

	add_data_point_to_plot('c-trie++', DataType.CONSTRUCTION, 0, fitlabel='n', skip_until=skip_until,markersize=1, label='\\texttt{c-trie++}')
	add_data_point_to_plot('ZT', DataType.CONSTRUCTION, 0, fitlabel='n', skip_until=skip_until,markersize=1)
	add_data_point_to_plot('MI-ZT', DataType.CONSTRUCTION, 0, fitlabel='n', skip_until=skip_until,markersize=1)
	add_data_point_to_plot('PZT', DataType.CONSTRUCTION, 0, cutoff=cutoff, fitlabel='n', skip_until=skip_until,markersize=1)
	add_data_point_to_plot('MI-PZT', DataType.CONSTRUCTION, 0, cutoff=cutoff, fitlabel='n', skip_until=skip_until,markersize=1)

	plot(figure_name='construction-time-num-keys-comparison', save=savefig, ylabel='Time (ns)', xlabel='Number of Keys ($n$)')


	plt.figure(num = 1, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	skip_until = 2**21

	add_data_point_to_plot('c-trie++', DataType.CONSTRUCTION, 1, fitlabel='N', skip_until=skip_until, label='\\texttt{c-trie++}')
	add_data_point_to_plot('ZT', DataType.CONSTRUCTION, 1, fitlabel='N', skip_until=skip_until)
	add_data_point_to_plot('MI-ZT', DataType.CONSTRUCTION, 1, fitlabel='N', skip_until=skip_until)
	add_data_point_to_plot('PZT', DataType.CONSTRUCTION, 1, cutoff=cutoff, fitlabel='N', skip_until=skip_until)
	add_data_point_to_plot('MI-PZT', DataType.CONSTRUCTION, 1, cutoff=cutoff, fitlabel='N', skip_until=skip_until)

	plot(figure_name='construction-time-total-key-length-comparison', save=savefig, ylabel='Time (ns)', xlabel='Total Key Length ($N$)')


	plt.figure(num = 2, figsize = (8, 5), dpi = get_dpi(savefig), facecolor = 'w', edgecolor = 'k') # type: ignore

	skip_until = 2**17

	add_data_point_to_plot('c-trie++', DataType.CONSTRUCTION, 2, fitlabel='L', skip_until=skip_until, label='\\texttt{c-trie++}')
	add_data_point_to_plot('ZT', DataType.CONSTRUCTION, 2, fitlabel='L', skip_until=skip_until)
	add_data_point_to_plot('MI-ZT', DataType.CONSTRUCTION, 2, fitlabel='L', skip_until=skip_until)
	add_data_point_to_plot('PZT', DataType.CONSTRUCTION, 2, cutoff=cutoff, fitlabel='L', skip_until=skip_until)
	add_data_point_to_plot('MI-PZT', DataType.CONSTRUCTION, 2, cutoff=cutoff, fitlabel='L', skip_until=skip_until)

	plot(figure_name='construction-time-total-LCP-length-comparison', save=savefig, ylabel='Time (ns)', xlabel='Total LCP Length ($L$)')

# def compare_removal_data(savefig: bool):
# 	plt.figure(num=200, figsize=(8, 5), dpi = get_dpi(savefig), facecolor='w', edgecolor='k') # type: ignore

# 	skip_until = 2**4
# 	cutoff = 2**13

# 	add_data_point_to_plot('c-trie++', DataType.REMOVAL, 0, fitlabel='n', skip_until=skip_until, markersize=1, label='\\texttt{c-trie++}')
# 	add_data_point_to_plot('ST', DataType.REMOVAL, 0, fitlabel='n', skip_until=skip_until, markersize=1)
# 	add_data_point_to_plot('PST', DataType.REMOVAL, 0, cutoff=cutoff, fitlabel='n', skip_until=skip_until, markersize=1)

# 	plot(figure_name='removal-time-n-comparison', save=savefig, ylabel='Time (ns)', xlabel='Number of Keys ($n$)')

# 	plt.figure(num=201, figsize=(8, 5), dpi = get_dpi(savefig), facecolor='w', edgecolor='k') # type: ignore

# 	skip_until = 2**18

# 	add_data_point_to_plot('c-trie++', DataType.REMOVAL, 1, fitlabel='N', skip_until=skip_until, label='\\texttt{c-trie++}')
# 	add_data_point_to_plot('ST', DataType.REMOVAL, 1, fitlabel='N', skip_until=skip_until)
# 	add_data_point_to_plot('PST', DataType.REMOVAL, 1, cutoff=cutoff, fitlabel='N', skip_until=skip_until)

# 	plot(figure_name='removal-time-N-comparison', save=savefig, ylabel='Time (ns)', xlabel='Total Key Length ($N$)')

def compare_search_data(savefig: bool):
	"""Create plots comparing search time for different implementations.

	Generates two plots:
	1. Search time vs. number of keys (n)
	2. Search time vs. LCP length (ℓ)

	Args:
		savefig: If True, save the figures to files; if False, display them
	"""
	plt.figure(num=300, figsize=(8, 5), dpi = get_dpi(savefig), facecolor='w', edgecolor='k') # type: ignore

	cutoff = 1

	add_data_point_to_plot('c-trie++', DataType.SEARCH, 0, draw_best_fit=False, markersize=4, label='\\texttt{c-trie++}')
	add_data_point_to_plot('ZT', DataType.SEARCH, 0, draw_best_fit=False, markersize=4)
	add_data_point_to_plot('MI-ZT', DataType.SEARCH, 0, draw_best_fit=False, markersize=4)
	add_data_point_to_plot('PZT', DataType.SEARCH, 0, cutoff=cutoff, draw_best_fit=False, markersize=4)
	add_data_point_to_plot('MI-PZT', DataType.SEARCH, 0, cutoff=cutoff, draw_best_fit=False, markersize=4)

	plot(figure_name='search-time-num-keys-comparison', save=savefig, ylabel='Time (ns)', xlabel='Number of Keys ($n$)')

	# plt.figure(num=301, figsize=(8, 5), dpi = get_dpi(savefig), facecolor='w', edgecolor='k') # type: ignore

	skip_until = 2**12

	# add_data_point_to_plot('c-trie++', DataType.SEARCH, 1, fitlabel='m', skip_until=skip_until)
	# add_data_point_to_plot('ZT', DataType.SEARCH, 1, fitlabel='m', skip_until=skip_until)
	# add_data_point_to_plot('MI-ZT', DataType.SEARCH, 1, fitlabel='m', skip_until=skip_until)
	# add_data_point_to_plot('PZT', DataType.SEARCH, 1, cutoff=cutoff, fitlabel='m', skip_until=skip_until)
	# add_data_point_to_plot('MI-PZT', DataType.SEARCH, 1, cutoff=cutoff, fitlabel='m', skip_until=skip_until)

	# plot(figure_name='search-time-m-comparison', save=savefig, ylabel='Time (ns)', xlabel='Key Length ($m$)')

	plt.figure(num=302, figsize=(8, 5), dpi = get_dpi(savefig), facecolor='w', edgecolor='k') # type: ignore

	add_data_point_to_plot('c-trie++', DataType.SEARCH, 2, fitlabel='\\ell', skip_until=skip_until, label='\\texttt{c-trie++}')
	add_data_point_to_plot('ZT', DataType.SEARCH, 2, fitlabel='\\ell', skip_until=skip_until)
	add_data_point_to_plot('MI-ZT', DataType.SEARCH, 2, fitlabel='\\ell', skip_until=skip_until)
	add_data_point_to_plot('PZT', DataType.SEARCH, 2, cutoff=cutoff, fitlabel='\\ell', skip_until=skip_until)
	add_data_point_to_plot('MI-PZT', DataType.SEARCH, 2, cutoff=cutoff, fitlabel='\\ell', skip_until=skip_until)

	plot(figure_name='search-time-lcp-length-comparison', save=savefig, ylabel='Time (ns)', xlabel='LCP Length ($\\ell$)')

# def compare_cutoffs(savefig: bool):
# 	plt.figure(num=400, figsize=(8, 5), dpi = get_dpi(savefig), facecolor='w', edgecolor='k') # type: ignore

# 	# add_data_point_to_plot('PST', DataType.SEARCH, 2, cutoff=1024, fitlabel='n', skip_until=2**4, markersize=1)
# 	# add_data_point_to_plot('PST', DataType.SEARCH, 2, cutoff=2048, fitlabel='n', skip_until=2**4, markersize=1)
# 	add_data_point_to_plot('ST', DataType.SEARCH, 2, fitlabel='\\ell', skip_until=2**15, markersize=1, label='$c=2^0$')
# 	add_data_point_to_plot('PST', DataType.SEARCH, 2, cutoff=2**13, fitlabel='\\ell', skip_until=2**12, markersize=1, label='$c=2^{13}$')
# 	add_data_point_to_plot('PST', DataType.SEARCH, 2, cutoff=2**14, fitlabel='\\ell', skip_until=2**12, markersize=1, label='$c=2^{14}$')
# 	add_data_point_to_plot('PST', DataType.SEARCH, 2, cutoff=2**15, fitlabel='\\ell', skip_until=2**12, markersize=1, label='$c=2^{15}$')
# 	add_data_point_to_plot('PST', DataType.SEARCH, 2, cutoff=2**16, fitlabel='\\ell', skip_until=2**12, markersize=1, label='$c=2^{16}$')

# 	plot(figure_name='construction-time-l-comparison-cutoff', save=savefig, ylabel='Time (ns)', xlabel='LCP Length ($\\ell$)')


if __name__ == '__main__':
	savefig = False
	# savefig = True

	compare_construction_data(savefig)
	# compare_removal_data(savefig) # Commented out removal benchmark
	compare_search_data(savefig)
	# compare_cutoffs(savefig)

