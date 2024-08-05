from enum import Enum, auto, unique
import matplotlib.pyplot as plt

from pathlib import Path
from collections import defaultdict
from functools import cache

from typing import NamedTuple

import numpy as np

import csv

DATA_DIRECTORY = Path('data')
HOSTNAME = ''
CONSTRUCTION_FILENAME = 'construction-data'
REMOVAL_FILENAME = 'removal-data'
SEARCH_FILENAME = 'search-data'

FIGURE_DIRECTORY = Path('figures')
FIGURE_DIRECTORY.mkdir(parents=True, exist_ok=True)

OKABE_COLORS = ['#000000', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7']
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=OKABE_COLORS) # type: ignore
plt.rcParams['text.usetex'] = True
DPI = 300

POINT_SIZE = 1

@unique
class DataType(Enum):
	CONSTRUCTION = auto()
	REMOVAL = auto()
	SEARCH = auto()

def get_file(method_name: str, data_type: DataType) -> Path:
	match data_type:
		case DataType.CONSTRUCTION:
			return DATA_DIRECTORY / f'{HOSTNAME}-{CONSTRUCTION_FILENAME}-{method_name}.csv'
		case DataType.REMOVAL:
			return DATA_DIRECTORY / f'{HOSTNAME}-{REMOVAL_FILENAME}-{method_name}.csv'
		case DataType.SEARCH:
			return DATA_DIRECTORY / f'{HOSTNAME}-{SEARCH_FILENAME}-{method_name}.csv'

class FitLine(NamedTuple):
	m: float
	b: float
	x: list[int]
	y: list[float]
	r2: float

def binary_search(x: list[int], key: int) -> int:
	low = 0
	high = len(x)

	while low < high:
		mid = (high + low) // 2

		if x[mid] < key:
			low = mid + 1
		else:
			high = mid

	return low

def fit_log_line(x: list[int], y: list[float], skip_until: int = 0, all_points: bool = False, add_next: list[int] = []) -> FitLine:
	if skip_until > 0:
		skip = binary_search(x, skip_until)
	else:
		skip = 0

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
		x = x[skip::] if all_points else x[skip::len(x) - skip - 1],
		y = 2 ** (expected_y if all_points else expected_y[::len(expected_y) - 1]), # type: ignore
		r2 = r2)
def plot(figure_name: str, save: bool, ylabel: str, xlabel: str = 'Num nodes (n)', legend_loc: str = 'best'):
	# plt.title(title) # type: ignore
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
	data = defaultdict[tuple[int, int, int], tuple[int, int]](lambda: (0, 0))

	with open(get_file(data_structure, data_type), newline='') as file:
		reader = csv.reader(file)

		for row in reader:
			if data_type == DataType.REMOVAL:
				n, m, time_ns, num_repetitions = map(int, row)
				l = m
			else:
				n, m, l, time_ns, num_repetitions = map(int, row)

			entry = data[(n, m, l)]
			data[(n, m, l)] = (entry[0] + time_ns, entry[1] + num_repetitions)

	return data

@cache
def load_data(data_structure: str, data_type: DataType, match: tuple[int, int, int]) -> tuple[list[int], list[float]]:
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

	return tuple(zip(*sorted(zip(x, y))))

@cache
def moving_average(x: list[int], y: list[float], window: float = 1.9) -> tuple[list[int], list[float]]:
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

def add_data_point_to_plot(method_name: str, data_type: DataType, match: tuple[int, int, int], skip_until: int = 0, fitlabel: str = '', maxvalue: int = -1, markersize: float = POINT_SIZE, draw_best_fit: bool = True, label: str = ''):
	x, y = load_data(method_name, data_type, match)

	if maxvalue > 0:
		x, y = zip(*[(x[i], y[i]) for i in range(len(x)) if x[i] <= maxvalue])

	x, y = moving_average(x, y) # type: ignore

	label = label if label else method_name

	if draw_best_fit:
		p = plt.loglog(x, y, base=2, marker = 'o', markersize = markersize, linestyle = 'None') # type: ignore
		pcolor = p[0].get_color()

		fit = fit_log_line(x, y, skip_until=skip_until)

		sign = '+' if fit.b >= 0 else '-'

		plt.loglog(fit.x, fit.y, base=2, color=pcolor, linestyle='--', # type: ignore
			label=f'{label}: $T({fitlabel}) = {fit.m:.1f} {fitlabel} {sign} {abs(fit.b):.1f}, R^2 = {fit.r2:.2f}$')
			# label=f'{label}: $T({fitlabel}) = {fit.m:.1f} \\cdot {fitlabel} + {fit.b:.1f}, R^2 = {fit.r2:.2f}$')
	else:
		plt.loglog(x, y, base=2, marker = 'o', markersize = markersize, linestyle = 'None', label=label) # type: ignore

def compare_construction_data(savefig: bool = False):
	# plt.figure(num = 0, figsize = (8, 5), dpi = DPI, facecolor = 'w', edgecolor = 'k') # type: ignore

	# add_data_point_to_plot('ctrie++', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='\\texttt{c-trie++}', draw_best_fit=False)
	# add_data_point_to_plot('skip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='ST', draw_best_fit=False)
	# add_data_point_to_plot('zip-zip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='ZZT', draw_best_fit=False)
	# add_data_point_to_plot('memory-efficient-zip-zip-trie', DataType.CONSTRUCTION, (-1, 2**10, 2**10), fitlabel='n', label='ME-ZZT', draw_best_fit=False)

	# plot(figure_name='construction-time-num-keys-comparison', save=savefig, ylabel='Time (ns)', xlabel='Number of Keys ($n$)')

	# plt.figure(num = 1, figsize = (8, 5), dpi = DPI, facecolor = 'w', edgecolor = 'k') # type: ignore

	# add_data_point_to_plot('ctrie++', DataType.CONSTRUCTION, (2**10, 2**22, -1), fitlabel='\\ell', label='\\texttt{c-trie++}', skip_until=2**14)
	# add_data_point_to_plot('skip-trie', DataType.CONSTRUCTION, (2**10, 2**22, -1), fitlabel='\\ell', label='ST', skip_until=2**14)
	# add_data_point_to_plot('zip-zip-trie', DataType.CONSTRUCTION, (2**10, 2**22, -1), fitlabel='\\ell', label='ZZT', skip_until=2**14)
	# add_data_point_to_plot('memory-efficient-zip-zip-trie', DataType.CONSTRUCTION, (2**10, 2**22, -1), fitlabel='\\ell', label='ME-ZZT', skip_until=2**14)

	# plot(figure_name='construction-time-lcp-comparison', save=savefig, ylabel='Time (ns)', xlabel='Mean LCP Length ($\\ell$)')

	plt.figure(num = 2, figsize = (8, 5), dpi = DPI, facecolor = 'w', edgecolor = 'k') # type: ignore

	add_data_point_to_plot('ctrie++', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='\\texttt{c-trie++}', skip_until=2**11)
	add_data_point_to_plot('skip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='ST', skip_until=2**11)
	add_data_point_to_plot('zip-zip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='ZZT', skip_until=2**11)
	add_data_point_to_plot('memory-efficient-zip-zip-trie', DataType.CONSTRUCTION, (2**10, -1, 2**10), fitlabel='m', label='ME-ZZT', skip_until=2**11)

	plot(figure_name='construction-time-key-length-comparison', save=savefig, ylabel='Time (ns)', xlabel='Key Length ($m$)')

def compare_search_data(savefig: bool = False):
	# plt.figure(num = 100, figsize = (8, 5), dpi = DPI, facecolor = 'w', edgecolor = 'k') # type: ignore

	# add_data_point_to_plot('ctrie++', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='\\texttt{c-trie++}', skip_until=2**13)
	# add_data_point_to_plot('skip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='ST', skip_until=2**13)
	# add_data_point_to_plot('zip-zip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='ZZT', skip_until=2**13)
	# add_data_point_to_plot('memory-efficient-zip-zip-trie', DataType.SEARCH, (-1, 2**10, 2**10), fitlabel='n', label='ME-ZZT', skip_until=2**13)

	# plot(figure_name='search-time-num-keys-comparison', save=savefig, ylabel='Time (ns)', xlabel='Number of Keys ($n$)')

	# plt.figure(num = 101, figsize = (8, 5), dpi = DPI, facecolor = 'w', edgecolor = 'k') # type: ignore

	# add_data_point_to_plot('ctrie++', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='\\texttt{c-trie++}', draw_best_fit=False)
	# add_data_point_to_plot('skip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='ST', draw_best_fit=False)
	# add_data_point_to_plot('zip-zip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='ZZT', draw_best_fit=False)
	# add_data_point_to_plot('memory-efficient-zip-zip-trie', DataType.SEARCH, (2**10, -1, 2**10), fitlabel='m', label='ME-ZZT', draw_best_fit=False)

	# plot(figure_name='search-time-key-length-comparison', save=savefig, ylabel='Time (ns)', xlabel='Key Length ($m$)')

	plt.figure(num = 102, figsize = (8, 5), dpi = DPI, facecolor = 'w', edgecolor = 'k') # type: ignore

	add_data_point_to_plot('ctrie++', DataType.SEARCH, (2**10, 2**22, -1), fitlabel='\\ell', label='\\texttt{c-trie++}', skip_until=2**14)
	add_data_point_to_plot('skip-trie', DataType.SEARCH, (2**10, 2**22, -1), fitlabel='\\ell', label='ST', skip_until=2**14)
	add_data_point_to_plot('zip-zip-trie', DataType.SEARCH, (2**10, 2**22, -1), fitlabel='\\ell', label='ZZT', skip_until=2**14)
	add_data_point_to_plot('memory-efficient-zip-zip-trie', DataType.SEARCH, (2**10, 2**22, -1), fitlabel='\\ell', label='ME-ZZT', skip_until=2**14)

	plot(figure_name='search-time-lcp-comparison', save=savefig, ylabel='Time (ns)', xlabel='Mean LCP Length ($\\ell$)')

if __name__ == '__main__':
	savefig = False
	# savefig = True

	# compare_construction_data(savefig)
	compare_search_data(savefig)
