// Stupid Fast Sudoku Solver.cpp
#include <iostream>
#include <array>
#include <algorithm>
#include <string_view>
#include "fmt/format.h"
#include "fmt/color.h"

using namespace std;
using namespace std::literals;

//https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
#include <limits>   

namespace Detail
{
	double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
	{
		return curr == prev
			? curr
			: sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
	}
}

/*
* Constexpr version of the square root
* Return value:
*   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
*   - Otherwise, returns NaN
*/
namespace patchwork {
	double constexpr sqrt(double x)
	{
		return x >= 0 && x < std::numeric_limits<double>::infinity()
			? Detail::sqrtNewtonRaphson(x, x, 0)
			: std::numeric_limits<double>::quiet_NaN();
	}
}

//https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
typedef struct RgbColor
{
	uint8_t r;
	uint8_t g;
	uint8_t b;
} RgbColor;

typedef struct HsvColor
{
	uint8_t h;
	uint8_t s;
	uint8_t v;
} HsvColor;

RgbColor HsvToRgb(HsvColor hsv)
{
	RgbColor rgb;
	unsigned char region, remainder, p, q, t;

	if (hsv.s == 0)
	{
		rgb.r = hsv.v;
		rgb.g = hsv.v;
		rgb.b = hsv.v;
		return rgb;
	}

	region = hsv.h / 43;
	remainder = (hsv.h - (region * 43)) * 6;

	p = (hsv.v * (255 - hsv.s)) >> 8;
	q = (hsv.v * (255 - ((hsv.s * remainder) >> 8))) >> 8;
	t = (hsv.v * (255 - ((hsv.s * (255 - remainder)) >> 8))) >> 8;

	switch (region)
	{
	case 0:
		rgb.r = hsv.v; rgb.g = t; rgb.b = p;
		break;
	case 1:
		rgb.r = q; rgb.g = hsv.v; rgb.b = p;
		break;
	case 2:
		rgb.r = p; rgb.g = hsv.v; rgb.b = t;
		break;
	case 3:
		rgb.r = p; rgb.g = q; rgb.b = hsv.v;
		break;
	case 4:
		rgb.r = t; rgb.g = p; rgb.b = hsv.v;
		break;
	default:
		rgb.r = hsv.v; rgb.g = p; rgb.b = q;
		break;
	}

	return rgb;
}

HsvColor RgbToHsv(RgbColor rgb)
{
	HsvColor hsv;
	unsigned char rgbMin, rgbMax;

	rgbMin = rgb.r < rgb.g ? (rgb.r < rgb.b ? rgb.r : rgb.b) : (rgb.g < rgb.b ? rgb.g : rgb.b);
	rgbMax = rgb.r > rgb.g ? (rgb.r > rgb.b ? rgb.r : rgb.b) : (rgb.g > rgb.b ? rgb.g : rgb.b);

	hsv.v = rgbMax;
	if (hsv.v == 0)
	{
		hsv.h = 0;
		hsv.s = 0;
		return hsv;
	}

	hsv.s = 255 * long(rgbMax - rgbMin) / hsv.v;
	if (hsv.s == 0)
	{
		hsv.h = 0;
		return hsv;
	}

	if (rgbMax == rgb.r)
		hsv.h = 0 + 43 * (rgb.g - rgb.b) / (rgbMax - rgbMin);
	else if (rgbMax == rgb.g)
		hsv.h = 85 + 43 * (rgb.b - rgb.r) / (rgbMax - rgbMin);
	else
		hsv.h = 171 + 43 * (rgb.r - rgb.g) / (rgbMax - rgbMin);

	return hsv;
}


template <size_t N> struct ansi_esc_sequence {
	uint8_t tail;
	uint8_t opts[N];
};

template <size_t N> struct fmt::formatter<ansi_esc_sequence<N>> : public fmt::formatter<string_view> {
	template <typename FormatContext>
	auto format(const ansi_esc_sequence<N>& p, FormatContext& ctx) noexcept {
		if constexpr (N > 0) {
			return fmt::format_to(ctx.out(), "\033[{}{}"sv, fmt::join(&p.opts[0], &p.opts[0] + N, ";"),
				(char)p.tail);
		}
		else {
			return fmt::format_to(ctx.out(), "\033[{}"sv, (char)p.tail);
		}
	}
};

using int_type = uint8_t;
using sudoku_mark = int_type;

struct sudoku_selection {
	int_type mark;
	int_type y;
	int_type x;
	//block coordinate is derived from x and y
	int_type b;
};

struct sudoku_selection_ordered {
	sudoku_selection selection;
	int_type idx;
};

/*
@param initial_grid one dimensional array of sudoku selections
@returns an array of selections which complete (or fail to complete if impossible) sudoku grid
*/
enum class overlap_error : uint8_t {
	no_overlap,
	row_overlap,
	col_overlap,
	block_overlap
};

void print_overlap_error(overlap_error err_type, int_type location, int_type mark) {
	std::array<std::string_view, 4> dim_str = { ""sv,"row"sv,"column"sv,"block"sv };
	fmt::print("error: unsolvable puzzle at least two cells in {} {} required to be {}\n"sv, dim_str[(uint8_t)err_type], location, mark);
}

template<size_t output_size, size_t input_size>
void print_sudoku(const std::array<sudoku_selection, input_size> initial_grid) {
	static_assert(output_size >= input_size, "output_size of sudoku solver must be >= than input_size");
	constexpr size_t N = patchwork::sqrt(output_size); //output_size should be perfect square

	static_assert(std::numeric_limits<int_type>::max() > (N + 2), "int_type needs to fit puzzle size N");
	constexpr int_type blk_size = patchwork::sqrt(N);

	static_assert(blk_size > 0, "the block size must be > 0");

	std::array<sudoku_selection_ordered, output_size> ret_puzzle = { 0 };

	int_type max_val = 0;
	for (size_t i = 0; i < input_size; i++) {
		const sudoku_selection& cell = initial_grid[i];
		size_t flat_index = cell.y * N + cell.x;
		//overwrite previous value
		sudoku_selection_ordered& ret_cell = ret_puzzle[flat_index];
		ret_cell.selection = cell;
		ret_cell.idx = i;
		max_val = std::max(max_val, cell.mark);
	}

	//force everything to be within a stack allocated buffer
	fmt::basic_memory_buffer<char, (output_size + (4 * N)) * 32> out;
	fmt::format_to(out, "{}"sv, max_val);
	const size_t width = out.size() + 2;
	out.clear();

	auto fmtit = fmt::format_to(out, "{:^{}}{:^{}}{:^{}}\n"sv, N, width, "x"sv, width, N, width);

	//left corner
	fmt::format_to(out, "{:^{}}"sv, " "sv, width);
	//x axis
	for (size_t x = 0; x < N; x++) {
		fmtit = fmt::format_to(fmtit, "{:^{}}"sv, x, width);
	}
	//horizontal seperators
	fmtit = fmt::format_to(fmtit, "\n{:^{}}{:_^{}}\n"sv, " "sv, width, ""sv, width * N);

	for (size_t y = 0; y < N; y++) {
		//y axis and seperator
		fmtit = fmt::format_to(fmtit, "{:^{}}|"sv, y, width);
		for (size_t x = 0; x < N; x++) {
			size_t flat_index = y * N + x;
			const sudoku_selection_ordered& cell = ret_puzzle[flat_index];
			/*
			fmt::format(fmt::emphasis::bold | fg(fmt::color::red),
                                      "The answer is {}", 42);
			*/
			//uint8_t{8} *uint8_t{cell.idx}
			if (cell.selection.mark) {
				//ansi_esc_sequence<5>color{ 'm', 38, 2, 220, 110 + ((std::numeric_limits<int_type>::max() - 110) / output_size) * cell.idx, (std::numeric_limits<int_type>::max() / output_size) * cell.idx };
				HsvColor _hsv{ ((std::numeric_limits<int_type>::max()) / N) * cell.selection.mark, 55 + ((std::numeric_limits<int_type>::max() - 55) / output_size) * cell.idx, 55 + ((std::numeric_limits<int_type>::max() - 55) / output_size) * cell.idx };
				RgbColor _rgb = HsvToRgb(_hsv);
				//ansi_esc_sequence<5>color{'m', 38, 2, 220, 110 + ((std::numeric_limits<int_type>::max() - 110) / N) * cell.selection.mark, ((std::numeric_limits<int_type>::max()) / output_size) * cell.idx };
				ansi_esc_sequence<5>color{ 'm', 38, 2, _rgb.r, _rgb.g, _rgb.b};

				ansi_esc_sequence<0>clear_ansi_line{ 'K' };
				fmtit = fmt::format_to(fmtit, "{}{}{:^{}}"sv, color, clear_ansi_line, cell.selection.mark, width);
			}
			else {
				//ansi_esc_sequence<5>color{ 'm', 38, 2, 220, 220 + (std::numeric_limits<int_type>::max() / output_size) * cell.idx, (std::numeric_limits<int_type>::max() / output_size) * cell.idx };
				ansi_esc_sequence<1>default_ansi_line{ 'm', 0 };
				ansi_esc_sequence<0>clear_ansi_line{ 'K' };
				fmtit = fmt::format_to(fmtit, "{}{}{:^{}}"sv, default_ansi_line, clear_ansi_line, "."sv, width);
			}
		}
		fmtit = fmt::format_to(fmtit, "{}\n"sv, ansi_esc_sequence<1>{'m', 0});
	}
	fmt::print("{}{}"sv, std::string_view{ out.data(), out.size() }, ansi_esc_sequence<1>{'m', 0});
}



//template<size_t N>
template<size_t output_size, size_t input_size>
std::array<sudoku_selection, output_size> solve_sudoku(const std::array<sudoku_selection, input_size> initial_grid) {
	static_assert(output_size >= input_size, "output_size of sudoku solver must be >= than input_size");

	constexpr size_t N = patchwork::sqrt(output_size); //output_size should be perfect square
	static_assert(std::numeric_limits<int_type>::max() > (N + 2), "int_type needs to fit puzzle size N and have extra room for two indexs");

	constexpr int_type blk_size = patchwork::sqrt(N);
	static_assert(blk_size > 0, "the block size must be > 0");

	//0 initialize
	std::array<std::array<sudoku_mark, N + 1>, N> row_mark = { 0 };
	std::array<std::array<sudoku_mark, N + 1>, N> col_mark = { 0 };
	std::array<std::array<sudoku_mark, N + 1>, N> blk_mark = { 0 };
	//per cell information
	std::array<std::array<int_type, N + 1>, N* N> cell_mark = { 0 };

	//ensure no duplicate indexes
	std::array<int_type, N* N> index_mark = { 0 };

	std::array<sudoku_selection, output_size> ret_puzzle = { 0 };

	const auto calc_blk = [](const int_type x, const int_type y, const int_type blk_width_and_height) -> int_type {
		int_type x_div = x / blk_width_and_height;
		int_type y_div = y / blk_width_and_height;
		int_type blk_ret = (y_div * blk_width_and_height + x_div);
		return blk_ret;
	};

	const auto calc_flat_index = [](const int_type x, const int_type y) -> int_type {
		return x * int_type{ N } + y;
	};
	/*
	const auto mark_cell = [&row_mark, &col_mark, &blk_mark](const sudoku_selection cell, int_type mark) -> void {
		sudoku_selection ret = cell;
		ret.mark = mark;

		row_mark[cell.y][mark] = mark;
		col_mark[cell.x][mark] = mark;
		blk_mark[cell.b][mark] = mark;
		return ret;
	};
	*/
	const auto unmark_cell = [&row_mark, &col_mark, &blk_mark](const sudoku_selection cell) -> sudoku_selection {
		sudoku_selection ret = cell;
		ret.mark = 0ULL;

		row_mark[cell.y][cell.mark] = 0ULL;
		col_mark[cell.x][cell.mark] = 0ULL;
		blk_mark[cell.b][cell.mark] = 0ULL;
		return ret;
	};

	const auto mark_cell_with = [&row_mark, &col_mark, &blk_mark](const sudoku_selection cell, int_type mark) -> sudoku_selection {
		sudoku_selection ret = cell;
		ret.mark = mark;
		//better access pattern
		row_mark[cell.y][cell.mark] = 0ULL;
		row_mark[cell.y][mark] = mark;

		col_mark[cell.x][cell.mark] = 0ULL;
		col_mark[cell.x][mark] = mark;

		blk_mark[cell.b][cell.mark] = 0ULL;
		blk_mark[cell.b][mark] = mark;
		return ret;
	};
	//initialize all positions in array
	{
		size_t idx = 0;
		for (size_t y = 0; y < N; y++) {
			for (size_t x = 0; x < N; x++) {
				sudoku_selection& cell = ret_puzzle[idx++];
				cell.x = x;
				cell.y = y;
				cell.b = calc_blk(cell.x, cell.y, blk_size);
			}
		}
	}

	size_t set_numbers = 0;
	//mark rows, columns and blocks as having used certain numbers from input
	for (size_t i = 0; i < input_size; i++) {
		const sudoku_selection& input_cell = initial_grid[i];
		int_type x_ = input_cell.x;
		int_type y_ = input_cell.y;
		int_type b_ = input_cell.b;
		int_type mark_ = input_cell.mark;

		size_t flat_index = x_ * N + y_;

		if (flat_index >= output_size) [[unlikely]] {
			fmt::print("error: badly formatted input, cell[{}] has an index {} which is >= {}\n"sv, i, flat_index, N);
			return ret_puzzle;
		}
		else if (index_mark[flat_index]) [[unlikely]] {
			fmt::print("error: badly formatted input, cell[{}] is a duplicate cell x:{} y:{}\n", i, x_, y_);
			return ret_puzzle;
		}
			//reorganize into grid in order
			//0 indexed for fma array access
		sudoku_selection& cell = ret_puzzle[i];
		cell.mark = mark_;

		//shift by one so we can test by checking > 0
		index_mark[flat_index] = flat_index + 1;

		if (mark_ > N) [[unlikely]] {
			fmt::print("error: badly formatted input, cell[{}] has a choice value ({}) > than the width of the puzzle ({})", i, mark_, N);
			return ret_puzzle;
		}
			if (y_ && row_mark[y_][mark_]) [[unlikely]] {
				fmt::print("error: unsolvable puzzle at least two cells in row {} are {}\n"sv, y_, mark_);
				return ret_puzzle;
			}
			else if (x_ && col_mark[x_][mark_]) [[unlikely]] {
				fmt::print("error: unsolvable puzzle, at least two cells in column {} are {}\n"sv, x_, mark_);
				return ret_puzzle;
			}
			else if (b_ && blk_mark[b_][mark_]) [[unlikely]] {
				fmt::print("error: unsolvable puzzle at least two cells in block {} are {}\n"sv, b_, mark_);
				return ret_puzzle;
			}

				//may perform dead writes to 0 initialized memory...do we care?
		cell = mark_cell_with(cell, mark_);
		set_numbers += mark_ ? int_type{ 1ULL } : int_type{ 0ULL };
	}

	//split into two arrays, an array where we've already made selections we will absolutely not undo and ones we will work with as a stack
	auto unsolved_it = std::partition(ret_puzzle.begin(), ret_puzzle.end(), [](sudoku_selection& cell) {
		return cell.mark;
		});
	size_t unsolved_index = std::distance(ret_puzzle.begin(), unsolved_it);

	//sort by least # of choices
	std::sort(ret_puzzle.begin() + unsolved_index, ret_puzzle.end(), [&row_mark, &col_mark, &blk_mark](sudoku_selection& lhs, sudoku_selection& rhs) {
		size_t lhs_marks = 0;
		size_t rhs_marks = 0;
		//we can ignore 0 b/c it's the throw away index
		for (size_t i = 1; i < (N + 1); i++) {
			lhs_marks += (row_mark[lhs.y][i] || col_mark[lhs.x][i] || blk_mark[lhs.b][i]);
		}
		for (size_t i = 1; i < (N + 1); i++) {
			rhs_marks += (row_mark[rhs.y][i] || col_mark[rhs.x][i] || blk_mark[rhs.b][i]);
		}
		return rhs_marks < lhs_marks;
		});

	//solver algorithm
	for (size_t x = unsolved_index; x >= unsolved_index && x < output_size;) {
		sudoku_selection& cell = ret_puzzle[x];
		//if we backtrack we start w/ this value +1 (ergo starting mark moves to 1 and N means failure)
#if _DEBUG
		for (size_t i = cell.mark; ++i <= N;) {
			if (row_mark[cell.y][i] || col_mark[cell.x][i] || blk_mark[cell.b][i]) {
				fmt::print("cell[{},{}] != {}\n"sv, cell.x, cell.y, i);
				continue;
			}
		}
#endif

		for (size_t i = cell.mark; ++i <= N;) {
			if (row_mark[cell.y][i] || col_mark[cell.x][i] || blk_mark[cell.b][i]) {
				//skip things we cannot be b/c of other selections
				continue;
			}
			//mark the cell (clears any previous value)
			cell = mark_cell_with(cell, i);
#if _DEBUG
			fmt::print("cell[{},{}] = {}\n"sv, cell.x, cell.y, i);
			print_sudoku<output_size>(ret_puzzle);
#endif
			//move on to next thing
			goto next_stack;
		}
		//could not fill the cell w/ a value!

		size_t last_x = x;
		--x;
		if (x >= unsolved_index && x < output_size) {
			//mark cell as empty in the event we revisit in future
			cell = unmark_cell(cell);
			continue;
		}
		else {
			//we've walked the entire tree back to the current state (like above, except we may not erase previous markers)
			fmt::print("{}"sv, "error: unsolveable sudoku"sv);
			break;
		}
	next_stack:
		//resort by the least# of choices left to increase speed
		auto it = std::min_element(ret_puzzle.begin() + (x + 1), ret_puzzle.end(), [&row_mark, &col_mark, &blk_mark](sudoku_selection& lhs, sudoku_selection& rhs) {
			size_t lhs_marks = 0;
			size_t rhs_marks = 0;
			//we can ignore 0 b/c it's the throw away index
			for (size_t i = 1; i < (N + 1); i++) {
				lhs_marks += (row_mark[lhs.y][i] || col_mark[lhs.x][i] || blk_mark[lhs.b][i]);
			}
			for (size_t i = 1; i < (N + 1); i++) {
				rhs_marks += (row_mark[rhs.y][i] || col_mark[rhs.x][i] || blk_mark[rhs.b][i]);
			}
			return rhs_marks < lhs_marks;
		});
		std::swap(ret_puzzle[x + 1], *it);
		++x;
		//forward progress we trample over stack
		//ret_puzzle[x] = unmark_cell(ret_puzzle[x]);
		continue;
	}

	return ret_puzzle;
}

int main()
{
	std::array<sudoku_selection, 9ULL * 9ULL> initial_sudoku_grid = { 0 };
	//setup coordinates
	for (size_t y = 0; y < 9ULL; y++) {
		for (size_t x = 0; x < 9UL; x++) {
			sudoku_selection& cell = initial_sudoku_grid[y * 9ULL + x];
			cell.y = y;
			cell.x = x;
		}
	}

	for (size_t y = 0; y < 1; y++) { //fill the 0'th row w/ 1 2 ... etc to 9
		for (size_t x = 0; x < 9UL; x++) {
			sudoku_selection& cell = initial_sudoku_grid[y * 9ULL + x];
			cell.mark = (x + 1);
		}
	}

	print_sudoku<9ULL * 9ULL>(initial_sudoku_grid);
	std::array<sudoku_selection, 9ULL * 9ULL> solved = solve_sudoku<9ULL * 9ULL>(initial_sudoku_grid);
	print_sudoku<9ULL * 9ULL>(solved);

	//duplicate the [1][0] cell w/ 1
	initial_sudoku_grid[1 * 9 + 0].mark = 4;
	initial_sudoku_grid[1 * 9 + 2].mark = 6;

	print_sudoku<9ULL * 9ULL>(initial_sudoku_grid);
	std::array<sudoku_selection, 9ULL * 9ULL> solved2 = solve_sudoku<9ULL * 9ULL>(initial_sudoku_grid);
	print_sudoku<9ULL * 9ULL>(solved2);

	return 0;
}
