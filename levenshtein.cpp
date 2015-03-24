#define _CRT_SECURE_NO_WARNINGS
/*
 * Levenshtein's Edit Distance
 */
#include "implementation.hpp"
#include <exception.hpp>
#include <stopwatch.hpp>
#include <interface.hpp>

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <cstdint>
#include <cstdio>


typedef std::uint32_t char_t;


void print_usage()
{
	std::cout << "Arguments: [ -debug ] <file1> <file2>" << std::endl;
	std::cout << "  -debug  - flag for debugging output" << std::endl;
	std::cout << "  <file>  - files with strings to be compared by edit distance" << std::endl;
}



/*
 * \bried Load an entire file into a vector of chars.
 */
void load_file(const std::string &fileName, std::vector<char_t> &res)
{
	// Open the file.
	std::FILE *fp = std::fopen(fileName.c_str(), "rb");
	if (fp == nullptr)
		throw (bpp::RuntimeError() << "File '" << fileName << "' cannot be opened for reading.");
	
	// Determine length of the file and 
	std::fseek(fp, 0, SEEK_END);
	std::size_t count = (std::size_t)(std::ftell(fp) / sizeof(char_t));
	std::fseek(fp, 0, SEEK_SET);
	res.resize(count);

	// Read the entire file.
	std::size_t offset = 0;
	while (offset < count) {
		std::size_t batch = std::min<std::size_t>(count - offset, 1024*1024);
		if (std::fread(&res[offset], sizeof(char_t), batch, fp) != batch)
			throw (bpp::RuntimeError() << "Error while reading from file '" << fileName << "'.");
		offset += batch;
	}
	
	std::fclose(fp);
}


// Main routine that performs the computation.
template<bool DEBUG>
std::size_t computeDistance(const std::vector<char_t> str1, const std::vector<char_t> &str2)
{
	// Initialize distance functor.
	EditDistance<char_t, std::size_t, DEBUG> distance;
	distance.init(str1.size(), str2.size());

	// Compute the distance.
	bpp::Stopwatch stopwatch(true);
	std::size_t res = distance.compute(str1, str2);
	stopwatch.stop();
	std::cout << "Time: " << stopwatch.getMiliseconds() << "\t";

	return res;
}


/*
 * Application Entry Point
 */
int main(int argc, char **argv)
{
	// Process arguments.
	--argc; ++argv;
	bool debug = false;
	if (argc == 3 && std::string(*argv) == std::string("-debug")) {
		--argc; ++argv;
		debug = true;
	}

	if (argc != 2) {
		print_usage();
		return 0;
	}


	// Load files.
	std::vector<char_t> str1, str2;
	try {
		load_file(argv[0], str1);
		load_file(argv[1], str2);
	}
	catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		print_usage();
		return 1;
	}

	std::size_t dist = (debug)
			? computeDistance<true>(str1, str2)
			: computeDistance<false>(str1, str2);
		std::cout << dist << std::endl;
	// Compute the distance.
	/*
	try {
		std::size_t dist = (debug)
			? computeDistance<true>(str1, str2)
			: computeDistance<false>(str1, str2);
		std::cout << dist << std::endl;
	}
	catch (std::exception &e) {
		std::cout << "FAILED" << std::endl;
		std::cerr << e.what() << std::endl;
		return 2;
	}
	*/
	return 0;
}
