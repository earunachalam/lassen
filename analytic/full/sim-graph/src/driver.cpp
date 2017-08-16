#include <array>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <random>

#include "Graph.tpp"
#include "types.h"


typedef std::array<unsigned int, 2> loc2d; // location of a point in the x,y plane

double dt = 0.001;
double alpha = 2.00;

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		std::cerr << "Error: Insufficient arguments provided." << std::endl;
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}
	
	char* fname = argv[1];
	FILE* fp = fopen(fname, "r");
	if (fp == nullptr)
	{
		std::cerr << "Error: File could not be opened." << std::endl;
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}

	char *endptr;
	double t = strtod(argv[2], &endptr);
	if (!*argv[2] || *endptr)
	{
		std::cerr << "Error: Invalid time provided: " << argv[2] << std::endl;
		std::cerr << __FILE__ << "\tL" << __LINE__ << std::endl;
		abort();
	}
	unsigned int ntimesteps = static_cast<unsigned int>(t/dt);

	char ignore[1000];

	// edge parameters
	unsigned int from_nseg, // 'origin' vertex for given edge: number of segments in configuration 
				 from_nmon,	// number of monomers in configuration
				 to_nseg,	// 'destination' vertex for given edge: number of segments in configuration
				 to_nmon,	// number of monomers in configuration
				 mult;		// multiplicity of given edge bundle (number of equivalent edges, or edges connecting degenerate origin and destination vertices
	
	fgets(ignore, 1000, fp);
	
	Graph<loc2d, summableEdge> g;

	while (!feof(fp))
	{
		unsigned int nread = fscanf(fp, "%u %u %u %u %u\n", &from_nseg, &from_nmon, &to_nseg, &to_nmon, &mult);
		if (nread != 5)
		{
			std::cerr << "Error: invalid row in datafile.\n" << std::endl;
			abort();
		}

		loc2d from_coords	 = {{from_nmon, from_nseg}};
		loc2d to_coords		 = {{to_nmon, to_nseg}};

		double rateconst_fwd = static_cast<double>(mult),
			   rateconst_rev = static_cast<double>(mult);

		if (from_nmon != to_nmon) rateconst_fwd *= alpha;
		summableEdge new_edge_fwd(rateconst_fwd);
		summableEdge new_edge_rev(rateconst_rev);
		g.add_edge_nodup_vertex(from_coords, to_coords, new_edge_fwd);
		g.add_edge_nodup_vertex(to_coords, from_coords, new_edge_rev);
	}

	std::mt19937 rng (time(0));
	std::uniform_real_distribution<double> unif(0.0, 1.0);

	for (unsigned int ct_i_vertex = 0; ct_i_vertex < g.g_n_vertices(); ++ct_i_vertex)
	{
		auto i_vertex_ptr = g.g_vertex_ptr(ct_i_vertex);
		i_vertex_ptr->s_population(1.00/g.g_n_vertices());//unif(rng)
	}

	unsigned int timestep = 0;

	std::cout << "NaN "; for (auto& elem: g.g_vertices()) { std::cout << elem->g_value()[0] << " "; }	std::cout << "\n";
	std::cout << "NaN "; for (auto& elem: g.g_vertices()) { std::cout << elem->g_value()[1] << " "; }	std::cout << "\n";
	
	std::cout << dt*static_cast<double>(timestep) << " ";
	for (auto& elem: g.g_vertices()) { std::cout << elem->g_population() << " "; }	std::cout << "\n";

	for (timestep = 1; timestep < ntimesteps; ++timestep)
	{
		std::vector<long double> dpdt = g.calculate_rates();
		for (auto& elem: dpdt) { elem *= dt; }
		g.update_populations(dpdt);

		std::cout << dt*static_cast<double>(timestep) << " ";
		for (auto& elem: g.g_vertices()) { std::cout << elem->g_population() << " "; }	std::cout << "\n";
	}
	return 0;
}
