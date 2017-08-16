#include <iostream>
#include <fstream>
#include <random>
#include <string>

#include "LayeredGraph.tpp"
#include "create_indiv_ustate_graph.h"
#include "types.h"

#include "debugtools.h"

unsigned int ntimesteps = 4000;
double dt = 0.001;

int main()
{
	std::mt19937 rng (time(0));
	std::uniform_real_distribution<double> unif(0.0, 1.0);

	unsigned int L = 4;

	LayeredGraph<ustate, activeEdge> cspace_raw;
	create_indiv_ustate_graph(cspace_raw, L, false);

	//cspace_raw.m_conn_manip();


	std::ofstream graph_structure; graph_structure.open("dat/data.json");
	cspace_raw.m_print_graph_json(graph_structure);
	graph_structure.close();
	std::ofstream vertex_legend; vertex_legend.open("dat/legend.dat");
	cspace_raw.m_print_vertex_legend(vertex_legend);
	vertex_legend.close();

	for (unsigned int ct_i_vertex = 0; ct_i_vertex < cspace_raw.g_n_vertices(); ++ct_i_vertex)
	{
		auto i_vertex_ptr = cspace_raw.g_vertex_ptr(ct_i_vertex);
		i_vertex_ptr->s_population(1.00/cspace_raw.g_n_vertices());//unif(rng)
	}

	unsigned int timestep = 0;
	std::cout << dt*static_cast<double>(timestep) << " ";
	for (auto& elem: cspace_raw.g_vertices()) { std::cout << elem->g_population() << " "; }	std::cout << "\n";

	for (timestep = 1; timestep < ntimesteps; ++timestep)
	{
		std::vector<long double> dpdt = cspace_raw.m_calculate_rates();
		for (auto& elem: dpdt) { elem *= dt; }
		cspace_raw.m_update_populations(dpdt);

		std::cout << dt*static_cast<double>(timestep) << " ";
		for (auto& elem: cspace_raw.g_vertices()) { std::cout << elem->g_population() << " "; }	std::cout << "\n";
	}

	return 0;
}
