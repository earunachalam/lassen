#ifndef CREATE_MICROSTATE_GRAPH_H
#define CREATE_MICROSTATE_GRAPH_H

#include <algorithm>
#include "LayeredGraph.tpp"
#include "types.h"

double k_fwd = 1.0;
double k_rev = 1.0;
double alpha = 1.0;

void create_indiv_ustate_graph(LayeredGraph<ustate, activeEdge>& cspace_raw, unsigned int L, bool ignore_multiplicity)
{
	// a base microstate, with all monomers unbound
	std::vector<unsigned int> seg_lengths(L);
	std::fill(seg_lengths.begin(), seg_lengths.end(), 1);
	auto new_vertex_ptr = cspace_raw.m_add_vertex(seg_lengths);
	cspace_raw.m_expand_layer(new_vertex_ptr, 0);

	unsigned int ct_loop = 0;
	bool not_terminated = true;
	while (not_terminated && /*just as a safeguard; not critical to function*/ (++ct_loop)<(L+5))
	{
		cspace_raw.m_insert_new_layer();
		auto next_layer_ptr = cspace_raw.g_layer_ptr(-1);
		auto current_layer_ptr = cspace_raw.g_layer_ptr(-2);
		for (unsigned int ct_i_vtx = 0; ct_i_vtx < current_layer_ptr->size(); ++ct_i_vtx)
		{
			auto i_vtx_ptr = current_layer_ptr->at(ct_i_vtx);
			unsigned int prev_j_el = -1;

			auto current_ustate = i_vtx_ptr->g_value();
			for (unsigned int ct_j_el = 0; ct_j_el < current_ustate.size(); ++ct_j_el)
			{
				unsigned int j_el = current_ustate.at(ct_j_el);
					
				if (ignore_multiplicity)
				{
					if (j_el == prev_j_el)	continue;
					else					prev_j_el = j_el;
				}

				unsigned int prev_k_el = -1;

				for (unsigned int ct_k_el = 1 + ct_j_el; ct_k_el < current_ustate.size(); ++ct_k_el)
				{
					unsigned int k_el = current_ustate.at(ct_k_el);
					if (ignore_multiplicity)
					{
						if (k_el == prev_k_el)	continue;
						else 					prev_k_el = k_el;
					}

					auto current_ustate_chg = current_ustate;
					unsigned int parent_idx_max = std::max(ct_j_el, ct_k_el);
					unsigned int parent_idx_min = std::min(ct_j_el, ct_k_el);
					auto b = current_ustate_chg.begin();
					current_ustate_chg.erase(b + parent_idx_max);
					current_ustate_chg.erase(b + parent_idx_min);
					current_ustate_chg.push_back(j_el + k_el);
					std::sort(current_ustate_chg.begin(), current_ustate_chg.end());

					double rate_fwd;
					if 		((j_el == 1) && (k_el == 1)) { rate_fwd = k_fwd + alpha; }
					else if	((j_el == 1) && (k_el != 1)) { rate_fwd = k_fwd + alpha; }
					else if	((j_el != 1) && (k_el == 1)) { rate_fwd = k_fwd + alpha; }
					else if	((j_el != 1) && (k_el != 1)) { rate_fwd = k_fwd; }

					std::string color_fwd;
                    if ((j_el == 1) || (k_el == 1)) { color_fwd = "#f00"; }
                    else { color_fwd = "#00f"; }

					std::string color_rev = "#00f";

					bool new_ustate_generate_fwdd = true;
					for (auto& existing_ustate_ptr : *next_layer_ptr)
					{
						if (existing_ustate_ptr->g_value() == current_ustate_chg)
						{
							new_ustate_generate_fwdd = false;
							auto ae_fwd = std::make_shared<activeEdge>(rate_fwd, color_fwd);
							cspace_raw.m_add_edge(i_vtx_ptr, existing_ustate_ptr, ae_fwd);
							auto ae_rev = std::make_shared<activeEdge>(k_rev, color_rev);
							cspace_raw.m_add_edge(existing_ustate_ptr, i_vtx_ptr, ae_rev);

							break;
						}
					}

					if (new_ustate_generate_fwdd)
					{
						unsigned int layer_idx = L - current_ustate_chg.size();
						new_vertex_ptr = cspace_raw.m_add_vertex(current_ustate_chg);
						cspace_raw.m_expand_layer(new_vertex_ptr, layer_idx);
						auto ae_fwd = std::make_shared<activeEdge>(rate_fwd, color_fwd);
						cspace_raw.m_add_edge(i_vtx_ptr, new_vertex_ptr, ae_fwd);
						auto ae_rev = std::make_shared<activeEdge>(k_rev, color_rev);
						cspace_raw.m_add_edge(new_vertex_ptr, i_vtx_ptr, ae_rev);
					}

					next_layer_ptr = cspace_raw.g_layer_ptr(-1);
					current_layer_ptr = cspace_raw.g_layer_ptr(-2);
				}
			}
		}

		auto last_layer_ptr = cspace_raw.g_layer_ptr(-1);
		if (last_layer_ptr->at(0)->g_value().size() == 1) not_terminated = false;
	}
}

#endif
