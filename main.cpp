
/**
* @mainpage Generic Physarum Solver
* @author Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @date 2012
*
* @section section_toc Table of Contents
*
* <ul>
*   <li> @ref intro
*   <li> @ref physarum
*   <li> @ref solver
*   <li> @ref experimental
*
* @subsection intro Introduction
* Here goes some introdoctory stuff about physarum. The minimum that is necessary to motivate the algorithm
*
* @subsection physarum A Model of Physarum Polycephalum
* Here goes the different model approaches as described in the papers.
* List strong points and weaknesses.
*
* @subsection solver Generic Physarum Solver
* Describe the idea of how to implement the algorithm
* Design choices and generic approaches. What are the algorithmic building blocks of the program. What are the choices for different core
* routines?
*
* @subsection experimental Testing and Experimental Results
* <ul>
* <li> Experimental setup (what rng + seed value, what high precision types, what types of input graphs, what random properties are in use, what hypotheses to test?)
* <li> Experimental execution (parameters, input sizes, what compiler, flags and libraries, special settings?)
* <li> Experimental analysis (error analysis, effects of quality of random numbers?, is the distribution of observables normal?, auto-correlation effects?
*                         some other artifacts?, systematic errors?, input bias )
* <li> Experimental results (dependence of the result on different instance size or input type, interpretation of errors, what do the observables imply? )
* </ul>
*
* General problems:
* <ul>
* <li> Whats the proper "scientific method" for CS experiments? How to bring together theory, observation (measurements), hypothesis testing and reproducibility???
* <li> How to I obtain meaninful and significant results concerning how my algorithm performs in real life? How do these results relate back to pure theory (i.e. the model)?
* <li> How do I distinguish between the stuff my algorithm does on paper and the (potentially) different thing that happens in the machine (e.g. systemeatic errors)?
* <li> How do I even guarantee that two runs are subject to the same systematic errors?
* <li> How to calibrate my experiment as to account for systematic influences?
* <li> How does an experiment capable of falsifying a hypothesis looks like in CS?
* </ul>
*
* I have no idea what it means to perform a proper CS experiment as of yet.
* If any of my experimental results shall ever find their way into some sort of publication, there MUST be a section on "Experimental methodology" including
* relevant points as mentiont above. In addition, programs and data have to be published alongside for reference and reproducibility.
*
* @todo What is the correct methodology for performing a reliable computer science experiment? How to relate experimental results to theory?
* @todo Make sure the main page is well rounded and ask Kurt how much precision on the error analysis is needed.
**/

#include "utility.h"
#include "develop.h"

#include <iostream>
#include <utility>

#include <boost/progress.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

#include "custom_graph.h"
#include "custom_randomize.h"
#include "config.h"

#include "preperations.h"


using namespace std;
using namespace boost;

int main() {

    GraphType graph;

    EdgeLengthMap length = get( length_t(), graph );
    EdgeDiameterMap diameter = get( diameter_t(), graph );
    EdgeCapacityMap capacity = get( edge_capacity, graph );                             //used for min-cut only
    EdgeReverseMap reverse_edge = get( edge_reverse, graph );
    EdgeTransposedMap transposed = get( transposed_t(), graph);
    EdgeNameMap name = get(edge_name, graph);
    EdgeFluxMap flux = get( flux_t(), graph);
    VertexIndexMap index = get( vertex_index, graph);
    VertexPotentialMap potential = get( potential_t(), graph);
    VertexPreviousPotentialMap previous_potential = get( previous_potential_t(), graph);

//    generate_rnd_graph( graph, n_vertices, n_edges, transposed, reverse_edge, rng);
    generate_wheatstone_graph( graph, transposed, reverse_edge);

    std::pair< vertex_descriptor, vertex_descriptor > source_sink = make_source_sink_pair( 0, 2 );
//    std::pair< vertex_descriptor, vertex_descriptor > source_sink = make_random_source_sink_pair( graph, rng );

    std::cout << "source_sink = ( " << source_sink.first << " , " << source_sink.second << " )\n";

    print_graph(graph, index );

    is_transposed_edge<EdgeTransposedMap> transposed_edge_filter( transposed );
    TransposedGraphType tgraph(graph, transposed_edge_filter);

    randomize_property< boost::property <       diameter_t,     BigFloat>::tag_type >( tgraph, random_edge_diameter );
    randomize_property< length_t >( tgraph, random_edge_length);
    randomize_property< edge_capacity_t >( tgraph, random_edge_length );

    BGL_FORALL_EDGES(v, tgraph, TransposedGraphType) {

        capacity[reverse_edge[v]] = capacity[v];
        length[reverse_edge[v]] = length[v];
        diameter[reverse_edge[v]] = diameter[v];
    }

    //making sure the edges are set up to mirror each other
    BOOST_ASSERT( mirror_edge_property_violated(graph, reverse_edge, capacity) == false );
    BOOST_ASSERT( mirror_edge_property_violated(graph, reverse_edge, length) == false );
    BOOST_ASSERT( mirror_edge_property_violated(graph, reverse_edge, diameter) == false );



    //This should be the algorithm

    unsigned int it_count = 0;

    compute_vertex_potentials(graph, source_sink.first, source_sink.second, potential, length, diameter, uncertainty_in_potentials);

    while ( calculate_maximal_potential_difference( graph, potential, previous_potential) > uncertainty_in_potentials ){

        compute_edge_flux( graph, potential, length, diameter, flux);
        store_current_vertex_potentials( graph, potential, previous_potential);
        compute_vertex_potentials( graph, source_sink.first, source_sink.second, potential, length, diameter, uncertainty_in_potentials);
        ++it_count;
    }

    //checking whether the edges still mirror each other
    BOOST_ASSERT( mirror_edge_property_violated(graph, reverse_edge, capacity) == false );
    BOOST_ASSERT( mirror_edge_property_violated(graph, reverse_edge, length) == false );
    BOOST_ASSERT( mirror_edge_property_violated(graph, reverse_edge, diameter) == false );

    std::cout << "diameter:\n";
    print_edge_property( tgraph, diameter, index);
    std::cout << "length:\n";
    print_edge_property( tgraph, length, index);
    std::cout << "diameter/length:\n";
    print_edge_property_quotient(tgraph, diameter, length, index);
    print_vertex_property( graph, potential, index);
    std::cout << "iteration count: " << it_count << std::endl;
}


