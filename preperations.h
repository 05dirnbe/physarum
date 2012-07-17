#ifndef PREPERATIONS_H_INCLUDED
#define PREPERATIONS_H_INCLUDED

/**
* @file    preperations.h
* @ingroup preprocessing
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @brief A collection of code used to prepare for the execution of the main algorithm.
**/
/*
* Author: Michael Dirnberger <mtd@mpi-inf.mpg.de>
*
* Copyright (c)  Max-Planck-Institute Saarbruecken (Germany)
*
* This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
* WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*
*/
/**
*   @defgroup preprocessing Preprocessing And Preperations
*   @brief Code used fpr preproceessing and preparing stuff.
*   @{
**/

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/mpl/if.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/random.hpp>
#include <boost/concept/assert.hpp>
#include <boost/graph/iteration_macros.hpp>

/**
* @brief Generates a strongly connected, directed graph.
*
* Adds @p V vertices and @p E edges, to @p graph. Source and target vertices of each edge are randomly choosen.
* If @p self_edges is false, then no edge will have the same source and targets.
* When an edge is added the reverse edge is added as well, effectively doubling the
* number of edges in the graph. The reversed/transposed edges are marked.
* I.e. they become part of the transposed edge set \f$ g_t \f$. The non-transposed edges
* are marked as well and become part of the edge set \f$ g \f$. The edge set of the
* whole graph is thus \f$ G_e = g \cup g_t \f$. Additionally the reversed edge set is recorded
* in the respective property map for subsequent use in algorithms.
*
* @tparam MutableGraph          A type that models boost::MutableGraphConcept
* @tparam RandNumGen            A type that models boost::PseudoRandomNumberGenerator
* @tparam TransposedEdgeMap     A type that models a Read/Write property map
* @tparam ReverseEdgeMap        A type that models a Read/Write property map
*
* @param[out]   graph               An empty mutable graph
* @param[in]    V                   Number of random vertices to be added
* @param[in]    E                   Half the number of random edges to be added
* @param[out]   edge_transposed     Describes the transposed edge property
* @param[out]   edge_reversed       Describes the reversed edge property
* @param[in]    gen                 A boost random number functor
* @param[in]    self_edges           Specifies if self edges are allowed
*
* @todo Make sure that it is clear why the reversed and the transposed are both in use instead of one.
*/
template < typename MutableGraph, typename RandNumGen, typename TransposedEdgeMap, typename ReverseEdgeMap  >
void generate_rnd_graph(       MutableGraph& graph,
                                typename boost::graph_traits< MutableGraph >::vertices_size_type V,
                                typename boost::graph_traits< MutableGraph >::vertices_size_type E,
                                TransposedEdgeMap edge_transposed,
                                ReverseEdgeMap edge_reverse,
                                RandNumGen& gen,
                                bool self_edges = false )   {

    typedef boost::graph_traits< MutableGraph >             Traits;
    typedef typename Traits::vertices_size_type             v_size_t;
    typedef typename Traits::edges_size_type                e_size_t;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_descriptor                edge_descriptor;


    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::MutableGraphConcept< MutableGraph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< ReverseEdgeMap, edge_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< TransposedEdgeMap, edge_descriptor > ));
    //Am I missing some other concept?

    //preconditions
    BOOST_ASSERT( num_vertices( graph ) == 0 );
    BOOST_ASSERT( num_edges( graph ) == 0 );

      for ( v_size_t i = 0; i < V; ++i )
        add_vertex( graph );

      for ( e_size_t j = 0; j < E; ++j ) {

        vertex_descriptor a = boost::random_vertex(graph, gen), b;
        std::pair< edge_descriptor, bool > edge, t_edge;

        do {
          b = boost::random_vertex(graph, gen);
        } while ( self_edges == false && a == b );

        edge = add_edge(a, b, graph);                    // add an edge A
        t_edge = add_edge(b, a, graph);         // add the transposed/reverse edge B


        edge_transposed[edge.first] = false;      //register A as the reverse edge of B
        edge_transposed[t_edge.first] = true;     //register B as the reverse edge of A

        edge_reverse[edge.first] = t_edge.first; //register A as the reverse edge of B
        edge_reverse[t_edge.first] = edge.first; //register B as the reverse edge of A

      }

    //postcondition
    BOOST_ASSERT( num_vertices( graph ) == V );
    BOOST_ASSERT( num_edges( graph ) == 2*E );

  }


/**
* @brief Generates a basic Wheatstone graph.
*
* When an edge is added the reverse edge is added as well, effectively doubling the
* number of edges in the graph. The reversed/transposed edges are marked.
* I.e. they become part of the transposed edge set \f$ g_t \f$. The non-transposed edges
* are marked as well and become part of the edge set \f$ g \f$. The edge set of the
* whole graph is thus \f$ G_e = g \cup g_t \f$. Additionally the reversed edge set is recorded
* in the respective property map for subsequent use in algorithms.
*
* @tparam MutableGraph          A type that models boost::MutableGraphConcept
* @tparam TransposedEdgeMap     A type that models a Read/Write property map
* @tparam ReverseEdgeMap        A type that models a Read/Write property map
*
* @param[out]   graph               An empty mutable graph
* @param[out]   edge_transposed     Describes the transposed edge property
* @param[out]   edge_reversed       Describes the reversed edge property
*
* @todo Make sure that it is clear why the reversed and the transposed are both in use instead of one.
*/
template < typename MutableGraph, typename TransposedEdgeMap, typename ReverseEdgeMap  >
void generate_wheatstone_graph( MutableGraph& graph,
                                TransposedEdgeMap edge_transposed,
                                ReverseEdgeMap edge_reverse )   {

    typedef boost::graph_traits< MutableGraph >             Traits;
    typedef typename Traits::vertices_size_type             v_size_t;
    typedef typename Traits::edges_size_type                e_size_t;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_descriptor                edge_descriptor;


    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::MutableGraphConcept< MutableGraph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< ReverseEdgeMap, edge_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< TransposedEdgeMap, edge_descriptor > ));
    //Am I missing some other concept?

    //preconditions
    BOOST_ASSERT( num_vertices( graph ) == 0 );
    BOOST_ASSERT( num_edges( graph ) == 0 );


    const v_size_t V = 4;
    const e_size_t E = V + 1;

      for ( v_size_t i = 0; i < V ; ++i )
        add_vertex( graph );

        std::pair< edge_descriptor, bool > edge, t_edge;

        edge = add_edge(0, 1, graph);
        t_edge = add_edge(1, 0, graph);

        edge_transposed[edge.first] = false;
        edge_transposed[t_edge.first] = true;
        edge_reverse[edge.first] = t_edge.first;
        edge_reverse[t_edge.first] = edge.first;


        edge = add_edge(1, 2, graph);
        t_edge = add_edge(2, 1, graph);

        edge_transposed[edge.first] = false;
        edge_transposed[t_edge.first] = true;
        edge_reverse[edge.first] = t_edge.first;
        edge_reverse[t_edge.first] = edge.first;

        edge = add_edge(2, 3, graph);
        t_edge = add_edge(3, 2, graph);

        edge_transposed[edge.first] = false;
        edge_transposed[t_edge.first] = true;
        edge_reverse[edge.first] = t_edge.first;
        edge_reverse[t_edge.first] = edge.first;

        edge = add_edge(3, 0, graph);
        t_edge = add_edge(0, 3, graph);

        edge_transposed[edge.first] = false;
        edge_transposed[t_edge.first] = true;
        edge_reverse[edge.first] = t_edge.first;
        edge_reverse[t_edge.first] = edge.first;

        edge = add_edge(1, 3, graph);
        t_edge = add_edge(3, 1, graph);

        edge_transposed[edge.first] = false;
        edge_transposed[t_edge.first] = true;
        edge_reverse[edge.first] = t_edge.first;
        edge_reverse[t_edge.first] = edge.first;

    //postcondition
    BOOST_ASSERT( num_vertices( graph ) == V );
    BOOST_ASSERT( num_edges( graph ) == 2*E );

}


  /** @} */

#endif // PREPERATIONS_H_INCLUDED
