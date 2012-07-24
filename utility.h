#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED

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
* @file    utility.h
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @ingroup utility
*
* Groups the little helpers.
*
**/

#include <iostream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/concept/assert.hpp>
#include <boost/graph/iteration_macros.hpp>

/**
*   @defgroup utility Helper Functions and Utility
*   @brief Helper functions and utility
*   @{
**/


/**
* @brief Prints the quotient of two edge roperties for every edge
*
* Prints the value of edge property A divided by edge property B for each edge.
* In addition, the connected vertices are listed by name.
*
* @tparam Graph                             A type that models boost::GraphConcept
* @tparam AnyEdgePropertyMapA               A type that models a Readable property map
* @tparam AnyEdgePropertyMapB               A type that models a Readable property map
* @tparam AnyVertexNamePropertyMap          A type that models a Readable property map
*
* @param[in]    graph                   A graph
* @param[in]    edge_property_A         Describes an edge property A
* @param[in]    edge_property_B         Describes an edge property B
* @param[in]    edge_name_property      Describes a vertex name property
*/
template < typename Graph, typename AnyEdgePropertyMapA, typename AnyEdgePropertyMapB, typename AnyVertexNamePropertyMap >
inline void print_edge_property_quotient(   const Graph& graph,
                                            const AnyEdgePropertyMapA edge_property_A,
                                            const AnyEdgePropertyMapB edge_property_B,
                                            const AnyVertexNamePropertyMap vertex_name_property )   {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::edge_descriptor                edge_descriptor;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_iterator                  edge_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyEdgePropertyMapA, edge_descriptor> ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyEdgePropertyMapB, edge_descriptor> ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyVertexNamePropertyMap, vertex_descriptor> ));
    //Am I missing some other concept?

    edge_iterator e_it, e_it_end;

    for (boost::tie(e_it, e_it_end) = edges(graph); e_it != e_it_end; ++e_it)
        std::cout   << "edge[ (" << vertex_name_property[source(*e_it, graph)] << "," << vertex_name_property[target(*e_it, graph)] << ") ] = "
                    << edge_property_A[*e_it] / edge_property_B[*e_it] << "\n";

}



/**
* @brief Print a given edge property for all edges in the graph
*
* In addition, the connected vertices are listed by name.
*
* @tparam Graph                             A type that models boost::GraphConcept
* @tparam AnyEdgePropertyMap                A type that models a Readable property map
* @tparam AnyVertexNamePropertyMap          A type that models a Readable property map
*
* @param[in]    graph                   A graph
* @param[in]    edge_property           Describes an edge property
* @param[in]    edge_name_property      Describes a vertex name property
*/
template < typename Graph, typename AnyEdgePropertyMap, typename AnyVertexNamePropertyMap >
inline void print_edge_property(  const Graph& graph, const AnyEdgePropertyMap edge_property, const AnyVertexNamePropertyMap vertex_name_property )   {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::edge_descriptor                edge_descriptor;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_iterator                  edge_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyEdgePropertyMap, edge_descriptor> ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyVertexNamePropertyMap, vertex_descriptor> ));
    //Am I missing some other concept?

    edge_iterator e_it, e_it_end;

    for (boost::tie(e_it, e_it_end) = edges(graph); e_it != e_it_end; ++e_it)
        std::cout   << "edge[ (" << vertex_name_property[source(*e_it, graph)] << "," << vertex_name_property[target(*e_it, graph)] << ") ] = "
                    << edge_property[*e_it] << "\n";

}





/**
* @brief Print a given vertex property for all vertices in the graph
*
* @tparam Graph                             A type that models boost::GraphConcept
* @tparam AnyVertexPropertyMap              A type that models a Readable property map
* @tparam AnyVertexNamePropertyMap          A type that models a Readable property map
*
* @param[in]    graph                   A graph
* @param[in]    vertex_property         Describes an edge property
* @param[in]    vertex_name_property      Describes a vertex name property
*/
template < typename Graph, typename AnyVertexPropertyMap, typename AnyVertexNamePropertyMap >
inline void print_vertex_property(  const Graph& graph, const AnyVertexPropertyMap vertex_property, const AnyVertexNamePropertyMap vertex_name_property )   {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::vertex_iterator                vertex_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyVertexPropertyMap, vertex_descriptor> ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyVertexNamePropertyMap, vertex_descriptor> ));
    //Am I missing some other concept?

    vertex_iterator v_it, v_it_end;

    for (boost::tie(v_it, v_it_end) = vertices(graph); v_it != v_it_end; ++v_it)
        std::cout << "vertex[ " << vertex_name_property[*v_it ] << " ] = "<< vertex_property[*v_it] << "\n";

}



/**
* @brief Constructs and returns a source-sink std::pair
*
* Note that the function cannot check wether the input vertices belong to the same graph.
* Im not sure whether type deduction can pick these things up or not.
* Use vertices and descriptors of the same graph to be on the safe side.
* Switch on optimizations (return value optimization) and have the compiler eliminate unecessary copies.
*
* @tparam VertexDesciptor                   A type that models boost::GraphConcept< Graph >::vertex_descriptor
*
* @param[in]    s                   The source vertex
* @param[in]    t                   The sink
*
* @return       tmp                 A s-t pair pair<vertex_descriptor,vertex_descriptor>
*
* @todo Check wether the compile can detect if s and t stem from different graphs. If the error is not detected, make up some Concept for it.
*/
template < typename VertexDesciptor >
inline std::pair < VertexDesciptor, VertexDesciptor> make_source_sink_pair( VertexDesciptor s, VertexDesciptor t ) {

    std::pair < VertexDesciptor, VertexDesciptor> tmp(s,t);
 return tmp;
}



/**
* @brief Constructs and returns a random source-sink std::pair
*
* The returned pair always satisfies s != 0.
* Note that the function cannot check whether the input vertices belong to the same graph.
* Im not sure whether type deduction can pick these things up or not.
* Use vertices and descriptors of the same graph to be on the safe side.
* Switch on optimizations (return value optimization) and have the compiler eliminate unecessary copies.
*
* @tparam Graph                     A type that models boost::GraphConcept
* @tparam RandNumGen                A type that models boost::PseudoRandomNumberGenerator
*
* @param[in]    graph               A graph
* @param[in]    gen              A boost random number functor
*
* @return       tmp                 A s-t pair<vertex_descriptor,vertex_descriptor>
*/
template < typename Graph, typename RandNumGen >
inline std::pair < typename boost::graph_traits< Graph >::vertex_descriptor, typename boost::graph_traits< Graph >::vertex_descriptor>
make_random_source_sink_pair( const Graph& graph, RandNumGen gen ) {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));

    vertex_descriptor s = random_vertex(graph, gen);
    vertex_descriptor t;

    do {

     t = random_vertex(graph, gen);

    } while (s == t );

 return make_source_sink_pair(s,t);
}


/**
* @brief Checks for each edge if the value of a given edge property is equal to the value of the reversed edge.
*
* @tparam Graph                             A type that models boost::GraphConcept
* @tparam ReverseEdgePropertyMap            A type that models a Readable property map
* @tparam AnyEdgePropertyMap                A type that models a Readable property map
*
* @param[in]    graph                   A graph
* @param[in]    edge_reverse            Describes the graphs reversed_edge property map
* @param[in]    edge_property           Describes an any edge property
*
* @return       violation_detected
*/
template < typename Graph, typename ReverseEdgeMap, typename AnyEdgePropertyMap >
inline bool mirror_edge_property_violated(     const Graph& graph,
                                            const ReverseEdgeMap edge_reverse,
                                            const AnyEdgePropertyMap edge_property )   {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::edge_descriptor                edge_descriptor;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_iterator                  edge_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<AnyEdgePropertyMap, edge_descriptor> ));
    //Am I missing some other concept?

    bool violation_detected = false;
    edge_iterator e_it, e_it_end;

    for (boost::tie(e_it, e_it_end) = edges(graph); e_it != e_it_end; ++e_it)
        if (edge_property[*e_it] != edge_property[ edge_reverse[*e_it] ])
            violation_detected = true;

    return violation_detected;
}




/** @} */
#endif // UTILITY_H_INCLUDED
