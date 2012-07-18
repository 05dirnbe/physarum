#ifndef DEVELOP_H_INCLUDED
#define DEVELOP_H_INCLUDED

/**
* @file    develop.h
* @ingroup develop
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @brief A collection of experimental code
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

#include <cmath>
#include <algorithm>
#include <iostream>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/mpl/if.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/random.hpp>
#include <boost/concept/assert.hpp>
#include <boost/graph/iteration_macros.hpp>


#include "config.h"

/**
*   @defgroup develop Development And Experimental Code
*   @brief Here we get things done!
*   @{
**/


/**
* @brief Computes the potentials of every vertex
*
* Takes a @p graph, a source vertex @p s, a sink vertex @p t and computes and stores
* the vertex potentials in @p vertex_potential. To this end it uses the edge lenght
* and edge diameter given by @p edge_length and @p edge_diameter respectively to
* solve a system of linear equations. The method in use is a variant of the Gauss-Seidel
* algorithm. The algorithm terminates once the maximum potential difference of a vertex between
* two iterations is below a fixed epsilon, named @p desired_precision.
*
*
* @tparam Graph                 A type that models boost::GraphConcept
* @tparam EdgeDiameterMap       A type that models a Readable property map
* @tparam EdgeLengthMap         A type that models a Readable property map
* @tparam VertexPotentialMap    A type that models a Read/Write property map
*
* @param[in]    graph                   An empty mutable graph
* @param[in]    s                       Source
* @param[in]    t                       Sink
* @param[out]   vertex_potential        Describes the vertex potential property
* @param[in]    edge_length             Describes the edge length property
* @param[in]    edge_diameter           Describes the edge diameter property
* @param[in]    desired_precision       Specifies the desired precision
*
* @todo Figure out how to put a description of the linear system in the comments.
* @todo This gauss-seidel type of method is strange. Maybe a more robust/transparent method should replace it. Check the blas part of boost.
* @todo Make up your mind about the stability/solvability of the linear systems produced. E.g. what about the conditioning of the resulting matric? Can it screw up?
* @todo Question: Is vertex_potential[ source ] always the largest possible potential? If so, there should be an assertion to check for it.
*/
template <class Graph, class EdgeDiameterMap, class EdgeLengthMap, class VertexPotentialMap >
void compute_vertex_potentials(     const Graph& graph,
                                    const typename boost::graph_traits< Graph >::vertex_descriptor s,
                                    const typename boost::graph_traits< Graph >::vertex_descriptor t,
                                    VertexPotentialMap vertex_potential,
                                    const EdgeLengthMap edge_length,
                                    const EdgeDiameterMap edge_diameter,
                                    const BigFloat desired_precision )  {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_descriptor                edge_descriptor;
    typedef typename Traits::vertex_iterator                vertex_iterator;
    typedef typename Traits::out_edge_iterator              out_edge_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< VertexPotentialMap, vertex_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept< EdgeLengthMap, edge_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept< EdgeDiameterMap, edge_descriptor > ));
    //Am I missing some other concept?



    BigFloat max_diff;

    do {

        max_diff = 0.0;
        vertex_iterator v_it, v_it_end;

        for ( boost::tie(v_it, v_it_end) = vertices(graph); v_it != v_it_end; ++v_it) {

                    if ( *v_it == t )
                        continue;

                    BigFloat R = 0.0;
                    BigFloat S = 0.0;
                    out_edge_iterator oe_it, oe_it_end;

                    for ( boost::tie(oe_it, oe_it_end) = out_edges(*v_it, graph); oe_it != oe_it_end; ++oe_it) {

                           vertex_descriptor adjacent = target(*oe_it, graph);
                           BigFloat edge_conductance = edge_diameter[*oe_it] / edge_length[*oe_it];

                           R += edge_conductance * vertex_potential[adjacent];
                           S += edge_conductance;
                    }


                    if (*v_it == s)
                        R += source_strength;


//
//                    if ( S < 0.00001 ) {
//                        continue;
//
//                    } else {

                        BigFloat previous_v_potential = vertex_potential[*v_it];
                        vertex_potential[*v_it] = R/S;

                        BigFloat current_diff = std::fabs(vertex_potential[*v_it] - previous_v_potential);

                        if ( current_diff > max_diff )
                            max_diff = current_diff;
//                    }
        }

//            std::cout << max_diff << " < " << desired_precision << "\n";

    } while ( max_diff > desired_precision );

    //postcondition
    BOOST_ASSERT( vertex_potential( t ) == 0.0 );
    // is vertex_potential( source ) the largest possible potential? If so, there should be an assertion to check for it

}


/**
* @brief Computes the flux through every edge
*
* Takes a @p graph, the potential of each vertex, @p vertex_potential, the diameter and the
* length of each edge, @p edge_diameter and @p edge_length respectively to compute
* the flux through each edge. In a second step the edge diameters are updated.
*
*
* @tparam Graph                 A type that models boost::GraphConcept
* @tparam EdgeDiameterMap       A type that models a Read/Write property map
* @tparam EdgeLengthMap         A type that models a Readable property map
* @tparam VertexPotentialMap    A type that models a Readable property map
* @tparam EdgeFluxMap           A type that models a Read/Write property map
*
* @param[in]    graph                   An empty mutable graph
* @param[in]    vertex_potential        Describes the vertex potential property
* @param[in]    edge_length             Describes the edge length property
* @param[out]   edge_diameter           Describes the edge diameter property
* @param[out]   edge_flux               Describes the edge diameter property
*
* @todo Figure out how to put a description of the physarum update equation in the documentation.
* @todo The way the core of the update works should be customizable from the outside. Maybe a compile time strategy pattern or some sort of policy design does the job?
* @todo Get rid of the magic numbers that appear in the code.
* @todo Future plans: All of this upddate business should be wrapped up using a clever functor and used as a plug-in to power a physarum solver object.
*/
template <class Graph, class EdgeDiameterMap, class EdgeLengthMap, class VertexPotentialMap, class EdgeFluxMap >
inline void compute_edge_flux(         const Graph& graph,
                                const VertexPotentialMap vertex_potential,
                                const EdgeLengthMap edge_length,
                                EdgeDiameterMap edge_diameter,
                                EdgeFluxMap edge_flux )  {


    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_descriptor                edge_descriptor;
    typedef typename Traits::vertex_iterator                vertex_iterator;
    typedef typename Traits::edge_iterator                  edge_iterator;
    typedef typename Traits::out_edge_iterator              out_edge_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept< VertexPotentialMap, vertex_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< EdgeFluxMap, edge_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept< EdgeLengthMap, edge_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< EdgeDiameterMap, edge_descriptor > ));
    //Am I missing some other concept?

    edge_iterator e_it, e_it_end;

    for ( boost::tie(e_it, e_it_end) = edges(graph); e_it != e_it_end; ++e_it) {

        BigFloat edge_conductance = edge_diameter[*e_it] / edge_length[*e_it] ;

        edge_flux[*e_it] = ( edge_conductance  ) * ( vertex_potential[ source(*e_it, graph) ] - vertex_potential[ target(*e_it, graph) ]  );
        edge_diameter[*e_it] = edge_diameter[*e_it] + 0.1 * ( std::fabs(edge_flux[*e_it]) - edge_diameter[*e_it]) ;

//        if ( edge_diameter[*e_it] <= 0.0000001 )
//            edge_diameter[*e_it] = 0.0;
    }

}


/**
* @brief Stores the current potential values for each vertex as the previous potential values
*
* Takes @p grap and iterates over all vertices in order to store the current potential
* values using the previous potentials property of each vertex.
*
*
* @tparam Graph                         A type that models boost::GraphConcept
* @tparam CurrentPotentialMap           A type that models a Readable property map
* @tparam PreviousPotentialMap          A type that models a Read/Write property map
*
* @param[in]    graph                   An empty mutable graph
* @param[in]    current_potential       Describes the current vertex potential property
* @param[out]   current_potential       Describes the previous vertex potential property
*/
template < typename Graph, typename CurrentPotentialMap, typename PreviousPotentialMap >
void store_current_vertex_potentials(   const Graph& graph,
                                        const CurrentPotentialMap current_potential,
                                        PreviousPotentialMap previous_potential )   {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::vertex_iterator                vertex_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<CurrentPotentialMap, vertex_descriptor> ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept<PreviousPotentialMap, vertex_descriptor> ));
    //Am I missing some other concept?

    vertex_iterator v_it, v_it_end;

    for (boost::tie(v_it, v_it_end) = vertices(graph); v_it != v_it_end; ++v_it)
       previous_potential[*v_it] = current_potential[*v_it];

}


/**
* @brief Computes the maximal potential difference between a vertex of iteration k and iteration k - 1
*
* Takes @p grap and iterates over all vertices in order to return the maximal potential difference
* between a vertex u of iteration k and the same vertex u of iteration k - 1. The vertex potential values
* corresponding to iteration k are stored in @p current_potential. @p previous_potential holds the
* potential values for each vertex as computed in iteration k - 1.
*
*
* @tparam Graph                         A type that models boost::GraphConcept
* @tparam CurrentPotentialMap           A type that models a Readable property map
* @tparam PreviousPotentialMap          A type that models a Readable property map
*
* @param[in]    graph                   An empty mutable graph
* @param[in]    current_potential       Describes the current vertex potential property
* @param[in]   previous_potential       Describes the previous vertex potential property
*
* @return max_diff                      The maximal potential difference of vertex u
*
* @todo Use proper notation to explain the history of vertices and iterations and stuff ...
*/
template < typename Graph, typename CurrentPotentialMap, typename PreviousPotentialMap >
BigFloat calculate_maximal_potential_difference(     const Graph& graph,
                                                     const CurrentPotentialMap current_potential,
                                                     const PreviousPotentialMap previous_potential)   {

    typedef boost::graph_traits< Graph >                Traits;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::vertex_iterator                vertex_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<CurrentPotentialMap, vertex_descriptor> ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept<PreviousPotentialMap, vertex_descriptor> ));
    //Am I missing some other concept?

    BigFloat max_diff = 0.0;
    BigFloat current_diff = 0.0;

    vertex_iterator v_it, v_it_end;

    for (boost::tie(v_it, v_it_end) = vertices(graph); v_it != v_it_end; ++v_it) {

      current_diff = std::fabs( previous_potential[*v_it] - current_potential[*v_it]);
            if ( current_diff  > max_diff)
                    max_diff = current_diff;
    }

    return max_diff;
}

/**
* @brief Computes the potentials of every vertex
*
* Takes a @p graph, a source vertex @p s, a sink vertex @p t and computes and stores
* the vertex potentials in @p vertex_potential. To this end it uses the edge lenght
* and edge diameter given by @p edge_length and @p edge_diameter respectively to
* solve a system of linear equations. The method in use is a variant of the Gauss-Seidel
* algorithm. The algorithm terminates once the maximum potential difference of a vertex between
* two iterations is below a fixed epsilon, named @p desired_precision.
*
*
* @tparam Graph                 A type that models boost::GraphConcept
* @tparam EdgeDiameterMap       A type that models a Readable property map
* @tparam EdgeLengthMap         A type that models a Readable property map
* @tparam VertexPotentialMap    A type that models a Read/Write property map
*
* @param[in]    graph                   An empty mutable graph
* @param[in]    s                       Source
* @param[in]    t                       Sink
* @param[out]   vertex_potential        Describes the vertex potential property
* @param[in]    edge_length             Describes the edge length property
* @param[in]    edge_diameter           Describes the edge diameter property
* @param[in]    desired_precision       Specifies the desired precision
*
* @todo Figure out how to put a description of the linear system in the comments.
* @todo This gauss-seidel type of method is strange. Maybe a more robust/transparent method should replace it. Check the blas part of boost.
* @todo Make up your mind about the stability/solvability of the linear systems produced. E.g. what about the conditioning of the resulting matric? Can it screw up?
* @todo Question: Is vertex_potential[ source ] always the largest possible potential? If so, there should be an assertion to check for it.
*/
template <class Graph, class EdgeDiameterMap, class EdgeLengthMap, class VertexPotentialMap >
void compute_vertex_potentials2(     const Graph& graph,
                                    const typename boost::graph_traits< Graph >::vertex_descriptor s,
                                    const typename boost::graph_traits< Graph >::vertex_descriptor t,
                                    VertexPotentialMap vertex_potential,
                                    const EdgeLengthMap edge_length,
                                    const EdgeDiameterMap edge_diameter,
                                    const BigFloat desired_precision )  {

    typedef boost::graph_traits< Graph >                    Traits;
    typedef typename Traits::vertex_descriptor              vertex_descriptor;
    typedef typename Traits::edge_descriptor                edge_descriptor;
    typedef typename Traits::vertex_iterator                vertex_iterator;
    typedef typename Traits::out_edge_iterator              out_edge_iterator;

    //concept checks
    BOOST_CONCEPT_ASSERT(( boost::GraphConcept< Graph > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadWritePropertyMapConcept< VertexPotentialMap, vertex_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept< EdgeLengthMap, edge_descriptor > ));
    BOOST_CONCEPT_ASSERT(( boost::ReadablePropertyMapConcept< EdgeDiameterMap, edge_descriptor > ));
    //Am I missing some other concept?



    BigFloat max_diff;

    do {

        max_diff = 0.0;
        vertex_iterator v_it, v_it_end;

        for ( boost::tie(v_it, v_it_end) = vertices(graph); v_it != v_it_end; ++v_it) {

                    if ( *v_it == t )
                        continue;

                    BigFloat R = 0.0;
                    BigFloat S = 0.0;
                    out_edge_iterator oe_it, oe_it_end;

                    for ( boost::tie(oe_it, oe_it_end) = out_edges(*v_it, graph); oe_it != oe_it_end; ++oe_it) {

                           vertex_descriptor adjacent = target(*oe_it, graph);
                           BigFloat edge_resistance = edge_length[*oe_it] / edge_diameter[*oe_it];

                            if (edge_resistance < desired_precision)
                                continue;

                           R += edge_resistance * vertex_potential[adjacent];
                           S += edge_resistance;
                    }


                    if (*v_it == s)
                        R += source_strength;



                    if ( S < 0.00001 ) {
                        continue;

                    } else {

                        BigFloat previous_v_potential = vertex_potential[*v_it];
                        vertex_potential[*v_it] = R/S;

                        BigFloat current_diff = std::fabs(vertex_potential[*v_it] - previous_v_potential);

                        if ( current_diff > max_diff )
                            max_diff = current_diff;
                    }
        }

//            std::cout << max_diff << " < " << desired_precision << "\n";

    } while ( max_diff > desired_precision );

    //postcondition
    BOOST_ASSERT( vertex_potential( t ) == 0.0 );
    // is vertex_potential( source ) the largest possible potential? If so, there should be an assertion to check for it

}


/** @} */

#endif // DEVELOP_H_INCLUDED
