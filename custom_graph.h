#ifndef CUSTOM_GRAPH_H_INCLUDED
#define CUSTOM_GRAPH_H_INCLUDED

/**
* @file    custom_graph.h
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @ingroup graph_config
* @brief Configuration file containing important graph specific controls
*
* This header file defines the type of graph that is used by the main program.
* Edge and vertex properties are included from their respective headers.
* In addition, a custom edge-filtered graph definition is provided.
* A number of often used typedefs intended for use in main.cpp is given.
*
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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/assert.hpp>
#include <boost/graph/transpose_graph.hpp>

#include "custom_edge_properties.h"
#include "custom_vertex_properties.h"

/**
*   @addtogroup graph_config
*   @todo Make sure to only define properties that are really needed. Not sure at the moment though.
*   @todo Create a proper documentation of the available graph options, i.e the parameters to adjacency_list.
*   @todo Link documentation of edge and vertex properties against each other.
*   @{
**/

// O( V + E) space requirement, neighbours can easily be accessed --> natural for exploring graphs
typedef boost::adjacency_list   <           boost::vecS,                               //container representing the edge-list for each of the vertices. listS == std::list
                                            boost::vecS,                               //container used to represent the vertex-list of the graph. vecS == std::vector
                                            boost::directedS,                          //class of the graph. Options are: directedS, undirectedS, and bidirectionalS
                                            VertexProperties,                          //custom vertex properties
                                            EdgeProperties,                            //custom edge properties
                                            boost::no_property,                        //for specifying property storage for the graph object
                                            boost::listS                               //container used to represent the edge-list for the graph. listS == std::list
                                >           GraphType;


/**
*   @defgroup edge_config Edge Controls
*   @ingroup graph_config
*   @brief Defines the edge properties
*   @{
**/
typedef boost::graph_traits< GraphType >::edge_descriptor       edge_descriptor;
///describes the number of edges in the graph
typedef boost::graph_traits< GraphType >::edges_size_type       edges_size_type;
typedef boost::graph_traits< GraphType >::edge_iterator     edge_iterator;
typedef boost::graph_traits< GraphType >::out_edge_iterator     out_edge_iterator;

/**
*   @name Internal edge properties
*   @brief The following property maps are used together with a mix of user-defined and boost predefined edge properties.
*   @{
**/
///Edge diameter
typedef boost::property_map< GraphType, diameter_t >::type                          EdgeDiameterMap;    //property maps are used to acces the respective edge or vertex properties
///Edge length
typedef boost::property_map< GraphType, length_t >::type                            EdgeLengthMap;
///Flux through an edge
typedef boost::property_map< GraphType, flux_t >::type                              EdgeFluxMap;
///Edge color (needed for s-t maxflow routine)
typedef boost::property_map< GraphType, edge_color_t >::type                        EdgeColorMap;
///Edge capacity (needed for s-t maxflow routine)
typedef boost::property_map< GraphType, boost::edge_capacity_t >::type              EdgeCapacityMap;
///residual edge capacity (needed for s-t maxflow routine)
typedef boost::property_map< GraphType, boost::edge_residual_capacity_t >::type     EdgeResidualCapacityMap;
///Reversed edge (needed for s-t maxflow routine)
typedef boost::property_map< GraphType, boost::edge_reverse_t >::type               EdgeReverseMap;
///Transposed edge. Used to define a filter distinguishing between the transposed == reversed edges and the normal edges
typedef boost::property_map< GraphType, transposed_t >::type                        EdgeTransposedMap;
///Edge name (needed for s-t maxflow routine)
typedef boost::property_map< GraphType, boost::edge_name_t >::type                  EdgeNameMap;
/** @} */


/** @} */


/**
*   @defgroup vertex_config Vertex Controls
*   @ingroup graph_config
*   @brief Defines the vertex properties
*   @{
**/

typedef boost::graph_traits< GraphType >::vertex_descriptor     vertex_descriptor;    //the type for vertex representative objects
///Describes the number of vertices
typedef boost::graph_traits< GraphType >::vertices_size_type    vertices_size_type;
typedef boost::graph_traits< GraphType >::vertex_iterator   vertex_iterator;

/**
*   @name Internal vertex properties
*   @brief The following property maps are used together with a mix of user-defined and boost predefined vertex properties.
*   @{
**/
///Vertex color. (needed for s-t maxflow routine)
typedef boost::property_map< GraphType, vertex_color_t >::type              VertexColorMap;
///Vertex potential
typedef boost::property_map< GraphType, potential_t >::type                 VertexPotentialMap;
///Previous vertex potential
typedef boost::property_map< GraphType, previous_potential_t >::type        VertexPreviousPotentialMap;
///Vertex index. Can only be used with vecS selected as vertex storage
typedef boost::property_map< GraphType, boost::vertex_index_t  >::type      VertexIndexMap;         //this property comes automatically with vecS as vertex representation

/** @} */

/** @} */

/**
*   @defgroup filter_config Filtered Graph Controls
*   @ingroup graph_config
*   @brief Defines the properties of a filtered graph.
*   @todo Add a filtered graph containing the non-transposed vertices.
*   @todo Find a proper naming system for the transposed view, the non-transposed view and the union (the whole graph). Make sure that the respective
*   properties are available for use in the other parts of the program.
*   @{
**/

typedef boost::filtered_graph<GraphType, is_transposed_edge<EdgeTransposedMap> > TransposedGraphType;       //defines a view of a graph according to a filter

/** @} */

/** @} */

#endif // CUSTOM_GRAPH_H_INCLUDED
