#ifndef CUSTOM_EDGE_PROPERTIES_H_INCLUDED
#define CUSTOM_EDGE_PROPERTIES_H_INCLUDED

/**
* @file    custom_edge_properties.h
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @ingroup edge_config
* @brief Configuration file containing important edge specific controls
*
* This header file defines the edge properties needed in the graph.
* Note, that a mix of user-defined and boost properties is used.
* Furthermore a filter-class is set up intended to be used in
* the filtered_graph type.
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

#include "config.h"

/**
*   @addtogroup edge_config
*   @{
**/

/**
*    @brief User defined edge diameter tag
**/
struct diameter_t {
  typedef boost::edge_property_tag kind;
};

/**
*    @brief User defined edge flux tag
**/
struct flux_t {
  typedef boost::edge_property_tag kind;
};

/**
*    @brief User defined edge length tag
**/
struct length_t {
  typedef boost::edge_property_tag kind;
};

/**
*    @brief User defined edge color tag
**/
struct edge_color_t {
  typedef boost::edge_property_tag kind;
};

/**
*    @brief User defined edge transposed tag
**/
struct transposed_t {
  typedef boost::edge_property_tag kind;
};


typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;

typedef boost::property <       diameter_t,     BigFloat,
        boost::property <       flux_t,         BigFloat,
        boost::property <       length_t,       BigFloat,
        boost::property <       transposed_t,   bool,
        boost::property <       boost::edge_capacity_t, BigInt,
        boost::property <       boost::edge_residual_capacity_t, BigInt,
        boost::property <       boost::edge_reverse_t, Traits::edge_descriptor >  >  >  >   >   >   > EdgeProperties;

//boost::property< tag, type, nextproperty >




/**
*   @brief Custom filter class that can be used as a filter to destinguish
*   Transposed edges from non-transposed edges.
*
*   This class can be used as a filter to create filtered sub-views of another graph.
*   The primary idea is to use filters to distinguish between the set of edges of the whole graph  \f$ G \f$ ,
*   the subgraph of transposed edges \f$ g_t \f$ and the subgraph of non-transposed edges \f$ g \f$.
*   Note that,
*
*    \f[
*    G = g \cup g_{t}
*    \f]
*
*   holds as an invariant troughout the program.
*   The use of filtered views helps to eliminate the edge set we are not interested to work with. A filtered graph can modify
*   the contents of the unfiltered graph via references. Iteration is realized via MultipassIterators.
*
*   @todo Fix up the filter class such as to allow for filtering in both directions. Maybe via partial specialisation.
*   @todo Come up with a solid explanation of why filters are useful. Explain, that the transposed edges are needed to compute the s-t mincut
*           but are not used otherwise.
*   @todo Figure out how to use proper mathematical notation in the documentation.
*/
template <typename EdgeTransposedMap>               //defines a filter that returns true if a filtered edge has the property transposed
struct is_transposed_edge {

    is_transposed_edge() : m_transposed() {}
    is_transposed_edge(EdgeTransposedMap transposed) : m_transposed(transposed) { }
    template <typename Edge>
    bool operator()(const Edge& e) const {
        return true == get(m_transposed, e);
    }

    EdgeTransposedMap m_transposed;
};

/** @} */

#endif // CUSTOM_EDGE_PROPERTIES_H_INCLUDED
