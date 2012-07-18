#ifndef CUSTOM_VERTEX_PROPERTIES_H_INCLUDED
#define CUSTOM_VERTEX_PROPERTIES_H_INCLUDED

/**
* @file    custom_vertex_properties.h
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @ingroup vertex_config
* @brief Configuration file containing important vertex specific controls
*
* This header file defines the vertex properties needed in the graph.
* Note, that a mix of user-defined and boost properties is used.
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
*   @addtogroup vertex_config
*   @{
**/

/**
*    @brief User defined vertex potential tag
**/
struct potential_t {
  typedef boost::vertex_property_tag kind;
};

/**
*    @brief User defined previous vertex potential tag
**/
struct previous_potential_t {
  typedef boost::vertex_property_tag kind;
};

/**
*    @brief User defined vertex color tag
**/
struct vertex_color_t {
  typedef boost::vertex_property_tag kind;
};


typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;

typedef boost::property <       potential_t,        BigFloat,
        boost::property <       previous_potential_t,        BigFloat,
        boost::property <       vertex_color_t,     boost::default_color_type,
        boost::property <       boost::vertex_distance_t, BigInt,
        boost::property <       boost::vertex_predecessor_t, Traits::edge_descriptor  >  > >  >   > VertexProperties;

//boost::property< tag, type, nextproperty >

/** @} */

#endif // CUSTOM_VERTEX_PROPERTIES_H_INCLUDED
