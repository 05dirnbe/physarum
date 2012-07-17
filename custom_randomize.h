#ifndef CUSTOM_RANDOMIZE_H_INCLUDED
#define CUSTOM_RANDOMIZE_H_INCLUDED

/**
* @file    custom_randomize.h
* @ingroup random_numbers
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

#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

#include "config.h"

/**
*   @addtogroup random_numbers
*   @{
**/

RandomNumberGeneratorType rng(random_seed);

EdgeLengthDistribution edge_length_distribution(min_edge_lenght,max_edge_lenght);
EdgeDiameterDistribution edge_diameter_distribution(min_edge_diameter,max_edge_diameter);

boost::variate_generator<RandomNumberGeneratorType&, EdgeLengthDistribution >       random_edge_length( rng, edge_length_distribution );
boost::variate_generator<RandomNumberGeneratorType&, EdgeDiameterDistribution >       random_edge_diameter( rng, edge_diameter_distribution );

/** @} */
#endif // CUSTOM_RANDOMIZE_H_INCLUDED
