#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED

/**
* @file    config.h
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @ingroup Configuration
* @brief Configuration file containing important program controls
*
* This header file is intended to bundle the control parameters of the program.
* The hope is that interesting settings can be explored more easily.
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

/**
*   @defgroup config Program Configuration
*   @brief Parameters and types that control the program
*
*   The purpose of this module is to group together parameters and types that define
*   the function of the program. Examples include arithmetic types, random number generation and
*   the definition of graph properties.
*
*   @{
**/

#include <ctime>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>

/**
*   @defgroup number_types High Precision Number Types
*   @brief Number types allowing for high precision arithmetic
*   @todo Replace @c long double and @c long int with something better. Options are: LEDA types or the gnu bigfloat library.
*   @warning Kurt says that high precision arithmetic is absolutely essential. Do not trust the output based on current
*   number types. However, I do not understand currently at which point in the program the need for more precision does arise.
*   @{
**/
typedef long double BigFloat;   ///<High precision @c float type
typedef long int BigInt;        ///<High precision @c int type
/** @} */

/**
*   @defgroup random_numbers Pseudo Random Number Generation
*   @brief This sub-module groups together the parts of the program
*   which control the generation of pseudo random numbers and the
*   the specific dsitributions used for actual sampling.
*   @note All boost generators sample closed intervalls.
*   @{
**/
///The type of random number generator
typedef boost::mt19937 RandomNumberGeneratorType;               //mt == Mersenne Twister
///Uniform real distribution used for random edge lengths
typedef boost::random::uniform_int_distribution<> EdgeLengthDistribution;
///Uniform real distribution used for random edge diameters
typedef boost::uniform_real<> EdgeDiameterDistribution;



/**
*   @name Random Seed
*   @todo Replace std::time(NULL) with a better seed. E.g. the output of a weaker rng seeded with std::time(NULL)
*   @{
**/
const unsigned int random_seed = 3;
//const unsigned int random_seed = std::time(NULL);
/** @} */

/**
*   @name Edge lenght range
*   @{
**/
///lower bound of range    (included)
const double min_edge_lenght = 1;
///upper bound of range    (included)
const double max_edge_lenght = 10;
/** @} */

/**
*   @name Edge diameter range
*   @{
**/
///lower bound of range   (included)
const double min_edge_diameter = 0;
///upper bound of range    (included)
const double max_edge_diameter = 1;
/** @} */
/** @} */

/**
*   @defgroup graph_config Graph Controls
*   @brief
*
*   This sub-module groups together all graph specific
*   program controls including a large number of useful
*   or at least necessary typedefs.
*   @{
**/

/** @name Basic controls
*   @{
**/
///Number of vertices in the graph
unsigned int n_vertices = 10;
///Number of edges in the graph
unsigned int n_edges = 2*n_vertices;
///Source strenght
unsigned int source_strength = 1;
///This specifies the error we allow for the values of the vertex potentials.
BigFloat uncertainty_in_potentials = 0.0000000001;
/** @} */

/** @} */
#endif // CONFIG_H_INCLUDED
