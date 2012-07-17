#ifndef KURTSMAIN_H_INCLUDED
#define KURTSMAIN_H_INCLUDED

/**=========================================================================
// Copyright (c)  Max-Planck-Institute Saarbruecken (Germany)
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Michael Dirnberger <mtd@mpi-inf.mpg.de>
//
//=========================================================================*/




#include "Physarum.h"

int main() {

    int number_of_graphs = 0, number_of_correct_graphs = 0, number_of_violations = 0;

    ///while (true)
    for (int k = 0; k < 100; k++) {

        NT::set_precision(256);
        newgraph:GRAPH<int,double> G;           //associates an int (lenght) and a double (diameter) with all nodes and edges of the graph

        /**
        G.read("example8-16.gw");
        int n = G.number_of_nodes();
        int m = G.number_of_edges();

        */

        //int n = 4; int m = 5;


        int n = 10;            //number or nodes
        int m = 3*n;            //number of edges

        random_source S;

        random_graph(G,                 // graph
                     n,                 // number of nodes
                     m,                 // number of edges
                     false,             // anti-parallel edges
                     true,              // loops
                     false);            // self-loops

        node s = G.first_node();
        node t = G.last_node();

        /** Physarum works on an undirected graph. I need to compute minimum cuts in
        this graph. LEDA's maximum flow (minimum cut) procedure works on a directed graph.
        I therefore make by graph bidirected by adding edge reversals. I make sure that both
        edges have the same diameter.

        When computing flows, I can treat the edges in a pair separately.
        */


        list<edge> NE;                          //edge list

        G.make_map(NE);                         //fill the edge list according to the graph
        m = G.number_of_edges();


        edge_array<NT>  L(G),                   //adds a bigfloat edge property L to G == length
                        D(G),                   //adds a bigfloat edge property D to G == diameter
                        Q(G),                   //adds a bigfloat edge property Q to G == flux
                        predicted_rate(G),      //adds a bigfloat edge property predicted_rate to G == ???
                        actual_rate(G),         //adds a bigfloat edge property actual_rate to G == ???
                        Q_prev(G);              //adds a bigfloat edge property Q_prev to G == ???

        edge_array<int> number_of_changes(G);      //adds a int edge property number_of_changes to G

        edge_array<bool>    e_in_G(G),              //adds a bool edge property e_in_G to G
                            in_SP(G);               //adds a bool edge property in_SP to G == ???

        node_array<NT>      p(G),                   //adds a bigfloat node property p to G == ???
                            predicted_p(G);         //adds a bigfloat node property predicted_p to G == ???

        node_array<bool>   v_in_G(G);               //adds a bool node property v_in_G to G



        int     Lmax = 100,             //fix maximum edge length
                Dmax = 100;             //fix maximum edge diameter

        edge e;                         //dummy edge that faciliates iteration over all edges in forall_edges
        node v;                         //dummy node that faciliates iteration over all nodes in forall_nodes


        forall_edges(e,G) {             //all edges shall be in the graph at the start
            e_in_G[e] = true;
        }

        forall_nodes(v,G) {             //all nodes shall be in the graph at the start
            v_in_G[v] = true;
        }

        forall_edges(e,G) {             //randomly set diameter and length of all edges
            L[e] = S(1,Lmax);                   //lenght between 1 and Lmax
            D[e]= S(1,Dmax);                    //diameter between 1 and Dmax
            number_of_changes[e] = 0;           //no changes are recorded at the start
        }

        forall_edges(e,G) {             //put all respective reversed edges in the graph
            edge f = G.reverse(e);
            L[e] = L[f];
            D[e] = D[f];
        }



        //G.write();

        NT      PH_D_st,    // some bigfloat == ??
                D_st;       // another bigfloat == ??

        float T = used_time();                  //start the timer

        int m_SP = predict_potentials(G, n, m, s, t, L, e_in_G, v_in_G, D, Q, p, predicted_p, in_SP, D_st, predicted_rate);

        std::cout << "\n\ntime to predict potentials = " << used_time(T);



        if ( m_SP == 0 ) {
            continue;
        }

        forall_edges(e,G) if (in_SP[e]) in_SP[G.reversal(e)] = true;

        bigfloat::set_output_precision(10);


        NT MAX_PREDICTED_RATE = -1;
        forall_edges(e,G)
        if ( !in_SP[e] && (predicted_rate[e] > MAX_PREDICTED_RATE) ) MAX_PREDICTED_RATE = predicted_rate[e];

        std::cout <<"\n\nmax predicted rate = " << (1 + MAX_PREDICTED_RATE);


        NT epsilon = 1.0*I0/(n*m*Lmax*Dmax);

        number_of_graphs++;

        NT uncertainty_in_potentials = 0.00001;

        initialize_potentials(G,
                              n,
                              m,
                              s,
                              t,
                              L,            //edge lenght property
                              e_in_G,       //edges in the graph
                              v_in_G,       //vertices in the graph
                              D,            //edge diameter property
                              Q,            //edge flux property
                              p,            //node potential property
                              uncertainty_in_potentials);

        std::cout << "\n\ntime to intialize potentials =  " << used_time(T);

        NT NC_PREV = 0;
        NT MC_PREV;
        bool flag = true;

        NT MC_initial;

        MC_initial = MinCut(G,D,in_SP,s,t);
        MC_PREV = MC_initial;

        int dir = ( MC_initial < I0/2.0 ? 1 : -1);

        std::cout << "\n\ntime to compute min-cut = " << used_time(T);

        std::cout << "\n\ndir is " << dir;

        NT M = 1;
        NT M_PREV = 1;

        NT DIFF, DIFF_PREV;

        int max_number_of_changes = 0;
        edge max_edge = G.first_edge();


        double i = 0;
        while (true) {
            std::cout << "\n\n\niteration " << i;
            if ( update(    G,
                            n,
                            m,
                            s,
                            t,
                            L,
                            e_in_G,
                            v_in_G,
                            D,
                            Q,
                            p,
                            uncertainty_in_potentials,
                            number_of_changes,
                            predicted_rate,
                            max_number_of_changes,
                            max_edge) < 0.00001 ) break;

            //statistics_on_changes(G,n,m,s,t,L,e_in_G,v_in_G,D,Q,p,uncertainty_in_potentials,number_of_changes,predicted_rate);

            NT MC = MinCut(G,D,in_SP,s,t);

            if ( dir > 0 && MC < MC_PREV - 0.01 ) {
                std::cout <<"\n\nexpected Mincut to increase " << MC;
                leda::wait(1.0);
            }

            if ( dir < 0 && MC > MC_PREV + 0.01 ) std::cout <<"\n\nexpected Mincut to decrease";

            if ( MC < MC_initial && MC < I0/4 ) {
                std::cout <<"\n\nMC dropped below I0/4";
            }


            //std::cout << "\n\niteration i " << i << "  MinCut = " << MC;

            NT NC = normalized_cost(G, s, t, L, e_in_G, v_in_G, D, p, MC );

            if (NC > NC_PREV + 0.01 && i > 1 ) {
                number_of_violations++;
                std::cout << "\n\nFALSE " << MC << " " << NC << " " << NC_PREV << " " << i << " " << number_of_violations;
                flag = false;
                int count = 1;
                forall_nodes(v,G) G[v] = count++;
                forall_edges(e,G) G[e] = L[e].to_double();

                G.write("example.gw");

                forall_edges(e,G) G[e] = D[e].to_double();

                G.write("example1.gw");

            }
            M = max_ratio(G, in_SP, e_in_G, D, MC);

            std::cout << "\n\nmax ratio off/on = " << M << " ratio to prev  " << M/M_PREV;



            std::cout << "\n\nnormalized cost (should converge to " << 2*D_st << ") = " << NC;

            DIFF = leda::abs(NC - 2*D_st);
            std::cout << "\n\ndistance to optimum " << DIFF << "   ratio of improvement = " << DIFF/DIFF_PREV;
            leda::wait(3.0);

            DIFF_PREV = DIFF;



            M_PREV = M;

            NC_PREV = NC;
            MC_PREV = MC;


            i++;
        }

        std::cout << "\n\n" << used_time(T);

        PH_D_st = analysis_of_network(G, n, s, t, L, in_SP, e_in_G, v_in_G, D, Q, p, D_st);

        if ( post_mortem_analysis(G,n,m,s,t,L,e_in_G,v_in_G,D,Q,p,predicted_p,predicted_rate,0.01) ) number_of_correct_graphs++;


        std::cout<<"\n\nThis took " << i << " iterations.\n\nD_st = " << D_st << " PH_D_st = " << PH_D_st << "\nmax number of changes in this run = " << max_number_of_changes << " predicted rate = " << predicted_rate[max_edge] << " predicted potential diff divided by length= " << (abs(predicted_p[G.source(max_edge)] - predicted_p[G.target(max_edge)]))/L[max_edge] - 1;

        std::cout<< "\n\nSo far we tried " << number_of_graphs << " graphs out of which Physarum did " << number_of_correct_graphs << "\n\n";



        std::cout<<"\n\nThis took " << i << " iterations.\n\nSo far we tried " << number_of_graphs << " graphs. Number of violations = " << number_of_violations;


    }
    std::cout << "\n\n";
    return 0;
}






#endif // KURTSMAIN_H_INCLUDED
