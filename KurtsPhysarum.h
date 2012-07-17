#ifndef KURTSPHYSARUM_H_INCLUDED
#define KURTSPHYSARUM_H_INCLUDED

#ifndef PHYSARUM_H_INCLUDED
#define PHYSARUM_H_INCLUDED

/**=========================================================================
// Copyright (c)  Max-Planck-Institute Saarbruecken (Germany)
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Michael Dirnberger <mtd@mpi-inf.mpg.de>
//
//=========================================================================*/

#include <LEDA/graph/graph.h>
#include <LEDA/numbers/matrix.h>
#include <LEDA/numbers/vector.h>
#include <LEDA/numbers/bigfloat.h>
#include <LEDA/numbers/integer.h>

#include <LEDA/core/array.h>
#include <LEDA/core/stack.h>
#include <LEDA/core/b_queue.h>
#include <LEDA/graph/node_pq.h>
#include <LEDA/system/assert.h>
#include <LEDA/core/random_source.h>
#include <LEDA/graph/templates/max_flow.h>


const double I0 = 400;

using namespace leda;

typedef three_tuple<bigfloat,node,list<edge> > simple_path;

int my_compare(const simple_path& p1, const simple_path& p2) {
    return compare(p1.first(),p2.first());
}

double abs(double x) {
    return fabs(x);
}

typedef bigfloat NT;

//kurts custom dijkstra routine
template<class NT>
NT DIJKSTRA_T(const graph& G,                           //a graph
              node s,                                   //source node
              node t,                                   //sink node
              const edge_array<bool>& e_in_G,           //bool edge property is e in G
              const edge_array<NT>& cost,               //bigfloat edge property cost of the edges
              node_array<NT>& dist,                     //bigfoat node property distance
              node_array<edge>& pred) {                 //node property edge to predecessor

    node_pq<NT> PQ(G);                      //define a node priority queue
    PQ.insert(s,0);                         //insert s

    dist[s] = 0;                            //set the distance from s to s to 0

    node v;
    edge e;
    forall_nodes(v,G) pred[v] = nil;        //set all edges to predecessors to nil

    assert(s != t);                         //check for null path

    while ( !PQ.empty() ) {

        node u = PQ.del_min();

        forall_inout_edges(e,u) {

            if ( ! e_in_G[e] )
                continue;

            node v = G.opposite(e,u);
            NT c = dist[u] + cost[e];

            if ( pred[v] == nil && v != s)
                PQ.insert(v,c); // first path to v

            else if (c < dist[v])
                PQ.decrease_p(v,c); // better path

            else
                continue;

            dist[v] = c;
            pred[v] = e;
        }
    }

    return dist[t];
}

/*
template<class NT>
NT PH_DIJKSTRA_T(const graph& G, node s, node t,
                 const edge_array<NT>& cost,
                 const edge_array<NT> D,
                 node_array<NT>& dist,
                 node_array<edge>& pred) {
    node_pq<NT> PQ(G);
    PQ.insert(s,0);

    dist[s] = 0;

    node v;
    edge e;
    forall_nodes(v,G) pred[v] = nil;

    assert(s != t);

    int n = G.number_of_nodes();

    while ( !PQ.empty() ) {
        node u = PQ.del_min();
        forall_inout_edges(e,u) {
            if (D[e] < 4/(n*n)) continue;
            node v = G.opposite(e,u);
            NT c = dist[u] + cost[e];
            if ( pred[v] == nil && v != s)
                PQ.insert(v,c); // first path to v
            else if (c < dist[v]) PQ.decrease_p(v,c); // better path
            else continue;
            dist[v] = c;
            pred[v] = e;
        }
    }
    assert(pred[t] != nil);

    return dist[t];

}
*/

template<class NT>
bool DECOMPOSITION_DIJKSTRA_T(const graph& G, node s, node t,
                              const edge_array<NT>& L,
                              const edge_array<NT> Q,
                              node_array<NT>& dist,
                              node_array<edge>& pred,
                              NT eps) {
    node_pq<NT> PQ(G);
    PQ.insert(s,0);

    dist[s] = 0;

    node v;
    edge e;
    forall_nodes(v,G) pred[v] = nil;


    while ( !PQ.empty() ) {
        node u = PQ.del_min();
        forall_inout_edges(e,u) {
            if ( abs(Q[e]) < eps ) continue;
            node v = G.opposite(e,u);
            if ( ( Q[e] > 0 && u == G.target(e) ) || ( Q[e] < 0 && u == G.source(e) ) ) continue;
            NT c = dist[u] + L[e];
            if ( pred[v] == nil && v != s)
                PQ.insert(v,c); // first path to v
            else if (c < dist[v]) PQ.decrease_p(v,c); // better path
            else continue;
            dist[v] = c;
            pred[v] = e;
        }
    }


    return ( pred[t] != nil );

}

/*
//I am removing all edges from the graph whose conductance and flow is below a certain epsilon.
//If all edges incident to a vertex are removed, I am also removing the vertex. If epsilon is less than
//1/n it should never happen that I remove s or t
//After pruning small edges, the graph my become disconnected. I should only keep the nodes that
//are reachable from s
void prune_graph(const graph& G, node s, node t, const edge_array<NT>& D, const edge_array<NT>& Q,
                 const edge_array<bool>& in_SP, edge_array<bool>& e_in_G,
                 node_array<bool>& v_in_G, int& n, int& m, NT epsilon)

{
    node v;
    edge e;
    bool did_prune = false;
    forall_edges(e,G) {
        if ( ! e_in_G[e]) continue;   // e was previously removed
        if ( (abs(Q[e]) < epsilon) && (D[e] < epsilon)) {
            e_in_G[e] = false;
            m--;
            did_prune = true;
            //assert ( !in_SP[e] );
            if ( in_SP[e] ) std::cout << "\nremoved edge of shortest path";
            std::cout << "\nremoved edge";
        }
    }
    forall_nodes(v,G) {
        if ( ! v_in_G[v] ) continue;   // v was previously removed
        bool flag = false;
        forall_inout_edges(e,v) {
            if ( e_in_G[e] ) flag = true;
        }
        if ( flag == false ) {
            v_in_G[v] = false;
            n--;
            did_prune = true;
            std::cout << "\nremoved node";
        }
    }

    assert(v_in_G[s] && v_in_G[t]);

    node_array<bool> reachable(G,false);
    forall_nodes(v,G) reachable[v] = false;
    stack<node> R;
    R.push(s);
    reachable[s] = true;
    while ( ! R.empty() ) {
        v = R.pop();
        forall_inout_edges(e,v) {
            if ( ! e_in_G[e] ) continue;   // only use edges that are still present
            node w = G.opposite(e,v);
            if ( ! reachable[w] ) {
                reachable[w] = true;
                R.push(w);
            }
        }
    }

    forall_nodes(v,G) {
        if ( reachable[v] || ! v_in_G[v] ) continue;
        // v is still in graph, but not reachable, remove it and all incident edges;
        forall_inout_edges(e,v) {
            if ( ! e_in_G[e]) continue;   // e was previously removed
            e_in_G[e] = false;
            m--;
        }
        v_in_G[v] = false;
        n--;
        did_prune = true;
    }
    assert(v_in_G[s] && v_in_G[t]);

    if ( did_prune ) std::cout << "\n\nI am down to " << n << " nodes and " << m << " edges.";
}
*/

/*
//I should only keep the nodes that are reachable from s
void prune_graph_modified(const graph& G, node s, node t, const edge_array<NT>& D, const edge_array<NT>& Q,
                          const edge_array<bool>& in_SP, edge_array<bool>& e_in_G,
                          node_array<bool>& v_in_G, int& n, int& m) {
    node v;
    edge e;
    bool did_prune = false;

    node_array<bool> reachable(G,false);
    forall_nodes(v,G) reachable[v] = false;
    stack<node> R;
    R.push(s);
    reachable[s] = true;
    while ( ! R.empty() ) {
        v = R.pop();
        forall_inout_edges(e,v) {
            if ( ! e_in_G[e] ) continue;   // only use edges that are still present
            node w = G.opposite(e,v);
            if ( ! reachable[w] ) {
                reachable[w] = true;
                R.push(w);
            }
        }
    }

    forall_nodes(v,G) {
        if ( reachable[v] || ! v_in_G[v] ) continue;
        // v is still in graph, but not reachable, remove it and all incident edges;
        forall_inout_edges(e,v) {
            if ( ! e_in_G[e]) continue;   // e was previously removed
            e_in_G[e] = false;
            m--;
        }
        v_in_G[v] = false;
        n--;
        did_prune = true;
    }
    assert(v_in_G[s] && v_in_G[t]);

    if ( did_prune ) std::cout << "\n\nI am down to " << n << " nodes and " << m << " edges.";
}
*/

/*
//I am computing the potentials. I am restricting
//myself to the network of edges and vertices that are still present. In order to connect vertices to indices,
//I number the nodes that are still present.
//I am setting up a dangerous system as all row sums are equal to zero, except for the row corresponding to t.
//Maybe the system becomes more stable, if I give node t the number one.
void initialize_potentials_old(const graph& G, int n, node s, node t, const edge_array<NT>& L,
                   const edge_array<bool>& e_in_G, const node_array<bool>& v_in_G,
                  edge_array<NT>& D, edge_array<NT>& Q, node_array<NT>& p, NT epsilon)  {

matrix M(n,n); vector b(n);
  node v; edge e;
  for (int i = 0; i < n; i++)
  { b[i] = 0;
    for (int j = 0; j < n; j++) M(i,j) = 0;
  }
  node_array<int> num(G);
  num[t] = 0;
  int count = 1;
  forall_nodes(v,G) { if ( v_in_G[v] && v != t ) num[v] = count++; }
  assert(count == n);

for every node v different from t, I am generating the equation corresponding to flow conservation.
The right hand side is zero, except for s, where it is I0, the desired flow from s to t.
For t, I am generating the equation, p[t] = 0
  forall_nodes(v,G)
  { if ( ! v_in_G[v] ) continue;
    int i = num[v];
    if ( v == t ) { M(i,i) = 1; b[i] = 0; }
    else
    { NT S = 0;
      forall_inout_edges(e,v)
      { if ( ! e_in_G[e] ) continue;
        node w = G.opposite(e,v);
        assert( ! (v == w)); // no self-loops
        M(i,num[w]) = -D[e]/L[e];
        S += D[e]/L[e];
      }
      M(i,i) = S;
      if ( v == s ) b[i] = I0;
    }
  }
  vector p0(n);
  p0 = M.solve(b);
  forall_nodes(v,G) if ( v_in_G[v] ) p[v] = p0[num[v]];
}

*/

/*
//I am computing the potentials and the flows and then update the conductances. I am restricting
//myself to the network of edges and vertices that are still present. In order to connect vertices to indices,
//I number the nodes that are still present.
//I am setting up a dangerous system as all row sums are equal to zero, except for the row corresponding to t.
//Maybe the system becomes more stable, if I give node t the number one.
void one_iteration(const graph& G, int n, node s, node t, const edge_array<NT>& L,
                   const edge_array<bool>& e_in_G, const node_array<bool>& v_in_G,
                  edge_array<NT>& D, edge_array<NT>& Q, node_array<NT>& p, NT epsilon){

                      matrix M(n,n); vector b(n);
  node v; edge e;
  for (int i = 0; i < n; i++)
  { b[i] = 0;
    for (int j = 0; j < n; j++) M(i,j) = 0;
  }
  node_array<int> num(G);
  num[t] = 0;
  int count = 1;
  forall_nodes(v,G) { if ( v_in_G[v] && v != t ) num[v] = count++; }
  assert(count == n);

for every node v different from t, I am generating the equation corresponding to flow conservation.
The right hand side is zero, except for s, where it is I0, the desired flow from s to t.
For t, I am generating the equation, p[t] = 0
  forall_nodes(v,G)
  { if ( ! v_in_G[v] ) continue;
    int i = num[v];
    if ( v == t ) { M(i,i) = 1; b[i] = 0; }
    else
    { NT S = 0;
      forall_inout_edges(e,v)
      { if ( ! e_in_G[e] ) continue;
        node w = G.opposite(e,v);
        assert( ! (v == w)); // no self-loops
        M(i,num[w]) = -D[e]/L[e];
        S += D[e]/L[e];
      }
      M(i,i) = S;
      if ( v == s ) b[i] = I0;
    }
  }
  vector p0(n);
  p0 = M.solve(b);
  forall_nodes(v,G) if ( v_in_G[v] ) p[v] = p0[num[v]];
  forall_edges(e,G) if ( e_in_G[e] ) Q[e] = (D[e]/L[e])*(p0[num[G.source(e)]] - p0[num[G.target(e)]]);
  forall_edges(e,G) if ( e_in_G[e] )
   { D[e] = D[e] + epsilon*(abs(Q[e]) - D[e]);  if (D[e] < 0 ) D[e] = 0; }
}
*/


//I am computing the potentials. I am restricting
//myself to the network of edges and vertices that are still present. In order to connect vertices to indices,
//I number the nodes that are still present.
//I am setting up a dangerous system as all row sums are equal to zero, except for the row corresponding to t.
//Maybe the system becomes more stable, if I give node t the number one.
void initialize_potentials(const graph& G, int& n, int& m, node s, node t, const edge_array<NT>& L,
                           edge_array<bool>& e_in_G, node_array<bool>& v_in_G,
                           edge_array<NT>& D, edge_array<NT>& Q, node_array<NT>& p, NT epsilon){
    node v;
    edge e;
    forall_nodes(v,G) p[v] = 0;

    while ( true ) {

        NT maxdiff = 0;
        forall_nodes(v,G) {

            if ( (! v_in_G[v]) || v == t ) continue;  // p[t] = 0 always

                NT S = 0;
                NT R = 0;

                forall_inout_edges(e,v) {

                    if ( ! e_in_G[e] ) continue;

                        node w = G.opposite(e,v);
                        assert( ! (v == w)); // no self-loops
                        NT g = D[e]/L[e];
                        if ( g < epsilon ) {
                            e_in_G[e] = false;
                            m--;
                            continue;
                        }

                    R += g * p[w];
                    S += g;
                }

                if ( v == s ) R += I0;
                if ( S == 0 ) continue;
                /*{ v_in_G[v] = false; n--;
                  forall_inout_edges(e,v)
                    if ( e_in_G[e] ) { e_in_G[e] = false; m--; }
                }
                */
                else {

                    NT oldp = p[v];
                    p[v] = R/S;
                    if ( abs(p[v] - oldp) > maxdiff ) maxdiff = abs(p[v] - oldp);
                }
        }
        if (maxdiff < epsilon) break;
    }
}

NT normalized_cost(const graph& G, node s, node t, const edge_array<NT>& L,
                   edge_array<bool>& e_in_G, node_array<bool>& v_in_G,
                   edge_array<NT>& D, node_array<NT>& p, NT MC) {
    node v;
    edge e;

    NT C = 0;
    NT SD = 0;

    forall_edges(e,G) {
        if ( ! e_in_G[e] ) continue;
        C += D[e] * L[e];
    }
    return C/MC;
}

//At first, I computed the maximum diameter of an edge not on the shortest path. This function does not always decrease and I should not have expected it. Even for the case of three
//links it does not. In the three-links demo, the diameter of the middle edge goes up first and
//the goes down.
//
//I therefore modified the function. I now compute the ratio of largest diameter off the shortest path divided by minimum diameter on the shortest path.
NT max_ratio(const graph& G, const edge_array<bool>& in_SP,
             const edge_array<bool>& e_in_G,
             const edge_array<NT>& D, NT MC)

{
    node v;
    edge e;

    NT MAX_NOT_IN_SP = 0;
    NT MIN_IN_SP = 50*I0;   // hack

    forall_edges(e,G) {
        assert(e_in_G[e]);
        if ( in_SP[e] && ( D[e] < MIN_IN_SP ) ) MIN_IN_SP = D[e];
        if ( !in_SP[e] && ( D[e] >  MAX_NOT_IN_SP ) ) MAX_NOT_IN_SP = D[e];
    }
    return MAX_NOT_IN_SP/MIN_IN_SP;
}

/*
NT MinCutWheatStone(const graph& G, const edge_array<NT>& D, edge e1, edge e2, edge e3, edge e4, edge e5) {
    NT C1 = D[e1] + D[e3];
    NT C2 = D[e2] + D[e5] + D[e3];
    NT C3 = D[e1] + D[e5] + D[e4];
    NT C4 = D[e2] + D[e4];
    return min(min(C1,C2),min(C3,C4));
}
*/

NT MinCut(const graph& G, const edge_array<NT>& D, const edge_array<bool>& in_SP, node s, node t) {


//    edge_array<integer> f(G);                         //this line makes the debugger cry

    edge_array<int> f(G);

    NT Large = 1;
    for (int i = 0; i < 20; i++) Large = 2*Large;
//    edge_array<integer> DI(G);                        //this line makes the debugger cry
    edge_array<int> DI(G);
    edge e;
    NT MIND = 0;
    forall_edges(e,G) {
        if ( in_SP[e] && ( MIND == 0 || D[e] < MIND )) MIND = D[e];
//        DI[e] = floor(D[e]*Large);
        DI[e] = floor(D[e]*Large).to_double() ;         //ugly fix: fist convert to double and then implicitely convert to int
    }
    NT F = bigfloat(MAX_FLOW_T(G,s,t,DI,f))/Large;
    if ( F < MIND - 1 ) {
        std::cout << "\n\n Oberfaul" << MIND << " " << F << " " << Large;
        forall_edges(e,G) std::cout << "\n" << D[e] << " " << DI[e];
        exit(1);
    }
    return F;
}

/*
//I decompose the flow; I assume that the node potentials correspond to the current D-values
//The flow induces an acyclic graph. I compose the flow into paths, by iteratively computing a shortest
//s-t path and then routing as much flow as possible on this path. I return the quantity
//                    sum_P  flow(P)/length(P).
//The conjecture is that this quantity goes down. Did I verify this for parallel links?
NT flow_decomposition(const graph& G, node s, node t, const edge_array<NT>& L,
                      edge_array<bool>& e_in_G, node_array<bool>& v_in_G,
                      edge_array<NT>& D, node_array<NT>& p, NT epsilon) {
    edge_array<NT> Q(G);

    node v;
    edge e;
    forall_edges(e,G) {
        if ( !e_in_G[e] ) continue;
        Q[e] = (D[e]/L[e])*(p[G.source(e)] - p[G.target(e)]);
    }

    NT NFS = 0;   // normalized flow sum

    node_array<NT> dist(G);
    node_array<edge> pred(G);

    while ( DECOMPOSITION_DIJKSTRA_T(G, s, t, L, Q, dist, pred, epsilon) ) {
        NT maxflow = I0;
        NT TotalResistance = 0;
        node v = t;
        while ( v != s ) {
            edge e = pred[v];
            if ( abs(Q[e]) < maxflow ) maxflow = abs(Q[e]);
            TotalResistance += L[e]*log(D[e].to_double());
            v = G.opposite(e,v);
        }

        NFS += maxflow * TotalResistance;

        v = t;
        while ( v != s ) {
            edge e = pred[v];
            Q[e] += ( Q[e] < 0 ? maxflow : -maxflow );
            v = G.opposite(e,v);
        }
    }
    return NFS;
}
*/


/*

*/
void statistics_on_changes(const graph& G, int& n, int& m, node s, node t, const edge_array<NT>& L,
                           edge_array<bool>& e_in_G, node_array<bool>& v_in_G,
                           edge_array<NT>& D, edge_array<NT>& Q, node_array<NT>& p, NT epsilon,
                           edge_array<int>& number_of_changes, const edge_array<NT>& predicted_rate)    {
    //dictionary<int,list<edge> > ordered_changes;
    edge e;
    forall_edges(e,G)
    if ( number_of_changes[e] > 0 ) {
        std::cout << "\n\nnumber of changes" << number_of_changes[e] << "  predicted_rate = " << predicted_rate[e];
    }
    return;
}


//I
//-- compute the flows
//-- update the conductivities
//-- update the potentials
NT update(const graph& G, int& n, int& m, node s, node t, const edge_array<NT>& L,
          edge_array<bool>& e_in_G, node_array<bool>& v_in_G,
          edge_array<NT>& D, edge_array<NT>& Q, node_array<NT>& p, NT epsilon,
          edge_array<int>& number_of_changes, const edge_array<NT>& predicted_rate, int& max_number_of_changes, edge& max_edge) {
    NT oldSum = 0;
    NT newSum = 0;

    node v;
    edge e;
    forall_edges(e,G) {

        if ( !e_in_G[e] ) continue;

        NT oldflow = Q[e];

        Q[e] = (D[e]/L[e])*(p[G.source(e)] - p[G.target(e)]);  //this is it


        if ( oldflow*Q[e] < 0 ) {

            number_of_changes[e]++;

            std::cout << "\n\nnumber of changes is now = " << number_of_changes[e] << " predicted rate = " << predicted_rate[e];

            if ( number_of_changes[e] > max_number_of_changes ) {

                max_number_of_changes++;
                max_edge = e;
                std::cout << "\n NEX MAX ";
            }
        }

        oldSum += abs(Q[e])*L[e];
        D[e] = D[e] + 0.1*(abs(Q[e]) - D[e]);
        if (D[e] < 0 ) D[e] = 0;
    }

    node_array<NT> old_p(G);
    forall_nodes(v,G) old_p[v] = p[v];

    while ( true ) {
        NT maxdiff = 0;
        forall_nodes(v,G) {
            if ( (! v_in_G[v]) || v == t ) continue;  // p[t] = 0 always
            NT S = 0;
            NT R = 0;
            forall_inout_edges(e,v) {
                if ( ! e_in_G[e] ) continue;
                node w = G.opposite(e,v);
                assert( ! (v == w)); // no self-loops
                NT g = D[e]/L[e];
                //if ( g < epsilon ) { e_in_G[e] = false; m--; continue; }
                R += g * p[w];
                S += g;
            }
            if ( v == s ) R += I0;
            if ( S == 0 ) continue;
            /*
            { v_in_G[v] = false; n--;
              forall_inout_edges(e,v)
                if ( e_in_G[e] ) { e_in_G[e] = false; m--; }
            }
            else
            */
            {
                NT oldp = p[v];
                p[v] = R/S;
                if ( abs(p[v] - oldp) > maxdiff ) maxdiff = abs(p[v] - oldp);
            }
        }
        if (maxdiff < epsilon) break;
    }


    /*
    forall_edges(e,G)
    { if ( !e_in_G[e] ) continue;
      Q[e] = (D[e]/L[e])*(p[G.source(e)] - p[G.target(e)]);
      newSum += abs(Q[e])*L[e];
    }
    if (newSum >= oldSum + 100*epsilon ) { std::cout << "\n\n\nCONJECTURE IS FALSE\nold sum and new sum =" << oldSum << " " << newSum; leda::wait(1.0); }
    */



    NT max_diff = 0;
    forall_nodes(v,G)
    if ( abs(old_p[v] - p[v]) > max_diff ) max_diff = abs(old_p[v] - p[v]);
    return max_diff;
}

//I
//-- compare the computed potential with the predicted potentials
bool post_mortem_analysis(const graph& G, int n, int m, node s, node t, const edge_array<NT>& L,
                          const edge_array<bool>& e_in_G, const node_array<bool>& v_in_G,
                          const edge_array<NT>& D, const edge_array<NT>& Q,
                          const node_array<NT>& p, const node_array<NT>& predicted_p,
                          const edge_array<NT>& predicted_rate, NT max_diff) {
    node v;
    edge e;
    int count = 0;
    forall_nodes(v,G) {
        std::cout << "\npredicted and actual potentials: " << predicted_p[v] << " " << p[v];
        if ( v == t || predicted_p[v] < 0 || (abs(p[v] - predicted_p[v])/p[v] < max_diff) ) count++;
    }

    std::cout << "\n\nnumber of nodes for which predicted and actual potential agree " << count << " out of " << n;

    int countm = 0;
    forall_edges(e,G) {
        NT Qcurrent = abs(D[e]/L[e] *(p[G.target(e)] - p[G.source(e)]));
        double measured_rate = (Qcurrent/D[e] - 1).to_double();
        double predicted_r = predicted_rate[e].to_double();
        std::cout << "\n\nconductances " << D[e] << " flow " << Qcurrent << " predicted rate " << predicted_r << " measured_rate " << measured_rate;
        if ( abs(predicted_r - measured_rate) < max_diff ) countm++;
    }
    std::cout << "\n\nnumber of edges for which predicted and actual agree " << countm << " out of " << m;

    return ( count == n && countm == m);
}


/*
bool predict_potentials_first_try(const graph& G, int n, int m, node s, node t, const edge_array<NT>& L,
                                  edge_array<bool>& e_in_G, node_array<bool>& v_in_G,
                                  const edge_array<NT>& D, const edge_array<NT>& Q,
                                  const node_array<NT>& p, node_array<NT>& predicted_p,
                                  edge_array<bool>& in_SP,
                                  node_array<NT>& dists, node_array<edge>& preds,
                                  node_array<NT>& distt, node_array<edge>& predt, NT& D_st)
// I predict the potentials
//-- If t is not reachable from s I return 0, otherwise, I return the number of edges in the shortest
//path from s to t
//-- first I compute shortest path distances from s and t, respectively
//-- nodes that are not reachable

{
    node v;
    edge e;
    forall_nodes(v,G) {
        v_in_G[v] = true;
        dists[v] = distt[v] = predicted_p[v] = 0;
    }
    forall_edges(e,G) {
        e_in_G[e] = true;
        in_SP[e] = false;
    }

    D_st = DIJKSTRA_T(G, s, t, e_in_G, L, dists, preds);

    if (preds[t] == nil) return 0;

//    I am checking that the shortest path from s to t is unique and I am marking all edges
//    on the shortest path.
//    The shortest path is unique if for every vertex on the shortest path there is only one incoming edge
//    that defines the distance of the node

    int m_SP = 0;

    v = t;
    while ( v != s) {
        edge e = preds[v];
        m_SP++;
        in_SP[e] = true;
        node w = G.opposite(e,v);
        edge f;
        forall_inout_edges(f,v) {
            if ( f == e) continue;
            node z = G.opposite(f,v);
            if ( dists[z] + L[f] == dists[v]) return 0;
        }
        v = w;
    }
// t is reachable from s and the shortest path from s to t is unique

    std::cout << "\n\nThe shortest s-t path consists of " << m_SP << " edges and has length " << D_st;

    DIJKSTRA_T(G, t, s, e_in_G, L, distt, predt);  // distances from $t$

    assert( distt[s] == D_st );

    forall_nodes(v,G) predicted_p[v] = 0;

// I first assign potentials to the nodes on the shortest path from s to t
    v = t;
    while ( v != s) {
        edge e = preds[v];
        assert(in_SP[e]);
        node w = G.opposite(e,v);
        predicted_p[w] = predicted_p[v] + L[e];
        v = w;
    }
    assert(predicted_p[s] == D_st);



    NT current_dist = D_st;

    while (true) {
        NT smallest_new_dist = 100000;
        node newnode = nil;

        forall_nodes(v,G) std::cout << "\nnode = " << v << " predicted potential " << predicted_p[v] << " sum of distances = " << distt[v] + dists[v];


        forall_nodes(v,G) {
            NT sum = distt[v] + dists[v];  // sum == 0 means not reachable
            if ( sum == 0 || sum < current_dist || v == t || predicted_p[v] != 0 ) continue;
            if ( sum < smallest_new_dist ) {
                smallest_new_dist = sum;
                newnode = v;
            }
        }

        current_dist = smallest_new_dist;

        std::cout << "\nnew current dist = " << current_dist << " new node = " << newnode;

// the nodes on the shortest s-t path via newnode are assigned their potentials.

        if ( newnode == nil ) return m_SP;

        NT Lsum = 0;

        v = newnode;

// walk towards s

        while ( predicted_p[v] == 0 && v != s ) {
            edge e = preds[v];
            node w = G.opposite(e,v);
            Lsum += L[e];
            v = w;
        }

        NT pot_high = predicted_p[v];

// walks towards t
        stack<edge> ES;
        v = newnode;
        while ( predicted_p[v] == 0 && v != t ) {
            edge e = predt[v];
            ES.push(e);
            node w = G.opposite(e,v);
            Lsum += L[e];
            v = w;
        }


        NT pot_low = predicted_p[v];

        NT pot_diff = pot_high - pot_low;

        std::cout << "\nnew potdiff  = " << pot_diff  << " length  = " << Lsum;

// Lsum is the length of the new path

// v has a predicted potential and v is on the path from newnode to t

        while ( ! ES.empty() ) {
            e = ES.pop();
            node w = G.opposite(e,v);
            predicted_p[w] = predicted_p[v] + pot_diff * L[e]/Lsum;
            std::cout << " L[e] = " << L[e];
            v = w;
        }

        assert(v == newnode);

        while (true) {
            edge e = preds[v];
            node w = G.opposite(e,v);
            if ( predicted_p[w] != 0 ) break;
            predicted_p[w] = predicted_p[v] + pot_diff * L[e]/Lsum;
            std::cout << " L[e] = " << L[e];
            v = w;
        }

    } // end while true

}
*/

/** I predict the potentials
-- If t is not reachable from s I return 0, otherwise, I return the number of edges in the shortest
path from s to t
-- first I compute shortest path distances from s and t, respectively
-- nodes that are not reachable
*/
int predict_potentials(const GRAPH<int,double>& G,      //a bidirected graph
                       int n,                           //number of nodes
                       int m,                           //number of edges
                       node s,                          //source of the path
                       node t,                          //sink of the path
                       const edge_array<NT>& L,         //edge property edge length
                       edge_array<bool>& e_in_G,        //edge property edge in the graph
                       node_array<bool>& v_in_G,        //node property node in the graph
                       const edge_array<NT>& D,         //edge property edge diameter
                       const edge_array<NT>& Q,         //edge property edge flux
                       const node_array<NT>& p,         //node property node potential
                       node_array<NT>& predicted_p,     //node property predicted node potential
                       edge_array<bool>& in_SP,         //edge property edge in_SP
                       NT& D_st,                        //bigfloat distance(s,t)
                       edge_array<NT>& predicted_rate)  //edge property edge predicted rate
{

    edge e;                         //dummy edge that faciliates iteration over all edges in forall_edges
    node v;                         //dummy node that faciliates iteration over all nodes in forall_nodes

    forall_nodes(v,G) {
        v_in_G[v] = true;           //reset all nodes to be part of the graph
        predicted_p[v] = -1;        //reset the predicted potentials to -1 for all nodes, i.e. not processed
    }

    forall_edges(e,G) {
        e_in_G[e] = true;               //reset all edges to be part of the graph
        in_SP[e] = false;               //reset all edges to not be part of the graph
        predicted_rate[e] = -1;         //reset the predicted potentials to -1 for all edges, i.e. not processed
    }

    node_array<NT> dist(G);                 //define a bigfloat node property distance for G
    node_array<edge>   pred(G,nil);         //define a node property of type edge initialized to nil for G
    edge_array<bool>   used(G,false);       //define a bool edge property initialized to false for G

    forall_edges(e,G) assert(used[e] == false);     //check wether all edges are unused
    forall_nodes(v,G) assert(pred[v] == nil);       //check wether all nodes do not have a pred edge


    // D_st == distance(s,t)
    D_st = DIJKSTRA_T(  G,              //graph
                        s,              //start of the path
                        t,              //end of the path
                        e_in_G,         //are the nodes in G?
                        L,              //length of the paths == cost of the path
                        dist,           //will hold the distances computed
                        pred);          //will hold a list of edges to the predecessor of a certain node



    if (pred[t] == nil) {               //if the sink has no edge leading to a predecessor, it cannot be reached from the source
        return 0;
    }


    /** I am checking that the shortest path from s to t is unique and I am marking all edges
    on the shortest path. The shortest path is unique if for every vertex on the shortest path there is only one incoming edge that defines the distance of the node.
    I predict the potential for the nodes on the shortest
    path. The potential is the length of the path to t.*/

    int m_SP = 0;               //number of edges in the shortes path

    v = t;
    predicted_p[t] = 0;         //set the sink potential

    while ( v != s) {           //as long as we havent reached the source

        edge e = pred[v];                       //fetch the edge leading to the predecessor node
        m_SP++;                                 //increase the number of edges in the shortest path

        in_SP[e] = true;                        //remember that we've visited this edge
        in_SP[G.reverse(e)] = true;             //and dont forget to mark its reversed edge

        predicted_rate[e] = 0;                  //since the edge is connected to the sink, the potential must be 0
        predicted_rate[G.reverse(e)] = 0;       //the same goes for the reversed edge

        used[e] = true;                         //e has been processed
        used[G.reverse(e)] = true;              //so has its reversed counterpart

        node w = G.opposite(e,v);                   //fetch the actual predecessor node
        predicted_p[w] = predicted_p[v] + L[e];     //here the potential changes according the length of the connecting edge


        edge f;                                     //iteration dummy
        forall_inout_edges(f,v) {

            if ( f == e || (G.reverse(e) != nil && f == G.reverse(e)) )
                continue;

            node z = G.opposite(f,v);

            if ( dist[z] + L[f] == dist[v]) {                   //if i can take a detour from the determined shortest path that does not exceed its length

                std::cout << "\nshortest path is not unique";
                return 0;
            }
        }

        v = w;                                      //proceed from this node which is in a shortest path
    }


    std::cout << "\n\nThe shortest s-t path consists of " << m_SP << " edges and has length " << D_st;

    /** I am determining the best-ratio-path (largest ratio of potential difference and length)
    connecting two nodes with an a predicted potential. A node v
    has a predicted potential if and only if predicted_[v] >= 0. I am only using edges that are not on
    a previous path

    I am starting Dijkstra at all nodes which have a predicted potential. */

    while (true) {
        // I am continuing until all nodes are predicted

        NT largest_ratio = 0;
        node    start_node = nil,
                goal_node = nil;

        forall_nodes(v,G) {

            if ( predicted_p[v] < 0 )
                continue;

            // v is predicted and I am starting a shortest path computation

            node_pq<NT> PQ(G);

            PQ.insert(v,0);
            dist[v] = 0;

            node u;
            edge e;
            forall_nodes(u,G) pred[u] = nil;

            while ( !PQ.empty() ) {

                node u = PQ.del_min();

                if ( u != v && predicted_p[u] >= 0 ) {

                    NT pot_diff = abs(predicted_p[u] - predicted_p[v]);

                    if ( pot_diff/dist[u] > largest_ratio ) {
                        start_node = v;
                        goal_node = u;
                        largest_ratio = pot_diff/dist[u];
                    }

                    continue;
                }

                forall_inout_edges(e,u) {

                    if ( used[e] )
                        continue;

                    node z = G.opposite(e,u);
                    NT c = dist[u] + L[e];

                    if ( pred[z] == nil && z != s)
                        PQ.insert(z,c);         // first path to v

                    else if (c < dist[z])
                        PQ.decrease_p(z,c);     // better path

                    else
                        continue;

                    dist[z] = c;
                    pred[z] = e;
                }
            }
        }

        /** start_node and goal_node are the predicted nodes with the largest potdiff/length ratio;
        I am repeating the shortest path computation, in order to determine the path and its length*/

        if ( start_node == nil ) {
            return m_SP;
            break;
        }

        {
            v = start_node;

            node_pq<NT> PQ(G);



            PQ.insert(v,0);
            dist[v] = 0;

            node u;
            edge e;
            forall_nodes(u,G) pred[u] = nil;

            while ( !PQ.empty() ) {

                node u = PQ.del_min();
                if ( u == goal_node )
                    break;

                forall_inout_edges(e,u) {

                    if ( used[e] )
                        continue;

                    node z = G.opposite(e,u);
                    NT c = dist[u] + L[e];

                    if ( pred[z] == nil && z != s)
                        PQ.insert(z,c);                         // first path to v

                    else if (c < dist[z])
                        PQ.decrease_p(z,c);   // better path

                    else
                        continue;

                    dist[z] = c;
                    pred[z] = e;
                }
            }
        }

        NT Lsum = dist[goal_node];

        /** I am now walking back on the path and predict potentials */

        NT pot_diff = predicted_p[start_node] - predicted_p[goal_node]; // may be positive or negative

        v = goal_node; //std::cout <<"\nnew path of length " << Lsum << " from " << G[goal_node];
        while (true) {
            edge e = pred[v];
            predicted_rate[e] = predicted_rate[G.reverse(e)] = abs(pot_diff)/Lsum - 1;
            used[e] = used[G.reverse(e)] = true;
            node w = G.opposite(e,v); //std::cout << " via " << G[w];
            if ( w == start_node ) break;
            predicted_p[w] = predicted_p[v] + pot_diff * L[e]/Lsum;
            v = w;
        }

    }
//std::cout << "\n\n";
}


bool is_on_path(const GRAPH<int,double>& G, node w, const list<edge>& path) {
    edge e;
    forall(e,path) {
        if ( w == G.source(e) || w == G.target(e) )
            return true;
    }
    return false;
}

void construct_all_simple_paths(const GRAPH<int,double>& G, const edge_array<NT>& L, node t, node end_point,
                                NT current_length,
                                const list<edge>& path, list<simple_path>& all_simple_paths) {
    if ( end_point == t ) all_simple_paths.append(simple_path(current_length,t,path));
    edge e;
    forall_inout_edges(e,end_point) {
        node w = G.opposite(e,end_point);
        if ( is_on_path(G,w,path) ) continue;
        list<edge> new_path = path;
        new_path.append(e);
        construct_all_simple_paths(G,L,t,w,current_length + L[e],new_path,all_simple_paths);
    }
}

/*
//I predict the potentials
//-- If t is not reachable from s I return 0, otherwise, I return the number of edges in the shortest
//path from s to t
//-- first I compute shortest path distances from s and t, respectively
//-- nodes that are not reachable
int predict_potentials_third_try(const GRAPH<int,double>& G, int n, int m, node s, node t, const edge_array<NT>& L,
                                 edge_array<bool>& e_in_G, node_array<bool>& v_in_G,
                                 const edge_array<NT>& D, const edge_array<NT>& Q,
                                 const node_array<NT>& p, node_array<NT>& predicted_p,
                                 edge_array<bool>& in_SP, NT& D_st) {
    node v;
    edge e;
    forall_nodes(v,G) {
        v_in_G[v] = true;
        predicted_p[v] = -1;
    }
    forall_edges(e,G) {
        e_in_G[e] = true;
        in_SP[e] = false;
    }

    node_array<NT> dist(G);
    node_array<edge>   pred(G,nil);
    edge_array<bool>   used(G,false);

    forall_edges(e,G) assert(used[e] == false);


    D_st = DIJKSTRA_T(G, s, t, e_in_G, L, dist, pred);

    if (pred[t] == nil) return 0;

//    I am checking that the shortest path from s to t is unique and I am marking all edges
//    on the shortest path. The shortest path is unique if for every vertex on the shortest path there is only one incoming edge that defines the distance of the node. I am predict the potential for the nodes on the shortest
//    path. The potential is the length of the path to t.

    int m_SP = 0;

    predicted_p[t] = 0;
    predicted_p[s] = D_st;

    v = t;
    while ( v != s) {
        edge e = pred[v];
        m_SP++;
        in_SP[e] = true;
        //used[e] = true;
        node w = G.opposite(e,v);
        //predicted_p[w] = predicted_p[v] + L[e];
        edge f;
        forall_inout_edges(f,v) {
            if ( f == e) continue;
            node z = G.opposite(f,v);
            if ( dist[z] + L[f] == dist[v]) return 0;
        }
        v = w;
    }
// t is reachable from s and the shortest path from s to t is unique

    std::cout << "\n\nThe shortest s-t path consists of " << m_SP << " edges and has length " << D_st;

    simple_path path(0,s,list<edge>());

    list<simple_path> all_st_paths;
    construct_all_simple_paths(G,L,t,s,0,list<edge>(),all_st_paths);

    all_st_paths.sort(my_compare);

    forall(path,all_st_paths) {
        list<edge> cur_path = path.third();

        bool nothing_new = true;
        forall(e,cur_path)
        if ( predicted_p[G.source(e)] < 0 || predicted_p[G.target(e)] < 0 ) nothing_new = false;

        if (nothing_new) continue;


        std::cout << "\n\nsimple path of length " << path.first() << "   node sequence: ";
        v = s;
        forall(e,cur_path) {
            std::cout << " " << G[v];
            v = G.opposite(e,v);
        }
        std::cout << " " << G[t] << "\nwe are next decomposing this path.";

        while( true ) {
            list_item it;
            while (it = cur_path.first() ) {
                edge f = cur_path[it];
                if ( predicted_p[G.source(f)] >= 0 && predicted_p[G.target(f)] >= 0 )
                    cur_path.del_item(it);
                else break;
            }
            // the first item is mixed or the current path is empty

            if ( cur_path.empty() ) break;
            e = cur_path[it];
            node start_node = ( predicted_p[G.source(e)] >= 0 ? G.source(e) : G.target(e));

            stack<edge> S;
            S.push(e);
            cur_path.del_item(it);

            NT Lsum = L[e];
            while (it = cur_path.first() ) {
                edge e = cur_path[it];
                if ( predicted_p[G.source(e)] < 0 && predicted_p[G.target(e)] < 0 ) {
                    S.push(e);
                    Lsum += L[e];
                    cur_path.del_item(it);
                } else break;
            }
            // the first edge of cur_path is mixed.

            it = cur_path.first();
            e = cur_path[it];
            S.push(e);
            Lsum += L[e];
            cur_path.del_item(it);

            node goal_node =  ( predicted_p[G.source(e)] >= 0 ? G.source(e) : G.target(e));

            NT pot_diff = predicted_p[start_node] - predicted_p[goal_node]; // may be positive or negative

            v = goal_node;
            std::cout <<"\nnew path of length " << Lsum << " from " << G[goal_node];
            while (true) {
                edge e = S.pop();
                node w = G.opposite(e,v);
                std::cout << " via " << G[w];
                if ( w == start_node ) break;
                predicted_p[w] = predicted_p[v] + pot_diff * L[e]/Lsum;
                v = w;
            }
        }
    }
    return m_SP;
}

*/


//I am analysing the network
NT analysis_of_network(const graph& G, int n, node s, node t, const edge_array<NT>& L,
                       const edge_array<bool>& in_SP,
                       const edge_array<bool>& e_in_G, const node_array<bool>& v_in_G,
                       const edge_array<NT>& D, const edge_array<NT>& Q, const node_array<NT>& p, NT D_st)  {

    std::cout << "\n\n\n\npotential of s = " << p[s] << " should converge to shortest path distance which is " << D_st;
    NT Ds = 0;
    node v;
    edge e;
    forall_inout_edges(e,s) Ds += D[e];
    std::cout << "\n\nconductivitiy of edges out of s = " << Ds << " should converge to I0 which is " << I0;
    node_array<NT> dist(G);
    node_array<edge> pred(G);
    NT PH_D_st = DIJKSTRA_T(G, s, t, e_in_G, L, dist, pred);
    assert (PH_D_st <= D_st);
    if ( PH_D_st > D_st ) std::cout<<"\n\nD_st = " << D_st << " PH_D_st = " << PH_D_st;
    NT minD_inSP = 1000;
    NT maxD_outsideSP = 0;
    forall_edges(e,G) {
        if ( in_SP[e] && D[e] < minD_inSP ) minD_inSP = D[e];
        if ( !in_SP[e] && (G.reverse(e) == nil || !in_SP[G.reverse(e)]) && D[e] > maxD_outsideSP ) maxD_outsideSP = D[e];
    }
    std::cout << "\n\nminD_inSP = " << minD_inSP << " should converge to I0; \nmaxD_outsideSP = (should converge to zero) " << maxD_outsideSP;
    std::cout.flush();
    return PH_D_st;
    //read_int("\ninput an integer to continue");
}





#endif // PHYSARUM_H_INCLUDED


#endif // KURTSPHYSARUM_H_INCLUDED
