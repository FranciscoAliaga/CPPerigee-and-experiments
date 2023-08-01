#ifndef GRAPH_H
#define GRAPH_H

#include "Constants.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <span>
#include <set>
#include <string>
#include <vector>


/* Graph: An adjacency-list based data structure for low-degree / sparse directed graphs
in this graph, a node can be reached from others via in_neighbors and can reach others via out_neighbors

This approach models the topology of a network of nodes in which every node has a preferential set of "out neighbors" 
that each node will change to optimize its interests, and can be reached through a set of "in neighbors" that may contact
each node without permission (other than to not exceed the maximum number of incoming nodes). In summary,
 - Each node is in charge of deciding whom to set as out_going node. 
 - Each node `u` who is `outgoing` from `v` has `v` as `incoming`

*/

class Graph{
    public:
        typedef unsigned int node;
        typedef node neighbor;
        typedef std::vector<neighbor> neighborList;

        /* adjacency list concentrating both incoming and outgoing edges */
        class adjacencyList { // not elegant to use lowercase 
            public: 
                neighborList in_neighbors {};
                neighborList out_neighbors {};

                adjacencyList(unsigned int dout, unsigned int din) : 
                    in_neighbors {neighborList(0)}, out_neighbors {neighborList(0)}
                {
                    in_neighbors.reserve(din);
                    out_neighbors.reserve(dout);
                }

                // size of all the neighbors
                size_t size() const {return in_neighbors.size()+out_neighbors.size();}
                // allows to fetch each of the neighbors (useful in a 'for loop') NOTE: both incoming and outgoing
                node operator[](const size_t& i) const {
                    return (i<out_neighbors.size()) ? out_neighbors[i] : in_neighbors[i-out_neighbors.size()];
                    }
        };
    protected:
        // name of the graph
        std::string name = "G";
        // number of nodes in the graph
        const unsigned int V = 1;
        // graph restrictions
        const unsigned int max_d_out = 1;
        const unsigned int max_d_in  = 1;

    private:
        // graph structure
        std::vector<adjacencyList> graphAdjacency;

    public:
        // constructor
        Graph( const std::string& name, unsigned int V, unsigned int max_d_out, unsigned int max_d_in) : 
            name {name},  V{V}, max_d_out {max_d_out}, max_d_in{max_d_in}, graphAdjacency { std::vector<adjacencyList>(V,adjacencyList{max_d_in,max_d_out}) }
        {};

    /* functions */

    private:
        /* these functions do the underground work of adding, searching, removing and mantaining the graph structure
        needed to implement the more abstract public functions */
        // code convention: always treat the edge (u,v) as the concept  "u->v" or "u is incoming seen from v, v is outgoing seen from u".

        /* edge addition */

        // private handler of add_edge
        // returns 1 if succesfull
        bool _add_edge(node u,node v){
            if(_can_add_edge(u,v)){
                graphAdjacency[v].in_neighbors.push_back(u);
                graphAdjacency[u].out_neighbors.push_back(v);
                return 1;
            }
            return 0;
        }
        /// private handler 
        bool _can_add_edge(node u,node v) const {
            // require: uv is not a loop
            bool cond1 = (u!= v); if(!cond1) return false;
            // require : u has enough space in out_going
            bool cond2 = (out_degree(u)  < max_d_out  ); if (!cond2) return false;
            // require : v has enough space in in_coming
            bool cond3 = (in_degree(v)  < max_d_in  ); if (!cond3) return false;
            // require: u is not already connected to v
            bool cond4 = !(is_edge(u,v)); if (!cond4) return false;

            // passed every condition.
            return true;
        }

        // handles edge removal. first checks if the edge is well-defined by searching the index of the target
        // vertices in their respective adjacency lists, and if they both exist, then executes the removal
        bool _remove_edge(node u,node v){
            // Can we remove this edge?
            // does this edge exists in both vertices?
            int iuv = search_in_edge(u,v);
            int ouv = search_out_edge(u,v);
            bool cond_existence1 = ( iuv >= 0 ); // guarantees iuv will be unsigned
            bool cond_existence2 = ( ouv >= 0 );
            if (!cond_existence1 || !cond_existence2){
                return 1; // exited because of condition
            }
            // perform the removal
            // these methods depend on the search behaviour of search_in/out_edge.
            // iuv, ouv
            remove_from_neighbors(graphAdjacency[v].in_neighbors , static_cast<unsigned int>(iuv));
            remove_from_neighbors(graphAdjacency[u].out_neighbors, static_cast<unsigned int>(ouv));
            return 0;
        }

        // removes a vertex of index `index` from a general `neighbors` adjacency list
        bool remove_from_neighbors(neighborList& neighbors,unsigned int index){
            // overwrite the last neighbor onto this one
            neighbors[index] = neighbors[neighbors.size()-1];
            // remove the last neighbor
            neighbors.pop_back();
            return 0;
        }

        /* edge search */


        // returns the index of the edge u->v in the in_neighbors adjacency of v.
        int search_in_edge(node u,node v) const {
            //rule(u->v)
            return neighbor_linear_search(graphAdjacency[v].in_neighbors,u);
        }

        // returns the index of the edge u->v in the out_neighbors adjacency of u.
        int search_out_edge(node u,node v) const {
            //rule(u->v)
            return neighbor_linear_search(graphAdjacency[u].out_neighbors,v);
        }

        // implements the search for a vertex index within an adjacency list
        int neighbor_linear_search(const neighborList& neighbors,node u) const {
            for(unsigned int i=0; i<neighbors.size();i++){
                if (neighbors[i]==u) return static_cast<int>(i);
            }
            return -1; // error: not found
        }
    
    public:
        unsigned int getV() const { return V; }
        unsigned int in_degree(node v) const { return graphAdjacency[v].in_neighbors.size();}
        unsigned int out_degree(node v) const { return graphAdjacency[v].out_neighbors.size();}
        unsigned int get_in_degree(node v) const { return in_degree(v);}
        unsigned int get_out_degree(node v) const {return out_degree(v);}

        // returns true if (u->v) is an edge inside the graph
        bool is_edge(node u,node v) const {return (search_out_edge(u,v)>=0);} // search_in_edge should be consistent with search_out_edge

        bool add_edge(node u,node v){
            return _add_edge(u,v);
        }

        bool remove_edge(node u,node v){
            return _remove_edge(u,v);
        }

        const adjacencyList& get_neighbors(node v) const {
            return (graphAdjacency[v]);
        }

        std::set<node> get_neighbors_set(node v) const { 
            std::set<node> neighbors_;
            neighbors_.insert(graphAdjacency[v].in_neighbors.begin(),graphAdjacency[v].in_neighbors.end());
            neighbors_.insert(graphAdjacency[v].out_neighbors.begin(),graphAdjacency[v].out_neighbors.end());
            return neighbors_;
        }

        const neighborList& get_out_neighbors(node v) const {
            return graphAdjacency[v].out_neighbors;
        }

        const neighborList& get_in_neighbors(node v) const {
            return graphAdjacency[v].in_neighbors;
        }

        bool can_add_edge(node u,node v) const {return _can_add_edge(u,v);}

        const neighborList out_neighbors(node v) const {return get_out_neighbors(v);}
        const neighborList in_neighbors(node v) const {return get_in_neighbors(v);}
        const adjacencyList& get_adjacencyList(node v) const {return graphAdjacency[v];}

        unsigned int get_max_d_in() const { return max_d_in; }
        unsigned int get_max_d_out() const { return max_d_out; }
};

#endif