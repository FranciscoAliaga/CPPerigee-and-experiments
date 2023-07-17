#include "Constants.hpp"
#include "Graph.hpp"
#include "InternetTopology.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <map>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <vector>

#ifndef GRAPH_ALGORITHMS_H
#define GRAPH_ALGORITHMS_H

// used to store graph algorithm predecessors in a safe way
class pred {
    public: 
        unsigned int node = 0;
        bool valid        = false;
    
    // constructors
    pred(unsigned int n,bool val) : node{n}, valid{val} {}
    pred(unsigned int n) : pred(n,true) {}
    pred() : pred(0,false) {}

};

inline bool operator==(const unsigned int u, const pred& p){ if (!p.valid){return false;} return p.node == u; }
inline bool operator==(const pred& p, const unsigned int u){ return u==p;}
inline bool operator==(const pred& p, const pred& q){return (p.node==q.node) && (p.valid == q.valid) && p.valid;}
inline bool operator>(const pred& p, const pred& q){return (p.node>q.node) && (p.valid) && (q.valid);}

class GraphAlgorithms{
    public:

        /* To store vectors of distances and path predecessors to a source for single source shortest paths for a Graph
        source:= source node
        distances:= distances from node to source
        predecessors:= predecessors to form the shortest path between node and source*/
        struct distanceAndPaths {
            unsigned int source = 0; 
            std::vector<double> distance {};
            std::vector<pred> predecessor {};
            // NOTE: predecessor[v]=v iff v is source.
            // NOTE: invalid predecessors are chequed from v.valid
        };

    public:
        /* calculates directed shortest paths, implementing Dikjstra's algorithm for directed graphs.*/
        static distanceAndPaths DistanceAndPaths(const Graph& G, const InternetTopology& internet, unsigned int source){
            // what is to be returned
            distanceAndPaths data;

            // input
            unsigned int V = G.getV();
            unsigned int s = source;

            // fetched to simplify notation
            std::vector<double>& dist  = data.distance;
            std::vector<pred>&   predecessor = data.predecessor;
            data.source = source;

            pred invalid {0,false};

            // initializations
            dist.resize(V);
            predecessor.resize(V);
            std::fill(dist.begin(),dist.end(), std::numeric_limits<double>::infinity() ); 
            std::fill(predecessor.begin(),predecessor.end(),invalid);

            // define priority queue
            using par = std::pair<double,pred>;
            auto comparator  = [](const auto& lhs, const auto& rhs) -> bool {
                if ( lhs.first > rhs.first ) return true;
                //else if ( lhs.first == rhs.first ) return lhs.second > rhs.second;
                return false; };

            std::priority_queue<par,std::vector<par>,decltype(comparator)> Q(comparator);

            predecessor[s] = pred(s);
            dist[s] = 0.;
            Q.push(par{0.0,{s}});

            while(Q.size()>0){
                par v_pair = Q.top(); Q.pop();
                auto [d_v,v] = v_pair;

                unsigned int iv = v.node;

                if(d_v > dist[iv]) continue;

                std::set<unsigned int> v_neighbors = G.get_neighbors_set(iv);
                for(auto iu : v_neighbors){
                    auto u = pred(iu);
                    double lvu = internet.latency(iv,iu);
                    if(dist[iv] + lvu < dist[iu]){
                        dist[iu] = dist[iv] + lvu;
                        predecessor[iu] = v;
                        Q.push({dist[iu],u});
                    }
                }
            } 
            return data;
        }

        // Calculates distance and paths, while skipping nodes that are not active (given by activeNodes)
        static distanceAndPaths DistanceAndPaths(
            const Graph& G,
            const InternetTopology& internet,
            const std::vector<bool>& activeNodes, 
            unsigned int source){
            // what is to be returned
            distanceAndPaths data;

            // input
            unsigned int V = G.getV();
            unsigned int s = source;

            // fetched to simplify notation
            std::vector<double>& dist  = data.distance;
            std::vector<pred>&   predecessor = data.predecessor;
            data.source = source;

            pred invalid {0,false};

            // initializations
            dist.resize(V);
            predecessor.resize(V);
            std::fill(dist.begin(),dist.end(), std::numeric_limits<double>::infinity() ); 
            std::fill(predecessor.begin(),predecessor.end(),invalid);

            // define priority queue
            using par = std::pair<double,pred>;
            auto comparator  = [](const auto& lhs, const auto& rhs) -> bool {
                if ( lhs.first > rhs.first ) return true;
                //else if ( lhs.first == rhs.first ) return lhs.second > rhs.second;
                return false; };

            std::priority_queue<par,std::vector<par>,decltype(comparator)> Q(comparator);

            predecessor[s] = pred(s);
            dist[s] = 0.;
            Q.push(par{0.0,{s}});

            while(Q.size()>0){
                par v_pair = Q.top(); Q.pop();
                auto [d_v,v] = v_pair;

                unsigned int iv = v.node;

                if(d_v > dist[iv]) continue;

                std::set<unsigned int> v_neighbors = G.get_neighbors_set(iv);
                for(auto iu : v_neighbors){
                    if (!activeNodes[iu]) continue; // NOTE < inactive node skipping is implemented here
                    auto u = pred{iu};
                    double lvu = internet.latency(iv,iu);
                    if(dist[iv] + lvu < dist[iu]){
                        dist[iu] = dist[iv] + lvu;
                        predecessor[iu] = v;
                        Q.push({dist[iu],u});
                    }
                }
            } 
            return data;
        }
    
    // gets each node under radious k from source (using jump distance)
    // V(source, k) = {v in V(G): D_G(source,v) <= k }\{source}
    static std::vector<unsigned int> get_neighborhood(const Graph& G,unsigned int source,unsigned int k){

        using par = std::pair<unsigned int, unsigned int >;
        std::priority_queue<par,std::vector<par>,std::greater<par>> Q;

        // here we use a modified version of dikjstra to calculate d(v,u)
        // that skips processing vertices u that have no chance of having d(v,u)<=k

        std::map<unsigned int, unsigned int> distance;
        distance[source] = 0;
        Q.push(par{0,source});

        while(Q.size()>0){
            par v_par = Q.top(); Q.pop();
            auto [d_v, v] = v_par;

            if (d_v >= k){ // skip any vertex at the border
                continue;  // since, by the priority queue invariant,
            }              // their neighbors must be at distance at least k+1 from v

            std::set<unsigned int> v_neighbors = G.get_neighbors_set(v);
            for(auto u: v_neighbors){
                // if not seen, then distance is d_v+1
                // note: not seen <==> not found on distance map/dictionary <==> expression inside if
                if(distance.find(u)==distance.end()){
                    distance[u] = d_v + 1;
                    Q.push(par{d_v+1,u});
                }
                else
                {
                    // seen, but might have higher distance
                    if (distance[u] > d_v + 1){
                        distance[u] = d_v + 1;
                        Q.push(par{d_v+1,u});
                    }
                }
            }
        }
        // at algorithm end, distance[] map is only populated with keys having
        // jump distance to v lesser than or equal to k.

        // report those in a vector. 
        std::vector<unsigned int> ans(distance.size()-1); // removing source from report
        unsigned int index = 0;
        for(auto it = distance.begin();it!=distance.end(); it++){
            auto n = it->first; if(n==source) continue;
            ans[index] = n; // first = dicrionary key = node
            index+=1;
        }
        return ans;
    }

};

#endif

