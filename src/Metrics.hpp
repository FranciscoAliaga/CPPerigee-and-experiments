#include "Constants.hpp"
#include "Graph.hpp"
#include "GraphAlgorithms.hpp"
#include "HashPower.hpp"
#include "InternetTopology.hpp"

#include <random>
#include <algorithm>
#include <limits>
#include <vector>

#ifndef METRICS_H
#define METRICS_H

class propagationLatencyEdges {
    public:
        // used to report calculation result of a single node propagating message
        struct answer_propagation {
            unsigned int     source{};
            double           weight{}; // the full sum of every edge latency involved in the message propagation
            bool          succesful{}; // true if the message was succesfully broadcasted to the entierety of the network
            unsigned int peers_lost{}; // if the flag succesful is false, reports the number of peers missing the message.
        };

        static answer_propagation calculate_single (
            const Graph& G,
            const std::vector<bool>& activeNodes,
            const InternetTopology& internet,
            unsigned int source
            ){
                GraphAlgorithms::distanceAndPaths dp = GraphAlgorithms::DistanceAndPaths(G,internet,activeNodes,source);
                auto V = dp.predecessor.size();
                double weight = 0.;
                bool node_was_invalid_flag = false;
                unsigned int peers_lost    = 0;
                for(unsigned int iv = 0; iv < V; iv++){
                    auto pred = dp.predecessor[iv];
                    // if a node didn't get the message, write flag and continue (skip adding weight)
                    if (!pred.valid){
                        node_was_invalid_flag = true;
                        peers_lost++;
                        continue; }
                    weight += internet.latency(iv,pred.node);
                }
            return {source,weight,!node_was_invalid_flag,peers_lost};
        }

};

class timeToHashPowerPercentile {
    public:
        static std::vector<double> calculate(
                const Graph& G,
                const std::vector<double>& activeHashPower,
                const std::vector<bool>& activeNodes,
                const InternetTopology& internet,
                const std::vector<double> percentiles,
                unsigned int v,
                bool bound = false) {

            // first, depending on if calculating actual distances or bounds, 
            // we fill and sort a temporary array.
            unsigned int V = G.getV();
            std::vector<std::pair<double,unsigned int>> temp(V); 
            if(!bound){
                put_distance(G,activeNodes,internet,v,temp);
            }else{
                put_latencies(internet,v,temp);
            }
            std::sort(temp.begin(),temp.end());

            // prepare answer
            std::vector<double> answer(percentiles.size(),std::numeric_limits<double>::infinity());
            // search for the first node index exceding target hash , annotate it's distance/latency, and search for the next target...
            double accumulated_hash = 0.0;
            double target = percentiles[0];
            unsigned int    ans_index = 0;
            // iterate
            for(unsigned int iu=0;iu<V;iu++){
                unsigned int node = temp[iu].second;
                accumulated_hash += activeHashPower[node];
                if (accumulated_hash >= target){ // answer for target percentile reached
                    answer[ans_index] = temp[iu].first; // distance to u
                    ans_index+=1;
                    if (ans_index>=percentiles.size()) break;
                    target = percentiles[ans_index];
                }
                if (temp[iu].first==std::numeric_limits<double>::infinity()) break; // case we find that distances start to be infinity...
            }
            // else, target hash was never found... (unlikely, but we cover the border case with infinity)
            return answer;
        }

        static void put_distance(
            const Graph& G,
            const std::vector<bool>& activeNodes,
            const InternetTopology& internet,
            unsigned int v, std::vector<std::pair<double,unsigned int>>& temp){
                GraphAlgorithms::distanceAndPaths dp = GraphAlgorithms::DistanceAndPaths(G,internet,activeNodes,v);
                unsigned int V = temp.size();
                    for(unsigned int iu = 0; iu<V; iu++){
                        temp[iu]= std::move(std::pair<double,unsigned int>(dp.distance[iu],iu));
                    }
            }

        static void put_latencies(
            const InternetTopology& internet,
            unsigned int v, std::vector<std::pair<double,unsigned int>>& temp){
            unsigned int V = temp.size();
            for(unsigned int iu = 0; iu<V; iu++){
                temp[iu]=std::pair<double,unsigned int>(internet.latency(v,iu), iu);
            }
        }

        static double calculate(const Graph& G, const std::vector<double>& activeHashPower, const std::vector<bool>& activeNodes, const InternetTopology& internet, double percentile,unsigned int v){
            std::vector<double> percentiles(1,percentile);
            return calculate(G,activeHashPower,activeNodes,internet, percentiles, v)[0];
        }

};

#endif
