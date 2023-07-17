#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include "GraphAlgorithms.hpp"
#include "Parallelization.hpp"
#include "HashPower.hpp"
#include "InternetTopology.hpp"
#include "Network.hpp"
#include "Protocol.hpp"

using nodeRef = Node*;

int main(){
    // preliminary: hardware multithreading capabilities
    unsigned int threads; GET_HARDWARE_THREADS(&threads);

    unsigned int seed = 4165 ; // sampled at random from RANDOM.ORG
    const double link_validation_scale = 0.5;

    // general network characteristics
    unsigned int V = 10000 ;
    unsigned int d = 10    ; // hypercube [0,1]^d dimension

    // TCP graph constraints
    unsigned int dout       =  20;
    unsigned int scored_out =  16;
    unsigned int din        =  80;

    // rounds
    unsigned int blocks_per_round  = 100;
    unsigned int visibility_radius = 2;
    unsigned int visiblity_radius_updates = 200;

    // node characteristics
    double blockGenerationHashPower   = 1.    ; // [block hashes/s] (relative)
    double blockValidationHashPower   = 1000. ; // [MB/s]  (absolute)
    double blockSize                  = 50.   ; // [MB]
    double nodeDownloadBandwidth      = 1.    ; // [MB/s]
    double nodeUploadBandwidth        = 1.    ; // [MB/s]

    // internet topology construction
    FlatTopology internetTopology{V,d,2000.,500.*link_validation_scale};

    // protocol construction
    EdgePriority scoringProtocol{"Edge Priority",scored_out,blocks_per_round};  // NOTE edge priority

    // connection construction
    ConnectionProtocol connectionProtocol("Basic random connection protocol",dout,din,visibility_radius,visiblity_radius_updates);

    // node construction
    std::vector<nodeRef> nodes(V,nullptr);
    for(unsigned int i = 1; i <= V; i++){
        nodeRef node = new Node(
            scoringProtocol,
            connectionProtocol,
            true,
            blockGenerationHashPower,
            blockValidationHashPower,
            nodeDownloadBandwidth,
            nodeUploadBandwidth);
        nodes[i-1] = node;
    }

    // network construction
    Network N(nodes,internetTopology,blockSize,dout,din,threads,seed);

    //////////////////////////////////////////////////////////////////////////////////////////////////// SIMULATION
    metrics::parameters pars {};
    pars.percentiles = {90.};
    pars.randomized = true;
    pars.randomSamples = static_cast<unsigned int>(10000 * 0.05);

    // NOTE: we measure metrics at startup

    N.writeMetricsCSV(
        "percentile_90_bound",
        metrics::kind::TimeToHashPowerPercentBound,
        pars);

    N.writeMetricsCSV(
        "percentile_90_start",
        metrics::kind::TimeToHashPowerPercent,
        pars);

    N.writeMetricsCSV(
        "edge_weight_start",
        metrics::kind::EdgePropagationWeight,
        pars);

    // also save a copy of the graph at startup
    N.saveGraphCSV("randomized_connections_startup");

    for(unsigned int i = 1; i <= 10; i++){
        // run the simulation
        N.run(100*10,100,0,0,1);

        // save both metrics
        std::string filename_p = "edge_priority_percentile_after_" + std::to_string(i*1000) + "_blocks";
        N.writeMetricsCSV(filename_p,metrics::kind::TimeToHashPowerPercent,pars);

        std::string filename_w = "edge_priority_edge_weight_after_" + std::to_string(i*1000) + "_blocks";
        N.writeMetricsCSV(filename_w,metrics::kind::EdgePropagationWeight,pars);

        // save a copy of the graph
        std::string filename_g = "edge_priority_connections_after_" + std::to_string(i*1000) + "_blocks";
        N.saveGraphCSV(filename_g);

    }

    return 0;
}
