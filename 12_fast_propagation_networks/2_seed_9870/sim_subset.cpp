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

    unsigned int seed = 9870 ; // sampled at random from RANDOM.ORG

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

    // internet topology construction               // NOTE fast distribution network!
    FastFlatTopology internetTopology{V,d,0.10,0.1,2000.,500.};

    // protocol construction
    PerigeeSubset scoringProtocol{"Perigee Subset",scored_out,90.,blocks_per_round}; // NOTE perigee subset

    // connection construction
    ConnectionProtocol connectionProtocol("Basic random connection protocol",dout,din,visibility_radius,visiblity_radius_updates);

    // NOTE: here we add power to the first 10%
    unsigned int high_power_nodes = static_cast<unsigned int>(V*0.10);

    // node construction
    std::vector<nodeRef> nodes(V,nullptr);
    for(unsigned int i = 1; i <= V; i++){
        nodeRef node = new Node(
            scoringProtocol,
            connectionProtocol,
            true,
            blockGenerationHashPower * ((i < high_power_nodes)? 9. : 1./9. ), // NOTE first 10% get 9 times as much hash power as the rest.
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

    // NOTE: we dont measure metrics at startup nor save the connections. already done at edge priority simulation

    for(unsigned int i = 1; i <= 10; i++){
        // run the simulation
        N.run(100*10,100,1,0,1);

        // save both metrics
        std::string filename_p = "subset_percentile_after_" + std::to_string(i*1000) + "_blocks";
        N.writeMetricsCSV(filename_p,metrics::kind::TimeToHashPowerPercent,pars);

        std::string filename_w = "subset_edge_weight_after_" + std::to_string(i*1000) + "_blocks";
        N.writeMetricsCSV(filename_w,metrics::kind::EdgePropagationWeight,pars);

        // save a copy of the graph
        std::string filename_g = "subset_connections_after_" + std::to_string(i*1000) + "_blocks";
        N.saveGraphCSV(filename_g);
    }

    return 0;
}
