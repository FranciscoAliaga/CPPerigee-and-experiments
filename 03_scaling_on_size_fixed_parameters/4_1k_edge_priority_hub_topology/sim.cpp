
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

    // general network characteristics
    unsigned int V = 50000 ;
    unsigned int d = 10    ; // hypercube [0,1]^d dimension

    // TCP graph constraints
    unsigned int dout       = 25;
    unsigned int scored_out = 20;
    unsigned int din        = 125;

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
    HubTopology internetTopology{V,d,1000.,50.}; // TODO verify this is what we want

    // protocol construction
    EdgePriority scoringProtocol{"Edge Priority",scored_out,blocks_per_round};

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
    Network N(nodes,internetTopology,blockSize,dout,din,threads,113322);

    //////////////////////////////////////////////////////////////////////////////////////////////////// SIMULATION
    metrics::parameters pars {};
    pars.percentiles = {90.};
    pars.randomized = true;
    pars.randomSamples = static_cast<unsigned int>( 50.000 * 0.02 );

    N.writeMetricsCSV(
        "percentile_90_bound",
        metrics::kind::TimeToHashPowerPercentBound,
        pars);

    N.writeMetricsCSV(
        "percentile_90_start",
        metrics::kind::TimeToHashPowerPercent,
        pars);

    N.writeMetricsCSV(
        "edge_weights_start",
        metrics::kind::EdgePropagationWeight,
        pars);

    for(unsigned int i = 1; i <= 10; i++){
        N.run(100*10,100,0,1,1);
        // percentile
        std::string p_filename = "percentile_after_" + std::to_string(i*1000) + "_blocks";
        N.writeMetricsCSV(p_filename,
        metrics::kind::TimeToHashPowerPercent,
        pars);
        // edge weights
        std::string q_filename = "edge_weights_after_" + std::to_string(i*1000) + "_blocks";
        N.writeMetricsCSV(q_filename,
        metrics::kind::EdgePropagationWeight,
        pars);
    }

    return 0;
}
