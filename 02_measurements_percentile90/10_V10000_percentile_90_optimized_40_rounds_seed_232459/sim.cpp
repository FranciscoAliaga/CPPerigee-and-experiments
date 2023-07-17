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
    unsigned int V = 10000 ;
    unsigned int d = 10    ; // hypercube [0,1]^d dimension

    // TCP graph constraints
    unsigned int dout       = 15;
    unsigned int scored_out = 12;
    unsigned int din        = 25;

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
    FlatTopology internetTopology{V,d,1.,0.};

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
    Network N(nodes,internetTopology,blockSize,dout,din,threads,232459);

    //////////////////////////////////////////////////////////////////////////////////////////////////// SIMULATION
    metrics::parameters pars {};
    pars.percentiles = {90.};
    pars.randomized = false;

    N.run(40*100,1000,0,1,1);

    N.writeMetricsCSV(
        "percentile_90_bound",
        metrics::kind::TimeToHashPowerPercentBound,
        pars);

    N.writeMetricsCSV(
        "percentile_90",
        metrics::kind::TimeToHashPowerPercent,
        pars);

    return 0;
}
