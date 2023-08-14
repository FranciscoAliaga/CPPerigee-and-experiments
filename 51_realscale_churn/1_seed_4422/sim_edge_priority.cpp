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

#include "Churn.hpp"

using nodeRef = Node*;

int main(){
    // preliminary: hardware multithreading capabilities
    unsigned int threads; GET_HARDWARE_THREADS(&threads);

    unsigned int seed = 4422 ; // sampled at random from RANDOM.ORG

    // general network characteristics
    unsigned int V = 50000 ;
    unsigned int d = 20    ; // hypercube [0,1]^d dimension

    // TCP graph constraints
    unsigned int dout       =  20;
    unsigned int scored_out =  15;
    unsigned int din        =  80;

    // rounds
    unsigned int blocks_per_round  = 150;
    unsigned int visibility_radius = 2;
    unsigned int visiblity_radius_updates = 200;

    // node characteristics
    double blockGenerationHashPower   = 1.    ; // [block hashes/s] (relative)
    double blockValidationHashPower   = 1000. ; // [MB/s]  (absolute)
    double blockSize                  = 50.   ; // [MB]
    double nodeDownloadBandwidth      = 1.    ; // [MB/s]
    double nodeUploadBandwidth        = 1.    ; // [MB/s]

    // internet topology construction
    FlatTopology internetTopology{V,d,10000.,100.};

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

    // churn simulator construction
    long double churn_per_block = static_cast<long double> ( 200. / blocks_per_round ) ; // from bitnodes percentage...
    long double inactive_node_limit_percentage = 0.50;
    unsigned int churn_out_at_start = static_cast<unsigned int> (V * 0.15 ); // five percent

    ChurnSimulation Simulation{N,nodes,churn_per_block,inactive_node_limit_percentage,seed*seed};
    // churn some...
    Simulation.churn_out(churn_out_at_start);

    //////////////////////////////////////////////////////////////////////////////////////////////////// SIMULATION
    metrics::parameters pars {};
    pars.percentiles = {90.};
    pars.randomized = true;
    pars.randomSamples = static_cast<unsigned int>(10000 * 0.03);

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

    // 20000 = 5000 * 4
    for(unsigned int i = 1; i <= 4; i++){
        // NOTE run the simulation with churn.
        Simulation.run(5000,0,1);

        // save both metrics
        std::string filename_p = "edge_priority_percentile_after_" + std::to_string(i*5000) + "_blocks";
        N.writeMetricsCSV(filename_p,metrics::kind::TimeToHashPowerPercent,pars);

        std::string filename_w = "edge_priority_edge_weight_after_" + std::to_string(i*5000) + "_blocks";
        N.writeMetricsCSV(filename_w,metrics::kind::EdgePropagationWeight,pars);

        // save a copy of the graph
        std::string filename_g = "edge_priority_connections_after_" + std::to_string(i*5000) + "_blocks";
        N.saveGraphCSV(filename_g);

    }

    // log churns
    Simulation.write_history_log("edge_priority_churn_in_out");

    return 0;
}
