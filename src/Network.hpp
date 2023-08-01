#ifndef NETWORK_H
#define NETWORK_H

#include "Constants.hpp"
#include "Graph.hpp"
#include "GraphAlgorithms.hpp"
#include "HashPower.hpp"
#include "InternetTopology.hpp"
#include "Protocol.hpp"
#include "Metrics.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <thread>  // parallelization
#include <unordered_map>
#include <vector>


void Network::run(unsigned int total_blocks, unsigned int blocks_per_loop, int full_message, bool full_visibility, bool update_visibility_flag){
    auto rounds = total_blocks / blocks_per_loop ;
    auto remainder_blocks = total_blocks % blocks_per_loop;

    // simulate rounds
    if(full_message==1)std::cout << "blocks broadcasted :" << std::endl;
    if(full_message==0) std::cout << "rounds of " << blocks_per_loop << " blocks: " << std::endl;

    for(unsigned int r = 0; r < rounds ; r++){
        if(full_message==0){std::cout << "\r round: " << (r+1) << " / " << rounds;}
        for(unsigned int ib = 0 ; ib < blocks_per_loop; ib ++){
            if(full_message==1){
                std::cout <<
                "\r" << blocks_per_loop << " * " << r <<
                " + " << ib <<
                " = " << blocks_per_loop*r + ib <<
                " / " << total_blocks << "      ";
            }
            /***********************************************************/
            /******/Loop(full_visibility,update_visibility_flag);/******/
            /***********************************************************/
        }
    } // if there are blocks remaining, finish them
    if (remainder_blocks>0){
        auto progress = total_blocks - remainder_blocks;
        for(unsigned int ib = 0 ; ib < remainder_blocks; ib ++){
            if(full_message==1) std::cout << "\r" << progress << " / " << total_blocks << "                               ";
            Loop(full_visibility);
        }
    }
    std::cout << std::endl << "Simulation of " << total_blocks << " blocks : DONE. " << std::endl;
    return;
}

// MainLoop
void Network::Loop(bool full_visibility,bool update_visibility_flag){
    // MAIN LOOP ///////////////////////////////////////////////////////////////////////////////////

    // event queue and visibility procedures  ////////////////////////////////
    // process the eventQueue
    processEventQueue(); 
    // checks and performs if nodes should update their visibility
    if (update_visibility_flag){ // NOTE: severely affects performance
        LoopVisibilityProcessing();
    }

    //// MAIN COMPONENTS //////////////
    // pool (new connections) processing ///////////////////////////////////////////////////////////
    //     processes the nodes requests for connection among peers
    ProcessPool(full_visibility);
    // generate block +  propagation ///////////////////////////////////////////////////////////////
    //    generates a block and broadcasts it through the network,
    //    updating each active node observation (only for nodes that see the block)
    LoopBroadcastBlock();
    // node processing /////////////////////////////////////////////////////////////////////////////
    //    for each node, performs the internal protocol decision process,
    //    for example, deciding if it's time to score and remove neighbors based on observations
    LoopNodeProcessing();
    // pool (new connections) processing ///////////////////////////////////////////////////////////
    //     processes the nodes requests for connection among peers
    ProcessPool(full_visibility);

    ////////////////////////////////////////////////////////////////////////////////////////////////
}

void Network::LoopBroadcastBlock(){
    //    first, generates a block, then calculates routes
    unsigned int block_origin = hashPower.generateBlock();
    GraphAlgorithms::distanceAndPaths blockPropagationData = 
    GraphAlgorithms::DistanceAndPaths(G,internet,activeNodes,block_origin);

    //    second, uses the information to update node observations.
    //    It's parallelized via threads.

    auto fun = [this,blockPropagationData](thread_work_range range){
        ProcessRangeObservations(range,blockPropagationData);
    };

    std::vector<std::thread> blockPropagationThreads {};
    for(thread_work_range range : thread_ranges){
        std::thread t {fun , range};
        blockPropagationThreads.push_back(std::move(t));
    } for(auto& t : blockPropagationThreads) t.join(); 
}

void Network::LoopNodeProcessing(){
    // first  : build the array with which the nodes will drop the neighbors. (parallel)
    // second : perform node disconnection (interacts with graph G, so it must be sequential)

    //      do first
    unsigned int P = thread_ranges.size();
    std::vector<std::vector<Node::nodeRemoval>> removalList(P,std::vector<Node::nodeRemoval>(0));
    std::vector<std::thread> nodeSelectionThreads {};
    std::vector<rangeSelect> selectionTasks{};
    auto process = [this](rangeSelect arg){ ProcessRangeSelectNeighbors(arg); };
    for(unsigned int p = 0; p < P; p++){
        selectionTasks.push_back({thread_ranges[p],&(removalList[p])});
    }
    for( auto& task : selectionTasks){
        std::thread t {process, task};
        nodeSelectionThreads.push_back(std::move(t));
    } for(auto& t : nodeSelectionThreads) t.join(); 

    //      do second
    // iterate over removalList making the removals
    for(const auto& list : removalList){
        for(const auto& request : list){
            disconnect(request.callerId,request.targetId);
            // NOTE this must be done sequentially
            //   or through a mutex, but it's not extremely computationally heavy.
            // O(|V|*(d_out-d_scored))
        }
    }
    // updates the each node connection requests
    for(unsigned int iv = 0; iv<V; iv++) {
        if (!activeNodes[iv]) continue;
        auto n = nodes[iv];
        n->processConnectionRequest();
        // NOTE this must be done sequentially,
        // it's not extremely computationally heavy.
    }
}


void Node::updateObservations(const GraphAlgorithms::distanceAndPaths& newBlockData){
    // check if the block was broadcasted to this node.
    // a node has unvalid predecessor iff the block was not delivered to it.
    bool delivered = (newBlockData.predecessor[id]).valid;
    if(!delivered){ return; }

    // get outgoing neighbors for scoring mechanismns
    std::vector<nodeId> neighbors = network->getOutNeighbors(id);

    // update the neighbormap in case a new outgoing neighbor has appeared
    for(unsigned int i = 0; i < neighbors.size();i++){
        nodeId neighbor = neighbors[i];
        if (neighborMapForward.find(neighbor) == neighborMapForward.end()) {
            // update the neighbormap to add this new neighbor
            unsigned int new_index = neighborMapForward.size();  // remember the neighbormap is 0,1...
            neighborMapForward[neighbor] = new_index;
            neighborMapBackward[new_index] = neighbor;
            if (new_index >= maxNeighbors) continue; // NOTE what? this means the number of outgoing neighbors is complete
        }
    }
    // update number of neighbors seen during this round
    observations.neighbors = neighborMapForward.size();

    ////////// now separate in cases depending on the scoring protocol flags ////////////////////////////

    /// case: protocol scores using the predecessor for the blocks
    // add newly seen predecessor
    if (scoringProtocol.get_predecessor_flag){
        // we update the predecessor for the current block.
        // by correctness of the dikjstra algorithm, the predecessor MUST be an incoming neighbor,
        // however, it may not be an outgoing neighbor!
        unsigned int predecessor = (newBlockData.predecessor[id]).node; // by requirements of the function, predecessor is valid
        auto neighbor = neighborMapForward.find(predecessor); // search among outgoing neighbors

        if (neighbor!=neighborMapForward.end()){ // case that IT IS an outgoing neighbor.
            observations.predecessors[currentScoringBlock] = static_cast<int>(neighbor->second);
        }                                 // case it is not, then it is not scored.
        else{observations.predecessors[currentScoringBlock] = -1;}
    }

    /// case: protocol scores using the first time the block has been seen
    // add newly seen times 
    if (scoringProtocol.get_first_time_flag){
        // we update the first time a block is seen by the network distance between the block and the vertex.
        // check GraphAlgorithms (implemented using dikjstra)
        // since predecesor is valid, this distance must be finite.
        observations.firstTime[currentScoringBlock] = newBlockData.distance[id];
    }

    // cases: if we need latency between neighbors or the times at which the neighbors broadcast blocks
    if (scoringProtocol.get_latency_flag || scoringProtocol.get_times_flag){
        for(unsigned int i = 0; i<neighbors.size();i++){
            nodeId neighbor = neighbors[i];
            unsigned int neighborIndex  = neighborMapForward[neighbor];
            // latency from u to v := l(u,v) (internet latency) if u is connected to v, else, infinity.

            // house cleaning: if a new neighbor has been seen, then resize.
            if (neighborIndex >= observations.latency.size()){
                observations.latency.resize(neighborIndex+1);
                observations.latency[neighborIndex] = (
                    std::numeric_limits<double>::infinity());
            }

            // add newly seen latencies.
            if (network->isConnected(neighbor,id)){ // query u<->v 
                observations.latency[neighborIndex] = network->latency(neighbor,id);
            }else{
                observations.latency[neighborIndex] = std::numeric_limits<double>::infinity(); // case not connected...
            }

    // case: if protocol uses times the neighbors delivered a block
            // add newly seen times 
            if (scoringProtocol.get_times_flag){
                // add the times for the neighbor

                // house-cleaning. make space for newly seen neighbor
                if (neighborIndex >= observations.times.size()){
                    observations.times.resize(neighborIndex+1);
                    observations.times[neighborIndex] = (
                        std::vector<double>(
                            scoringProtocol.totalBlocks,
                            std::numeric_limits<double>::infinity())
                    );
                }

                // add times
                // T(b,v) = latency(u,v) + T(b,u) .
                // remember this can be infinite if either latency(u,v) is infinite or T(b,u) is infinite.
                // by protocol, T(b,u) should be finite if u is among this node neighbors.

                bool notRelayed = (id == (newBlockData.predecessor[neighbor]).node);
                
                if (!notRelayed){ // (if relayed)
                observations.times[neighborIndex][currentScoringBlock] = (
                    newBlockData.distance[neighbor] +   // time to neighbor
                    observations.latency[neighborIndex]);  // time from neighbor to self (might be infinite)
                }// else, this is infinity
            }
        }
    }
    currentScoringBlock++;
    currentVisibilityBlock++;
}

std::vector<Node::nodeRemoval> Node::processSelectNeighbors(){
    std::vector<Node::nodeRemoval> ans {};
    if(currentScoringBlock>=scoringProtocol.totalBlocks){
        ans = performNeighborSelection();
        if (ans.size()>0) requestNodeFlag = true;
        restartRound();
    }
    return ans;
}

void Node::processConnectionRequest(){
    if (!requestNodeFlag) return;

    int neighbors_needed =
        static_cast<int>(connectionProtocol.maxOutgoing)
      - static_cast<int>(network->getOutDegree(id)) ;

    neighbors_needed = neighbors_needed >= 0 ? neighbors_needed : 0;
    network->nodePoolRequest(id,static_cast<unsigned int>(neighbors_needed));

    requestNodeFlag = false;
}


void Node::restartRound(){
    // restart the first times
    if (scoringProtocol.get_first_time_flag){
        unsigned int n_blocks = scoringProtocol.totalBlocks;
        if (observations.firstTime.size()<n_blocks){
            observations.firstTime.resize(n_blocks);
        }
        std::fill(observations.firstTime.begin(),observations.firstTime.end(),std::numeric_limits<double>::infinity());
    }
    // restart the times
    if (scoringProtocol.get_times_flag){
        unsigned int n_blocks = scoringProtocol.totalBlocks;
        for(unsigned int iu = 0; iu<observations.times.size();iu++){
            if (observations.times[iu].size()<n_blocks){
                observations.times[iu].resize(n_blocks);
            }
            std::fill(observations.times[iu].begin(),observations.times[iu].end(),std::numeric_limits<double>::infinity());
        }
    }
    // restart the predecessors
    if (scoringProtocol.get_predecessor_flag){
        unsigned int n_blocks = scoringProtocol.totalBlocks;
        if (observations.predecessors.size()<n_blocks){
            observations.predecessors.resize(n_blocks);
        }
        observations.predecessors = std::vector<int>(n_blocks,-1);
    }
    // restart the latencies
    if (scoringProtocol.get_latency_flag){
        observations.latency = std::vector<double>(maxNeighbors,std::numeric_limits<double>::infinity());
    }

    // restart the neighbor map and block
    neighborMapForward.clear();
    neighborMapBackward.clear();
    currentScoringBlock=0;
}

std::vector<Node::nodeRemoval> Node::performNeighborSelection(){
    std::vector<Node::nodeRemoval> ans{};
    // this list is false on the neighbor index that will be dropped
    std::vector<bool> selectNeighbors = scoringProtocol.selectNeighbors(observations);
    // we run through each neighbor seen and decide if it is time to drop it
    for(auto it = neighborMapForward.begin();it!=neighborMapForward.end();it++){
        unsigned int index = it->second;
        if(!selectNeighbors[index]){
            nodeId dropped = it->first;
            ans.push_back({id,dropped});
        }
    } 
    return ans;
}

void Network::ProcessRangeObservations(const thread_work_range range, const GraphAlgorithms::distanceAndPaths& data){
    for(unsigned int iv = range.first; iv < range.second ; iv++){
        if (activeNodes[iv]){ // if is active, update observation
            nodes[iv]->updateObservations(data);
        }
    }
}

void Network::ProcessRangeSelectNeighbors(Network::rangeSelect arg){
    auto range = arg.r;
    auto ans = arg.output;
    for(unsigned int iv = range.first; iv < range.second ; iv++){
        auto v = nodes[iv];
        if(v->getActive()){ // if is active, perform
            auto nodeAns = v->processSelectNeighbors();
            ans->insert(ans->end(),nodeAns.begin(),nodeAns.end()); // append
        }
    }
    return;
}

void Network::LoopVisibilityProcessing(){
    updateVisibility();
}

void Network::updateVisibility(){
    auto task = [this](thread_work_range range){ProcessVisibilityRange(range);}; // uses helper function (below)
    std::vector<std::thread> visibilityWorkThreads {};
    for(thread_work_range r : thread_ranges){
        std::thread t {task,r};
        visibilityWorkThreads.push_back(std::move(t));
    } for(auto& t : visibilityWorkThreads) t.join();
}

void Network::ProcessVisibilityRange(thread_work_range range){
    auto [start,end] = range;
    for(unsigned int iv = start; iv < end ; iv++){
        if (!activeNodes[iv]) continue;
            auto v = nodes[iv]; 
            if(v->firstVisibilityUpdateFlag){ // if it is necessary to perform first visibility update
                auto k = v->connectionProtocol.visibleRadius;
                auto neighborhood = GraphAlgorithms::get_neighborhood(G,iv,k);
                v->setVisibleNodes(neighborhood);
                v->firstVisibilityUpdateFlag = false;
            } else
            if(v->currentVisibilityBlock >= v->connectionProtocol.visibilityUpdateBlocks){ // if sufficient blocks have been spotted
                auto k = v->connectionProtocol.visibleRadius;
                auto neighborhood = GraphAlgorithms::get_neighborhood(G,iv,k);
                v->setVisibleNodes(neighborhood);
                v->currentVisibilityBlock=0;
            }
        }
    return;
}

void Network::BootstrapNetwork(){
    std::cout << "Bootstrapping network...";
    // generate requests for nodes to connect at random
    for(unsigned int i = 0; i < V; i++){
        if (!activeNodes[i]) continue; // skip if node inactive
        auto n = nodes[i];
        auto request_amount = n->connectionProtocol.maxOutgoing;
        nodePool.updateRequest(n->id,request_amount);
    }
    // generate the connections
    ProcessPool(true);
    // update the visibility of each node
    updateVisibility();
    std::cout << " DONE. " << std::endl;
}

//// Node Pool Implementation ///////////////////////////////////////////////////////////////////////////

void Network::ProcessPool(bool bootstrap){
    unsigned int startRequests = nodePool.getUnsolvedRequests();
    if(startRequests==0) return;
    int tries = static_cast<int>(startRequests * nodePool.repetition + DEFAULT_NODEPOOL_MINIMUM_TRIES);
    while((tries-->0) && nodePool.requests.size() > 0 ){ 
        // get random node request /////////////////////////////////////
        auto request = nodePool.getRandomRequest(getRandomNumber());
        nodeId u = request.node;
        nodeId v = 0;
        if(bootstrap){
            v = getRandomNode();
        }else{
            v = nodes[u]->getRandomVisibleNode(getRandomNumber());
        }
        // require: requester and other must be distinct
        if (u==v){tries--;continue;}
        // require: target is active
        if(!activeNodes[v]){tries--;continue;}
        // require: request is sane 
        if(request.amount==0){tries--;continue;}
        // require: the requester hasn't connected to the subscriber
        if(isConnected(u,v)){tries--;continue;}
        
        // request is acceptable, perform connection //////////////////
        bool success = connect(u,v);
        //// if successfull, requester should be updated
        if (success){
            nodePool.updateRequest(u,request.amount-1);
        }
    }
}

/* Network primitives */
bool Network::connect(nodeId u,nodeId v){
    auto nu = nodes[u];
    auto nv = nodes[v];
    // require: both u and v agree that connection is okay with their protocols
    bool u_agrees_out = (G.get_out_degree(u) < nu->connectionProtocol.maxOutgoing);
    if (!u_agrees_out){return false;}
    bool v_agrees_in  = (G.get_in_degree(v) < nv->connectionProtocol.maxIncoming);
    if (!v_agrees_in){return false;}
    // attempt to connect them via Graph
    // implicit require: G can recieve edge
    if(G.add_edge(u,v)) return true;
    // else,
    return false;
}

void Network::disconnect(nodeId u, nodeId v){
    G.remove_edge(u,v);
}

bool Network::isConnected(nodeId u,nodeId v) const {
    const Graph::adjacencyList neighbors = G.get_neighbors(u);
    for(unsigned int i = 0; i<neighbors.size(); i++){
        if(neighbors[i]==v) return true;
    } // else
    return false;
}

bool Network::isConnectedDirected(nodeId u, nodeId v) const {
    std::vector<Graph::node> neighbors = G.get_out_neighbors(u);
    for(unsigned int i=0; i < neighbors.size(); i++){
        if(neighbors[i]==v) return true;
    } // else
    return false;
}
std::vector<nodeId> Network::getOutNeighbors(nodeId askingNode) const {
    return G.get_out_neighbors(askingNode);
}

std::vector<nodeId> Network::getInNeighbors(nodeId askingNode) const  {
    return G.get_in_neighbors(askingNode);
}

unsigned int Network::getOutDegree( nodeId caller) const {
    return G.get_out_degree(caller);
}

/* Node pool primitives */

NodePool::nodeRequest NodePool::getRandomRequest(unsigned int randomNumber){
    return requests[randomNumber%requests.size()];
}

// network manager
void Network::nodePoolRequest(nodeId caller,unsigned int request_amount){
    nodePool.updateRequest(caller, request_amount);
};

void NodePool::updateRequest(nodeId client, unsigned int request_amount){
    // 1. check if client is already present by performing hash map lookup
    auto mapElement = requestMap.find(client);
    if(mapElement != requestMap.end()){// the element is in the map
        // update the amount requested on the unsolvedRequest counter
        unsigned int index          = mapElement->second;

        unsigned old_request_amount = requests[index].amount;
        unsolvedRequests += (request_amount - old_request_amount);
        // update or remove depending on the amount requested.
        if ( request_amount == 0 ) {
            // remove the request.
            removeRequest(client); }
        else {
            // request > 0, update the request data.
            unsigned int index = requestMap[client]; 
            requests[index].amount = request_amount; }
        return;
    }
    // else, the client is not present
    // create the request and put it at the end of the request vector
    unsigned int new_index = requests.size();
    nodeRequest request = nodeRequest{client,request_amount};
    requests.push_back(request);
    // ad this index to the request map
    requestMap[client] = new_index;
    // finish by adding up to the unsolved requests counter
    unsolvedRequests += request_amount;
    return;
}

void NodePool::removeRequest(nodeId client){
    // if the client index is not found , then return
    if(requestMap.find(client) == requestMap.end()){return;}

    // remove client from the map
    unsigned int index = requestMap.at(client);
    //      separate by cases
    //          index is at the end (easy case)
    if(index == requests.size()-1){
        requestMap.erase(client);
        requests.pop_back(); 
        return;
    }
    //          index is inside (hard case, needs swap)
    requestMap.erase(client); // erase from map
    unsigned int last = requests.size()-1; // last index used by requests
    std::swap(requests[index],requests[last]); // swap last to the interior
    requests.pop_back(); // get rid of the request at the back
    // update the requestMap
    requestMap.at(requests.at(index).node) = index; 
    return;
}

// Writing functions 

// This function saves the network graph onto GRAPH_DIRECTORY/folder_name
void Network::saveGraphCSV(std::string folder_name){
    std::cout << "Saving graph...";
    ////
    //    this function performs several steps:
    //  1. makes sure the directories are well setup
    //  2. creates the writing streams
    //  3. writes
    //

    // 1. make sure directories are well setup /////////////////////////////////////////////////////

    // if graph directory does not exist, create it
    std::filesystem::path graph_directory_path(GRAPH_DIRECTORY);
    if(!std::filesystem::exists(graph_directory_path)){
        std::filesystem::create_directories(graph_directory_path);
    }

    // if specific graph folder does not exist, create it
    std::filesystem::path graph_folder_path(GRAPH_DIRECTORY + folder_name);
    if(!std::filesystem::exists(graph_folder_path)){
        std::filesystem::create_directories(graph_folder_path);
    }

    // 2. stream setup /////////////////////////////////////////////////////////////////////////////

    std::string graph_output = GRAPH_DIRECTORY + folder_name + GRAPH_FILENAME_END  ;
    std::string metadata     = GRAPH_DIRECTORY + folder_name + GRAPH_METADATA_FILENAME_END  ;
    std::string hashpower_file =  GRAPH_DIRECTORY + folder_name + GRAPH_HASHPOWER_FILENAME_END  ;

    std::ofstream myGraphFile(graph_output);
    std::ofstream myGraphMetadata(metadata);

    // 3. stream writing //////////////////////////////////////////////////////////////////////////

    // add metadata
    myGraphMetadata << G.getV() << "\n";
    myGraphMetadata << G.get_max_d_out() << "\n";
    myGraphMetadata << G.get_max_d_in() << "\n";
    myGraphMetadata.close();

    // save hashpower
    hashPower.saveToCsv(hashpower_file);

    // set double precision
    myGraphFile << std::setprecision(LATENCY_DOUBLE_PRECISION_OUTPUT) << std::fixed;
    unsigned int max_d_out = G.get_max_d_out();
    unsigned int max_d_in  = G.get_max_d_in();

    for(unsigned int iv = 0; iv<V; iv++){
        Graph::adjacencyList v_adj = G.get_adjacencyList(iv);
        // first the outgoing neighbors
        myGraphFile << iv;
        for(unsigned int d=0;d<max_d_out;d++){
            if(d<v_adj.out_neighbors.size()){
                // v->u
                // v:  (u,w_vu, w_uv)
                unsigned int v = iv;
                unsigned int u = v_adj.out_neighbors[d];
                double w_uv = internet.latency(u,v);
                double w_vu = internet.latency(v,u);
                myGraphFile << "," << u << "," << w_vu << "," << w_uv;
            }else{
                myGraphFile << "," << -1 << "," << "Inf" << "," << "Inf" ;
            }
        }
        myGraphFile << "\n";
        // incoming neighbors
        myGraphFile << iv;
        for(unsigned int d=0;d<max_d_in;d++){
            if(d<v_adj.in_neighbors.size()){
                // u->v
                // v:  (u,w_uv, w_vu)
                unsigned int v = iv;
                unsigned int u = v_adj.in_neighbors[d];
                double w_uv = internet.latency(u,v);
                double w_vu = internet.latency(v,u);
                myGraphFile << "," << u << "," << w_uv << "," << w_vu;
            }else{
                myGraphFile << "," << -1 << "," << "Inf" << "," << "Inf" ;
            }
        }
        myGraphFile << "\n";
    }
    myGraphFile.close();
    std::cout << " DONE. " << std::endl;
}

void formatPercentiles(std::vector<double>& percentiles){
    // prepare proper formatting for the percentile array (sorted and with values between [0,1])
    std::transform(percentiles.begin(),percentiles.end(),percentiles.begin(),[](double &p){
        double target = p/100.0;
        if (target<0) target = 0.0;
        if (target>1) target = 1.0;
        return target;
        });
    std::sort(percentiles.begin(),percentiles.end());
}


void Network::writeMetricsCSV(
    std::string filename_suffix,
    const metrics::kind kind,
    metrics::parameters parameters) {
    std::cout << "performing calculations to save metrics to csv...";
    // setup ///////////////////////////////////////////////////////////////////////////////////////

    // result : 
    struct results {
        std::vector<unsigned int> indices {};
        std::vector<std::vector<double>> values {};
    };
    // construct results. number of measures depends on randomization
    unsigned int N = parameters.randomized ? parameters.randomSamples : V;
    results result; result.indices.resize(N); result.values.resize(N);
    size_t P = parameters.percentiles.size();

    if(kind == metrics::kind::TimeToHashPowerPercent){
        formatPercentiles(parameters.percentiles);
        for(unsigned int i = 0; i<N; i++) result.values[i].reserve(P);
    }
    if(kind == metrics::kind::TimeToHashPowerPercentBound){
        formatPercentiles(parameters.percentiles);
        for(unsigned int i = 0; i<N; i++) result.values[i].reserve(P);
    }
    if(kind == metrics::kind::EdgePropagationWeight){
        for(unsigned int i = 0; i<N; i++) result.values[i] = {-1.};
    }

    //
    if(parameters.randomized){ // samples are random
        for(unsigned int i = 0; i< N; i++){result.indices[i]=getRandomNode();}
    }else{ // samples are the entire network
    for(unsigned int i = 0; i< N; i++){ result.indices[i]=i;}}

    // thread setup
    const size_t N_threads = thread_ranges.size();
    thread_work_ranges workload = getWorkRanges(N,N_threads); // separate workload of N nodes in N_threads
    struct thread_args {
        thread_work_range range;
        results* result;
    };

    // Calculation ///////////////////////////////////////////////////////////////////////////////////

    // specialzied definitions (for each kind of metric)
    auto edge_weight_calculation = [this,parameters](thread_args arg){
            auto range  = arg.range;
            auto result = arg.result;
            for(unsigned int iv = range.first; iv<range.second ; iv++){
                auto v = result->indices[iv];
                auto res = propagationLatencyEdges::calculate_single(
                    G,
                    activeNodes,
                    internet,
                    v
                );
                result->values[iv] = {res.weight};
            }
        return;
    };

    auto timetohash_calculation = [this,parameters](thread_args arg){
            auto range  = arg.range;
            auto result = arg.result;
            for(unsigned int iv = range.first; iv<range.second ; iv++){
                auto v = result->indices[iv];
                result->values[iv] =
                timeToHashPowerPercentile::calculate(
                    G,
                    hashPower.getActiveHashPower(),
                    activeNodes,
                    internet,
                    parameters.percentiles,
                    v,false);
                }
    };

    auto timetohash_bound_calculation = [this,parameters](thread_args arg){
            auto range  = arg.range;
            auto result = arg.result;
            for(unsigned int iv = range.first; iv<range.second ; iv++){
                auto v = result->indices[iv];
                result->values[iv] =
                timeToHashPowerPercentile::calculate(
                    G,
                    hashPower.getActiveHashPower(),
                    activeNodes,
                    internet,
                    parameters.percentiles,
                    v,true);
                }
    };

    // calculation threads launch and join
    std::vector<std::thread> threads{};
    for(auto range : workload){
        auto argument = thread_args{range,&result};
        if(kind == metrics::kind::TimeToHashPowerPercent){
                std::thread t{timetohash_calculation,argument};
                threads.push_back(std::move(t));
        }
        if(kind == metrics::kind::TimeToHashPowerPercentBound){
                std::thread t{timetohash_bound_calculation,argument};
                threads.push_back(std::move(t));
        }
        if(kind == metrics::kind::EdgePropagationWeight){
                std::thread t{edge_weight_calculation,argument};
                threads.push_back(std::move(t));
        }
    } for( auto& t : threads) t.join();

    // writing ///////////////////////////////////////////////////////////////////////////////////////

    std::string input = filename_suffix+".csv";
    std::ofstream myFile(input);
    if(!myFile){
        return;
    }
    // headers: [NODE] | p1 | p2 | p3 | ...
    // write header
    myFile << "v,"; for(unsigned int m=0;m<P;m++){
        if(kind!=metrics::kind::EdgePropagationWeight)
            myFile << parameters.percentiles[m];
        if (m!=P-1) myFile << ",";
    }
    myFile << "\n";

    myFile << std::setprecision(LATENCY_DOUBLE_PRECISION_OUTPUT) << std::fixed;

    // write each line
    for(unsigned int i = 0;i < N; i++){
        myFile << result.indices[i] << ","; for(unsigned int m=0;m<P;m++){
            myFile << result.values[i][m];
            if (m!=P-1) myFile << ",";
        }
        if (i!=N-1) myFile << "\n";
    }
    myFile.close();
    std::cout << " DONE. " << std::endl;
    return;
}

// requestable node management

//std::vector<nodeId> requestable {};                        // these are used to provide an O(1) updateable and
//std::unordered_map<nodeId,unsigned int> requestableMap{};  // O(1) random sample among {v: activeNodes[v]==true}

void Network::updateActiveRequestable(nodeId v, bool value){
    if(!value){ // value = false means wants to make it unreachable
        removeActiveRequestable(v);
        return;
    }
    // else, cases
    auto mapElement = requestableMap.find(v);

    // if v is present, and we dont want to set it
    // to false, then we do nothing
    if(mapElement!=requestableMap.end()) return;

    // if v it is not present, we add it at the end
    unsigned int new_index = requestable.size();
    requestable.push_back(v);
    requestableMap[v] = new_index;
}

void Network::removeActiveRequestable(nodeId v){
    auto mapElement = requestableMap.find(v);
    // if the client index is not found, then return
    if(mapElement==requestableMap.end()) return;
    // remove from map
    unsigned int index = mapElement->second;
    // separate by cases
    //     index is at the end (easy)
    if(index == requestable.size()-1){
        requestableMap.erase(v);
        requestable.pop_back();
        return;
    }
    // else, index is inside.
    // what we do is swap it to the end and pop from there
    requestableMap.erase(v);
    unsigned int last = requestable.size()-1;
    std::swap(requestable[index],requestable[last]);
    requestable.pop_back();  
    // update the requestmap
    requestableMap.at(requestable.at(index)) = index;
}

///// ProcessEventQueue /////

void Network::processEventQueue(){
    bool work_left = (nodeEventQueue.size() > 0);
    if (!work_left) return;

    // process. must be done sequentially
    std::vector<nodeEvent> activity{};
    std::vector<nodeEvent> blockValidationHash{};
    std::vector<nodeEvent> blockGenerationHash{};
    std::vector<nodeEvent> downloadBandwidth{};
    std::vector<nodeEvent> uploadBandwidth{};

    for(const auto req : nodeEventQueue){
        if(req.event==nodeEventType::activeChange) activity.push_back(req);
        if(req.event==nodeEventType::blockValidationHashPowerChange) blockValidationHash.push_back(req);
        if(req.event==nodeEventType::blockGenerationHashPowerChange) blockGenerationHash.push_back(req);
        if(req.event==nodeEventType::downloadBandwidthChange) downloadBandwidth.push_back(req);
        if(req.event==nodeEventType::uploadBandwidthChange) uploadBandwidth.push_back(req);
    }

    // activity update //////////////
    for(const auto req : activity){
        nodeId id       = req.node;
        bool is_active  = (bool) req.value;
        // perform actions
        activeNodes[id] = is_active;            // 
        updateActiveRequestable(id,is_active);  // be reachable / unreachable

        // Connection request Pool
        if(is_active){
            // node went active
            // send a request for connections
            int neighbors_needed =
                static_cast<int>(nodes[id]->connectionProtocol.maxOutgoing)
              - static_cast<int>(getOutDegree(id)) ;

            neighbors_needed = neighbors_needed >= 0 ? neighbors_needed : 0;
            nodePoolRequest(id,static_cast<unsigned int>(neighbors_needed));
        } else { // node went inactive
            //// two things to do:
            // 1. send a request to be unreachable within nodePool
            nodePoolRequest(id,0);
            // 2. disconnect each of those that have this node as preferred
            //    first get them (by copy, not reference! because the adjacency list will mutate afterwards...)
            std::vector<nodeId> nodes_that_preferred_this_one = G.get_in_neighbors(id);
            //    for each one
            for(auto u: nodes_that_preferred_this_one){
                // disconnect from graph
                disconnect(u,id);
                // update request
                int neighbors_needed =
                    static_cast<int>(nodes[u]->connectionProtocol.maxOutgoing)
                  - static_cast<int>(getOutDegree(u)) ;
                neighbors_needed = neighbors_needed >= 0 ? neighbors_needed : 0;
                nodePoolRequest(u,static_cast<unsigned int>(neighbors_needed));
            }
        }

    }
    // block generation hash update /////////////
    for(const auto req : blockGenerationHash){hashCapabilities[req.node] = req.value;}
    // finally, update hashpower if needed
    bool block_generation_hash_power_affected = 
        (activity.size()>0) || (blockGenerationHash.size()>0);

    if (block_generation_hash_power_affected)
    hashPower.updateHashDistribution(hashCapabilities, activeNodes);

    // latency related updates ////////////////
    for(const auto req : blockValidationHash) internet.setValidationHash(req.node,req.value);
    for(const auto req : downloadBandwidth) internet.setDownloadBandwidth(req.node,req.value);
    for(const auto req : uploadBandwidth) internet.setUploadBandwidth(req.node,req.value);

    // NOTE network.setBlockSize is handled directly...
}

#endif