#ifndef NETWORK_OBJECTS_H
#define NETWORK_OBJECTS_H

#include "Constants.hpp"
#include "Graph.hpp"
#include "GraphAlgorithms.hpp"
#include "HashPower.hpp"
#include "InternetTopology.hpp"
#include "Parallelization.hpp"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>


// class declarations

class Node;
class Network;
class NodePool;

using nodeId  = unsigned int;
using integer = unsigned int;

class ScoringProtocol;        // Specifies how the node selects which nodes to keep connection with.
class ConnectionProtocol;     // Specifies how the node searches and establishes connection with other nodes.

class ScoringProtocol {
    public:
        struct observation {
            unsigned int neighbors = 0 ; 
            std::vector<double> firstTime {};
            std::vector<double> latency   {};          // rememeber to initialize these values to infinity
            std::vector<std::vector<double>> times {};
            std::vector<int> predecessors {};
        };

        std::string name;
        const unsigned int scoredNeighbors  ;
        const unsigned int totalBlocks;

        const bool get_first_time_flag   = false;
        const bool get_times_flag        = false;
        const bool get_predecessor_flag  = false;
        const bool get_latency_flag      = false;
        
        // constructor
        ScoringProtocol(
            const std::string& name,
            const unsigned int scoredNeighbors = 5,
            const unsigned int totalBlocks = 40,
            const bool get_first_time_flag   = false,
            const bool get_times_flag        = false,
            const bool get_predecessor_flag  = false,
            const bool get_latency_flag      = false) :
            name{name} , scoredNeighbors{scoredNeighbors}, totalBlocks{totalBlocks},
            get_first_time_flag{get_first_time_flag},
            get_times_flag{get_times_flag},
            get_predecessor_flag{get_predecessor_flag},
            get_latency_flag{get_latency_flag} {}

        virtual std::vector<bool> selectNeighbors(const observation& data) const = 0;
        virtual ~ScoringProtocol() = default;
};

std::vector<bool> ScoringProtocol::selectNeighbors(const ScoringProtocol::observation& data) const {
    return std::vector<bool>(data.neighbors,true);
}

class ConnectionProtocol {

    public:
        std::string name;
        
        // constructor
        ConnectionProtocol(
            const std::string name,
            unsigned int maxOutgoing = 8,
            unsigned int maxIncoming = 25,
            unsigned int visibleRadius = 2,
            unsigned int visibilityUpdateBlocks = BIG_NUMBER_OF_BLOCKS // by default, big number means to not perform visibility update
            ) : 
            name{name},
            maxIncoming{maxIncoming},
            maxOutgoing{maxOutgoing},
            visibleRadius{visibleRadius},
            visibilityUpdateBlocks{visibilityUpdateBlocks} {}

    public:
        const unsigned int maxIncoming;
        const unsigned int maxOutgoing;
        const unsigned int visibleRadius;
        const unsigned int visibilityUpdateBlocks;
};

class NodePool {
    public:
        struct nodeRequest {
            nodeId node;
            unsigned int amount;
        }; // node node wants amount neighbors

        std::unordered_map<nodeId, unsigned int> requestMap {}; // provides index in the request vector
        std::vector<nodeId>                     subscribers {};
        std::vector<nodeRequest>                   requests {};
        unsigned int unsolvedRequests = 0;

    private:
        Network& network;
    
    public:
        const unsigned int repetition = 2;
    
    public:
        NodePool(Network& N, unsigned int rep = 2) : network{N}, repetition{rep} {}

        unsigned int getUnsolvedRequests() const {return unsolvedRequests;}
        void updateRequest(nodeId client, unsigned int request_amount = 0);
        void removeRequest(nodeId client);
        nodeRequest getRandomRequest(unsigned int randomNumber); 

// the data structure used for random sampling requests is based on:
//  https://www.geeksforgeeks.org/design-a-data-structure-that-supports-insert-delete-search-and-getrandom-in-constant-time/
};

// nodeEventType and nodeEvent are used to update information to network from node.

enum class nodeEventType {
    null,
    activeChange,
    blockValidationHashPowerChange,
    blockGenerationHashPowerChange,
    downloadBandwidthChange,
    uploadBandwidthChange
};

struct nodeEvent {
    nodeEventType event = nodeEventType::null;
    nodeId node = 0;
    double value = 0.;
};

// node

class Node {
    public:
        // Static variables
        static nodeId nodeCount;
        // Member variables
    public:
        ScoringProtocol& scoringProtocol;
        ConnectionProtocol& connectionProtocol;
        integer maxNeighbors = 0;
        const integer id = 0;

    private:
        bool   active = false;             
        double blockGenerationHashPower  = 1.  ;
        double blockValidationHashPower  = 50. ;
        double downloadBandwidth   = 1.  ;
        double uploadBandwidth     = 1.  ;

        //      network interface
        Network* network{nullptr};

        //      Observations for round
        ScoringProtocol::observation observations {0,{},{},{},{}};
        std::map<nodeId,unsigned int> neighborMapForward  {};
        std::map<unsigned int,nodeId> neighborMapBackward {};

        unsigned int currentScoringBlock = 0; 

        // for updating visibility
    public:
        bool firstVisibilityUpdateFlag = true;
        unsigned int currentVisibilityBlock = 0; 

    private:
        //      relationship to other nodes
        std::vector<nodeId> visibleNodes {};
        bool requestNodeFlag = false;

        // Member functions
    public:
        // set/getters. setters have a nodeEvent dispatch
        void setActive     (bool   val)
        { dispatchNodeEvent(nodeEventType::activeChange                   , val) ; active = val;}
        void setBlockValidationHashPower (double val)
        { dispatchNodeEvent(nodeEventType::blockValidationHashPowerChange , val) ; blockValidationHashPower = val;}
        void setBlockGenerationHashPower (double val)
        { dispatchNodeEvent(nodeEventType::blockGenerationHashPowerChange , val) ; blockGenerationHashPower  = val;}
        void setDownloadBandwidth   (double val)
        { dispatchNodeEvent(nodeEventType::downloadBandwidthChange        , val) ; downloadBandwidth   = val;}
        void setUploadBandwidth     (double val)
        { dispatchNodeEvent(nodeEventType::uploadBandwidthChange           , val) ; uploadBandwidth     = val;}

        inline bool   getActive                   () { return active;                    }
        inline double getBlockValidationHashPower () { return blockValidationHashPower;  }
        inline double getBlockGenerationHashPower () { return blockGenerationHashPower;  }
        inline double getDownloadBandwidth        () { return downloadBandwidth;         }
        inline double getUploadBandwidth          () { return uploadBandwidth;           }

    public:
        // constructor
        Node(
            ScoringProtocol& SP,
            ConnectionProtocol& CP,
            bool  active,
            double blockGenerationHashPower = DEFAULT_BLOCK_GENERATION_HASH_POWER ,
            double blockValidationHashPower = DEFAULT_BLOCK_VALIDATION_HASH_POWER ,
            double downloadBandwidth        = DEFAULT_DOWNLOAD_BANDWIDTH ,
            double uploadBandwidth          = DEFAULT_UPLOAD_BANDWIDTH 
            ) :
                scoringProtocol{SP},
                connectionProtocol{CP},
                maxNeighbors{CP.maxOutgoing},
                id{nodeCount},
                active{active},
                blockGenerationHashPower{blockGenerationHashPower},
                blockValidationHashPower{blockValidationHashPower},
                downloadBandwidth{downloadBandwidth},
                uploadBandwidth{uploadBandwidth}
                {
                    nodeCount++ ;
                    restartRound();
                };

        void setNetwork(Network* N){
            network = N;
        }

        void updateObservations(
            const GraphAlgorithms::distanceAndPaths& blockPropagationData);

        void setVisibleNodes(const std::vector<nodeId>& nodes){
            visibleNodes = nodes; // copy assignment !
        }
    
    private:
        void dispatchNodeEvent(nodeEventType type, double val);
        void restartRound(); // cleanup between rounds
    
    public:
        struct nodeRemoval {
            nodeId callerId;
            nodeId targetId;
        };

        std::vector<nodeRemoval> processSelectNeighbors();
        void processConnectionRequest();
    
    private:
        std::vector<nodeRemoval> performNeighborSelection();
    
    public:
        nodeId getRandomVisibleNode(unsigned int randomValue) 
        {return visibleNodes[randomValue%visibleNodes.size()];}
    
    public:
        // since Node has pointer data members, we need to
        // delete these functions:
        Node(const Node&) = delete;
        Node(Node&&) = delete;
        Node& operator=(const Node&) = delete;
        Node& operator=(Node&&) = delete;
        ~Node() {}
};

nodeId Node::nodeCount = 0;

using nodeRef = Node*;

namespace metrics {
    enum class kind {
        TimeToHashPowerPercent,
        TimeToHashPowerPercentBound,
        EdgePropagationWeight
    };

    struct parameters {
        std::vector<double> percentiles {{0.}};
        bool randomized = false;
        unsigned int randomSamples = 0;
        parameters() : percentiles{{0.}}, randomized{false}, randomSamples {0} {};
    };
}


class Network {
    public:
        const unsigned int V = 0;

    private:
        // theoretically important
        Graph G;
        NodePool nodePool;
        std::vector<nodeRef>& nodes;
        InternetTopology& internet;
        HashPower hashPower;

        std::vector<bool> activeNodes {};
        std::vector<double> hashCapabilities {};

        double blockSize = 0.;

        // computing
        std::vector<nodeEvent> nodeEventQueue {};
        const thread_work_ranges thread_ranges;


    public:
        Network(std::vector<nodeRef>& nodes, InternetTopology& internet, double blocksize = DEFAULT_BLOCK_SIZE,
        unsigned int maxDout  = DEFAULT_D_OUT,
        unsigned int maxDin   = DEFAULT_D_IN,
        unsigned int nThreads = DEFAULT_PARALLELIZATION_THREADS,
        unsigned int seed     = 0
        );

        double getBlockSize(){return blockSize;}
        void setBlockSize(double bSize){internet.setBlockSize(bSize);}
    
        void Loop(bool full_visibility,bool update_visibility_flag = false);
        void ProcessPool(bool full_visibility);

        // Sub loop
        void LoopVisibilityProcessing();
        void LoopBroadcastBlock();
        void LoopNodeProcessing();

        // multithread structs

        struct rangeSelect {
            thread_work_range r;
            std::vector<Node::nodeRemoval>* output;
        };

    private:
        // multithreaded range functions
        void ProcessRangeObservations(const thread_work_range range, const GraphAlgorithms::distanceAndPaths& data);
        void ProcessRangeSelectNeighbors(rangeSelect arg);

        // random node picking (for bootstrap and eventQueue)
        std::uniform_int_distribution<unsigned int> uniformRandom {};     // standard randomness support
        std::mt19937 generator{};                                         //
        unsigned int getRandomNumber(){return uniformRandom(generator);}  //

        std::vector<nodeId> requestable {};                        // these are used to provide an O(1) updateable and
        std::unordered_map<nodeId,unsigned int> requestableMap{};  // O(1) random sample among {v: activeNodes[v]==true}
        nodeId getRandomNode(){ 
            if (requestable.size()==0){
             throw std::runtime_error("tried to get random active node when all the nodes were inactive"); }
            unsigned int random_index = getRandomNumber() % requestable.size();
            return requestable[random_index];}

        // EventQueue
        void processEventQueue();

    
    public:
        // public methods
        bool isConnected(nodeId u,nodeId v) const; 
        bool isConnectedDirected(nodeId u, nodeId v) const;
        bool connect(nodeId u, nodeId v);  // attempts to connect. (u->v) return false if failed.
        void disconnect(nodeId u, nodeId v);  // deletes connection. (u->v)

        // NodePool management
        void BootstrapNetwork(); 
        void nodePoolRequest(nodeId caller,unsigned int request_amount);
        void updateVisibility(); 
        void ProcessVisibilityRange(thread_work_range range); 


        std::vector<nodeId> getOutNeighbors(nodeId askingNode) const;
        std::vector<nodeId> getInNeighbors(nodeId askingNode) const;
        unsigned int getOutDegree(nodeId caller) const;
        double latency (unsigned int u, unsigned int v) const { return internet.latency(u,v);} 

        // active node requestable management
        void updateActiveRequestable(nodeId v, bool value);
        void removeActiveRequestable(nodeId v);

        // Writing/Saving 
        void saveGraphCSV(std::string filename_suffix);
        void writeMetricsCSV(
            std::string filename_suffix,
            const metrics::kind kind,
            metrics::parameters parameters);
        
        // eventQueue
        void pushNodeEvent(nodeEvent event){ nodeEventQueue.push_back(event); }

    void run(
        unsigned int total_blocks,               // total blocks
        unsigned int blocks_per_loop,            //    batch sizes of blocks
        int full_message = 1,                // message format option
        bool full_visibility = false,        // true if each node can communicate with every other node, false if only bounded to its visibility radius
        bool update_visibility_flag = false  // true if nodes can update their visibility periodically. NOTE severely affects performance.
        );

};

// constructors ////////////////

Network::Network(
    std::vector<nodeRef>& nodes_,
    InternetTopology& internet,
    double blockSize,
    unsigned int dout,
    unsigned int din,
    unsigned int nThreads,
    unsigned int seed
    ) :
        V{static_cast<unsigned int>(nodes_.size())},
        G{Graph("Network Graph",V,dout,din)},
        nodePool{NodePool(*this)},
        nodes {nodes_},
        internet{internet},
        hashPower {HashPower(V,true)},
        blockSize{blockSize},
        nodeEventQueue{std::vector<nodeEvent>(0)},
        thread_ranges{getWorkRanges(V,nThreads)}
        {
            auto network_internet_topology_size_match = (internet.V == V);
            if (!network_internet_topology_size_match){
                throw std::runtime_error("Error: at Network Constructor: internet V and network V does not match");
            }
            generator.seed(seed); // sets seed for randomizer (affects randomization and randomized topologies)

            // we can assume V==internet.V
            std::vector<double> blockValidationHash(V,0.);
            std::vector<double> downLoadBandwidth(V,0.);
            std::vector<double> uploadBandwidth(V,0.);

            std::vector<double> BlockGenerationHashPower(V,0.);

            activeNodes.resize(V);

            for(unsigned int iv = 0; iv<V; iv++ ){
                nodeRef N = nodes[iv];
                bool is_active = N->getActive();
                activeNodes[iv] = is_active;
                updateActiveRequestable(N->id,N->getActive());

                BlockGenerationHashPower[iv] = N->getBlockGenerationHashPower();
                blockValidationHash[iv] = N->getBlockValidationHashPower();
                downLoadBandwidth[iv]    = N->getDownloadBandwidth();
                uploadBandwidth[iv]      = N->getUploadBandwidth();

                N->setNetwork(this);
            }

            // internet mangling
            internet.setValidationHash(blockValidationHash);
            internet.setDownloadBandwidth(downLoadBandwidth);
            internet.setUploadBandwidth(uploadBandwidth);

            // manual hash power initialization
            hashCapabilities = BlockGenerationHashPower;
            hashPower.updateHashDistribution(BlockGenerationHashPower, activeNodes);

            // give some room to nodeEventQueue
            nodeEventQueue.reserve(1000);

            // Bootstrap
            BootstrapNetwork();
        }

// utilities

void Node::dispatchNodeEvent(nodeEventType type, double val) {
    auto event = nodeEvent{type, this->id, val};
    if(network){ // only works if network is not nullptr, i.e. node forms part of a network
        network->pushNodeEvent(event);
    }
}

#endif