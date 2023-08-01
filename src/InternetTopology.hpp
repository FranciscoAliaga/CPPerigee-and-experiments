#ifndef INTERNET_TOPOLOGY_H
#define INTERNET_TOPOLOGY_H

#include "Constants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <tuple>
#include <utility>
#include <vector>


constexpr double pi    = 3.141592653;
constexpr double sqrt2 = 1.41412;

// pairwiseInternetLatency(u,v) := physical communication latency = ??? has a defined name in literature
// bandwidthLatency(u,v)        := time it takes to pass a block through the node's bandwidth, :=  block_size / min(upload bandwidth u, download bandwidth v)
// block_validation_time        := latency incurred by validation of a block, := block_size / block_validation_capacity
// latency(u,v)                 := pairwise_link_latency(u,v) + bandwidthLatency(u,v) + block_validation_time(u)

class InternetTopology {
    public:
        unsigned int V;

    private:
        double blockSize  = DEFAULT_BLOCK_SIZE;
        std::vector<double> blockValidationHash {};
        std::vector<double> downloadBandwidth   {};
        std::vector<double>   uploadBandwidth   {};
    
    public:
        InternetTopology(unsigned int V) :
        V {V},
        blockValidationHash{std::vector<double>(V,1.)},
        downloadBandwidth{std::vector<double>(V,1.)},
        uploadBandwidth{std::vector<double>(V,1.)}
        {} ;
        InternetTopology() : InternetTopology(0) {};
        virtual ~InternetTopology() = default;

        // returns the pairwise internet latency between two nodes in the network,
        // due to the time a message needs to travel through the internet 
        // behaviour should be overriden 
        virtual double pairwiseInternetLatency([[maybe_unused]]unsigned int u, [[maybe_unused]] unsigned int v) const {
            if (u<V) [[likely]] { return 10.; } 
            return std::numeric_limits<double>::infinity();
        }

        // returns the bandwidth internet latency between two nodes in the network,
        // due to the time a block of a certain size passes through the bandwidth of each node
        // behaviour should be overriden 
        virtual double bandwidthInternetLatency([[maybe_unused]]unsigned int u, [[maybe_unused]] unsigned int v) const {
            if (u<V) [[likely]] {
                return blockSize/std::min(downloadBandwidth[u],uploadBandwidth[v]);
                } 
            return std::numeric_limits<double>::infinity();
        }

        // returns the validation latency from a node,
        // due to the time a node spends validating the block
        // behaviour should be overriden */
        virtual double blockValidationTime([[maybe_unused]]unsigned int u) const {
            if (u<V) [[likely]] {
                double time = blockSize/blockValidationHash[u];
                return time; }
            return std::numeric_limits<double>::infinity();
        }

        /* returns the total latency from propagating a block from a node towards node v .
        behaviour could be overriden */
        virtual double latency(unsigned int u, unsigned int v) const {
            return pairwiseInternetLatency(u,v) + bandwidthInternetLatency(u,v) + blockValidationTime(u);
        }
        ////////// setters /////////////////////
        void setValidationHash(const std::vector<double>& t){
            blockValidationHash = t;
        }
        void setDownloadBandwidth(const std::vector<double>& t){
            downloadBandwidth = t;
        }
        void setUploadBandwidth(const std::vector<double>& t){
            uploadBandwidth = t;
        }

        void setValidationHash(unsigned int u, double val){
            blockValidationHash[u] = val;
        }
        void setDownloadBandwidth(unsigned int u, double val){
            downloadBandwidth[u] = val;
        }
        void setUploadBandwidth(unsigned int u, double val){
            uploadBandwidth[u] = val;
        }

        void setBlockSize(double bSize){blockSize = bSize;}
};

class FlatTopology : public InternetTopology{
    private:
        unsigned int d;
        double latencyScaleFactor;
        double validationLatency;
        std::vector<std::vector<double>> locations {};
    public:
        double p {2.};
        double invp {0.5};
    
    public:
        unsigned int get_d() const {return d;}

    public:
        // constructor
        FlatTopology(
            unsigned int V,
            unsigned int d,
            double maxPairwise,
            double validationLatency,
            unsigned int seed = 0,
            double p = 2.
            ) :
            InternetTopology(V),
            d {d>0 ? d: 1}, // d should be greater than 0
            latencyScaleFactor{maxPairwise/std::pow(static_cast<double>(d),1./p)}, // factor for correcting latency scale
            validationLatency{validationLatency},
            locations{std::vector<std::vector<double>>(V,std::vector<double>(d,0.0))},
            p{p},
            invp{1./p}
        {
            // distribute points along [0,1]^d
            std::mt19937 generator(seed);
            std::uniform_real_distribution<double> uniform(0.0,1.0);
            for(unsigned int i = 0; i<V; i++){
                for(unsigned int j = 0; j < d; j++){
                    locations[i][j]= uniform(generator);
                }
            }
            
        }
        virtual double normDistance(const std::vector<double>& X, const std::vector<double>& Y) const {
            double res = 0.;
            for(unsigned int j=0; j<d ;j++){
                res += std::pow(X[j]-Y[j],p);
            }
            double x = std::pow(res,invp);
            return x;
        }

        /* returns the pairwise internet latency between two nodes in the network. */
        virtual double pairwiseInternetLatency(unsigned int u, unsigned int v) const override {
            return normDistance(locations[u],locations[v])*latencyScaleFactor;
        }

        /* returns the validation latency from a node. */
        virtual double blockValidationTime([[maybe_unused]]unsigned int u) const override {
            return validationLatency;
        }

        // returns the total latency for propagating a block
        virtual double latency(unsigned int u, unsigned int v) const override {
            return blockValidationTime(u) + pairwiseInternetLatency(u,v);
        };

};

class FastFlatTopology : public FlatTopology {

    private:
        const double high_power_percent {10.};
        const double reduction {0.1};

    public: 
        FastFlatTopology(unsigned int V,
            unsigned int d,
            double percent_high_power,
            double reduction,
            double maxPairwise,
            double validationLatency,
            unsigned int seed = 0,
            double p = 2.
            ) : FlatTopology(V,d,maxPairwise,validationLatency,seed=0,p=2.), high_power_percent{percent_high_power}, reduction{reduction} {}
    
        // returns the total latency for propagating a block, with connection between high power neighbors faster
    double latency(unsigned int u, unsigned int v) const override {
            bool u_is_fast = u < (V*high_power_percent);
            bool v_is_fast = v < (V*high_power_percent);
            if(u_is_fast && v_is_fast){
                return blockValidationTime(u) + reduction*pairwiseInternetLatency(u,v);
            }

            return blockValidationTime(u) + pairwiseInternetLatency(u,v);
        };


};


class HubTopology : public FlatTopology {
    public:

        HubTopology( unsigned int V,
            unsigned int d,
            double maxPairwise,
            double validationLatency,
            unsigned int seed = 0,
            double p = 2.
            ) : FlatTopology(V,d,maxPairwise,validationLatency,seed=0,p=2.) {}

        // overrides
        double normDistance(const std::vector<double>& X, const std::vector<double>& Y) const override {
            double res = 0.;
            for(unsigned int j=0; j<get_d() ;j++){
                res += std::pow(X[j]-Y[j],p);
            }
            double x = std::pow(res,invp);
            return x;
            double sx = std::sqrt(x);
            double fx = sx*sx*(3.-2.*sx);
            double gx = 0.5 + 4.*std::pow(fx-0.5,3.);
            gx = std::min({1.,gx});
            gx = std::max({0.0,gx});
            return gx;
        }

};

class UniformTopology : public InternetTopology{
    public:
        UniformTopology(unsigned int V_) : InternetTopology(V_) {};

        /* returns the pairwise internet latency between two nodes in the network.
        behaviour should be overriden */
        double pairwiseInternetLatency([[maybe_unused]]unsigned int u,[[maybe_unused]]unsigned int v) const override {
            return 1.0;
        }

        double latency(unsigned int u, unsigned int v) const override {
            return pairwiseInternetLatency(u,v);
        }
};


#endif
