#include "Constants.hpp"

#include <string>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>


#ifndef HASHPOWER_H
#define HASHPOWER_H

/* Takes care of producing the source of a block in the network, given a distribution of hashpower from the nodes.
Can be updated to support varying hash power capabilities as well as deactivating unactive nodes.
The default class implements uniform hash power. */
class HashPower {

    private:
        std::uniform_real_distribution<double> uniformRandom = std::uniform_real_distribution<double>(0.0,1.0); // in order to sample using the generalized inverse of uniform variable theorem.
        std::mt19937 generator = std::mt19937(0); // pseudorandom number generator using seed 0.

    private:
        unsigned int V;
        std::vector<double> distribution {};      // distribution (in the probability sense) used for the no sampling method
        std::vector<unsigned int> nodeMap {};     // used for the node sampling method
        std::vector<double> hashCapabilities {}; // used to store last hashcapabilities.
        std::vector<double> activeHashPower {};   // holds the percentage of hash power each node holds (0 if innactive)

        /* constructor */

    public:
        HashPower(unsigned int V, bool manual_initialization) :
        V{V}
        {  
            if(!manual_initialization){
                std::vector<double> UniformHashPowerCapabilities(V,1.0);
                std::vector<bool>   EveryoneIsActive(V,true);
                updateHashDistribution(UniformHashPowerCapabilities,EveryoneIsActive);
            }
        }
        HashPower(unsigned int V) : HashPower(V, false) {}

        /* main functions */
    
    public:
        /* produce a list of `NumberOfBlocks` valid blocks sources in the network via random sampling */
        std::vector<unsigned int> generateBlocks(unsigned int numberOfBlocks){
            std::vector<unsigned int> results(numberOfBlocks,0);
            for (unsigned int b = 0; b<numberOfBlocks; b++){
                results[b] = generateBlock();
            }
            return results;
        }

        inline unsigned int generateBlock(){
                double random_val = uniformRandom(generator);
                return generalizedInverse(random_val);
        }

        /* Updates the sampling distribution by using the values of hashing capabilities `hashCapabilities` of the nodes in the network,
        considering if they are active via the `activeNodes` list.`` */
        void updateHashDistribution(const std::vector<double>& _hashCapabilities, const std::vector<bool>& activeNodes){
            /* NOTE this is the most important critical part of the sampling method implementation */
            if ( &_hashCapabilities!= &hashCapabilities){
                hashCapabilities = _hashCapabilities; // by copy
            }

            if ((_hashCapabilities.size()!= V) || (activeNodes.size()!=V)){
                std::cout << "Error on hash power vector size. HashDistribution Error." << std::endl;
                return;
            }

            // note: this can be done more efficiently, but for now, it is expected to be an operation
            // not performed so often (just in a few node churning phases) therefore, we favour readability
            // over efficiency. Maybe it is the case that the compiler optimizes.

            // calculate total hashing
            double total_hashing = 0;
            for(unsigned int iv=0;iv<V;iv++){
                if (activeNodes[iv]){
                total_hashing += _hashCapabilities[iv];
                }
            }

            // construct a map from the first M integers to active nodes, where M is the number of active nodes.
            unsigned int currentNode = 0 ; 
            nodeMap.clear();nodeMap.reserve(V);
            for(unsigned int iv=0;iv<V;iv++){
                if (activeNodes[iv]){
                    nodeMap.push_back(iv);
                    currentNode += 1; 
                }
            }
            unsigned int lastNode = currentNode; // currentNode ends up counting the number M of distinct active nodes.

            // now we construct the accumulated distribution of hash power.

            // cleaning
            distribution.resize(lastNode);
            for(unsigned int i=0; i<lastNode;i++){
                distribution[i] = 0.0;
            }
            
            // this 'if' is to take care of a border case
            if (lastNode>=1){
                distribution[0] = _hashCapabilities[nodeMap[0]]/total_hashing;
                for(unsigned int i=1;i<lastNode;i++){
                    distribution[i] = distribution[i-1] + _hashCapabilities[nodeMap[i]]/total_hashing;
                }
            }
            else{
                distribution[0]=1.0;
            }

            // active hash power
            if (activeHashPower.size()!=V){
                activeHashPower.resize(V);
            }
            for(unsigned int i = 0; i<V ; i++){
                activeHashPower[i] = _hashCapabilities[i]/total_hashing * (double)(activeNodes[i]);
            }

        }

        void updateHashDistribution(const std::vector<bool>& activeNodes){
            updateHashDistribution(hashCapabilities,activeNodes);
        }

        /* sets the seed of the mersenne twister pseudorandom number generator */
        void setSeed(unsigned int seed){ generator.seed(seed); }

        /* returns reference to the vector activeHashPower which holds the percentage of current hash power of each
        node in the network (0 if not active)  */
        const std::vector<double>& getActiveHashPower() const{
            return activeHashPower;
        }
    
    public:
        void saveToCsv(std::string& filename){
            // assumes file can be written. 
            std::ofstream myHashPowerFile(filename);

            // set double precision
            myHashPowerFile << std::setprecision(HASHPOWER_DOUBLE_PRECISION_OUTPUT) << std::fixed;

            // write headers
            myHashPowerFile << "Node, HashPower" << "\n";

            // write data
            for(unsigned int iv = 0; iv < V ; iv++){
                myHashPowerFile << iv << "," << activeHashPower[iv] << "\n";
            }

            myHashPowerFile.close();
        }

            
    
    private:

        // using binary search, returns the least x such that distribution(x) >= target 
        unsigned int generalizedInverse(double target) const {
            // note: current implementation seems stable, but could be made safer using std::upperbound from the <algorithm> header.
            int L = 0;
            int R = distribution.size()-1;
            // the search space begins as [L .... R]
            // each iteration, it is cut by [L ..cut...R]
            unsigned int res = static_cast<unsigned int>(R);
            while(R-L>0){
                unsigned int cut = static_cast<unsigned int>(L + ((R-L)/2)); // cuts the search space in half
                if (distribution[cut]>=target){
                    // the cut is equal or greater than x, therefore, there is no need to keep searching at the right of the cut.
                    // the remaining search space is then [L ... cut]
                    res = cut;
                    R   = static_cast<int>(cut)-1;
                }
                else{
                    // distribution[cut]<target
                    // the cut is is less than x, ...
                    L = static_cast<int>(cut)+1;
                }
            }
            // border case
            if (distribution[static_cast<unsigned int>(L)]>=target){
                res= static_cast<unsigned int>(L);
            }

            return nodeMap[res];
        }
};

#endif