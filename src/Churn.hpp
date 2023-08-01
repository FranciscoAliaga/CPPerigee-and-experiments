#ifndef CHURN
#define CHURN

#include <algorithm>
#include <random>
#include <functional>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "Network.hpp"


class ChurnSimulation {
    using ii = unsigned int;
    using nodeRef = Node*;

    private:
        ii block_height {0};
        Network* N{nullptr};
        std::vector<nodeRef>* Nodes {nullptr};

        std::mt19937 generator{0};
        std::uniform_int_distribution<ii> uniformRandom {};
        std::binomial_distribution<ii> binomialRandom {};

        std::vector<unsigned int> inactive_nodes {};
        std::unordered_set<unsigned int> inactive_node_set {};
        ii inactive_nodes_limit   {0};
        long double churn_per_block{0.};

    public:
        // function to turn in random disconnected nodes 
        void churn_in(){
            ii amount = get_random_amount_of_churn();
            for(ii i = 1;i<=amount; i++){
                if (inactive_nodes.size()==0) break;
                pop();
            }
        }

        // function to remove random nodes from the network
        void churn_out(ii how_many = 0){
            // if there are too many nodes out, exit.
            if(inactive_nodes.size()>= inactive_nodes_limit) return;
            
            auto churning_out = (how_many > 0) ? how_many : get_random_amount_of_churn();

            ii successes = 0;
            while(successes < churning_out){
                auto random = uniformRandom(generator);
                successes += churn_one_out(random); // >_<
            }
        }

    private:
        struct log {
            ii height{0};
            ii   whom{0};
            bool  out{0};
        };
        std::vector<log> history {};

        ii churn_one_out(ii random){
            // two cases: inactive is empty or not
            if(inactive_node_set.size()>0){
                // look if it is inside otherwise add
                auto check = std::find(inactive_node_set.begin(),inactive_node_set.end(),random);
                if (check==inactive_node_set.end()){
                    // not in, add
                    push(random);
                    return 1; // success
                } else {
                    // do nothing
                    return 0; // unsuccesful
                }
            } else {
                push(random);
                return 1; // success
            }
        }

        void push(ii random){

            inactive_nodes.push_back(random);
            inactive_node_set.insert(random);
            // network set inactive
            Nodes->at(random)->setActive(false);
            // log a node went out
            history.push_back(
                {block_height,random,true}
            );
        }

        void pop(){
            ii S = inactive_nodes.size();
            if(S==0) return;
            ii random_id = uniformRandom(generator)%S;
            // swap to last place and pop
            std::swap(inactive_nodes[random_id],inactive_nodes[S-1]);
            ii random = inactive_nodes[S-1];
            inactive_nodes.pop_back();
            inactive_node_set.erase(random);
            // network set active
            Nodes->at(random)->setActive(true);
            // log that a node went in
            history.push_back(
                {block_height,random,false}
            );
        }

        ii get_random_amount_of_churn(){
            // MODEL: each node of the network might activate/deactivate randomly and independent of each other
            // then, amount = #{ v : v activated/deactivated at random during this block } ~ Binomial(V,p)
            // where p = churn at each block / nodes in the network.

            ii amount = binomialRandom(generator); // complexity: amortized constant.
            
            // NOTE: p is set at constructor
            return amount;
        }

    public:
        void write_history_log(std::string name){
            std::string input = name+".csv";
            std::ofstream myFile(input);
            if(!myFile){
                return;
            }
            // write header
            myFile << "block height, v, in/out\n";
            // start writing
            for(auto log : history){
                myFile << log.height << ',';
                myFile << log.whom << ',';
                myFile << (log.out ? "out" : "in") << '\n';
            }
            myFile.close();
            return;
        }


    public:
        ChurnSimulation(Network& N, std::vector<nodeRef>& Nodes,
        long double churn_per_block,
        long double inactive_node_percentage_limit,
        ii seed = 0) :
            N{&N}, Nodes{&Nodes},
            inactive_nodes_limit{static_cast<ii>(N.V * inactive_node_percentage_limit)},
            churn_per_block{churn_per_block}
        {
            generator.seed(seed);
            const long double p = static_cast<long double>(
                churn_per_block / N.V
            );
            // set binomial
            binomialRandom.param(
                std::binomial_distribution<ii>::param_type(N.V,p));
            // set uniform random
            uniformRandom.param(
                std::uniform_int_distribution<ii>::param_type(0,N.V-1)
            );
        }

        ChurnSimulation operator=(const ChurnSimulation&) = delete;
        ChurnSimulation(const ChurnSimulation&) = delete;

        void run(ii blocks, bool full_visibility = 0, bool update_visibility = 1){
            std::cout << "Simulation of " << blocks << " blocks... \n";
            for(ii b=1; b <= blocks ; b++){
                // message
                std::cout << 
                "\r" << b << " / " << blocks << "        ";
                // the actual stuff...
                N->Loop(full_visibility,update_visibility);
                churn_out();
                churn_in();
                block_height+=1;
            }
            std::cout << std::endl << "Simulation of " << blocks << " blocks : DONE. " << std::endl;
        }

};

#endif