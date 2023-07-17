#include "Constants.hpp"

#include <algorithm>
#include <functional> // std::greater
#include <iostream>
#include <limits>
#include <queue> // std::priority_queue
#include <vector>

#ifndef PERCENTILE_H
#define PERCENTILE_H

class Percentile{
    public:
        /* returns the p-th percentile of the inputVector */
        static double calculate(const std::vector<double>& inputVector,double p){ //empirically faster on instances of V=1k, 10k, 50k
            // copy into temporary array
            std::vector<double> temp(inputVector.size());
            for(unsigned int i=0;i<temp.size();i++){
                temp[i]=inputVector[i];
            }
            // sort array
            std::sort(temp.begin(),temp.end());

            // get pth value index
            int index = (int)((p/100.0)*((double)(temp.size()-1)));

            // bound safeguards
            index = (index >= 0)          ? index : 0;
            index = (index < static_cast<int>(temp.size())) ? index : static_cast<int>(temp.size())-1;

            // return the value
            return temp[static_cast<unsigned int>(index)];
        }

        // empirically proven to be slower than the other method ...
        static double calculateSlow(const std::vector<double>& inputVector,double p){
            // number of elements in queue
            //   = |{values greater than percentile p}|
            //   =  inputVector.size() * (1-p)
            size_t queue_length = static_cast<size_t>((1-(p/100.))*inputVector.size()) + 1;
            std::priority_queue<double,std::vector<double>,std::greater<double>> Q;
            for(size_t i = 0; i < inputVector.size(); i++){
                auto x = inputVector[i];
                if (Q.size() < queue_length){
                    Q.push(x);
                }else if(x <= Q.top()){
                    continue;
                }else{
                    Q.pop();Q.push(x);
                }
            }
            return Q.top();
        }
};

#endif
