#include "Constants.hpp"

#include "PercentileTools.hpp"
#include "NetworkObjects.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#ifndef PROTOCOL_H
#define PROTOCOL_H


/* Each protocol should be initialized and then referenced by nodes who use the reference to perform their tasks.
The protocol constructor should contain its parameters (ex: d_out,d_in, etc.). */


/*                   Scoring Protocol                          
Each scoring protocol has a function to select or discard neighbors. in concrete terms,
SelectNeighbors returns an array of neighbors to select (true or false for given index)
false means the neighbors is to be dropped. by default, base class returns true for all neighbors,
and therefore, should be overriden to prescribe each protocol : perigee vanilla, perigee subset, EdgePriority etc..

base behaviour:
std::vector<bool> ScoringProtocol::selectNeighbors(const observation& data){
    return std::vector<bool>(data.neighbors,true);
}

*/

class EdgePriority : public ScoringProtocol {
    public: 
        EdgePriority(
            const std::string name,
            unsigned int scoredNeighbors_,
            unsigned int blocks=50
            ) : ScoringProtocol(
                name,
                scoredNeighbors_,
                blocks,false,false,true,false) {}

        virtual std::vector<bool> selectNeighbors(const observation& data) const override {
            if(data.neighbors < scoredNeighbors || scoredNeighbors == 0){ return std::vector<bool>(data.neighbors,true); }
            std::vector<bool>  answer(data.neighbors,false);
            std::vector<double> score(data.neighbors,0.0);

            for(unsigned int ib = 0; ib<totalBlocks; ib++){
                if (data.predecessors[ib]!=-1){ // importante
                    score[static_cast<unsigned int>(data.predecessors[ib])]+=1.0;
                }
            }
            
            std::vector<double> score_copy = score;
            std::sort(score_copy.begin(),score_copy.end());
            double cutoff = score_copy[score_copy.size() - scoredNeighbors];
            
            unsigned int selected = 0; // used to solve cutoff ties.
            for(unsigned int u = 0; u<data.neighbors;u++){
                if(score[u]>=cutoff){answer[u]=true; selected+=1;}
                if(selected>=scoredNeighbors) break;
            }

            return answer;
        }
};

class EdgePriorityChurn : public ScoringProtocol {
    public:
        const double base_importance   = 1.;
        const double target_importance = 1.;

    public: 
        EdgePriorityChurn(
            const std::string name,
            unsigned int scoredNeighbors_,
            unsigned int blocks=50,
            double base_importance   = 1.,
            double target_importance = 3.
            ) : ScoringProtocol(
                name,
                scoredNeighbors_,
                blocks,false,false,true,false),
                base_importance{base_importance},
                target_importance{target_importance} {}

        virtual std::vector<bool> selectNeighbors(const observation& data) const override {
            if(data.neighbors < scoredNeighbors || scoredNeighbors == 0){ return std::vector<bool>(data.neighbors,true); }
            std::vector<bool>  answer(data.neighbors,false);
            std::vector<double> score(data.neighbors,0.0);

            // importance = base_value + (target_value-base_value)*(b+1/B)

            for(unsigned int ib = 0; ib<totalBlocks; ib++){
                if (data.predecessors[ib]!=-1){
                    unsigned int whom = static_cast<unsigned int>(data.predecessors[ib]);
                    double round_progress = (double)(ib+1) / (double)(totalBlocks);
                    double importance 
                    = base_importance + round_progress*(target_importance-base_importance);
                    score[whom] += importance;
                }
            }
            
            std::vector<double> score_copy = score;
            std::sort(score_copy.begin(),score_copy.end());
            double cutoff = score_copy[score_copy.size() - scoredNeighbors];
            
            unsigned int selected = 0; // used to solve cutoff ties.
            for(unsigned int u = 0; u<data.neighbors;u++){
                if(score[u]>=cutoff){answer[u]=true; selected+=1;}
                if(selected>=scoredNeighbors) break;
            }

            return answer;
        }
};



// perigee protocol implementation
class PerigeeVanilla : public ScoringProtocol {
    public:
        double alpha = 90.; // target percentile

        // constructor
        PerigeeVanilla(
            const std::string& name,
            unsigned int scoredNeighbors_,
            double alpha = 90.,
            unsigned int blocks=100) :
            ScoringProtocol(name,scoredNeighbors_,blocks, true, true, false, true), alpha{alpha} {}

        // neighbors selection
        virtual std::vector<bool> selectNeighbors(const observation& data) const override {
            // default
            if(data.neighbors<scoredNeighbors){
                return std::vector<bool>(data.neighbors,false);
            }
            // score each neighbor individually
            std::vector<double> scores(data.neighbors);

            // first, time-normalize the observations
            std::vector<unsigned int> blocksIgnored;
            auto times_copy = data.times;
            for(unsigned int b=0;b<totalBlocks;b++){
                for(unsigned int u=0;u<data.neighbors;u++){
                        times_copy[u][b] = data.times[u][b]-data.firstTime[b];
                    }
            }
            // score them according to the time-normalized observations
            for(unsigned int u = 0; u<data.neighbors;u++){
                scores[u]=Percentile::calculate(times_copy[u],alpha);
            }

            std::vector<double> score_copy = scores;
            std::sort(score_copy.begin(),score_copy.end());
            double cutoff = score_copy[scoredNeighbors-1];
            std::vector<bool> selected(data.neighbors,false);

            unsigned int selected_count = 0;
            for(unsigned int u=0;u<data.neighbors;u++){
                if (scores[u]<=cutoff){
                    selected[u] = true;selected_count+=1;
                }
                if(selected_count>=scoredNeighbors) break;
            }

            return selected;
        }
};

class PerigeeSubset : public ScoringProtocol {
    public:
        double alpha; // target percentile
        // constructor
        PerigeeSubset(
            const std::string& name,
            unsigned int scoredNeighbors_,
            double alpha = 90,
            unsigned int blocks=100) :
            ScoringProtocol(name,scoredNeighbors_,blocks,true,true,false,true),
            alpha{alpha} {}
        
        // nota breve sobre la función:
        // - que la función tome observation& se usa para ahorrar memoria
        // - override se usa porque este es una función que "reemplaza" otra "selectNeighbors" de la clase ScoringProtocol 
        //     este tipo de funciones se les llama "virtuales" porque no necesariamente tienen forma concreta en la clase original.

        // entrada: observaciones ; salida: una lista de demensión [vecinos] de valores true/false que determina si el vecino se queda o se descarta en la ronda
        virtual std::vector<bool> selectNeighbors(const observation& data) const override { 
            // data es de tipo observation y sigue aproximadamente esta estructura: 
            // struct observation{                                  para un nodo v
            //     int neighbors = 0                              ; número de vecinos salientes de v 
            //     std::vector<double>                firstTime   ; lista de dimension [bloques] con los tiempos en que bloques b fueron visto por primera vez (en paper de Mao, t_v^b), o bien D(b,v)
                                                                    // sección 4.2.1 Paper
            //     std::vector<std::vector<double>>   times       ; matriz de dimensión [vecinos x bloques] con los tiempos en que los vecinos propagaron el bloque b a v.
            // };

            // un comportamiento default: si hay menos vecinos salientes que vecinos a evaluar, se asume que todos estos vecinos se conservan
            //      (y el nodo generará una request para lleanr todos los vecinos salientes)
            // en la práctica esto no sucede mucho, dado que el protocolo de reconexión, la gran mayoría de las veces, obtiene los 8 vecinos
            // y eso supera scoredNeighbors (=5).
            if(data.neighbors < scoredNeighbors){
                return std::vector<bool>(data.neighbors,true);
            }

            // Obtener los tiempos normalizados para cada bloque y nodo.  //////////////////////
            auto normalized_time = data.times; // copiamos el arreglo de tiempos
            for(unsigned int b=0;b<totalBlocks;b++){
                // if block has been seen, time normalize: time_normalized = time - first time seen (Ecuación (2) del paper, en 4.2.1)
                for(unsigned int u=0;u<data.neighbors;u++){
                        normalized_time[u][b] = data.times[u][b]-data.firstTime[b];
                    }
            }

            // Sección scoring ////////////////////////////////////////////////////////////////
            std::vector<bool> selected(data.neighbors,false);   // selected[nodo] == true <=> se queda en la ronda

            // find minimum among time normalized time observations

            // algoritmo argmin lineal clásico para encontrar el mejor puntaje
            int    best_u = -1;                                         // candidato
            double best_s = std::numeric_limits<double>::infinity();    // puntaje candidato
            for(unsigned int u = 0; u < data.neighbors ; u++ ){
                double score_u = Percentile::calculate(normalized_time[u],alpha); // calcula el percentil alpha (90) entre data.times[u] (tiempos normalizados)
                if (score_u<best_s){best_u = static_cast<int>(u);best_s=score_u;} 
            }
            if (best_u==-1){ return std::vector<bool>(data.neighbors,false); }

            // best_u has been selected, add it to the solution.
            selected[static_cast<unsigned int>(best_u)]=true;

            // from k best selected vertices, calculate k+1.
            unsigned int last_best = static_cast<unsigned int>(best_u); // to remember the last bests
            for(unsigned int k = 1; k<scoredNeighbors; k++){
                // update the set of transformed observations
                // for each t_ub, we will use min(t_ub, min(t_vb : v selected)) (sección 4.3 Paper)
                for(unsigned int u=0;u<data.neighbors;u++){
                    if (selected[u]){continue;} // ignore the nodes already selected. 

                    for(unsigned int b=0;b<totalBlocks;b++){
                        normalized_time[u][b] = std::min(normalized_time[last_best][b],normalized_time[u][b]);
                        // note that this works inductively:
                        // the first pass, we redefine t_ub as min(t_ub,t_vb) for v the best one.
                        // the next  pass, we redefine t_ub as min(t_ub,t_wb) for w the best last one.
                        // but in this case, t_wb = min(t_wb,t_vb), therefore
                        // t_ub = min(t_ub,min(t_wb, t_vb)) and so its the minimum of the two selected (v,w), and u.
                        // inductively the formula follows...
                    }
                }

                // find the best vertex among the transformed observations
                int    best_u = -1;
                double best_s = std::numeric_limits<double>::infinity();
                for(unsigned int u = 0; u<data.neighbors;u++){
                    if (selected[u]){continue;}
                    double score_u = Percentile::calculate(data.times[u],alpha);
                    if (score_u<best_s){best_u = static_cast<int>(u);best_s=score_u;}
                }
                // take care of a border case: if best_u is still -1, that means these vertices have infinite score.
                // its better to dispose of every remaining node...
                if(best_u == -1){return selected;} // this keeps the selected nodes 
                else{ selected[static_cast<unsigned int>(best_u)]=true; } // if not, then the best u is added to the selected set.
                last_best = static_cast<unsigned int>(best_u);
                // iterate...
            }

            return selected;
        }

};

class Naive : public ScoringProtocol { // this protocol does absolutely nothing
    public:
        // constructor sets flags to do nothing
        Naive(const std::string name) :
        ScoringProtocol(name,0,
        BIG_NUMBER_OF_BLOCKS // infinity here means to never process.
        ) {}

        virtual std::vector<bool> selectNeighbors(const observation& data) const override {
            return std::vector<bool>(data.neighbors,true); // everyone is to not be dropped.
        }

};

#endif