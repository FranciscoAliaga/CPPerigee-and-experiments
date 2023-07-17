#ifndef SIMULATION_CONSTANTS_H
#define SIMULATION_CONSTANTS_H


const double DEFAULT_D_OUT                = 8  ; // 25
const double DEFAULT_D_SCORED             = 5  ; // 10 
const double DEFAULT_D_IN                 = 25 ; // 125

const double DEFAULT_VISIBILITY_RADIUS    = 2  ;

const double DEFAULT_BLOCK_SIZE                  = 50.0 ;
const double DEFAULT_BLOCK_GENERATION_HASH_POWER = 1.   ; // [block hashes/s] (relative)
const double DEFAULT_BLOCK_VALIDATION_HASH_POWER = 1.   ; // [MB/s]  
const double DEFAULT_DOWNLOAD_BANDWIDTH          = 1.   ; // [MB/s]
const double DEFAULT_UPLOAD_BANDWIDTH            = 1.   ; // [MB/s]
                                                         
const unsigned int HASHPOWER_DOUBLE_PRECISION_OUTPUT      = 32;
const unsigned int LATENCY_DOUBLE_PRECISION_OUTPUT        = 32;

const unsigned int DEFAULT_NODEPOOL_MINIMUM_TRIES = 1000;

const unsigned int BIG_NUMBER_OF_BLOCKS = 100000000;
// it is used in some protocolos that depend on number of block seen, to represent "never perform this task".

#define GRAPH_DIRECTORY "./graphs/"
#define GRAPH_FILENAME_END "/graph.csv"
#define GRAPH_METADATA_FILENAME_END "/metadata.csv"
#define GRAPH_HASHPOWER_FILENAME_END "/hashpower.csv"
#define GRAPH_PERCENTILE_FILENAME_END "/percentile_metric.csv"

#define CURVE_DIRECTORY "./curves/"
#define CURVE_FILENAME_END "/metric_curve.csv"
#define CURVE_TIME_FILENAME_END  "/metric_time_spent.csv"

#endif