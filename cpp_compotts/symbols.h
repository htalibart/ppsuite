#ifndef SYMBOLS_1
#define SYMBOLS_1
using namespace std;

//Real special values
#undef INFINITY
#define INFINITY 1000000000

//Stats of an object
#define NOT_INITIALIZED 0
#define INITIALIZED     1
#define WORKING         2
#define FINISHED        3

//Stats of a solution
#define OPTIMAL         100
#define EPSILON         101
#define APPROXIMATE     102
#define MAX_ITER        103
#define NOT_SIMILAR     104

//Reason of an approximation
#define TIME_LIMIT      200
#define NODE_NUMBER     201
#define ERROR           202

//Error codes
#define MEM_ALLOCATION  300
#define UNKNOW          301

#endif
