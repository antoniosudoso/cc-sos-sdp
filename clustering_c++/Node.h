#ifndef CLUSTERING_NODE_H
#define CLUSTERING_NODE_H

#include <armadillo>
#include <map>
#include <set>
#include <vector>
#include <MatlabDataArray/TypedArray.hpp>

class Node {

public:

    // must link constraints for kmeans
    std::vector<std::pair<int, int>> global_ml_pairs;
    // cannot link constraints for kmeans
    std::vector<std::pair<int, int>> global_cl_pairs;

    // lower bound
    double lb;
    // upper bound
    double ub;
    // node id
    int id;

};


class SDPNode : public Node {


public:

    std::vector<std::vector<arma::sp_mat>> B_vector;

};

typedef struct NodeData {

    SDPNode *node;
    int i;
    int j;

} NodeData;

typedef struct JobData {

    int type;
    NodeData *node_data;

} JobData;

typedef struct SDPResult {

    int n; // current problem size
    int sdp_flag;
    double lb;
    double ub;
    arma::sp_mat X_assignment;
    int n_ineq;
    int cp_iter;
    int cp_flag;
    std::pair<int, int> branching;
    std::vector<std::vector<arma::sp_mat>> B_vector;

} SDPResult;


#endif //CLUSTERING_NODE_H
