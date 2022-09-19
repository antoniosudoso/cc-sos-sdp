#ifndef CLUSTERING_SDP_BRANCH_AND_BOUND_H
#define CLUSTERING_SDP_BRANCH_AND_BOUND_H

#include <armadillo>
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include "JobQueue.h"

typedef struct UserConstraints {

	double gamma, delta;
	std::vector<std::pair<int,int>> ml_pairs, cl_pairs;

} UserConstraints;

typedef struct MatlabStruct {

    std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;
    matlab::data::ArrayFactory factory;

} MatlabStruct;


typedef struct SharedData {

    // Between workers and main
    std::condition_variable mainConditionVariable;
    std::vector<bool> threadStates;

    // Queue of requests waiting to be processed
    JobAbstractQueue *queue;
    // This condition variable is used for the threads to wait until there is work to do
    std::condition_variable queueConditionVariable;
    // Mutex to protect queue
    std::mutex queueMutex;

    double global_ub;
    arma::sp_mat global_X;
    double gap;
    int n_nodes;

} SharedData;

typedef struct InputData {

    arma::mat Ws;
    arma::mat C; // cluster sizes

} InputData;


arma::sp_mat sdp_branch_and_bound(arma::mat &Ws, arma::mat &C);
std::pair<JobData *, JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data);
std::pair<JobData *, JobData *> build_cl_problem(NodeData *job_data, InputData *input_data, SharedData *shared_data);
std::pair<JobData *, JobData *> build_ml_problem(NodeData *job_data, InputData *input_data, SharedData *shared_data);
#endif //CLUSTERING_SDP_BRANCH_AND_BOUND_H
