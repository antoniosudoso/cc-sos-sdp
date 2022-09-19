#include <thread>
#include "matlab_util.h"
#include "sdp_branch_and_bound.h"
#include "JobQueue.h"
#include "util.h"
#include "config_params.h"
#include "Node.h"
#include "ThreadPool.h"

// root
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr,
                    matlab::data::ArrayFactory &factory, arma::mat &Ws, arma::mat &C) {

    // convert data
    matlab::data::TypedArray<double> matlab_Ws = arma_to_matlab_matrix(factory, Ws);
    matlab::data::TypedArray<double> matlab_C = arma_to_matlab_matrix(factory, C);

    // Create StructArray
    std::vector<std::string> f = {"n_threads", "bb_tol", "sdp_verbose", "sdp_tol",
                                  "sdp_maxiter", "sdp_maxtime", "cp_maxiter", "cp_tol",
                                  "cp_maxineq", "cp_percineq", "cp_epsineq", "cp_activeineq",
                                  "cp_inheritineq", "gurobi_verbose"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1},f);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(matlab_session_threads_root);
    struct_matlab[0]["bb_tol"] = factory.createScalar<double>(branch_and_bound_tol);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(sdp_solver_tol);
    struct_matlab[0]["sdp_maxiter"] = factory.createScalar<double>(sdp_solver_max_iter);
    struct_matlab[0]["sdp_maxtime"] = factory.createScalar<double>(sdp_solver_max_time);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(cp_max_iter_root);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(cp_max_ineq);
    struct_matlab[0]["cp_percineq"] = factory.createScalar<double>(cp_perc_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(cp_eps_active);
    struct_matlab[0]["cp_inheritineq"] = factory.createScalar<int>(cp_inherit_ineq);
    struct_matlab[0]["gurobi_verbose"] = factory.createScalar<double>(heuristic_verbose);

    std::vector<matlab::data::Array> args({matlab_Ws, matlab_C, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    matlabPtr->eval(u"clear");
    std::vector<matlab::data::Array> result;
    if (sdp_relaxation == MATRIX_LIFTING)
        result = matlabPtr->feval(u"call_solve_matrix_lifting_root",n_return, args);
    else if (sdp_relaxation == VECTOR_LIFTING)
        result = matlabPtr->feval(u"call_solve_vector_lifting_root",n_return, args);
    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::SparseArray<double> field_best_Xass = structArray[0]["best_Xass"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::TypedArray<double> field_termcode_list = structArray[0]["termcode_list"];
    matlab::data::TypedArray<double> field_ineq_list = structArray[0]["ineq_list"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_B_cell"];
    // matlab::data::TypedArray<double> field_n = structArray[0]["n"];
    matlab::data::TypedArray<double> field_i_idx = structArray[0]["i_idx"];
    matlab::data::TypedArray<double> field_j_idx = structArray[0]["j_idx"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::sp_mat X_assignment = matlab_to_arma_sparse(field_best_Xass);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    int sdp_flag = (int) field_termcode_list[cp_iter];
    int n_ineq = (int) field_ineq_list[cp_iter];
    std::vector<std::vector<arma::sp_mat>> B_vector;
    if (sdp_relaxation == MATRIX_LIFTING) {
        B_vector.reserve(1);
        std::vector<arma::sp_mat> B_inner = matlab_to_arma_vector_sp_mat(field_best_B_cell);
        B_vector.push_back(B_inner);
    } else if (sdp_relaxation == VECTOR_LIFTING) {
        B_vector = matlab_to_arma_B_cell(field_best_B_cell);
    }
    int i_idx = (int) field_i_idx[0];
    int j_idx = (int) field_j_idx[0];
    std::pair<int, int> branching = std::make_pair(i_idx, j_idx);
    // int n = (int) field_n[0];
    int n = (int) Ws.n_rows;

    return SDPResult{n, sdp_flag, lb, ub, X_assignment, n_ineq, cp_iter, cp_flag, branching, B_vector};
}

// lower bound must link and cannot link
SDPResult solve_sdp(std::unique_ptr<matlab::engine::MATLABEngine> &matlabPtr, matlab::data::ArrayFactory &factory,
        arma::mat &Ws, arma::mat &C, std::vector<std::vector<arma::sp_mat>> &parent_B_vector, double global_ub, arma::sp_mat &global_X,
        std::vector<std::pair<int, int>> &global_ml_pairs, std::vector<std::pair<int, int>> &global_cl_pairs) {

    // convert data
    matlab::data::TypedArray<double> matlab_Ws = arma_to_matlab_matrix(factory, Ws);
    matlab::data::TypedArray<double> matlab_C = arma_to_matlab_matrix(factory, C);
    matlab::data::TypedArray<double> matlab_ml = vector_pair_to_matlab_matrix(factory, global_ml_pairs);
    matlab::data::TypedArray<double> matlab_cl = vector_pair_to_matlab_matrix(factory, global_cl_pairs);
    matlab::data::CellArray matlab_B_cell = factory.createCellArray({0, 0});
    if (sdp_relaxation == MATRIX_LIFTING)
        matlab_B_cell = arma_to_matlab_cell(factory, parent_B_vector[0]);
    else if (sdp_relaxation == VECTOR_LIFTING)
        matlab_B_cell = arma_to_matlab_vector_cell(factory, parent_B_vector);
    matlab::data::SparseArray<double> matlab_global_X = arma_to_matlab_sparse(factory, global_X);
    matlab::data::TypedArray<double> matlab_gub = factory.createScalar<double>(global_ub);

    // Create StructArray
    std::vector<std::string> f = {"n_threads", "bb_tol", "sdp_verbose", "sdp_tol",
                                  "sdp_maxiter", "sdp_maxtime", "cp_maxiter", "cp_tol",
                                  "cp_maxineq", "cp_percineq", "cp_epsineq", "cp_activeineq",
                                  "cp_inheritineq", "gurobi_verbose"};

    matlab::data::StructArray struct_matlab = factory.createStructArray({1, 1},f);
    struct_matlab[0]["n_threads"] = factory.createScalar<int>(matlab_session_threads_child);
    struct_matlab[0]["bb_tol"] = factory.createScalar<double>(branch_and_bound_tol);
    struct_matlab[0]["sdp_verbose"] = factory.createScalar<int>(sdp_solver_verbose);
    struct_matlab[0]["sdp_tol"] = factory.createScalar<double>(sdp_solver_tol);
    struct_matlab[0]["sdp_maxiter"] = factory.createScalar<double>(sdp_solver_max_iter);
    struct_matlab[0]["sdp_maxtime"] = factory.createScalar<double>(sdp_solver_max_time);
    struct_matlab[0]["cp_maxiter"] = factory.createScalar<int>(cp_max_iter_child);
    struct_matlab[0]["cp_tol"] = factory.createScalar<double>(cp_tol);
    struct_matlab[0]["cp_maxineq"] = factory.createScalar<double>(cp_max_ineq);
    struct_matlab[0]["cp_percineq"] = factory.createScalar<double>(cp_perc_ineq);
    struct_matlab[0]["cp_epsineq"] = factory.createScalar<double>(cp_eps_ineq);
    struct_matlab[0]["cp_activeineq"] = factory.createScalar<double>(cp_eps_active);
    struct_matlab[0]["cp_inheritineq"] = factory.createScalar<int>(cp_inherit_ineq);
    struct_matlab[0]["gurobi_verbose"] = factory.createScalar<double>(heuristic_verbose);

    std::vector<matlab::data::Array> args({matlab_Ws, matlab_C, matlab_ml, matlab_cl, matlab_B_cell,
                                           matlab_gub, matlab_global_X, struct_matlab});

    // Call MATLAB function and return result
    const size_t n_return = 1;
    matlabPtr->eval(u"clear");
    std::vector<matlab::data::Array> result;
    if (sdp_relaxation == MATRIX_LIFTING)
        result = matlabPtr->feval(u"call_solve_matrix_lifting_child",n_return, args);
    else if (sdp_relaxation == VECTOR_LIFTING)
        result = matlabPtr->feval(u"call_solve_vector_lifting_child",n_return, args);

    matlab::data::StructArray structArray = result[0];
    matlab::data::TypedArray<double> field_best_lb = structArray[0]["best_lb"];
    matlab::data::TypedArray<double> field_best_ub = structArray[0]["best_ub"];
    matlab::data::SparseArray<double> field_best_Xass = structArray[0]["best_Xass"];
    matlab::data::TypedArray<double> field_cp_iter = structArray[0]["cp_iter"];
    matlab::data::TypedArray<double> field_cp_flag = structArray[0]["cp_flag"];
    matlab::data::TypedArray<double> field_termcode_list = structArray[0]["termcode_list"];
    matlab::data::TypedArray<double> field_ineq_list = structArray[0]["ineq_list"];
    matlab::data::CellArray field_best_B_cell = structArray[0]["best_B_cell"];
    matlab::data::TypedArray<double> field_n = structArray[0]["n"];
    matlab::data::TypedArray<double> field_i_idx = structArray[0]["i_idx"];
    matlab::data::TypedArray<double> field_j_idx = structArray[0]["j_idx"];

    double lb = field_best_lb[0];
    double ub = field_best_ub[0];
    arma::sp_mat X_assignment = matlab_to_arma_sparse(field_best_Xass);
    int cp_iter = (int) field_cp_iter[0];
    int cp_flag = (int) field_cp_flag[0];
    int sdp_flag = (int) field_termcode_list[cp_iter];
    int n_ineq = (int) field_ineq_list[cp_iter];
    std::vector<std::vector<arma::sp_mat>> B_vector;
    if (sdp_relaxation == MATRIX_LIFTING) {
        B_vector.reserve(1);
        std::vector<arma::sp_mat> B_inner = matlab_to_arma_vector_sp_mat(field_best_B_cell);
        B_vector.push_back(B_inner);
    } else if (sdp_relaxation == VECTOR_LIFTING) {
        B_vector = matlab_to_arma_B_cell(field_best_B_cell);
    }    int i_idx = (int) field_i_idx[0];
    int j_idx = (int) field_j_idx[0];
    std::pair<int, int> branching = std::make_pair(i_idx, j_idx);
    int n = (int) field_n[0];

    return SDPResult{n, sdp_flag, lb, ub, X_assignment, n_ineq, cp_iter, cp_flag, branching, B_vector};
}



std::pair<JobData *, JobData *> create_cl_ml_jobs(double node_gap, SDPNode *node, std::pair<int, int> branching_pair,
                                                  NodeData *parent, SharedData *shared_data) {

	if (std::isinf(node->lb) || node_gap <= branch_and_bound_tol) {
        delete(node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::make_pair(nullptr, nullptr);
    }

    int i = branching_pair.first;
    int j = branching_pair.second;

    if (i == -1 && j == -1) {

        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        log_file << "PRUNING BY OPTIMALITY " << node->id << "\n";
        delete (node);
        if (parent != nullptr) {
            delete (parent->node);
            delete (parent);
        }
        return std::make_pair(nullptr, nullptr);

        // mutex is automatically released when lock goes out of scope
    }

    auto *cl_data = new NodeData();
    cl_data->node = new SDPNode(*node);
    cl_data->i = i;
    cl_data->j = j;

    auto *ml_data = new NodeData();
    ml_data->node = new SDPNode(*node);
    ml_data->i = i;
    ml_data->j = j;

    auto *cl_job_data = new JobData();
    cl_job_data->type = CANNOT_LINK;
    cl_job_data->node_data = cl_data;

    auto *ml_job_data = new JobData();
    ml_job_data->type = MUST_LINK;
    ml_job_data->node_data = ml_data;

    if (parent != nullptr) {
        delete (parent->node);
        delete (parent);
    }

    delete (node);

    return std::make_pair(cl_job_data, ml_job_data);

}

std::pair<JobData *, JobData *> build_cl_problem(NodeData *node_data, InputData *input_data, SharedData  *shared_data) {

    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(sdp_solver_folder, gurobi_folder);

    // generate cannot link child
    auto cl_node = new SDPNode();
	auto parent = node_data->node;

	double parent_gap = (shared_data->global_ub - parent->lb) / shared_data->global_ub;
	if (parent_gap <= branch_and_bound_tol)
        return std::make_pair(nullptr, nullptr);

    cl_node->global_ml_pairs = parent->global_ml_pairs;
    cl_node->global_cl_pairs = parent->global_cl_pairs;
    cl_node->global_cl_pairs.emplace_back(node_data->i, node_data->j);

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult  sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
              input_data->Ws, input_data->C, parent->B_vector,
              shared_data->global_ub, shared_data->global_X,
              cl_node->global_ml_pairs, cl_node->global_cl_pairs);

    int n = sdp_result.n;
    int sdp_flag = sdp_result.sdp_flag;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    cl_node->lb = std::max(sdp_result.lb, parent->lb);
    cl_node->ub = sdp_result.ub;
    cl_node->B_vector = sdp_result.B_vector;
    arma::sp_mat X_assignment = sdp_result.X_assignment;
    std::pair<int, int> branching_pair = sdp_result.branching;

    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        bool ub_updated = false;
        if (cl_node->ub - shared_data->global_ub <= -branch_and_bound_tol) {
            // update global upper bound
            shared_data->global_ub = cl_node->ub;
            shared_data->global_X = X_assignment;
            ub_updated = true;
        }

        cl_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;
        int open = shared_data->queue->getSize();

        node_gap = (shared_data->global_ub - cl_node->lb) / shared_data->global_ub;
        double gap = node_gap;
        Node *min_lb_node = shared_data->queue->getMinLb();
        if (min_lb_node != nullptr)
            gap = (shared_data->global_ub - min_lb_node->lb) / shared_data->global_ub;

        shared_data->gap = gap;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        auto time = (double) duration.count();

        print_log_sdp(log_file, n, parent->id, cl_node->id,
                      parent->lb, cl_node->lb,sdp_flag, time, cp_iter, cp_flag, n_ineq,
                      cl_node->ub,shared_data->global_ub, node_data->i, node_data->j, node_gap,
                      shared_data->gap, open, ub_updated);

    }

    delete(matlab_struct);

    return create_cl_ml_jobs(node_gap, cl_node, branching_pair, node_data, shared_data);
}

std::pair<JobData *, JobData *> build_ml_problem(NodeData *node_data, InputData *input_data, SharedData *shared_data) {

    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(sdp_solver_folder, gurobi_folder);

    // generate must link child
    SDPNode *ml_node;
	SDPNode *parent = node_data->node;

	double parent_gap = (shared_data->global_ub - parent->lb) / shared_data->global_ub;
	if (parent_gap <= branch_and_bound_tol)
        return std::make_pair(nullptr, nullptr);

    ml_node = new SDPNode();
    ml_node->global_ml_pairs = parent->global_ml_pairs;
    ml_node->global_ml_pairs.emplace_back(node_data->i, node_data->j);
    ml_node->global_cl_pairs = parent->global_cl_pairs;

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult  sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
                                      input_data->Ws, input_data->C, parent->B_vector,
                                      shared_data->global_ub,  shared_data->global_X,
                                      ml_node->global_ml_pairs, ml_node->global_cl_pairs);

    int n = sdp_result.n;
    int sdp_flag = sdp_result.sdp_flag;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    ml_node->lb = std::max(sdp_result.lb, parent->lb);
    ml_node->ub = sdp_result.ub;
    ml_node->B_vector = sdp_result.B_vector;
    arma::sp_mat X_assignment = sdp_result.X_assignment;
    std::pair<int, int> branching_pair = sdp_result.branching;

    double node_gap;

    {
        const std::lock_guard<std::mutex> lock(shared_data->queueMutex);

        bool ub_updated = false;
        if (ml_node->ub - shared_data->global_ub <= -branch_and_bound_tol) {
            // update global upper bound
            shared_data->global_ub = ml_node->ub;
            shared_data->global_X = X_assignment;
            ub_updated = true;
        }

        ml_node->id = shared_data->n_nodes;
        shared_data->n_nodes++;

        int open = shared_data->queue->getSize();

        node_gap = (shared_data->global_ub - ml_node->lb) / shared_data->global_ub;

        double gap = node_gap;
        Node *min_lb_node = shared_data->queue->getMinLb();
        if (min_lb_node != nullptr)
            gap = (shared_data->global_ub - min_lb_node->lb) / shared_data->global_ub;

        shared_data->gap = gap;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
        auto time = (double) duration.count();

        print_log_sdp(log_file, n, parent->id, ml_node->id, parent->lb, ml_node->lb,
                      sdp_flag, time, cp_iter, cp_flag, n_ineq,ml_node->ub,
					  shared_data->global_ub, node_data->i, node_data->j, node_gap,
                      shared_data->gap, open, ub_updated);

    }

    // mutex is automatically released when lock goes out of scope

    delete(matlab_struct);

    return create_cl_ml_jobs(node_gap, ml_node, branching_pair, node_data, shared_data);
}

std::pair<JobData *, JobData *> build_root_problem(MatlabStruct *matlab_struct, InputData *input_data, SharedData *shared_data) {

    // init root
    SDPNode *root;
    root = new SDPNode();
    root->id = shared_data->n_nodes;
    root->global_ml_pairs = {};
    root->global_cl_pairs = {};

    auto start_time = std::chrono::high_resolution_clock::now();

    SDPResult sdp_result = solve_sdp(matlab_struct->matlabPtr, matlab_struct->factory,
                                     input_data->Ws, input_data->C);

    int n = sdp_result.n;
    int sdp_flag = sdp_result.sdp_flag;
    int cp_iter = sdp_result.cp_iter;
    int cp_flag = sdp_result.cp_flag;
    int n_ineq = sdp_result.n_ineq;
    root->lb = sdp_result.lb;
    root->ub = sdp_result.ub;
    root->B_vector = sdp_result.B_vector;
    arma::sp_mat X_assignment = sdp_result.X_assignment;
    std::pair<int, int> branching_pair = sdp_result.branching;

    // update shared data
    shared_data->global_ub = root->ub;
    shared_data->global_X = X_assignment;
    shared_data->n_nodes++;

    int open = shared_data->queue->getSize();

    double node_gap = (shared_data->global_ub - root->lb) / shared_data->global_ub;
    shared_data->gap = node_gap;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    auto time = (double) duration.count();

    print_log_sdp(log_file, n, -1, root->id, -std::numeric_limits<double>::infinity(),
                  root->lb,sdp_flag, time, cp_iter, cp_flag, n_ineq, root->ub,
				  shared_data->global_ub, -1, -1, node_gap, node_gap, open, true);

    //return std::make_pair(nullptr, nullptr);
    return create_cl_ml_jobs(node_gap, root, branching_pair, nullptr, shared_data);
}

bool is_thread_pool_working(std::vector<bool> &thread_state) {
    int count = 0;
    for (auto && i : thread_state) {
        if (i)
            count++;
    }
    if (count == 0)
        return false;
    return true;
}

void save_X_to_file(arma::sp_mat &X){

    std::ofstream f;
    f.open(result_path);
    for (size_t i = 0; i < X.n_rows; i++){
        int val = X(i,0);
        f << val;
        for (size_t j = 1; j < X.n_cols; j++){
            val = X(i,j);
            f << " " << val;
        }
        f << "\n";
    }
    f.close();
}


arma::sp_mat sdp_branch_and_bound(arma::mat &Ws, arma::mat &C) {

    int n_thread = branch_and_bound_parallel;

    JobAbstractQueue *queue;
    switch (branch_and_bound_visiting_strategy) {
        case DEPTH_FIRST:
            queue = new JobStack();
            break;
        case BEST_FIRST:
            queue = new JobPriorityQueue();
            break;
        case BREADTH_FIRST:
            queue = new JobQueue();
            break;
        default:
            queue = new JobPriorityQueue();
    }

    auto *shared_data = new SharedData();
    shared_data->global_ub = std::numeric_limits<double>::infinity();
    shared_data->n_nodes = 0;
    shared_data->queue = queue;

    shared_data->threadStates.reserve(n_thread);
    for (int i = 0; i < n_thread; i++) {
        shared_data->threadStates.push_back(false);
    }

    auto *input_data = new InputData();
    input_data->Ws = Ws; // data matrix
    input_data->C = C; // cluster sizes

    ThreadPool pool(shared_data, input_data, n_thread);
    
    print_header_sdp(log_file);

    auto start_all = std::chrono::high_resolution_clock::now();
    
    auto *matlab_struct = new MatlabStruct();
    matlab_struct->matlabPtr = start_matlab(sdp_solver_folder, gurobi_folder);

    std::pair<JobData *, JobData *> jobs = build_root_problem(matlab_struct, input_data, shared_data);

    delete (matlab_struct);
    
    double root_gap = shared_data->gap;

    JobData *cl_job = jobs.first;
    JobData *ml_job = jobs.second;
    if (cl_job != nullptr && ml_job != nullptr) {
        pool.addJob(cl_job);
        pool.addJob(ml_job);
    }

    while (true) {

        {
            std::unique_lock<std::mutex> l(shared_data->queueMutex);
            while (is_thread_pool_working(shared_data->threadStates) && shared_data->n_nodes < branch_and_bound_max_nodes) {
                shared_data->mainConditionVariable.wait(l);
            }

            if (shared_data->queue->empty() || shared_data->n_nodes >= branch_and_bound_max_nodes)
                break;
        }

    }

    auto end_all = std::chrono::high_resolution_clock::now();
    auto duration_all = std::chrono::duration_cast<std::chrono::seconds>(end_all - start_all);

    pool.quitPool();

    if (queue->empty())
        shared_data->gap = 0.0;

    log_file << "\n";
    log_file << "TIME: " << duration_all.count() << " sec\n";
    log_file << "NODES: " << shared_data->n_nodes << "\n";
    log_file << "ROOT_GAP: " << std::max(0.0, root_gap) << "\n";
    log_file << "GAP: " << std::max(0.0, shared_data->gap) << "\n";
    log_file << "OPT: " << shared_data->global_ub << "\n\n";

    arma::sp_mat result = shared_data->global_X;
    save_X_to_file(result);

    // free memory

    delete (input_data);
    delete (queue);
    delete (shared_data);

    return result;

}
