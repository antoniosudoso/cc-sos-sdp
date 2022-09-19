#ifndef CLUSTERING_CONFIG_PARAMS_H
#define CLUSTERING_CONFIG_PARAMS_H

#define CANNOT_LINK 0
#define MUST_LINK 1

#define VECTOR_LIFTING 0
#define MATRIX_LIFTING 1

#define BEST_FIRST 0
#define DEPTH_FIRST 1
#define BREADTH_FIRST 2

// cp_flag values
#define CP_FLAG_WORST -3
#define CP_FLAG_NO_SUCCESS -2
#define CP_FLAG_MAX_ITER -1
#define CP_FLAG_NO_VIOL 0
#define CP_FLAG_MAX_INEQ 1
#define CP_FLAG_PRUNING 2
#define CP_FLAG_CP_TOL 3
#define CP_FLAG_INFEAS 4
#define CP_FLAG_SDP_INFEAS 5

// data full path
extern const char *data_path;
extern const char *log_path;
extern const char *result_path;
extern std::ofstream log_file;

// branch and bound
extern double branch_and_bound_tol;
extern int branch_and_bound_parallel;
extern int branch_and_bound_max_nodes;
extern int branch_and_bound_visiting_strategy;

// matlab
extern int matlab_session_threads_root;
extern int matlab_session_threads_child;

// sdpnal
extern const char *sdp_solver_folder;
extern int sdp_relaxation;
extern double sdp_solver_tol;
extern int sdp_solver_max_iter;
extern int sdp_solver_max_time;
extern int sdp_solver_verbose;

// cutting plane
extern int cp_max_iter_root;
extern int cp_max_iter_child;
extern double cp_tol;
extern int cp_max_ineq;
extern double cp_perc_ineq;
extern int cp_inherit_ineq;
extern double cp_eps_ineq;
extern double cp_eps_active;

// heuristic
extern const char *gurobi_folder;
extern int heuristic_verbose;

#endif //CLUSTERING_CONFIG_PARAMS_H
