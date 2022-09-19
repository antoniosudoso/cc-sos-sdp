#include <iostream>
#include <map>
#include <algorithm>
#include <armadillo>
#include "sdp_branch_and_bound.h"
#include "config_params.h"

// data full path
const char *data_path;
const char *log_path;
const char *result_path;
std::ofstream log_file;

// branch and bound
double branch_and_bound_tol;
int branch_and_bound_parallel;
int branch_and_bound_max_nodes;
int branch_and_bound_visiting_strategy;

// matlab
int matlab_session_threads_root;
int matlab_session_threads_child;

// sdpnal
const char *sdp_solver_folder;

double sdp_solver_tol;
int sdp_relaxation;
int sdp_solver_max_iter;
int sdp_solver_max_time;
int sdp_solver_verbose;

// cutting plane
int cp_max_iter_root;
int cp_max_iter_child;
double cp_tol;
int cp_max_ineq;
double cp_perc_ineq;
int cp_inherit_ineq;
double cp_eps_ineq;
double cp_eps_active;

// heuristic
const char *gurobi_folder;
int heuristic_verbose;


std::map<std::string, std::string> read_params(std::string &config_file) {

    std::map<std::string, std::string> config_map = {};

    std::ifstream cFile (config_file);
    if (cFile.is_open()) {
        std::string line;
        while (getline(cFile, line)){
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find('=');
            auto key = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);
            config_map.insert(std::pair<std::string, std::string>(key, value));
        }

    }
    else {
        std::cerr << "Couldn't open config file for reading.\n";
    }

    return config_map;
}

// read data Ws
arma::mat read_data(const char *filename, int &n, int &d) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << strerror(errno) << "\n";
        exit(EXIT_FAILURE);
    }
    // read the header n, d
    file >> n >> d;
    arma::mat Ws(n, d);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < d; j++) {
            file >> Ws(i, j);
        }
    }
    return Ws;
}


void run(int argc, char **argv) {

    std::string config_file = "config.txt";
    std::map<std::string, std::string> config_map = read_params(config_file);

    if (argc < 4) {
        std::cerr << "Input: <DATA_PATH> <K> <SIZE_1> ... <SIZE_K> <LOG_PATH> <RESULT_PATH>" << std::endl;
        exit(EXIT_FAILURE);
    }

    data_path = argv[1];
    int n, d, k;
    arma::mat Ws = read_data(data_path, n, d);
    k = std::stoi(argv[2]);
    arma::vec cluster_sizes(k);
    for (int i=0; i < k; i++) {
        cluster_sizes(i) = std::stoi(argv[i+3]);
    }
    log_path = argv[3+k];
    log_file.open(log_path);
    result_path = argv[4+k];

    // branch and bound
    branch_and_bound_tol = std::stod(config_map["BRANCH_AND_BOUND_TOL"]);
    branch_and_bound_parallel = std::stoi(config_map["BRANCH_AND_BOUND_PARALLEL"]);
    branch_and_bound_max_nodes = std::stoi(config_map["BRANCH_AND_BOUND_MAX_NODES"]);
    branch_and_bound_visiting_strategy = std::stoi(config_map["BRANCH_AND_BOUND_VISITING_STRATEGY"]);

    // matlab
    matlab_session_threads_root = std::stoi(config_map["MATLAB_SESSION_THREADS_ROOT"]);
    matlab_session_threads_child = std::stoi(config_map["MATLAB_SESSION_THREADS_CHILD"]);

    // sdp solver
    sdp_solver_folder = config_map["SDP_SOLVER_FOLDER"].c_str();
    sdp_relaxation = std::stoi(config_map["SDP_RELAXATION"]);
    sdp_solver_tol = std::stod(config_map["SDP_SOLVER_TOL"]);
    sdp_solver_max_iter = std::stoi(config_map["SDP_SOLVER_MAX_ITER"]);
    sdp_solver_max_time = std::stoi(config_map["SDP_SOLVER_MAX_TIME"]);
    sdp_solver_verbose = std::stoi(config_map["SDP_SOLVER_VERBOSE"]);

    // cutting plane
    cp_max_iter_root = std::stoi(config_map["CP_MAX_ITER_ROOT"]);
    cp_max_iter_child = std::stoi(config_map["CP_MAX_ITER_CHILD"]);
    cp_tol = std::stod(config_map["CP_TOL"]);
    cp_max_ineq = std::stoi(config_map["CP_MAX_INEQ"]);
    cp_perc_ineq = std::stod(config_map["CP_PERC_INEQ"]);
    cp_inherit_ineq = std::stoi(config_map["CP_INHERIT_INEQ"]);
    cp_eps_ineq = std::stod(config_map["CP_EPS_INEQ"]);
    cp_eps_active = std::stod(config_map["CP_EPS_ACTIVE"]);

    // heuristic
    gurobi_folder = config_map["GUROBI_FOLDER"].c_str();
    heuristic_verbose = std::stoi(config_map["HEURISTIC_VERBOSE"]);

    log_file << "\n" << "DATA_PATH: " << data_path << " " << n << " " << d << " " << k <<"\n";
    log_file << "CLUSTER_SIZES: ";
    for (int i=0; i < k; i++)
        log_file << cluster_sizes(i) << " ";
    log_file << "\n";
    log_file << "LOG_PATH: " << log_path << "\n\n";

    log_file << "BRANCH_AND_BOUND_TOL: " << branch_and_bound_tol << "\n";
    log_file << "BRANCH_AND_BOUND_PARALLEL: " << branch_and_bound_parallel << "\n";
    log_file << "BRANCH_AND_BOUND_MAX_NODES: " <<  branch_and_bound_max_nodes << "\n";
    log_file << "BRANCH_AND_BOUND_VISITING_STRATEGY: " << branch_and_bound_visiting_strategy << "\n\n";

    log_file << "MATLAB_SESSION_THREADS_ROOT: " << matlab_session_threads_root << "\n";
    log_file << "MATLAB_SESSION_THREADS_CHILD: " << matlab_session_threads_child << "\n\n";

    log_file << "SDP_SOLVER_FOLDER: " << sdp_solver_folder << "\n";
    log_file << "SDP_RELAXATION: ";
    switch (sdp_relaxation) {
        case VECTOR_LIFTING:
            log_file << sdp_relaxation << " (VECTOR LIFTING)" << "\n";
            break;
        case MATRIX_LIFTING:
            log_file << sdp_relaxation << " (MATRIX LIFTING)" << "\n";
            break;
        default:
            log_file << "UNKNOWN RELAXATION" << "\n";
            exit(EXIT_FAILURE);
    }
    log_file << "SDP_SOLVER_TOL: " << sdp_solver_tol << "\n";
    log_file << "SDP_SOLVER_MAX_ITER: " << sdp_solver_max_iter << "\n";
    log_file << "SDP_SOLVER_MAX_TIME: " << sdp_solver_max_time << "\n";
    log_file << "SDP_SOLVER_VERBOSE: " << sdp_solver_verbose << "\n\n";

    log_file << "CP_MAX_ITER_ROOT: " << cp_max_iter_root << "\n";
    log_file << "CP_MAX_ITER_CHILD: " << cp_max_iter_child << "\n";
    log_file << "CP_TOL: " << cp_tol << "\n";
    log_file << "CP_MAX_INEQ: " << cp_max_ineq << "\n";
    log_file << "CP_PERC_INEQ: " << cp_perc_ineq << "\n";
    log_file << "CP_INHERIT_INEQ: " << cp_inherit_ineq << "\n";
    log_file << "CP_EPS_INEQ: " << cp_eps_ineq << "\n";
    log_file << "CP_EPS_ACTIVE: " << cp_eps_active << "\n\n";

    log_file << "GUROBI_FOLDER: " << gurobi_folder << "\n";
    log_file << "HEURISTIC_VERBOSE: " << heuristic_verbose << "\n\n";

    arma::mat C = arma::diagmat(cluster_sizes);
    sdp_branch_and_bound(Ws, C);

}

int main(int argc, char **argv) {

    run(argc, argv);

    return EXIT_SUCCESS;
}
