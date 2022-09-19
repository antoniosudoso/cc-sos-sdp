#include "matlab_util.h"

/*
// connect to existing matlab session
std::unique_ptr<matlab::engine::MATLABEngine> init_matlab(std::string &session_name,
                                                          const char *sdp_solver_folder) {

    std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;
    std::cout << "Connecting Matlab..." << "\n";
    // Connect MATLAB engine synchronously
    std::u16string session(session_name.begin(), session_name.end());
    matlabPtr = matlab::engine::connectMATLAB(session);
    // Set path
    std::string folder(sdp_solver_folder);
    std::string matlab_path = std::string("addpath(genpath('" + folder + "'));");
    std::u16string path_command(matlab_path.begin(), matlab_path.end());
    matlabPtr->eval(path_command);
    std::cout << "Matlab Connected!" << "\n";
    return matlabPtr;
}
*/

std::u16string generate_path_command(const char *folder_char) {
    std::string folder(folder_char);
    std::string matlab_path = std::string("addpath(genpath('" + folder + "'));");
    std::u16string path_command(matlab_path.begin(), matlab_path.end());
    return path_command;
}

// start new matlab session
std::unique_ptr<matlab::engine::MATLABEngine> start_matlab(const char *sdp_solver_folder,
                                                           const char *gurobi_folder) {

    // Start MATLAB with -r option
    std::unique_ptr<matlab::engine::MATLABEngine> matlabPtr;
    //std::cout << "Starting Matlab..." << "\n";
    std::vector<std::u16string> optionVec;
    optionVec.emplace_back(u"-nojvm");
    matlabPtr = matlab::engine::startMATLAB(optionVec);
    // Set SDPNAL+ path
    std::u16string sdp_path_command = generate_path_command(sdp_solver_folder);
    matlabPtr->eval(sdp_path_command);
    // Set gurobi path
    std::u16string gurobi_path_command = generate_path_command(gurobi_folder);
    matlabPtr->eval(gurobi_path_command);
    return matlabPtr;
}

/*
matlab::data::TypedArray<double> arma_to_matlab_vector(matlab::data::ArrayFactory &factory, arma::vec &v) {

    std::vector<double> vec_v = arma::conv_to<std::vector<double>>::from(v);
    auto matlab_vector = factory.createArray({v.size()}, vec_v.begin(), vec_v.end());
    return matlab_vector;

}
*/

matlab::data::TypedArray<double> arma_to_matlab_matrix(matlab::data::ArrayFactory &factory, arma::mat &X) {

    std::vector<double> vec_X = arma::conv_to<std::vector<double>>::from(arma::vectorise(X));
    auto matlab_matrix = factory.createArray({X.n_rows, X.n_cols}, vec_X.begin(), vec_X.end());
    return matlab_matrix;

}

matlab::data::SparseArray<double> arma_to_matlab_sparse(matlab::data::ArrayFactory &factory, arma::sp_mat &X) {

    arma::sp_mat::const_iterator start = X.begin();
    arma::sp_mat::const_iterator end = X.end();

    size_t nnz = X.n_nonzero;

    std::vector<double> data;
    data.reserve(nnz);
    std::vector<size_t> rows;
    rows.reserve(nnz);
    std::vector<size_t> cols;
    cols.reserve(nnz);

    for (arma::sp_mat::const_iterator it = start; it != end; ++it) {
        data.push_back(*it);
        rows.push_back(it.row());
        cols.push_back(it.col());
    }

    auto data_p = factory.createBuffer<double>(nnz);
    auto rows_p = factory.createBuffer<size_t>(nnz);
    auto cols_p = factory.createBuffer<size_t>(nnz);

    double* dataPtr = data_p.get();
    size_t* rowsPtr = rows_p.get();
    size_t* colsPtr = cols_p.get();

    std::for_each(data.begin(), data.end(), [&](const double& e) { *(dataPtr++) = e; });
    std::for_each(rows.begin(), rows.end(), [&](const size_t& e) { *(rowsPtr++) = e; });
    std::for_each(cols.begin(), cols.end(), [&](const size_t& e) { *(colsPtr++) = e; });

    auto sparse_matrix = factory.createSparseArray<double>({X.n_rows, X.n_cols}, nnz,
                                                           std::move(data_p), std::move(rows_p), std::move(cols_p));

    return sparse_matrix;

}

arma::sp_mat matlab_to_arma_sparse(matlab::data::SparseArray<double> &X_matlab) {

    const size_t n = X_matlab.getDimensions()[0];
    const size_t k = X_matlab.getDimensions()[1];
    arma::sp_mat X(n, k);

    for (auto it = X_matlab.begin(); it != X_matlab.end(); it++) {
        std::pair<size_t, size_t> pair = X_matlab.getIndex(it);
        X(pair.first, pair.second) = *it;
    }

    return X;
}

/*
// A is m x n * n
matlab::data::CellArray arma_to_matlab_cell(matlab::data::ArrayFactory &factory, arma::sp_mat &A) {

    size_t n = std::sqrt(A.n_cols);
    size_t m = A.n_rows;
    matlab::data::CellArray Acell = factory.createCellArray({1, m});
    for (size_t i = 0; i < m; i++) {
        arma::sp_mat A_i = arma::reshape(A.row(i), n, n);
        Acell[i] = arma_to_matlab_sparse(factory, A_i);
    }

    return Acell;
}
*/

// B operator
matlab::data::CellArray arma_to_matlab_cell(matlab::data::ArrayFactory &factory, std::vector<arma::sp_mat> &B_vector) {

    size_t m = B_vector.size();
    matlab::data::CellArray Bcell = factory.createCellArray({1, m});
    for (size_t i = 0; i < m; i++) {
        Bcell[i] = arma_to_matlab_sparse(factory, B_vector[i]);
    }

    return Bcell;
}

matlab::data::CellArray arma_to_matlab_vector_cell(matlab::data::ArrayFactory &factory,
                                                   std::vector<std::vector<arma::sp_mat>> &B_vector) {

    size_t k = B_vector.size();
    matlab::data::CellArray Bcell = factory.createCellArray({1, k});
    for (size_t j = 0; j < k; j++) {
        Bcell[j] = arma_to_matlab_cell(factory, B_vector[j]);
    }

    return Bcell;
}


matlab::data::TypedArray<double> vector_pair_to_matlab_matrix(matlab::data::ArrayFactory &factory,
                                                              std::vector<std::pair<int, int>> &pairs) {
    std::size_t n_constraints = pairs.size();
    auto M = factory.createArray<double>({n_constraints, 2});
    for (std::size_t c = 0; c < n_constraints; c++) {
        M[c][0] = pairs[c].first;
        M[c][1] = pairs[c].second;
    }
    return M;
}

/*
arma::mat matlab_to_arma_matrix(matlab::data::TypedArray<double> &X_matlab) {
    const size_t n = X_matlab.getDimensions()[0];
    const size_t k = X_matlab.getDimensions()[1];
    arma::mat X(n, k);
    for (size_t i = 0; i < n ; i++) {
        for (size_t j = 0; j < k ; j++) {
            X(i, j) = X_matlab[i][j];
        }
    }
    return X;
}
*/

/*
arma::vec matlab_to_arma_vector(matlab::data::TypedArray<double> &v_matlab) {
    const size_t n = v_matlab.getDimensions()[0];
    arma::vec v(n);
    for (size_t i = 0; i < n; i++) {
        v(i) = v_matlab[i];
    }
    return v;
}
*/

// B_matlab is a cell of m sparse matrices i.e. {sparse_1, ..., sparse_m}
std::vector<arma::sp_mat> matlab_to_arma_vector_sp_mat(matlab::data::CellArray &B_matlab) {

    std::vector<arma::sp_mat> B_vector;
    const size_t m = B_matlab.getDimensions()[1];
    B_vector.reserve(m);

    for (size_t i = 0; i < m; i++) {
        matlab::data::SparseArray<double> B_i = B_matlab[i];
        arma::sp_mat B_sparse = matlab_to_arma_sparse(B_i);
        B_vector.push_back(B_sparse);
    }

    return B_vector;
}

// B_matlab is a cell of k cell arrays i.e. {{1xm_1 cell}, {1xm_2 cell}, ... {1xm_k cell}}
// each cell {1xm_j} contains m_j sparse matrices
std::vector<std::vector<arma::sp_mat>> matlab_to_arma_B_cell(matlab::data::CellArray &B_matlab) {

    const size_t k = B_matlab.getDimensions()[1];
    std::vector<std::vector<arma::sp_mat>> B_outer_vector;
    B_outer_vector.reserve(k);

    for (size_t j = 0; j < k; j++) {
        matlab::data::CellArray B_j = B_matlab[j];
        std::vector<arma::sp_mat> B_inner_vector_j = matlab_to_arma_vector_sp_mat(B_j);
        B_outer_vector.push_back(B_inner_vector_j);
    }
    return B_outer_vector;
}


