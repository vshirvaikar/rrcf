/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (rrcf).

  rrcf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  rrcf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with rrcf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include "InstrumentalGLM.h"

namespace rrcf {

InstrumentalGLM::InstrumentalGLM(size_t dummy){ // constructor
    this->counter = new size_t[dummy];
}

InstrumentalGLM::~InstrumentalGLM(){ // destructor
    if (counter != nullptr) {
        delete[] counter;
    }
}

double InstrumentalGLM::dummy(){
    double output = (double) rand() / RAND_MAX;
    return output;
}

Eigen::VectorXd InstrumentalGLM::cwiseExp(Eigen::VectorXd input){
    return input.array().exp().matrix();
}

Eigen::VectorXd InstrumentalGLM::variance(std::string family, Eigen::VectorXd mu){
    if(family == "logistic"){
        return mu.cwiseProduct((1 - mu.array()).matrix());
    } else if(family == "poisson"){
        return mu;
    }
    return Eigen::VectorXd::Ones(mu.size());
}

Eigen::VectorXd InstrumentalGLM::invlink(std::string family, Eigen::VectorXd eta){
    if(family == "logistic"){
        Eigen::VectorXd denom = (eta.array().exp() + 1).matrix();
        return cwiseExp(eta).cwiseQuotient(denom);
    } else if(family == "poisson"){
        return cwiseExp(eta);
    }
    return eta;
}

Eigen::VectorXd InstrumentalGLM::invlink_prime(std::string family, Eigen::VectorXd eta){
    if(family == "logistic"){
        Eigen::VectorXd denom = (eta.array().exp() + 1).pow(2).matrix();
        return cwiseExp(eta).cwiseQuotient(denom);
    } else if(family == "poisson"){
        return cwiseExp(eta);
    }
    return Eigen::VectorXd::Ones(eta.size());
}

double InstrumentalGLM::glm_fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y,
                                std::string family, size_t maxit, double tol) {

    int n_cols = X.cols();
    int n_rows = X.rows();
    Eigen::VectorXd s = Eigen::VectorXd::Zero(n_cols);
    Eigen::VectorXd s_old;
    Eigen::VectorXd eta = Eigen::VectorXd::Ones(n_rows);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(X);
    Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n_rows, n_cols);
    Eigen::MatrixXd fullR = qr.matrixQR().triangularView<Eigen::Upper>();
    Eigen::MatrixXd R = Eigen::MatrixXd::Identity(n_cols, n_rows) * fullR;
    Eigen::MatrixXd V;
    bool converged = false;

    for (size_t i = 0; i < maxit; i++) {
        s_old = s;
        Eigen::VectorXd mu = invlink(family, eta);
        Eigen::VectorXd mu_p = invlink_prime(family, eta);
        Eigen::VectorXd z = eta + (y - mu).cwiseQuotient(mu_p);
        Eigen::VectorXd W = mu_p.array().pow(2).matrix().cwiseQuotient(variance(family, mu));
        Eigen::MatrixXd WQ = (Q.array().colwise() * W.array()).matrix();
        Eigen::LLT<Eigen::MatrixXd> cholesky(Q.transpose() * WQ);
        Eigen::VectorXd s1 = cholesky.matrixL().solve(Q.transpose() * W.cwiseProduct(z));
        s = cholesky.matrixU().solve(s1);
        eta = Q * s;

        double error = (s - s_old).norm();
        if (std::isnan(error)) {
            break;
        }
        if (error < tol) {
            V = (X.transpose() * W.asDiagonal() * X).inverse();
            converged = true;
            break;
        }
    }

    if(!converged){
        return 0;
    }
    Eigen::VectorXd coeffs = R.householderQr().solve(Q.transpose() * eta);
    Eigen::VectorXd stderrs = V.diagonal().cwiseSqrt();
    Eigen::VectorXd tstats = coeffs.cwiseQuotient(stderrs);
    return abs(tstats(tstats.rows() - 1));
}

} // namespace rrcf
