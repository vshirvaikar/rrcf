/*-------------------------------------------------------------------------------
  This file is part of generalized-random-forest.

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

#ifndef GRF_INSTRUMENTALGLM_H
#define GRF_INSTRUMENTALGLM_H

#include "commons/Data.h"
#include "Eigen/Dense"

namespace rrcf {

    class InstrumentalGLM final{
    public:

        InstrumentalGLM(size_t dummy);
        ~InstrumentalGLM();

        double dummy();

    double calculate_p_value(double z);

	Eigen::VectorXd cwiseExp(Eigen::VectorXd input);

	Eigen::VectorXd variance(std::string family, Eigen::VectorXd mu);

	Eigen::VectorXd invlink(std::string family, Eigen::VectorXd mu);

	Eigen::VectorXd invlink_prime(std::string family, Eigen::VectorXd eta);

	double glm_fit(const Eigen::MatrixXd& X, const Eigen::VectorXd& y,
		       std::string family, size_t maxit, double tol);

    private:

        size_t* counter;

        DISALLOW_COPY_AND_ASSIGN(InstrumentalGLM);
    };

} // namespace rrcf

#endif //GRF_INSTRUMENTALGLM_H
