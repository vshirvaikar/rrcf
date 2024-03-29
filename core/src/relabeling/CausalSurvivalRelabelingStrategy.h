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

#ifndef GRF_CAUSALSURVIVALRELABELINGSTRATEGY_H
#define GRF_CAUSALSURVIVALRELABELINGSTRATEGY_H

#include <vector>

#include "commons/Data.h"
#include "relabeling/RelabelingStrategy.h"

namespace rrcf {

class CausalSurvivalRelabelingStrategy final: public RelabelingStrategy {
public:
  bool relabel(
      const std::vector<size_t>& samples,
      const Data& data,
      Eigen::ArrayXXd& responses_by_sample) const;

};

} // namespace rrcf

#endif //GRF_CAUSALSURVIVALRELABELINGSTRATEGY_H
