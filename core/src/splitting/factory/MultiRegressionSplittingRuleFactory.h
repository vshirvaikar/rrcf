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

#ifndef GRF_MULTIREGRESSIONSPLITTINGRULEFACTORY_H
#define GRF_MULTIREGRESSIONSPLITTINGRULEFACTORY_H


#include "splitting/factory/SplittingRuleFactory.h"

namespace grf {

/**
 * A factory that produces standard regression splitting rules.
 *
 * In addition to performing standard regression splits, this rule applies
 * a penalty to avoid splits too close to the edge of the node's data.
 */
class MultiRegressionSplittingRuleFactory final: public SplittingRuleFactory {
public:
  MultiRegressionSplittingRuleFactory(size_t num_outcomes);

  std::unique_ptr<SplittingRule> create(size_t max_num_unique_values,
                                        const TreeOptions& options) const;
private:
  size_t num_outcomes;
  
  DISALLOW_COPY_AND_ASSIGN(MultiRegressionSplittingRuleFactory);
};

} // namespace rrcf

#endif //GRF_MULTIREGRESSIONSPLITTINGRULEFACTORY_H
