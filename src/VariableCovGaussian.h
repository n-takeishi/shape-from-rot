#pragma once

#include <gtsam/linear/NoiseModel.h>

namespace gtsam
{

namespace noiseModel
{

class GTSAM_EXPORT VariableCovGaussian: public Gaussian {
  void setCovariance(const Matrix& covariance) {
    size_t m = covariance.rows(), n = covariance.cols();
    if (m != n)
        throw std::invalid_argument("VariableCovGaussian::setCovariance: covariance not square");
    sqrt_information_ = inverse_square_root(covariance);
  };
};

} // namespace NoiseModel

} // namespace gtsam
