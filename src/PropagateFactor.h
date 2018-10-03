#pragma once

//#include "VariableCovGaussian.h"

#include <gtsam/nonlinear/NonlinearFactor.h>

#include <Eigen/Core>
#include <unsupported/Eigen/NumericalDiff>

namespace gtsam
{

/*******************************************************************/
class PropagateFactor : public NoiseModelFactor5<Vector3, Vector3, Vector3, Vector3, Vector2>
{

protected:
  double tf_;

public:
  typedef boost::shared_ptr<PropagateFactor> shared_ptr;

  PropagateFactor(
    double tf,
    Key j1, Key j2, Key j3, Key j4, Key j5,
    const SharedNoiseModel &noiseModel
    ) : NoiseModelFactor5(noiseModel, j1, j2, j3, j4, j5), tf_(tf)
  {};

  virtual ~PropagateFactor() {}

  Vector evaluateError(
    const Vector3 &alpha1, const Vector3 &alpha2,
    const Vector3 &omega1, const Vector3 &omega2,
    const Vector2 &k,
    boost::optional< Matrix & > H1 = boost::none,
    boost::optional< Matrix & > H2 = boost::none,
    boost::optional< Matrix & > H3 = boost::none,
    boost::optional< Matrix & > H4 = boost::none,
    boost::optional< Matrix & > H5 = boost::none
    ) const;

  static Vector6 computeError(
    double tf,
    const Vector3 &alpha1, const Vector3 &alpha2,
    const Vector3 &omega1, const Vector3 &omega2,
    const Vector2 &k
    );

  static Vector6 propagate(
    double tf,
    const Vector3 &alpha1,
    const Vector3 &omega1,
    const Vector2 &k
    );

  static Vector6 timeDerivative(
    const Vector6 &alpha_and_omega,
    const Vector2 &k
    );

};

/*******************************************************************/
template<typename _Scalar, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
struct BaseFunctor
{
  typedef _Scalar Scalar;
  enum {
      InputsAtCompileTime = NX,
      ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;
  const double tf_base;
  const Vector3 alpha1_base;
  const Vector3 alpha2_base;
  const Vector3 omega1_base;
  const Vector3 omega2_base;
  const Vector2 k_base;

  BaseFunctor(
    int inputs,
    double tf,
    const Vector3 &alpha1, const Vector3 &alpha2,
    const Vector3 &omega1, const Vector3 &omega2,
    const Vector2 &k
    ) : m_inputs(inputs), m_values(6),
        tf_base(tf),
        alpha1_base(alpha1), alpha2_base(alpha2),
        omega1_base(omega1), omega2_base(omega2),
        k_base(k) {}
  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
};

/*******************************************************************/
struct Functor_alpha1 : BaseFunctor<double>
{
  Functor_alpha1(
    double tf,
    const Vector3 &alpha1, const Vector3 &alpha2,
    const Vector3 &omega1, const Vector3 &omega2,
    const Vector2 &k
    ) : BaseFunctor(3, tf, alpha1, alpha2, omega1, omega2, k) {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int operator()(const Vector &alpha1_now, Vector &fvec) const {
    fvec = PropagateFactor::computeError(tf_base, alpha1_now, alpha2_base, omega1_base, omega2_base, k_base);
    return 0;
  }
};

/*******************************************************************/
struct Functor_omega1 : BaseFunctor<double>
{
  Functor_omega1(
    double tf,
    const Vector3 &alpha1, const Vector3 &alpha2,
    const Vector3 &omega1, const Vector3 &omega2,
    const Vector2 &k
    ) : BaseFunctor(3, tf, alpha1, alpha2, omega1, omega2, k) {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int operator()(const Vector &omega1_now, Vector &fvec) const {
    fvec = PropagateFactor::computeError(tf_base, alpha1_base, alpha2_base, omega1_now, omega2_base, k_base);
    return 0;
  }
};

/*******************************************************************/
struct Functor_k : BaseFunctor<double>
{
  Functor_k(
    double tf,
    const Vector3 &alpha1, const Vector3 &alpha2,
    const Vector3 &omega1, const Vector3 &omega2,
    const Vector2 &k
    ) : BaseFunctor(3, tf, alpha1, alpha2, omega1, omega2, k) {}

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int operator()(const Vector &k_now, Vector &fvec) const {
    fvec = PropagateFactor::computeError(tf_base, alpha1_base, alpha2_base, omega1_base, omega2_base, k_now);
    return 0;
  }
};

/*******************************************************************/
} // namespace gtsam
