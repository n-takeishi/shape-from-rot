#include "PropagateFactor.h"

namespace gtsam
{

/*******************************************************************/
Vector PropagateFactor::evaluateError(
  const Vector3 &alpha1, const Vector3 &alpha2,
  const Vector3 &omega1, const Vector3 &omega2,
  const Vector2 &k,
  boost::optional< Matrix & > H1,
  boost::optional< Matrix & > H2,
  boost::optional< Matrix & > H3,
  boost::optional< Matrix & > H4,
  boost::optional< Matrix & > H5
  ) const
{

  // compute errors
  Vector err(6);
  err = computeError(tf_, alpha1, alpha2, omega1, omega2, k);

  // change noise covariance
  // noiseModel_->setCovariance(...);

  // numerical differentitations
  if (H1) {
    // d(err)/d(alpha1)
    Matrix jac(6,3);
    Functor_alpha1 fun(tf_, alpha1, alpha2, omega1, omega2, k);
    Eigen::NumericalDiff<Functor_alpha1> numDiff(fun);
    numDiff.df(alpha1, jac);
    *H1 = jac;
  }

  if (H2) {
    // d(err)/d(alpha2)
    *H2 = zeros(6,3);
    H2->block<3,3>(0,0).setIdentity();
  }

  if (H3) {
    // d(err)/d(omega1)
    Matrix jac(6,3);
    Functor_omega1 fun(tf_, alpha1, alpha2, omega1, omega2, k);
    Eigen::NumericalDiff<Functor_omega1> numDiff(fun);
    numDiff.df(omega1, jac);
    *H3 = jac;
  }

  if (H4) {
    // d(err)/d(omega2)
    *H4 = zeros(6,3);
    H4->block<3,3>(3,0).setIdentity();
  }

  if (H5) {
    // d(err)/d(k)
    Matrix jac(6,2);
    Functor_k fun(tf_, alpha1, alpha2, omega1, omega2, k);
    Eigen::NumericalDiff<Functor_k> numDiff(fun);
    numDiff.df(k, jac);
    *H5 = jac;
  }

  return err;
}

/*******************************************************************/
Vector6 PropagateFactor::computeError(
    double tf,
    const Vector3 &alpha1, const Vector3 &alpha2,
    const Vector3 &omega1, const Vector3 &omega2,
    const Vector2 &k)
{
  Vector6 alpha2_and_omega2;
  alpha2_and_omega2.segment<3>(0) = alpha2;
  alpha2_and_omega2.segment<3>(3) = omega2;
  Vector6 alpha2p_and_omega2p = propagate(tf, alpha1, omega1, k);

  return alpha2_and_omega2 - alpha2p_and_omega2p;
}

/*******************************************************************/
Vector6 PropagateFactor::propagate(
  double tf,
  const Vector3 &alpha1,
  const Vector3 &omega1,
  const Vector2 &k
  )
{
  Vector6 alpha_and_omega;
  alpha_and_omega.segment<3>(0) = alpha1;
  alpha_and_omega.segment<3>(3) = omega1;

  // Runge-Kutta
  double t = 0.0, dt = 0.0, h = tf/100.0;
  bool done = false;
  Vector6 RK_1, RK_2, RK_3, RK_4;
  while (!done) {
    if(tf-h-t > 0) {
      dt = h;
    } else {
      dt = tf-t;
      done = true;
    }
    RK_1 = dt * timeDerivative(alpha_and_omega, k);
    RK_2 = dt * timeDerivative(alpha_and_omega+0.5*RK_1, k);
    RK_3 = dt * timeDerivative(alpha_and_omega+0.5*RK_2, k);
    RK_4 = dt * timeDerivative(alpha_and_omega+RK_3, k);
    alpha_and_omega = alpha_and_omega + (RK_1 + 2.0*RK_2 + 2.0*RK_3 + RK_4)/6.0;
    t += dt;
  }

  return alpha_and_omega;
}

/*******************************************************************/

Vector6 PropagateFactor::timeDerivative(
  const Vector6 &alpha_and_omega,
  const Vector2 &k
  )
{
  Vector3 alpha = alpha_and_omega.segment<3>(0);
  Vector3 omega = alpha_and_omega.segment<3>(3);

  // inertia moment matrix
  Matrix3 J = Matrix3::Zero();
  Matrix3 invJ = Matrix3::Zero();
  J(0,0) = exp(k(0)); J(1,1) = 1.0; J(2,2) = exp(-k(1));
  invJ(0,0) = exp(-k(0)); invJ(1,1) = 1.0; invJ(2,2) = exp(k(1));

  // cross product matrix
  Matrix3 omegaCpMat = skewSymmetric(omega(0), omega(1), omega(2));

  // coefficients
  double theta = alpha.norm();
  double gamma, eta;
  if (theta < 1e-8) {
      gamma = (12.0 - theta*theta) / 6.0;
      eta = omega.dot(alpha) * (60.0 + theta*theta) / 360.0;
  } else {
      double cot_half_theta = cos(0.5*theta) / sin(0.5*theta); // cot=cos/sin
      gamma = theta * cot_half_theta;
      eta = omega.dot(alpha) / theta * (cot_half_theta - 2.0/theta);
  }

  // derivatives
  Vector3 d_alpha = 0.5 * (gamma*omega + omegaCpMat*alpha - eta*alpha);
  Vector3 d_omega = -invJ * omegaCpMat * J * omega;

  Vector6 d_alpha_and_omega;
  d_alpha_and_omega.segment<3>(0) = d_alpha;
  d_alpha_and_omega.segment<3>(3) = d_omega;
  return d_alpha_and_omega;
}

} // namespace gtsam