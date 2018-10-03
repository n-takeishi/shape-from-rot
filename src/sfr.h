/*
 * Shape From Rotation
 * (c) Naoya Takeishi, 2018.
 */

#pragma once

#define TINYSTD          1e-6
#define ALPHA0_PRIOR_STD 3.0
#define OMEGA0_PRIOR_STD 0.1
#define XT_PRIOR_STD     1.0
#define K_PRIOR_STD      1.0

#define L_ABSPRIOR_STD   1.0  // if not using L abs-prior, do not define this
//#define ZT_ABSPRIOR_STD  1.0  // if not using Zt abs-prior, do not define this
#define HUBERPARAM       1.0  // if not using Huber loss, do not define this
//#define RANSAC_THRESHOLD 0.1 // if not using RANSAC (instead use LMEDS), do not define this

// original
#include "coordtrans.h"
#include "PropagateFactor.h"
//#include "VariableCovGaussian.h"

// gtsam
#include <gtsam/geometry/Point2.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/geometry/Cal3_S2.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/linear/NoiseModel.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>
#include <gtsam/slam/expressions.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>

// gtsam optimizers
#include <gtsam/nonlinear/DoglegOptimizer.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>

// other libs
#include <opencv2/opencv.hpp>
#include <boost/shared_ptr.hpp>

// C standard
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

namespace gtsam
{

class SFR
{

private:

  // camera intrinsic matrix
  const Cal3_S2::shared_ptr IntrMat_;

  // observation noise models for Y, dZr, & dZt
  // maintained to be shared among multiple factors
  #ifdef HUBERPARAM
    noiseModel::Robust::shared_ptr Y_noise_;
  #else
    noiseModel::Isotropic::shared_ptr Y_noise_;
  #endif
  noiseModel::Diagonal::shared_ptr dZr_noise_;
  noiseModel::Diagonal::shared_ptr dZt_noise_;

  // graph and value container for optimizer
  ExpressionFactorGraph graph_;
  Values values_;

  // obsevations
  vector< vector<Point2> > Y_; // landmark observation, n_X x n_L : Y_[i][...] in <D_i>
  vector< vector<bool> > M_;   // landmark mask
  vector<Rot3> dZr_;           // dZr_[i] = inv(Zr[i-1]) * Zr[i]  : in <I>, i.e., Zr[i] in <Zr[i-1]>
  vector<Point3> dZt_;         // dZt_[i] = Zt[i] - Zt[i-1]       : in <I>, i.e., Zt[i] in <Zt[i-1]>
  vector<double> Ts_;          // timestamps

  // estimated unknown variables container
  vector<Rot3> Zr_;       // rotation of spacecraft (r) : in <I>
  vector<Point3> Zt_;     // position of spacecraft (t) : in <I>
  vector<Vector3> alpha_; // log(rotation) of asteroid in each frame (a) : in <I>
  vector<Vector3> omega_; // angular veloc of asteroid in each frame (o) : in <I>
  Point3 Xt_;             // position of asteroid in the first frame (s) : in <I>
  Vector2 k_;             // inertial ratio (k)
  vector<Point3> L_;      // landmark locations (l) : in <B>

  // landmark activation flags (true = already initialized)
  vector<bool> A_;

public:

  // default constructor
  SFR(
    const Cal3_S2::shared_ptr IntrMat,
    double Y_noise_std,
    Vector3 &dZr_noise_std,
    Vector3 &dZt_noise_std
    ) : IntrMat_(IntrMat)
  {
    noiseModel::Isotropic::shared_ptr Base = noiseModel::Isotropic::Sigma(2, Y_noise_std);
    #ifdef HUBERPARAM
        noiseModel::mEstimator::Huber::shared_ptr mEstimator
          = noiseModel::mEstimator::Huber::Create(HUBERPARAM);
        Y_noise_ = noiseModel::Robust::Create(mEstimator, Base);
    #else
        Y_noise_ = Base;
    #endif
    dZr_noise_ = noiseModel::Diagonal::Sigmas(dZr_noise_std);
    dZt_noise_ = noiseModel::Diagonal::Sigmas(dZt_noise_std);
  };

  // initialization
  void initialize(
    const Rot3 &Xr_init,
    const Point3 &Xt_init
    );

  // update
  void update(
    const vector< vector<Point2> > &Y_new,
    const vector< vector<bool> > &M_new,
    const vector<Rot3> &dZr_new,
    const vector<Point3> &dZt_new,
    const vector<double> &Ts_new
    );

  // triangulation
  void triangulation(int i1, int i2, const vector<int> &N_now);

  // store current estimation
  void get_estimated_values(int n_frames_get, const Values v);

  // remove landmark from graph, values, and data
  //void remove_landmark(int landmark_idx);

  // get current number of frames
  int n_frames() const { return static_cast<int>(Y_.size()); };

  // get current number of landmarks
  int n_landmarks() const { return static_cast<int>(A_.size()); };

  // get estimated values
  inline Rot3 Zr(int i) const { return Zr_[i]; };
  inline Point3 Zt(int i) const { return Zt_[i]; };
  inline Vector3 alpha(int i) const { return alpha_[i]; };
  inline Vector3 omega(int i) const { return omega_[i]; };
  inline Point3 Xt() const { return Xt_; };
  inline Vector2 k() const { return k_; };
  inline Point3 L(int j) const { return L_[j]; };

};

} // namespace gtsam