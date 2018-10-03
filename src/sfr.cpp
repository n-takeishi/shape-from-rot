/*
 * Shape From Rotation
 * (c) Naoya Takeishi, 2018.
 */

#include "sfr.h"

using namespace std;

namespace gtsam
{

/*******************************************************************/
void SFR::initialize(
  const Rot3 &Xr_init,
  const Point3 &Xt_init
  )
{
  // NOTE: why using Constrained noise makes optimization much slow???

  // initialization and prior of Zr[0] (r)
  Zr_.resize(1); Zr_[0] = Rot3(); // origin
  values_.insert(Symbol('r',0), Zr_[0]);
  graph_.add(PriorFactor<Rot3>(Symbol('r',0), Zr_[0],
    noiseModel::Isotropic::Sigma(3, TINYSTD)));

  // initialization and prior of Zt[0] (t)
  Zt_.resize(1); Zt_[0] = Point3(); // origin
  values_.insert(Symbol('t',0), Zt_[0]);
  graph_.add(PriorFactor<Point3>(Symbol('t',0), Zt_[0],
    noiseModel::Isotropic::Sigma(3, TINYSTD)));

  // initialization and prior of alpha[0] (a)
  alpha_.resize(1); alpha_[0] = Rot3::Logmap(Xr_init);
  values_.insert(Symbol('a',0), alpha_[0]);
  graph_.add(PriorFactor<Vector3>(Symbol('a',0), alpha_[0],
    noiseModel::Isotropic::Sigma(3, ALPHA0_PRIOR_STD)));

  // initialization and prior of omega[0] (o)
  omega_.resize(1); omega_[0] << 0.0, 0.0, 0.0;
  values_.insert(Symbol('o',0), omega_[0]);
  graph_.add(PriorFactor<Vector3>(Symbol('o',0), omega_[0],
    noiseModel::Isotropic::Sigma(3, OMEGA0_PRIOR_STD)));

  // initialization and prior of Xt (s)
  Xt_ = Point3(Xt_init.x(), Xt_init.y(), Xt_init.z());
  values_.insert(Symbol('s',0), Xt_);
  graph_.add(PriorFactor<Point3>(Symbol('s',0), Xt_,
    noiseModel::Diagonal::Sigmas((Vector3() << XT_PRIOR_STD,XT_PRIOR_STD,TINYSTD).finished())));

  // initialization and prior of k (k)
  k_ << 0.0, 0.0;
  values_.insert(Symbol('k',0), k_);
  graph_.add(PriorFactor<Vector2>(Symbol('k',0), k_,
    noiseModel::Isotropic::Sigma(2, K_PRIOR_STD)));
}

/*******************************************************************/
void SFR::update(
  const vector< vector<Point2> > &Y_new,
  const vector< vector<bool> > &M_new,
  const vector<Rot3> &dZr_new,
  const vector<Point3> &dZt_new,
  const vector<double> &Ts_new
  )
{
  bool firstupdate = n_frames() < 1 ? true : false;

  // check inputs size
  size_t n_newframes = Y_new.size();
  assert(n_frames() + n_newframes > 1);
  assert(M_new.size() == n_newframes);
  assert(dZr_new.size() == n_newframes);
  assert(dZt_new.size() == n_newframes);
  assert(Ts_new.size() == n_newframes);
  for (int i=0; i<n_newframes; ++i)
    assert(M_new[i].size() == Y_new[i].size());

  // check timestamps order
  if (!firstupdate) assert(Ts_[n_frames()-1] < Ts_new[0]);
  if (n_newframes > 1) {
    for (int i=1; i<n_newframes; ++i) assert(Ts_new[i-1] < Ts_new[i]);
  }

  if (firstupdate) {
    // at first update, the number of landmarks must be the same for all new frames
    int n_landmarks = Y_new[0].size();
    for (int i=1; i<n_newframes; ++i) assert(Y_new[i].size() == n_landmarks);
    // memory allocation and initialization of activation flags
    L_.resize(n_landmarks);
    A_.resize(n_landmarks);
    for (int j=0; j<n_landmarks; ++j) A_[j] = false;
  } else {
    for (int i=0; i<n_newframes; ++i)
      // TODO: currently, the total number of landmarks cannot be increased or decreased in updates.
      //       instead, make it possible to increase the number of landmarks.
      assert(Y_new[i].size() == n_landmarks());
  }

  // add new observations
  Y_.insert(end(Y_), begin(Y_new), end(Y_new));
  M_.insert(end(M_), begin(M_new), end(M_new));
  dZr_.insert(end(dZr_), begin(dZr_new), end(dZr_new));
  dZt_.insert(end(dZt_), begin(dZt_new), end(dZt_new));
  Ts_.insert(end(Ts_), begin(Ts_new), end(Ts_new));

  // allocate memory
  Zr_.resize(n_frames());
  Zt_.resize(n_frames());
  alpha_.resize(n_frames());
  omega_.resize(n_frames());

  // main loop
  for (int i=firstupdate?1:n_frames()-n_newframes; i<n_frames(); ++i) {
    cout << "Frame #" << i << " start :::" << endl;

    /*******************
     INITIALIZATION
    *******************/

    // initialize Zr[i] (r)
    Zr_[i] = Zr_[i-1] * dZr_[i];
    values_.insert(Symbol('r',i), Zr_[i]);

    // initialize Zt[i] (t)
    Zt_[i] = Zt_[i-1] + dZt_[i];
    values_.insert(Symbol('t',i), Zt_[i]);

    if (i==1) {
      // convert vector<Point2> to vector<cv::Point2d>
      vector<Point2> Y1=Y_[0], Y2=Y_[1];
      vector<bool> M1=M_[0], M2=M_[1];
      vector<cv::Point2d> Y1_cv, Y2_cv;
      for (int j=0; j<n_landmarks(); ++j) {
        if (M1[j] && M2[j]) {
          Y1_cv.push_back( cv::Point2d(Y1[j].x(), Y1[j].y()) );
          Y2_cv.push_back( cv::Point2d(Y2[j].x(), Y2[j].y()) );
        }
      }

      // 5pt-method
      double f = 0.5 * (IntrMat_->fx()+IntrMat_->fy()); // focal length for OpenCV
      cv::Point2d pp(IntrMat_->px(), IntrMat_->py());   // image center for OpenCV
      #ifdef RANSAC_THRESHOLD
        cv::Mat E12_cv = cv::findEssentialMat(Y1_cv, Y2_cv, f, pp, cv::RANSAC, 0.999, RANSAC_THRESHOLD);
      #else
        cv::Mat E12_cv = cv::findEssentialMat(Y1_cv, Y2_cv, f, pp, cv::LMEDS);
      #endif
      cv::Mat R12_cv, t12_cv;
      cv::recoverPose(E12_cv, Y1_cv, Y2_cv, R12_cv, t12_cv, f, pp);

      // relative pose between camera#2 and camera#1 in World Coordinate
      Matrix3 R12_mat;
      R12_mat << R12_cv.at<double>(0,0), R12_cv.at<double>(0,1), R12_cv.at<double>(0,2),
                 R12_cv.at<double>(1,0), R12_cv.at<double>(1,1), R12_cv.at<double>(1,2),
                 R12_cv.at<double>(2,0), R12_cv.at<double>(2,1), R12_cv.at<double>(2,2);

      // initialize alpha[1] (a)
      alpha_[1] = Rot3::Logmap(Zr_[1] * Rot3(R12_mat) * Zr_[0].inverse() * Rot3::Expmap(alpha_[0]));
      values_.insert(Symbol('a',1), alpha_[1]);

      // initialize omega[1] (o)
      omega_[1] << 0.0, 0.0, 0.0;
      values_.insert(Symbol('o',1), omega_[1]);
    } else {
      // initialize alpha[i] (a)
      Vector6 alpha_and_omega = PropagateFactor::propagate(Ts_[i]-Ts_[i-1], alpha_[i-1], omega_[i-1], k_);
      alpha_[i] = alpha_and_omega.segment<3>(0);
      values_.insert(Symbol('a',i), alpha_[i]);

      // initialize omega[i] (o)
      omega_[i] = alpha_and_omega.segment<3>(3);
      values_.insert(Symbol('o',i), omega_[i]);
    }

    // investivate N, which is the set of indices of landmark that appeared twice at first in frame #i
    vector<int> N, N_inv;
    N.resize(0);
    N_inv.resize(0);
    for (int j=0; j<n_landmarks(); ++j) {
      // skip landmarks already activated or not observed now
      if (A_[j] || !M_[i][j]) continue;

      // survey frames from #i-1 to #0
      for (int ii=i-1; ii>=0; --ii) {
        if (M_[ii][j]) {
          A_[j] = true; // activate
          N.push_back(j);
          N_inv.push_back(ii);
          break;
        }
      }
    }

    // initialize L[N]
    if (N.size() > 0) {
      // do triangulations between frames #ii and #i at once, ii=min(N_inv):max(N_inv)
      int N_inv_max = *max_element(N_inv.begin(), N_inv.end());
      int N_inv_min = *min_element(N_inv.begin(), N_inv.end());
      for (int ii=N_inv_min; ii<=N_inv_max; ++ii) {
        vector<int> N_now(0);
        for (int jj=0; jj<N_inv.size(); ++jj) {
          if (N_inv[jj]==ii) N_now.push_back(N[jj]);
        }
        if (N_now.size() > 0) triangulation(ii, i, N_now);
      }

      // finally add initialized L's to values
      for (int jj=0; jj<N.size(); ++jj) values_.insert(Symbol('l',N[jj]), L_[N[jj]]);

      // cout << "  initialize L_new by triangulation" << endl;
    }

    /*******************
     SETTING FACTORS
    *******************/

    // add factor between Zr[i-1] and Zr[i] (r)
    graph_.add(BetweenFactor<Rot3>(Symbol('r',i-1), Symbol('r',i), dZr_[i], dZr_noise_));

    // add factor between Zt[i-1] and Zt[i] (t)
    graph_.add(BetweenFactor<Point3>(Symbol('t',i-1), Symbol('t',i), dZt_[i], dZt_noise_));

    // add factor between alpha[i-1], alpha[i] (a), omega[i-1], omega[i] (o), and k (k)
    graph_.add(
      PropagateFactor(
        Ts_[i]-Ts_[i-1],
        Symbol('a',i-1), Symbol('a',i), Symbol('o',i-1), Symbol('o',i), Symbol('k',0),
        noiseModel::Isotropic::Sigma(6, TINYSTD)
      )
    );

    // add absolute factor for {L[N]}
    #ifdef L_ABSPRIOR_STD
      for (size_t jj=0; jj<N.size(); ++jj) {
        graph_.add(PriorFactor<Point3>(Symbol('l',N[jj]), Point3(),
          noiseModel::Isotropic::Sigma(3, L_ABSPRIOR_STD)));
      }
    #endif

    // add absolute factor for Zt[i]
    #ifdef ZT_ABSPRIOR_STD
      graph_.add(PriorFactor<Point3>(Symbol('t',i), Point3(),
        noiseModel::Isotropic::Sigma(3, ZT_ABSPRIOR_STD)));
    #endif

    // add camera factors among L[N[jj]], alpha, Xt, dalpha, Zr[N_inv[jj]], & Zt[N_inv[jj]]
    //   N_inv[jj] is the index of frame where landmark #N[jj] appeared in the PAST (from frames #0 to #i-1)
    for (size_t jj=0; jj<N.size(); ++jj) {
      Point3_ L_in_I = BtoI(Vector3_('a',N_inv[jj]), Point3_('s',0), Point3_('l',N[jj]));
      Point3_ L_in_C = ItoC(Rot3_('r',N_inv[jj]), Point3_('t',N_inv[jj]), L_in_I);
      Point2_ Y_pred = CtoD(L_in_C, Cal3_S2_(*IntrMat_));
      graph_.addExpressionFactor(Y_pred, Y_[N_inv[jj]][N[jj]], Y_noise_);
    }

    // add camera factors among L[E_i], alpha, Xt, dalpha, Zr[i], & Zt[i]
    //   E_i is the set of indices of landmarks that EVER (from frames #0 to #i) appeared more than once
    for (int j=0; j<n_landmarks(); ++j) {
      if (A_[j] && M_[i][j]) {
        Point3_ L_in_I = BtoI(Vector3_('a',i), Point3_('s',0), Point3_('l',j));
        Point3_ L_in_C = ItoC(Rot3_('r',i), Point3_('t',i), L_in_I);
        Point2_ Y_pred = CtoD(L_in_C, Cal3_S2_(*IntrMat_));
        graph_.addExpressionFactor(Y_pred, Y_[i][j], Y_noise_);
      }
    }

    /*******************
     UPDATE
    *******************/

    // optimization by other optimizers
    // (DoglegOptimizer, LevenbergMarquardtOptimizer, etc.)
    values_ = LevenbergMarquardtOptimizer(graph_, values_).optimize();
    get_estimated_values(i+1, values_);

    // display
    cout << "::: Frame #" << i << " end" << endl;

    for (int ii=i; ii<=i; ++ii) {
      Vector3 r = Rot3::Logmap(Zr_[ii]);
      cout << "  Log(Zr_" << ii << ") = [" << r[0] << ", " << r[1] << ", " << r[2] << "]" << endl;
      cout << "  Zt_" << ii << "      = [" << Zt_[ii].x() << ", " << Zt_[ii].y() << ", " << Zt_[ii].z() << "]" << endl;
      cout << "  alpha_" << ii << "   = [" << alpha_[ii][0] << ", " << alpha_[ii][1] << ", " << alpha_[ii][2] << "]" << endl;
      cout << "  omega_" << ii << "   = [" << omega_[ii][0] << ", " << omega_[ii][1] << ", " << omega_[ii][2] << "]" << endl;
    }
  }

  return;
}

/*******************************************************************/
void SFR::triangulation(int i1, int i2, const vector<int> &N_now)
{
  size_t n = N_now.size();

  // convert vector<Point2> to vector<cv::Point2d>
  vector<Point2> Y1=Y_[i1], Y2=Y_[i2];
  vector<cv::Point2f> Y1_cv(n), Y2_cv(n);
  Point2 tmp;
  for (size_t jj=0; jj<n; ++jj) {
    int j = N_now[jj];
    tmp = IntrMat_->calibrate(Y1[j]); Y1_cv[jj] = cv::Point2f(tmp.x(), tmp.y());
    tmp = IntrMat_->calibrate(Y2[j]); Y2_cv[jj] = cv::Point2f(tmp.x(), tmp.y());
  }

  // Xr1, Xr2
  Rot3 Xr1=Rot3::Expmap(alpha_[i1]), Xr2=Rot3::Expmap(alpha_[i2]);

  // X1, X2
  Pose3 X1(Xr1, Xt_), X2(Xr2, Xt_);

  // Z1, Z2
  Pose3 Z1(Zr_[i1], Zt_[i1]), Z2(Zr_[i2], Zt_[i2]);

  // Compute projection matrices.
  Matrix P1_mat
    = Z1.inverse().compose(X1).matrix().block<3,4>(0,0);
  Matrix P2_mat
    = Z2.inverse().compose(X2).matrix().block<3,4>(0,0);
  cv::Mat P1_cv = (cv::Mat_<float>(3,4) <<
    P1_mat(0,0), P1_mat(0,1), P1_mat(0,2), P1_mat(0,3),
    P1_mat(1,0), P1_mat(1,1), P1_mat(1,2), P1_mat(1,3),
    P1_mat(2,0), P1_mat(2,1), P1_mat(2,2), P1_mat(2,3));
  cv::Mat P2_cv = (cv::Mat_<float>(3,4) <<
    P2_mat(0,0), P2_mat(0,1), P2_mat(0,2), P2_mat(0,3),
    P2_mat(1,0), P2_mat(1,1), P2_mat(1,2), P2_mat(1,3),
    P2_mat(2,0), P2_mat(2,1), P2_mat(2,2), P2_mat(2,3));

  // Do triangulation.
  cv::Mat L_now_cv;
  cv::triangulatePoints(P1_cv, P2_cv, Y1_cv, Y2_cv, L_now_cv);

  float Lx, Ly, Lz, Lw;
  for (size_t jj=0; jj<n; ++jj) {
    Lx = L_now_cv.at<float>(0,jj);
    Ly = L_now_cv.at<float>(1,jj);
    Lz = L_now_cv.at<float>(2,jj);
    Lw = L_now_cv.at<float>(3,jj);
    L_[N_now[jj]] = Point3(Lx/Lw, Ly/Lw, Lz/Lw);
  }

  return;
}

/*******************************************************************/
void SFR::get_estimated_values(int n_frames_get, const Values v)
{
  for (int i=0; i<n_frames_get; ++i) {
    Zr_[i] = v.at<Rot3>(Symbol('r',i));
    Zt_[i] = v.at<Point3>(Symbol('t',i));
    alpha_[i] = v.at<Vector3>(Symbol('a',i));
    omega_[i] = v.at<Vector3>(Symbol('o',i));
  }
  Xt_ = v.at<Point3>(Symbol('s',0));
  k_ = v.at<Vector2>(Symbol('k',0));
  for (int j=0; j<n_landmarks(); ++j) {
    if (A_[j]) L_[j] = v.at<Point3>(Symbol('l',j));
  }

  return;
};

/*******************************************************************/
// void SFR::remove_landmark(int landmark_idx)
// {
// }

/*******************************************************************/
} // namespace gtsam
