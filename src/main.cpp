/*
 * Shape From Rotation
 * (c) Naoya Takeishi, 2018.
 */

#include <string>
#include <fstream>
#include <boost/program_options.hpp>
#include "sfr.h"

using namespace std;
using namespace gtsam;

/*******************************************************************/
/*
  numX:(number of images)
  numL:(number of landmarks)
  X:(image#)
  L:(landmark#) u v
  L:(landmark#) u v
  ...
  X:(image#)
  ...
*/
int load_Y(
  const string filename,
  vector< vector<Point2> > &Y,
  vector< vector<bool> > &M)
{
  ifstream ifs(filename);
  string str;

  if (ifs.fail()) return 0;

  // number of frames and landmarks
  int num_X = 0, num_L = 0;
  while (getline(ifs, str)) {
    if (str[0] == '#' || str.size() < 3) {
      continue;
    } else if (str.substr(0,5) == "numX:") {
      num_X = stoi(str.substr(5));
    } else if (str.substr(0,5) == "numL:") {
      num_L = stoi(str.substr(5));
    }
    if (num_X > 0 && num_L > 0) break;
  }

  // nemory allocation
  Y.resize(num_X);
  M.resize(num_X);
  for (int i=0; i<num_X; ++i) {
    Y[i].resize(num_L);
    M[i].resize(num_L);
    for (int j=0; j<num_L; ++j) M[i][j] = false;
  }

  // read main data
  int X = 0, L = 0;
  double u = 0.0, v = 0.0;
  while (getline(ifs, str)) {
    if (str[0] == '#' || str.size() < 3) {
      continue;
    } else if (str.substr(0,2) == "X:") {
      X = stoi(str.substr(2));
    } else if (str.substr(0,2) == "L:") {
      stringstream ss(str.substr(2));
      string buf;
      getline(ss, buf, ' '); L = stoi(buf);
      getline(ss, buf, ' '); u = stod(buf);
      getline(ss, buf, ' '); v = stod(buf);
      Y[X][L] = Point2(u, v);
      M[X][L] = true;
    } else {
      return 0;
    }
  }

  return 1;
}

/*******************************************************************/
/*
  qw qx qy qz
  ...
*/
int load_dZr(
  const string filename,
  const int n_frames,
  vector<Rot3> &dZr)
{
  ifstream ifs(filename);
  string str;

  if (ifs.fail()) return 0;

  int i = 0;
  double qw = 1.0, qx = 0.0, qy = 0.0, qz = 0.0;
  dZr.resize(n_frames);
  while (getline(ifs, str)) {
    stringstream ss(str);
    string buf;
    getline(ss, buf, ' '); qw = stod(buf);
    getline(ss, buf, ' '); qx = stod(buf);
    getline(ss, buf, ' '); qy = stod(buf);
    getline(ss, buf, ' '); qz = stod(buf);
    dZr[i] = Rot3(Quaternion(qw,qx,qy,qz));
    if (i==n_frames-1) {
      break;
    } else {
      i++;
    }
  }

  return 1;
}

/*******************************************************************/
// int load_dZt(
//   const string filename,
//   const int n_frames,
//   vector<Point3> &dZt)
// {
// }

/*******************************************************************/
// int load_Ts(
//   const string filename,
//   const int n_frames,
//   vector<double> &Ts)
// {
// }

/*******************************************************************/
int main(int argc, char** argv)
{
  using namespace boost::program_options;

  // define comand-line options
  options_description desc_cmd("Command-line Options");
  desc_cmd.add_options()
    ("rootdir,d",     value<string>()->required(),         "Input/output root directory")
    ("outputdir",     value<string>()->default_value("."), "Output directory (under rootdir)")
    ("configfile",    value<string>()->default_value("config.txt"),   "Configuration file")
    ("cameradata",    value<string>()->default_value("camera.txt"),   "Camera observation")
    ("scmoverdata",   value<string>()->default_value("scmover.txt"),  "Spacecraft rotation observation (optional)")
    ("scmovetdata",   value<string>()->default_value("scmovet.txt"),  "Spacecraft translation observation (optional)")
    ("timestepdata",  value<string>()->default_value("timestep.txt"), "Timesteps (optional)")
  ;

  // define config file options
  options_description desc_cfg("Configuration File Options");
  desc_cfg.add_options()
    ("USE_SCMOVER",           value<bool>()->default_value(false), "")
    ("USE_SCMOVET",           value<bool>()->default_value(false), "")
    ("USE_TIMESTEP",          value<bool>()->default_value(false), "")
    ("CAMERA_SKEW",           value<double>()->default_value(0.0), "")
    ("CAMERA_FOCAL_X",        value<double>(), "")
    ("CAMERA_FOCAL_Y",        value<double>(), "")
    ("CAMERA_CENTER_X",       value<double>(), "")
    ("CAMERA_CENTER_Y",       value<double>(), "")
    ("NOISESTD_CAMERA",       value<double>(), "")
    ("NOISESTD_SCMOVER_X",    value<double>(), "")
    ("NOISESTD_SCMOVER_Y",    value<double>(), "")
    ("NOISESTD_SCMOVER_Z",    value<double>(), "")
    ("NOISESTD_SCMOVET_X",    value<double>(), "")
    ("NOISESTD_SCMOVET_Y",    value<double>(), "")
    ("NOISESTD_SCMOVET_Z",    value<double>(), "")
    ("PRIORMEAN_XT_X",        value<double>(), "")
    ("PRIORMEAN_XT_Y",        value<double>(), "")
    ("PRIORMEAN_XT_Z",        value<double>(), "")
  ;

  variables_map vm;
  string rootdir;

  // parse command-line argument (and specify config file)
  store(parse_command_line(argc, argv, desc_cmd), vm);
  notify(vm);
  rootdir = vm["rootdir"].as<string>();
  cout << "Root directory: " << rootdir << endl;

  // parse config file
  string configfile = rootdir+"/"+vm["configfile"].as<string>();
  ifstream ifs(configfile);
  store(parse_config_file(ifs, desc_cfg), vm);
  notify(vm);
  cout << "Configuration file: " << configfile << endl;

  // prepare IntrMat
  double fx, fy, s, cx, cy;
  fx = vm["CAMERA_FOCAL_X"].as<double>();
  fy = vm["CAMERA_FOCAL_Y"].as<double>();
  s  = vm["CAMERA_SKEW"].as<double>();
  cx = vm["CAMERA_CENTER_X"].as<double>();
  cy = vm["CAMERA_CENTER_Y"].as<double>();
  Cal3_S2::shared_ptr IntrMat = Cal3_S2::shared_ptr(new Cal3_S2(fx, fy, s, cx, cy));
  cout << "Set camera intrinsic parameters: "
    << "fx=" << fx << ", fy=" << fy << ", s=" << s << ", cx=" << cx << ", cy=" << cy << endl;

  // prepare Y
  vector< vector<Point2> > Y;
  vector< vector<bool> > M;
  load_Y(rootdir+"/"+vm["cameradata"].as<string>(), Y, M);
  int num_X = Y.size();
  printf("Load Y. (num_X=%d)\n", num_X);

  // prepare dZr if any
  vector<Rot3> dZr;
  if (vm["USE_SCMOVER"].as<bool>()) {
    load_dZr(rootdir+"/"+vm["scmoverdata"].as<string>(), num_X, dZr);
    printf("Load dZr. (%lu)\n", dZr.size());
  } else {
    dZr.resize(num_X);
    for (int i=0; i<num_X; ++i) dZr[i] = Rot3();
    printf("Set dZr to be zeros. (%lu)\n", dZr.size());
  }

  // prepare dZt if any
  vector<Point3> dZt;
  if (vm["USE_SCMOVET"].as<bool>()) {
    // load_dZt(rootdir+"/"+vm["scmovetdata"].as<string>(), num_X, dZt);
    // printf("Load dZt. (%lu)\n", dZt.size());
  } else {
    dZt.resize(num_X);
    for (int i=0; i<num_X; ++i) dZt[i] = Point3();
    printf("Set dZt to be zeros. (%lu)\n", dZt.size());
  }

  // prepare Ts if any
  vector<double> Ts;
  if (vm["USE_TIMESTEP"].as<bool>()) {
    // load_Ts(rootdir+"/"+vm["timestepdata"].as<string>(), num_X, Ts);
    // printf("Load Ts. (%lu)\n", Ts.size());
  } else {
    Ts.resize(num_X);
    for (int i=0; i<num_X; ++i) Ts[i] = (double)i;
    printf("Set Ts to be 0:%lu.\n", Ts.size()-1);
  }

  // make SFR instance
  double Y_noise_std = vm["NOISESTD_CAMERA"].as<double>();
  Vector3 dZr_noise_std(
    vm["NOISESTD_SCMOVER_X"].as<double>(),
    vm["NOISESTD_SCMOVER_Y"].as<double>(),
    vm["NOISESTD_SCMOVER_Z"].as<double>());
  Vector3 dZt_noise_std(
    vm["NOISESTD_SCMOVET_X"].as<double>(),
    vm["NOISESTD_SCMOVET_Y"].as<double>(),
    vm["NOISESTD_SCMOVET_Z"].as<double>());
  SFR S(IntrMat, Y_noise_std, dZr_noise_std, dZt_noise_std);
  cout << "Prepare S." << endl;

  // initialize SFR instance
  Point3 Xt_init(
    vm["PRIORMEAN_XT_X"].as<double>(),
    vm["PRIORMEAN_XT_Y"].as<double>(),
    vm["PRIORMEAN_XT_Z"].as<double>());
  S.initialize(Rot3(), Xt_init);
  cout << "Solver initialized." << endl;

  // main estimation
  cout << "Start estimation." << endl;
  S.update(Y, M, dZr, dZt, Ts);
  cout << "Finish estimation." << endl;

  // output k, Xt, Xr, omega
  string outfilename_X = rootdir+"/"+vm["outputdir"].as<string>()+"/output_k_Xt_Xr_omega.txt";
  ofstream ofs_X(outfilename_X.c_str());
  ofs_X << S.k()[0] << " " << S.k()[1] << endl;
  ofs_X << S.Xt().x() << " " << S.Xt().y() << " " << S.Xt().z() << endl;
  for (int i=0; i<S.n_frames(); ++i) {
    Vector4 Xr_q = Rot3::Expmap(S.alpha(i)).quaternion();
    Vector3 omega = S.omega(i);
    ofs_X << Xr_q[0] << " " << Xr_q[1] << " " << Xr_q[2] << " " << Xr_q[3] << " "
          << omega[0] << " " << omega[1] << " " << omega[2] << endl;
  }
  ofs_X.close();
  cout << "Output k, Xt, Xr, and omega." << endl;

  // output Zr, Zt
  string outfilename_Z = rootdir+"/"+vm["outputdir"].as<string>()+"/output_Zr_Zt.txt";
  ofstream ofs_Z(outfilename_Z.c_str());
  for (int i=0; i<S.n_frames(); ++i) {
    Vector4 Zr_q = S.Zr(i).quaternion();
    ofs_Z << Zr_q[0] << " " << Zr_q[1] << " " << Zr_q[2] << " " << Zr_q[3] << " "
          << S.Zt(i).x() << " " << S.Zt(i).y() << " " << S.Zt(i).z() << endl;
  }
  ofs_Z.close();
  cout << "Output Zr and Zt." << endl;

  // output L
  string outfilename_L = rootdir+"/"+vm["outputdir"].as<string>()+"/output_L.txt";
  ofstream ofs_L(outfilename_L.c_str());
  for (int j=0; j<S.n_landmarks(); ++j) {
    ofs_L << S.L(j).x() << " " <<  S.L(j).y() << " " << S.L(j).z() << endl;
  }
  ofs_L.close();
  cout << "Output L." << endl;

  return 0;
}
