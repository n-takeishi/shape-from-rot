#pragma once

#include <gtsam/slam/expressions.h>

namespace gtsam
{

typedef Expression<Vector3> Vector3_;

inline Rot3_ RotExpmap(const Vector3_ &alpha) {
  return Rot3_(Rot3::Expmap, alpha);
}

/*
  The following coordinate transformations are implemented as expressions of gtsam:
    <I> is defined by Zr[0] & Zt[0]
    Xr, Xt, & dXr defines trans. <B_i> --> <I>  : (p in <I>)   = (Xr*dXr^dt) * (p in <B_i>) + Xt
    Zr[i] & Zt[i] defines trans. <I> --> <C_i>  : (p in <C_i>) = inv(Zr[i]) * {(p in <I>) - Zt[i]}
    IntrMat       defines trans. <C_i> --> <D_i>
*/

inline Point3_ BtoI(const Vector3_ &alpha, const Point3_ &Xt, const Point3_ &p) {
  return rotate(RotExpmap(alpha), p) + Xt; // Xr = RotExpmap(alpha)
}

inline Point3_ ItoC(const Rot3_ &Zr, const Point3_ &Zt, const Point3_ &p) {
  return unrotate(Zr, p - Zt);
};

inline Point2_ CtoD(const Point3_ &p, const Cal3_S2_ &IntrMat) {
  return project3<Cal3_S2, Point3>(Pose3_(Pose3()), p, IntrMat);
};


} // namespace gtsam