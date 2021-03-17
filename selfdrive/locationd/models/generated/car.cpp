
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7116916163456468248) {
   out_7116916163456468248[0] = delta_x[0] + nom_x[0];
   out_7116916163456468248[1] = delta_x[1] + nom_x[1];
   out_7116916163456468248[2] = delta_x[2] + nom_x[2];
   out_7116916163456468248[3] = delta_x[3] + nom_x[3];
   out_7116916163456468248[4] = delta_x[4] + nom_x[4];
   out_7116916163456468248[5] = delta_x[5] + nom_x[5];
   out_7116916163456468248[6] = delta_x[6] + nom_x[6];
   out_7116916163456468248[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4779835383543979127) {
   out_4779835383543979127[0] = -nom_x[0] + true_x[0];
   out_4779835383543979127[1] = -nom_x[1] + true_x[1];
   out_4779835383543979127[2] = -nom_x[2] + true_x[2];
   out_4779835383543979127[3] = -nom_x[3] + true_x[3];
   out_4779835383543979127[4] = -nom_x[4] + true_x[4];
   out_4779835383543979127[5] = -nom_x[5] + true_x[5];
   out_4779835383543979127[6] = -nom_x[6] + true_x[6];
   out_4779835383543979127[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_1747493067526112467) {
   out_1747493067526112467[0] = 1.0;
   out_1747493067526112467[1] = 0.0;
   out_1747493067526112467[2] = 0.0;
   out_1747493067526112467[3] = 0.0;
   out_1747493067526112467[4] = 0.0;
   out_1747493067526112467[5] = 0.0;
   out_1747493067526112467[6] = 0.0;
   out_1747493067526112467[7] = 0.0;
   out_1747493067526112467[8] = 0.0;
   out_1747493067526112467[9] = 1.0;
   out_1747493067526112467[10] = 0.0;
   out_1747493067526112467[11] = 0.0;
   out_1747493067526112467[12] = 0.0;
   out_1747493067526112467[13] = 0.0;
   out_1747493067526112467[14] = 0.0;
   out_1747493067526112467[15] = 0.0;
   out_1747493067526112467[16] = 0.0;
   out_1747493067526112467[17] = 0.0;
   out_1747493067526112467[18] = 1.0;
   out_1747493067526112467[19] = 0.0;
   out_1747493067526112467[20] = 0.0;
   out_1747493067526112467[21] = 0.0;
   out_1747493067526112467[22] = 0.0;
   out_1747493067526112467[23] = 0.0;
   out_1747493067526112467[24] = 0.0;
   out_1747493067526112467[25] = 0.0;
   out_1747493067526112467[26] = 0.0;
   out_1747493067526112467[27] = 1.0;
   out_1747493067526112467[28] = 0.0;
   out_1747493067526112467[29] = 0.0;
   out_1747493067526112467[30] = 0.0;
   out_1747493067526112467[31] = 0.0;
   out_1747493067526112467[32] = 0.0;
   out_1747493067526112467[33] = 0.0;
   out_1747493067526112467[34] = 0.0;
   out_1747493067526112467[35] = 0.0;
   out_1747493067526112467[36] = 1.0;
   out_1747493067526112467[37] = 0.0;
   out_1747493067526112467[38] = 0.0;
   out_1747493067526112467[39] = 0.0;
   out_1747493067526112467[40] = 0.0;
   out_1747493067526112467[41] = 0.0;
   out_1747493067526112467[42] = 0.0;
   out_1747493067526112467[43] = 0.0;
   out_1747493067526112467[44] = 0.0;
   out_1747493067526112467[45] = 1.0;
   out_1747493067526112467[46] = 0.0;
   out_1747493067526112467[47] = 0.0;
   out_1747493067526112467[48] = 0.0;
   out_1747493067526112467[49] = 0.0;
   out_1747493067526112467[50] = 0.0;
   out_1747493067526112467[51] = 0.0;
   out_1747493067526112467[52] = 0.0;
   out_1747493067526112467[53] = 0.0;
   out_1747493067526112467[54] = 1.0;
   out_1747493067526112467[55] = 0.0;
   out_1747493067526112467[56] = 0.0;
   out_1747493067526112467[57] = 0.0;
   out_1747493067526112467[58] = 0.0;
   out_1747493067526112467[59] = 0.0;
   out_1747493067526112467[60] = 0.0;
   out_1747493067526112467[61] = 0.0;
   out_1747493067526112467[62] = 0.0;
   out_1747493067526112467[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_1748629533821697514) {
   out_1748629533821697514[0] = state[0];
   out_1748629533821697514[1] = state[1];
   out_1748629533821697514[2] = state[2];
   out_1748629533821697514[3] = state[3];
   out_1748629533821697514[4] = state[4];
   out_1748629533821697514[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_1748629533821697514[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_1748629533821697514[7] = state[7];
}
void F_fun(double *state, double dt, double *out_6428121403010973063) {
   out_6428121403010973063[0] = 1;
   out_6428121403010973063[1] = 0;
   out_6428121403010973063[2] = 0;
   out_6428121403010973063[3] = 0;
   out_6428121403010973063[4] = 0;
   out_6428121403010973063[5] = 0;
   out_6428121403010973063[6] = 0;
   out_6428121403010973063[7] = 0;
   out_6428121403010973063[8] = 0;
   out_6428121403010973063[9] = 1;
   out_6428121403010973063[10] = 0;
   out_6428121403010973063[11] = 0;
   out_6428121403010973063[12] = 0;
   out_6428121403010973063[13] = 0;
   out_6428121403010973063[14] = 0;
   out_6428121403010973063[15] = 0;
   out_6428121403010973063[16] = 0;
   out_6428121403010973063[17] = 0;
   out_6428121403010973063[18] = 1;
   out_6428121403010973063[19] = 0;
   out_6428121403010973063[20] = 0;
   out_6428121403010973063[21] = 0;
   out_6428121403010973063[22] = 0;
   out_6428121403010973063[23] = 0;
   out_6428121403010973063[24] = 0;
   out_6428121403010973063[25] = 0;
   out_6428121403010973063[26] = 0;
   out_6428121403010973063[27] = 1;
   out_6428121403010973063[28] = 0;
   out_6428121403010973063[29] = 0;
   out_6428121403010973063[30] = 0;
   out_6428121403010973063[31] = 0;
   out_6428121403010973063[32] = 0;
   out_6428121403010973063[33] = 0;
   out_6428121403010973063[34] = 0;
   out_6428121403010973063[35] = 0;
   out_6428121403010973063[36] = 1;
   out_6428121403010973063[37] = 0;
   out_6428121403010973063[38] = 0;
   out_6428121403010973063[39] = 0;
   out_6428121403010973063[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6428121403010973063[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6428121403010973063[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6428121403010973063[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6428121403010973063[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6428121403010973063[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6428121403010973063[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6428121403010973063[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6428121403010973063[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6428121403010973063[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6428121403010973063[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6428121403010973063[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6428121403010973063[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6428121403010973063[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6428121403010973063[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6428121403010973063[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6428121403010973063[56] = 0;
   out_6428121403010973063[57] = 0;
   out_6428121403010973063[58] = 0;
   out_6428121403010973063[59] = 0;
   out_6428121403010973063[60] = 0;
   out_6428121403010973063[61] = 0;
   out_6428121403010973063[62] = 0;
   out_6428121403010973063[63] = 1;
}
void h_25(double *state, double *unused, double *out_3724627254418807924) {
   out_3724627254418807924[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8719532070886801573) {
   out_8719532070886801573[0] = 0;
   out_8719532070886801573[1] = 0;
   out_8719532070886801573[2] = 0;
   out_8719532070886801573[3] = 0;
   out_8719532070886801573[4] = 0;
   out_8719532070886801573[5] = 0;
   out_8719532070886801573[6] = 1;
   out_8719532070886801573[7] = 0;
}
void h_24(double *state, double *unused, double *out_6978128924664945159) {
   out_6978128924664945159[0] = state[4];
   out_6978128924664945159[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6961217986172828315) {
   out_6961217986172828315[0] = 0;
   out_6961217986172828315[1] = 0;
   out_6961217986172828315[2] = 0;
   out_6961217986172828315[3] = 0;
   out_6961217986172828315[4] = 1;
   out_6961217986172828315[5] = 0;
   out_6961217986172828315[6] = 0;
   out_6961217986172828315[7] = 0;
   out_6961217986172828315[8] = 0;
   out_6961217986172828315[9] = 0;
   out_6961217986172828315[10] = 0;
   out_6961217986172828315[11] = 0;
   out_6961217986172828315[12] = 0;
   out_6961217986172828315[13] = 1;
   out_6961217986172828315[14] = 0;
   out_6961217986172828315[15] = 0;
}
void h_30(double *state, double *unused, double *out_3881813922769937069) {
   out_3881813922769937069[0] = state[4];
}
void H_30(double *state, double *unused, double *out_113313091873566763) {
   out_113313091873566763[0] = 0;
   out_113313091873566763[1] = 0;
   out_113313091873566763[2] = 0;
   out_113313091873566763[3] = 0;
   out_113313091873566763[4] = 1;
   out_113313091873566763[5] = 0;
   out_113313091873566763[6] = 0;
   out_113313091873566763[7] = 0;
}
void h_26(double *state, double *unused, double *out_845027871325604384) {
   out_845027871325604384[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8385775869127102491) {
   out_8385775869127102491[0] = 0;
   out_8385775869127102491[1] = 0;
   out_8385775869127102491[2] = 0;
   out_8385775869127102491[3] = 0;
   out_8385775869127102491[4] = 0;
   out_8385775869127102491[5] = 0;
   out_8385775869127102491[6] = 0;
   out_8385775869127102491[7] = 1;
}
void h_27(double *state, double *unused, double *out_4668137107320470359) {
   out_4668137107320470359[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1174268895963058549) {
   out_1174268895963058549[0] = 0;
   out_1174268895963058549[1] = 0;
   out_1174268895963058549[2] = 0;
   out_1174268895963058549[3] = 1;
   out_1174268895963058549[4] = 0;
   out_1174268895963058549[5] = 0;
   out_1174268895963058549[6] = 0;
   out_1174268895963058549[7] = 0;
}
void h_29(double *state, double *unused, double *out_4914495541496060071) {
   out_4914495541496060071[0] = state[1];
}
void H_29(double *state, double *unused, double *out_687871705729102938) {
   out_687871705729102938[0] = 0;
   out_687871705729102938[1] = 1;
   out_687871705729102938[2] = 0;
   out_687871705729102938[3] = 0;
   out_687871705729102938[4] = 0;
   out_687871705729102938[5] = 0;
   out_687871705729102938[6] = 0;
   out_687871705729102938[7] = 0;
}
void h_28(double *state, double *unused, double *out_6359701166335554340) {
   out_6359701166335554340[0] = state[5];
   out_6359701166335554340[1] = state[6];
}
void H_28(double *state, double *unused, double *out_2578335679334628390) {
   out_2578335679334628390[0] = 0;
   out_2578335679334628390[1] = 0;
   out_2578335679334628390[2] = 0;
   out_2578335679334628390[3] = 0;
   out_2578335679334628390[4] = 0;
   out_2578335679334628390[5] = 1;
   out_2578335679334628390[6] = 0;
   out_2578335679334628390[7] = 0;
   out_2578335679334628390[8] = 0;
   out_2578335679334628390[9] = 0;
   out_2578335679334628390[10] = 0;
   out_2578335679334628390[11] = 0;
   out_2578335679334628390[12] = 0;
   out_2578335679334628390[13] = 0;
   out_2578335679334628390[14] = 1;
   out_2578335679334628390[15] = 0;
}
}

extern "C"{
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
