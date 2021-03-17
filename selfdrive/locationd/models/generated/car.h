/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_7116916163456468248);
void inv_err_fun(double *nom_x, double *true_x, double *out_4779835383543979127);
void H_mod_fun(double *state, double *out_1747493067526112467);
void f_fun(double *state, double dt, double *out_1748629533821697514);
void F_fun(double *state, double dt, double *out_6428121403010973063);
void h_25(double *state, double *unused, double *out_3724627254418807924);
void H_25(double *state, double *unused, double *out_8719532070886801573);
void h_24(double *state, double *unused, double *out_6978128924664945159);
void H_24(double *state, double *unused, double *out_6961217986172828315);
void h_30(double *state, double *unused, double *out_3881813922769937069);
void H_30(double *state, double *unused, double *out_113313091873566763);
void h_26(double *state, double *unused, double *out_845027871325604384);
void H_26(double *state, double *unused, double *out_8385775869127102491);
void h_27(double *state, double *unused, double *out_4668137107320470359);
void H_27(double *state, double *unused, double *out_1174268895963058549);
void h_29(double *state, double *unused, double *out_4914495541496060071);
void H_29(double *state, double *unused, double *out_687871705729102938);
void h_28(double *state, double *unused, double *out_6359701166335554340);
void H_28(double *state, double *unused, double *out_2578335679334628390);
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
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
