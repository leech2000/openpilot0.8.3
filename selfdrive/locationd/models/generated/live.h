/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8097402282370743041);
void inv_err_fun(double *nom_x, double *true_x, double *out_257000985813042990);
void H_mod_fun(double *state, double *out_7673816066968292903);
void f_fun(double *state, double dt, double *out_4785678607724337276);
void F_fun(double *state, double dt, double *out_8750663048589219192);
void h_3(double *state, double *unused, double *out_5745668108571162435);
void H_3(double *state, double *unused, double *out_5714530975065264385);
void h_4(double *state, double *unused, double *out_2317849258972104220);
void H_4(double *state, double *unused, double *out_1645757970379999256);
void h_9(double *state, double *unused, double *out_2887333069678608230);
void H_9(double *state, double *unused, double *out_676747573490856714);
void h_10(double *state, double *unused, double *out_4194581018197861753);
void H_10(double *state, double *unused, double *out_4929995238301780605);
void h_12(double *state, double *unused, double *out_8352872033279664398);
void H_12(double *state, double *unused, double *out_460797541359210313);
void h_31(double *state, double *unused, double *out_1272238562537136905);
void H_31(double *state, double *unused, double *out_2214292233784209508);
void h_32(double *state, double *unused, double *out_6600884388431956993);
void H_32(double *state, double *unused, double *out_8690444182527386028);
void h_13(double *state, double *unused, double *out_9116012224960326151);
void H_13(double *state, double *unused, double *out_268871762386714667);
void h_14(double *state, double *unused, double *out_2887333069678608230);
void H_14(double *state, double *unused, double *out_676747573490856714);
void h_19(double *state, double *unused, double *out_238862175807085635);
void H_19(double *state, double *unused, double *out_5109427424593287895);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_31 = 7.814728;
void update_31(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_32 = 9.487729;
void update_32(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);