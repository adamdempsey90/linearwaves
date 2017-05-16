/* disk.c */
double scaleH(double x, Params params, Disk *disk);
double omegaK(double x, Params params, Disk *disk);
double nu_func(double x, Params params, Disk *disk);
double cs2_func(double x, Params params, Disk *disk);
double dlvr_func(double x, Params params, Disk *disk);
double vr_func(double x, Params params, Disk *disk);
double cs_func(double x, Params params, Disk *disk);
double omega2_func(double x, Params params, Disk *disk);
double kappa2_func(double x, Params params, Disk *disk);
double kappa_func(double x, Params params, Disk *disk);
double omega_func(double x, Params params, Disk *disk);
double k2om_func(double x, Params params, Disk *disk);
double sigma_func(double x, Params params, Disk *disk);
double pres_func(double x, Params params, Disk *disk);
double dsdr_func(double x, Params params, Disk *disk);
double dpdr_func(double x, Params params, Disk *disk);
double dc2dr_func(double x, Params params, Disk *disk);
double Dfunc(int i, double omp, int m, Params params, Disk *disk);
double complex sig(int i, double complex s, Params params, Disk *disk);
void alloc_disk(Disk *disk, int n);
void free_disk(Disk *disk);
void init_disk(char *fname, double *lr, Disk *disk, Params params);
/* viscosity.c */
void viscosity_coeffs_u(int indx, double r, double complex *res, int m, double fac, Params params, Disk *disk);
void viscosity_coeffs_v(int indx, double r, double complex *res, int m, double fac, Params params, Disk *disk);
void viscosity_dcoeffs_u(int indx, double r, double complex *res, int m, double fac, Params params, Disk *disk);
void viscosity_dcoeffs_v(int indx, double r, double complex *res, int m, double fac, Params params, Disk *disk);
void viscosity_d2coeffs_u(int indx, double r, double complex *res, double fac, Params params, Disk *disk);
void viscosity_d2coeffs_v(int indx, double r, double complex *res, double fac, Params params, Disk *disk);
/* coeffs.c */
void main_diag(int indx, double r, int m, double complex *res, Params params, Disk *disk);
void upper_diag(int indx, double r, int m, double complex *res, Params params, Disk *disk);
void lower_diag(int indx, double r, int m, double complex *res, Params params, Disk *disk);
void zero_inner_bc(double complex *md0, double complex *ud0);
void zero_outer_bc(double complex *mdn, double complex *ldn);
void lw_inner_bc(int indx, double r0, int m, int eps, double complex *md0, double complex *ud0, Params params, Disk *disk);
void lw_outer_bc(int indx, double rn, int m, int eps, double complex *mdn, double complex *ldn, Params params, Disk *disk);
void construct_matrix(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, double *dppot, double *drpot, int m, Params params, Disk *disk);
void construct_matrix_second(double *r, double complex *ld, double complex *md, double complex *ud, double complex *fd, double complex *Ru, double complex *Rv,double complex *Rs, int m, Params params, Disk *disk);
/* ctridiag.c */
void cmatmat(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta, int nA);
void cmatvec(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta, int nB);
void cmatmat3(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta);
void cmatvec3(double complex *A, double complex *B, double complex *C, double complex alpha, double complex beta);
void csolve(double complex *A, double complex *B, int nA, int nRHS);
void cthomas_alg(double complex *a, double complex *b, double complex *c, double complex *d, int n);
void cthomas_alg_block(double complex *a, double complex *b, double complex *c, double complex *d, int n, int m);
void cconstruct_total_matrix(double complex *ld, double complex *md, double complex *ud, double complex *mat, int n, int m);
/* fft.c */
void fft_potential(double *r, double *pot1, double *pot2, int num_modes, int nphi, int howmany, double eps2, int indirect);
/* interpolation.c */
void read_sigma(char *fname, double *lr, double *sigma, double *dlsdlr, double *d2lsdlr,double *omega, double *dlomdlr, double *d2lomdlr, int n,int readomega);
/* read_params.c */
void set_var(char *name, int int_val, double double_val, int bool_val, char *str_val, Params *params);
void parse_argument(int argc, char *argv[], Params *params);
void read_param_file(char *fname, int argc, char *argv[], Params *params);
/* torques.c */
void calc_torques(double *r, double complex *us, double complex *vs, double complex *ss, double *fw, double *drfw, double *lamex, double *lamdep, double complex *sol, double *dppot, double *TL, double *TR,int m, Params params, Disk *disk, int silent);
/* output.c */
void output_disk(char *fname, Params params, Disk *disk);
void output(double *r, double complex *sol, double *lamex, double *lamdep, double *drfw, double *fw, char *fname, Params params);
void output_torques(char *fname, Params params, Grid *grid);
void output_matrix(double complex *ld, double complex *md, double complex *ud, double complex *fd, Params params);
void output_second(double *r, double complex *sol, char *fname, Params params);
