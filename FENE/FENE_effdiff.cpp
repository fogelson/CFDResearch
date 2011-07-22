// added flag -turnfivept_on which will read in tr_t.dat and spit out tr_fivept_t.dat
// and will also turn on flag fivept_on

// working on flag -rand_tr which when activated
// will randomize the initial tracer points 


// added flag -fivept_on which will create a  5pt stencil around
// tracer points and update them for estimating lyap exponents
//
//
// added flag m and default falseod_off to turn off the mod in the tracer points
// use -mod_off to run tracer points without mod
//

// update trace with bicubic interpolation using FFT to get derivatives at corners

// 2/27/07 take away 2/3 de-alias  replace with Hou filter

// modify 1/31/07 to limit modes of u used

// 5/23/06 modify to have time dependent f

// 4/28/06

// fixed de-aliasing 4/28/06

// modified 4/25/06 to have dynamic filenames!

// agrees with matlab code (n=128, 1000 iterations differs O(10^-15) for S and O(10^-16) for U

// VERSION 5 updates Shat with ABCN with 2/3 de-aliasing

// VERSION 4 does all of below but updates Shat rather than S
// still with Runge Kutta

// VERSION 3 does 2/3 versus doubling for de-aliasing

// This is the fully working fftw C code as of 2/24/06
// for 2 d Oldroyd-B

//  modified for ibm fftw by estarose 1/26/06
//  modified by becca 1/3/06
//  rearranged the order to save after compute u,
//  now writing files with 20.16e precision
//  updating u in RK step!!

// Two dimension of Oldroyd-B.


#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fftw.h>

#include "FENETiming.h"

#include "FENEStress.h"

using namespace std;

#define EPS .0000001

// Fill in the following parameters: N, dt (optional), Total_iterations, Wi, NU, start_w_zero, pertid, timedep_F,
// lastsaved, time_lastsaved

int N;
double dt;
double Wi;
double NU;
int N2;
int KK;
double reciprocal_Wi;
double Beta;
int start_w_zero;
int pertid;
int timedep_F;
int tracer;
int restart_tr;
int mod_off;
int fivept_on;
int fivept;
int turnfivept_on;
int rand_tr;

// FENE Parameters
double deltaX, h, offset, D, H, Q0;

double lastsaved;
double finalTime;
double saveTime;
int dtOn;
int betaOn;
int tr_pts;

//  Claims of sub-functions.


void Hou_filter(double *Ahat_re, double *Ahat_im, int k_dim);

void Add_Shat(double *mid_Shat_Re, double *mid_shat_Im, double *S_hat_Re,
		double *S_hat_Im);

void Get_Frequency(int Freq[]);

void initialization_S(double *S);

void initialization_tr(double *tr);

void init_tr_rand(double *tr);

void initialization_Spert(double *S);

void Get_F_hat(double *F_hat_Re, double *F_hat_Im);

void Get_F_hat_t(double *F_hat_Re, double *F_hat_Im, double time,
		int iterations);

void Get_hat(double *A, double *A_hat_Re, double *A_hat_Im, int k_dim);

void Get_P_U_hat(double *P_hat_Re, double *P_hat_Im, double *U_hat_Re,
		double *U_hat_Im, double *S_hat_Re, double *S_hat_Im, double *F_hat_Re,
		double *F_hat_Im, int Freq[]);

void Get_gradUS_convolution(double *Re_gradUS_hat, double *Im_gradUS_hat,
		double *S_hat_Re, double *S_hat_Im, double *U_hat_Re, double *U_hat_Im,
		int Freq[]);

void Get_UgradS_convolution(double *Re_UgradS_hat, double *Im_UgradS_hat,
		double *S_hat_Re, double *S_hat_Im, double *U_hat_Re, double *U_hat_Im,
		int Freq[]);

void Get_U_S(double *U, double *U_hat_Re, double *U_hat_Im, int k_dim);

void update_S_hat(double *New_Shat_Re, double *New_Shat_Im, double *S_hat_Re,
		double *S_hat_Im, double *U_hat_Re, double *U_hat_Im, int Freq[]);

void Save_U_S(double *U, char filename[], int k_dim);

void Save_tr(double *tr, char filename[], int fivept);

void zeroboundary_hat(double *A_hat_Re, double *A_hat_Im, int k_dim);

void Readin_S(double *S, char filename[], int k_dim);

void Readin_tr(double *tr, char filename[], int fivept);

void Euler_Step(double *S_hat_Re, double *S_hat_Im, double *RHS_Re,
		double *RHS_Im, double nu, double ddt, int Freq[]);

void ABCN_Step(double *S_hat_Re, double *S_hat_Im, double *RHS_Re,
		double *RHS_Im, double *old_RHS_Re, double *old_RHS_Im, double nu,
		double ddt, int Freq[]);

void Euler_tr(double *tr, double *update, int fivept);

void AB_tr(double *tr, double *update, double *old_update, int fivept);

void get_fivept_tr(double *tr, double *tr_fivept);

void Get_update(double *tr, double *U, double *U_hat_Re, double *U_hat_Im,
		double *update, int Freq[], int fivept);

void bcuint(double *y, double *y1, double *y2, double *y12, double x1l,
		double x1u, double x2l, double x2u, double x1, double x2, double *ansy);

void readParameters(int argc, char **argv);

void Get_mixd2U(double *gradU_hat_re, double *gradU_hat_im, double *mixd2U,
		int Freq[]);

void Get_grad_hat(double *A_hat_re, double *A_hat_im, double *gradA_hat_re,
		double *gradA_hat_im, int Freq[], int k_dim);

//  END OF CLAIMS.

// define plans for fftw so that they are reusable

fftwnd_plan planN;
double PI2;

//#define BECCA_DEBUG

#ifdef BECCA_DEBUG
#define UXY(kx,ky) \
		({ \
	cerr << "kx=" << kx << " ky=" << ky << " " << 2*(kx*N+ky) << endl;	\
	2*(kx*N+ky); \
})
#else
#define UXY(kx,ky)   2*(kx*N+ky)
#endif

////////////////////////////////////////////////////////////////////////////
////////////////////////////---Main Part---/////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

	int Total_iterations;
	int Record_iterations;

	double *tr;
	double *update;
	double *old_update;
	double *tr_fivept;
	double *update_fivept;
	double *old_update_fivept;

	// double *tr_new;

	//      default parameters
	N = 512;
	Wi = 0.1;
	NU = 0.;
	start_w_zero = true;
	pertid = false;
	timedep_F = false;
	tracer = false;

	lastsaved = 1.;
	dtOn = true;
	betaOn = true;
	finalTime = 4.;
	saveTime = 1.;
	tr_pts = 64;
	restart_tr = false;
	mod_off = false;
	fivept_on = false;
	turnfivept_on = false;
	rand_tr = false;


	// ****************************************** //
	// FENE PARAMETERS
	Q0 = PI2 / 100.0;
	h = Q0 * (2.0 / 30.0);
	offset = h / 2;
	D = 1;
	H = 1;
	deltaX = PI2 / N;


	// END FENE PARAMETERS
	// ****************************************** //

	readParameters(argc, argv);


	CFD::FokkerPlanckSolver fps(N, deltaX, dt, h, offset, D, H, Q0);



	if (dtOn)
		dt = .01 / (pow(2, (log2(N) - 6)));
	Total_iterations = int(finalTime / dt + EPS);
	Record_iterations = int(saveTime / dt + EPS);
	//        Total_iterations = int(finalTime *100 + EPS);
	//       Record_iterations = int(saveTime * 50 + EPS);
	N2 = N * N;
	KK = N / 2;
	reciprocal_Wi = 1. / Wi;
	if (betaOn)
		Beta = 1. / (2. * Wi);
	if (turnfivept_on)
		fivept_on = true;

	if (fivept_on)
		fivept = 5;
	else
		fivept = 1;

	//  Allocate the matrices.

	//  tracer
	if (tracer) {
		tr = (double*) calloc((tr_pts * tr_pts * 2 * fivept), sizeof(double));
		update = (double*) calloc((tr_pts * tr_pts * 2 * fivept),
				sizeof(double));
		old_update = (double*) calloc((tr_pts * tr_pts * 2 * fivept),
				sizeof(double));

		if (fivept_on) {
			tr_fivept = (double*) calloc((tr_pts * tr_pts * 2 * fivept),
					sizeof(double));
			update_fivept = (double*) calloc((tr_pts * tr_pts * 2 * fivept),
					sizeof(double));
			old_update_fivept = (double*) calloc(
					(tr_pts * tr_pts * 2 * fivept), sizeof(double));
		}

	}

	//  Pressure.
	double *P = (double*) calloc(N2, sizeof(double));

	double *P_hat_Re = (double*) calloc(N2, sizeof(double));
	double *P_hat_Im = (double*) calloc(N2, sizeof(double));

	//  Velocity U=(u1, u2).
	double *U = (double*) calloc(N2 * 2, sizeof(double));

	double *U_hat_Re = (double*) calloc(N2 * 2, sizeof(double));
	double *U_hat_Im = (double*) calloc(N2 * 2, sizeof(double));

	//  Stress tensor S=(S1, S2; S2, S3);
	double *S = (double*) calloc(N2 * 3, sizeof(double));

	double *S_hat_Re = (double*) calloc(N2 * 3, sizeof(double));
	double *S_hat_Im = (double*) calloc(N2 * 3, sizeof(double));

	double *RHS_Re = (double*) calloc(N2 * 3, sizeof(double));
	double *RHS_Im = (double*) calloc(N2 * 3, sizeof(double));

	double *old_RHS_Re = (double*) calloc(N2 * 3, sizeof(double));
	double *old_RHS_Im = (double*) calloc(N2 * 3, sizeof(double));

	//  External force F=(f1, f2).
	double *F = (double*) calloc(N2 * 2, sizeof(double));

	double *F_hat_Re = (double*) calloc(N2 * 2, sizeof(double));
	double *F_hat_Im = (double*) calloc(N2 * 2, sizeof(double));

	//  Frequency.
	int *Freq = (int *) calloc(N, sizeof(int));
	int i;

	//  General variables.

	int iterations;
	double time;

	//  Define data files.

	char filename[132];
	char dirname[80];
	char command[132];

	/////////////////////////////////////////////////////////////////////////
	///////////////////--- STEP 1. Initialization ---////////////////////////
	/////////////////////////////////////////////////////////////////////////
	cout << " Wi = " << Wi << endl;
	cout << " N = " << N << endl;
	cout << " nu = " << NU << endl;
	cout << " dt = " << dt << endl;
	cout << " tracer = " << tracer << endl;
	cout << " pert id = " << pertid << endl;
	cout << " time dep f = " << timedep_F << endl;
	cout << " save time = " << saveTime << endl;
	cout << " final time = " << finalTime << endl;
	cout << " restart tr = " << restart_tr << endl;
	cout << " mod_off = " << mod_off << endl;
	cout << " fivept_on = " << fivept_on << endl;
	cout << " rand_tr = " << rand_tr << endl;

	iterations = 0;
	PI2 = 8. * atan((double) 1.);

	planN = fftw2d_create_plan(N, N, FFTW_FORWARD, FFTW_ESTIMATE
			| FFTW_IN_PLACE);

	sprintf(dirname, "./wi%.2f_n%d_nu%.6f_pi%d_td%d", Wi, N, NU, pertid,
			timedep_F);

	if (access(dirname, F_OK)) { /* data directory not present */
		sprintf(command, "mkdir %s", dirname);
		system(command);
	}

	if (start_w_zero) {
		time = 0;
		if (pertid)
			initialization_Spert(S);
		else
			initialization_S(S);
		if (tracer) {
			if (rand_tr)
				init_tr_rand(tr);
			else
				initialization_tr(tr);
			if (fivept_on)
				get_fivept_tr(tr, tr_fivept);

		}
	}
	else {
		time = lastsaved;
		sprintf(filename, "%s/S%3.3f.dat", dirname, lastsaved);

		Readin_S(S, filename, 3);

		sprintf(filename, "%s/RHS_re%3.3f.dat", dirname, lastsaved);
		Readin_S(RHS_Re, filename, 3);
		sprintf(filename, "%s/RHS_im%3.3f.dat", dirname, lastsaved);
		Readin_S(RHS_Im, filename, 3);
		if (tracer && !restart_tr) {

			if (turnfivept_on) {
				sprintf(filename, "%s/tr_%3.3f.dat", dirname, lastsaved);
				Readin_tr(tr, filename, 1);
				sprintf(filename, "%s/update%3.3f.dat", dirname, lastsaved);
				Readin_tr(update, filename, 1);
				get_fivept_tr(tr, tr_fivept);

			}
			else
				if (fivept_on && !turnfivept_on) {
					sprintf(filename, "%s/tr_fivept_%3.3f.dat", dirname,
							lastsaved);
					Readin_tr(tr_fivept, filename, fivept);
					sprintf(filename, "%s/update_fivept_%3.3f.dat", dirname,
							lastsaved);
					Readin_tr(update_fivept, filename, fivept);
				}
				else {
					sprintf(filename, "%s/tr_%3.3f.dat", dirname, lastsaved);
					Readin_tr(tr, filename, fivept);
					sprintf(filename, "%s/update%3.3f.dat", dirname, lastsaved);
					Readin_tr(update, filename, fivept);

				}
		}
		else {
			if (rand_tr)
				init_tr_rand(tr);
			else
				initialization_tr(tr);
			if (fivept_on)
				get_fivept_tr(tr, tr_fivept);
		}
	}

	Get_hat(S, S_hat_Re, S_hat_Im, 3);

	zeroboundary_hat(S_hat_Re, S_hat_Im, 3);

	Get_Frequency(Freq);

	if (timedep_F)
		Get_F_hat_t(F_hat_Re, F_hat_Im, time, iterations);
	else
		Get_F_hat(F_hat_Re, F_hat_Im);

	/////////////////////////////////////////////////////////////////////////
	///////////////////--- STEP 2. Updating S, U ---/////////////////////////
	/////////////////////////////////////////////////////////////////////////


	// for saving data.
	while (iterations <= Total_iterations) {

		//1.    Get P_hat and U_hat.
		Get_P_U_hat(P_hat_Re, P_hat_Im, U_hat_Re, U_hat_Im, S_hat_Re, S_hat_Im,
				F_hat_Re, F_hat_Im, Freq);

		//4.    Saving U and S.

		// get U and S and save
		if ((iterations % Record_iterations) == 0) {
			sprintf(filename, "%s/U%3.3f.dat", dirname, time);

			Get_U_S(U, U_hat_Re, U_hat_Im, 2);

			Save_U_S(U, filename, 2);
			sprintf(filename, "%s/S%3.3f.dat", dirname, time);

			Get_U_S(S, S_hat_Re, S_hat_Im, 3);

			Save_U_S(S, filename, 3);
			//cout << "bef if tracer " << endl;
			if (tracer) {
				if (fivept_on) {
					sprintf(filename, "%s/tr_fivept_%3.3f.dat", dirname, time);
					Save_tr(tr_fivept, filename, fivept);
				}
				else {
					sprintf(filename, "%s/tr_%3.3f.dat", dirname, time);
					Save_tr(tr, filename, fivept);
				}
			}

			if (iterations > 0 || start_w_zero == false) {
				sprintf(filename, "%s/RHS_re%3.3f.dat", dirname, time);
				Save_U_S(RHS_Re, filename, 3);
				sprintf(filename, "%s/RHS_im%3.3f.dat", dirname, time);
				Save_U_S(RHS_Im, filename, 3);
				if (tracer) {
					if (fivept_on) {
						sprintf(filename, "%s/update_fivept_%3.3f.dat",
								dirname, time);
						Save_tr(update_fivept, filename, fivept);
					}
					else {
						sprintf(filename, "%s/update%3.3f.dat", dirname, time);
						Save_tr(update, filename, fivept);
					}
				}
			}
		}

		//2.    Update S using Second-Order ABCN.

		bool doFENE = true;
		if(doFENE){
			cout << "Started updatePolymers()" << endl;
			fps.updatePolymersAndCalculateStressTensor(U,S);
			cout << "Finished updatePolymers()" << endl;
#ifdef FeneTiming
			cout << "Ran updatePolymersAndCalculateStressTensor" << endl;

			cout << "Calls to solver: " << CFD::Timing::callSolver;
			cout << ". Total time in method: " << CFD::Timing::callSolverTime << endl;

			cout << "Calls to stressAtPoint: " << CFD::Timing::stressAtPoint;
			cout << ". Total time in method: " << CFD::Timing::stressAtPointTime << endl;

			cout << "Total calls to updatePolymersAndCalculateStressTensor: " << CFD::Timing::callUpdatePolymers;
			cout << ". Total time in method: " << CFD::Timing::callUpdatePolymersTime << endl;
#endif
			Get_hat(S, S_hat_Re, S_hat_Im, 3);
		}
		else{
			if (iterations == 0 && start_w_zero) {
				update_S_hat(RHS_Re, RHS_Im, S_hat_Re, S_hat_Im, U_hat_Re,
						U_hat_Im, Freq);
				Euler_Step(S_hat_Re, S_hat_Im, RHS_Re, RHS_Im, NU, dt, Freq);
				if (tracer && !fivept_on) {
					Get_U_S(U, U_hat_Re, U_hat_Im, 2);
					Get_update(tr, U, U_hat_Re, U_hat_Im, update, Freq, fivept);
					Euler_tr(tr, update, fivept);
				}
				if (tracer && fivept_on) {
					Get_U_S(U, U_hat_Re, U_hat_Im, 2);
					Get_update(tr_fivept, U, U_hat_Re, U_hat_Im, update, Freq,
							fivept);
					Euler_tr(tr_fivept, update_fivept, fivept);
				}
			}
			else {
				for (i = 0; i < 3 * N2; i++) {
					old_RHS_Re[i] = RHS_Re[i];
					old_RHS_Im[i] = RHS_Im[i];
				}

				update_S_hat(RHS_Re, RHS_Im, S_hat_Re, S_hat_Im, U_hat_Re,
						U_hat_Im, Freq);

				ABCN_Step(S_hat_Re, S_hat_Im, RHS_Re, RHS_Im, old_RHS_Re,
						old_RHS_Im, NU, dt, Freq);
				if (tracer) {
					if (restart_tr) {
						if (fivept_on) {
							Get_U_S(U, U_hat_Re, U_hat_Im, 2);
							Get_update(tr_fivept, U, U_hat_Re, U_hat_Im, update,
									Freq, fivept);
							Euler_tr(tr_fivept, update_fivept, fivept);
						}
						else {
							Get_U_S(U, U_hat_Re, U_hat_Im, 2);
							Get_update(tr, U, U_hat_Re, U_hat_Im, update, Freq,
									fivept);
							Euler_tr(tr, update, fivept);
						}
						restart_tr = false;
					}
					else {
						if (fivept_on) {
							Get_U_S(U, U_hat_Re, U_hat_Im, 2);
							for (i = 0; i < (2 * tr_pts * tr_pts * fivept); i++) {
								old_update_fivept[i] = update_fivept[i];
							}

							Get_update(tr_fivept, U, U_hat_Re, U_hat_Im,
									update_fivept, Freq, fivept);

							AB_tr(tr_fivept, update_fivept, old_update_fivept,
									fivept);
						}
						else {
							Get_U_S(U, U_hat_Re, U_hat_Im, 2);
							for (i = 0; i < (2 * tr_pts * tr_pts * fivept); i++) {
								old_update[i] = update[i];
							}

							Get_update(tr, U, U_hat_Re, U_hat_Im, update, Freq,
									fivept);

							AB_tr(tr, update, old_update, fivept);
						}
					}
				}
			}
		}

		//cout << "bef if tracer3a " << endl;
		iterations += 1;
		time = time + dt;
		if (timedep_F)
			Get_F_hat_t(F_hat_Re, F_hat_Im, time, iterations);

		if ((iterations % 100) == 0)
			cout << " iterations === " << iterations << endl;
	}
	cout << "end reached" << endl;

	/////////////////////////////////////////////////////////////////////////
	///////////////////--- STEP 3. Saving U      ---/////////////////////////
	//////////////////////////////////////////////////////////////////////////

	//  FREE the assigned matrices.

	free(P);
	free(P_hat_Re);
	free(P_hat_Im);

	free(U);
	free(U_hat_Re);
	free(U_hat_Im);

	free(S);
	free(S_hat_Re);
	free(S_hat_Im);

	free(RHS_Re);
	free(RHS_Im);

	free(old_RHS_Re);
	free(old_RHS_Im);

	free(F);
	free(F_hat_Re);
	free(F_hat_Im);

	if (tracer) {
		free(tr);
		free(old_update);
		//free(tr_new);
		free(update);
	}

	fftwnd_destroy_plan(planN);
}

////////////////////////////////////////////////////////////////////////////

void Readin_S(double *S, char filename[], int k_dim) {
	FILE *fp;
	int i, j, k, index;
	double s;

	fp = fopen(filename, "r");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			index = (i * N + j) * k_dim;
			for (k = 0; k < k_dim; k++) {
				fscanf(fp, "%lf", &s);
				//fscanf(fp, "   ");
				S[index + k] = s;
			}
		}
		fscanf(fp, "\n");
	}
	fclose(fp);
}

void Readin_tr(double *tr, char filename[], int fivept) {
	FILE *fp;
	int i;
	double s;

	fp = fopen(filename, "r");
	for (i = 0; i < (tr_pts * tr_pts * 2 * fivept); i++) {
		fscanf(fp, "%lf", &s);
		tr[i] = s;
	}
	fscanf(fp, "\n");

	fclose(fp);
}

void initialization_S(double *S) {

	double mesh_size = PI2 / N;

	int i, j, index;

	double x, y;

	index = 0;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {

			x = i * mesh_size;
			y = j * mesh_size;

			S[index] = 1.0;
			S[index + 1] = 0.0;
			S[index + 2] = 1.0;

			index += 3;
		}
}

//void    initialization_Spert( double *S )
//{
//    
//    double  mesh_size = PI2/N;
//    
//    int     i, j, index;
//    
//    double x, y;
//    
//    index = 0;
//    
//    for( i = 0; i < N; i ++ )
//        for( j = 0; j < N; j ++ )
//    {
//        
//        x       = i*mesh_size;
//       y       = j*mesh_size;
//        
//        
//        S[index]        = 1. -(.05)*cos(y)*(-2*sin(x)-(3/2)*sin(2*x));
//        
//        S[index+1]      = (.02)*cos(2*y)*sin(x);
//        
//        S[index+2]      = 1. +(.06)*cos(x)*(-2*sin(y)-(3/2)*sin(2*y));
//        
//        index   += 3;
//        }
//}


void initialization_Spert(double *S) {

	double mesh_size = PI2 / N;

	int t, i, j, index, index1;

	int start_val;

	double wi_mult; // add a stronger perturbation if Wi is large because beta gets tiny so the effect on the
	// velocity gets smaller

	double x, y, rand_locx, rand_locy, nn1, nn2, scale, rr;
	fftw_complex *g1 = (fftw_complex *) calloc(N2, sizeof(fftw_complex));
	fftw_complex *g2 = (fftw_complex *) calloc(N2, sizeof(fftw_complex));

	fftw_complex *g3 = (fftw_complex *) calloc(N2, sizeof(fftw_complex));

	if (Wi < 10)
		wi_mult = 1;
	else
		wi_mult = Wi / 10.;

	start_val = 60;
	cout << "start_val = " << start_val << endl;
	srand(start_val);

	for (t = 1; t < 21; t++) {

		rand_locx = PI2 * ((double) rand() / ((double) (RAND_MAX)
				+ (double) (1)));
		rand_locy = PI2 * ((double) rand() / ((double) (RAND_MAX)
				+ (double) (1)));
		nn1 = 100. * ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		nn2 = 100. * ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		scale = wi_mult * 1. / (pow((nn1 + nn2), .5)) * (1 / (pow(2., (nn1
				+ nn2))));

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				index1 = i * N + j;
				x = (N + i) * mesh_size;
				y = j * mesh_size;
				g1[index1].re = g1[index1].re + scale * (pow((1. + sin(x
						- rand_locx)), nn1)) * (pow((1. + sin(y - rand_locy)),
								nn2));
				g1[index1].im = 0;
			}
		}
	}
	srand(start_val + 1);
	for (t = 1; t < 21; t++) {

		rand_locx = PI2 * ((double) rand() / ((double) (RAND_MAX)
				+ (double) (1)));
		rand_locy = PI2 * ((double) rand() / ((double) (RAND_MAX)
				+ (double) (1)));
		nn1 = 100. * ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		nn2 = 100. * ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		scale = wi_mult * 1. / (pow((nn1 + nn2), .5)) * (1 / (pow(2., (nn1
				+ nn2))));

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				index1 = i * N + j;
				x = (N + i) * mesh_size;
				y = j * mesh_size;
				g2[index1].re = g2[index1].re + scale * (pow((1. + sin(x
						- rand_locx)), nn1)) * (pow((1. + sin(y - rand_locy)),
								nn2));
				g2[index1].im = 0;
			}
		}
	}

	srand(start_val + 2);

	for (t = 1; t < 21; t++) {

		rand_locx = PI2 * ((double) rand() / ((double) (RAND_MAX)
				+ (double) (1)));
		rand_locy = PI2 * ((double) rand() / ((double) (RAND_MAX)
				+ (double) (1)));
		nn1 = 100. * ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		nn2 = 100. * ((double) rand() / ((double) (RAND_MAX) + (double) (1)));
		scale = wi_mult * 1. / (pow((nn1 + nn2), .5)) * (1 / (pow(2., (nn1
				+ nn2))));

		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				index1 = i * N + j;
				x = (N + i) * mesh_size;
				y = j * mesh_size;
				g3[index1].re = g3[index1].re + scale * (pow((1. + sin(x
						- rand_locx)), nn1)) * (pow((1. + sin(y - rand_locy)),
								nn2));
				g3[index1].im = 0;
			}
		}
	}

	index = 0;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {

			x = (N + i) * mesh_size;
			y = j * mesh_size;

			index1 = i * N + j;

			S[index] = 1. + g1[index1].re; //-(.05)*cos(y)*(-2*sin(x)-(3/2)*sin(2*x));
			S[index + 1] = g2[index1].re; //(.02)*cos(2*y)*sin(x);
			S[index + 2] = 1. + g3[index1].re; //+(.06)*cos(x)*(-2*sin(y)-(3/2)*sin(2*y));

			index += 3;
		}

	//      free the assigned arrays.
	free(g1);
	free(g2);
	free(g3);
}

void initialization_tr(double *tr) {
	double mesh_size = PI2 / tr_pts;

	int i, j, index;

	double x, y;

	index = 0;
	for (i = 0; i < tr_pts / 2; i++)
		for (j = 0; j < tr_pts / 2; j++) {
			x = i * mesh_size;
			y = j * mesh_size;

			tr[index] = x;
			tr[index + 1] = y;
			index += 2;
		}
	for (i = 0; i < tr_pts / 2; i++)
		for (j = tr_pts / 2; j < tr_pts; j++) {
			x = i * mesh_size;
			y = j * mesh_size;

			tr[index] = x;
			tr[index + 1] = y;
			index += 2;
		}

	for (i = tr_pts / 2; i < tr_pts; i++)
		for (j = tr_pts / 2; j < tr_pts; j++) {
			x = i * mesh_size;
			y = j * mesh_size;

			tr[index] = x;
			tr[index + 1] = y;
			index += 2;
		}
	for (i = tr_pts / 2; i < tr_pts; i++)
		for (j = 0; j < tr_pts / 2; j++) {
			x = i * mesh_size;
			y = j * mesh_size;

			tr[index] = x;
			tr[index + 1] = y;
			index += 2;
		}
}

void init_tr_rand(double *tr) {

	int i, j, index;
	double XMin, XMax, YMin, YMax;

	index = 0;
	srand(11);
	XMin = 1.;
	XMax = 1.05;
	YMin = 1.;
	YMax = 1.05;
	for (i = 0; i < tr_pts * tr_pts; i++) {
		tr[index] = XMin + rand() * (XMax - XMin) / RAND_MAX;
		//PI2*((double)rand()/((double)(RAND_MAX)+(double)(1)));
		tr[index + 1] = YMin + rand() * (YMax - YMin) / RAND_MAX;
		// PI2*((double)rand()/((double)(RAND_MAX)+(double)(1)));
		index += 2;

	}
}

void get_fivept_tr(double *tr, double *tr_fivept) {
	int i, index, index_5pt;
	double spacex, spacey;

	//spacex = .0001;
	//spacey = .0001;

	spacex = .00000001;
	spacey = .00000001;

	index_5pt = 0;
	index = 0;
	for (i = 0; i < tr_pts * tr_pts; i++) {
		index = 2 * i;
		tr_fivept[index_5pt] = tr[index]; // Cx
		tr_fivept[index_5pt + 1] = tr[index + 1]; // Cy
		tr_fivept[index_5pt + 2] = tr[index] - spacex; // Lx
		tr_fivept[index_5pt + 3] = tr[index + 1]; // Ly
		tr_fivept[index_5pt + 4] = tr[index] + spacex; // Rx
		tr_fivept[index_5pt + 5] = tr[index + 1]; // Ry
		tr_fivept[index_5pt + 6] = tr[index]; // Bx
		tr_fivept[index_5pt + 7] = tr[index + 1] - spacey; // By
		tr_fivept[index_5pt + 8] = tr[index]; // Tx
		tr_fivept[index_5pt + 9] = tr[index + 1] + spacey; // Ty


		index_5pt = index_5pt + 10;
	}

}
void Get_F_hat(double *F_hat_Re, double *F_hat_Im) {
	double mesh_size = PI2 / N;

	int i, j, k, index1, index2;
	double x, y;

	fftw_complex *F12 = (fftw_complex *) calloc(N2, sizeof(fftw_complex));

	for (k = 0; k < 2; k++) {
		//  the external force.
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				index1 = i * N + j;
				x = i * mesh_size;
				y = j * mesh_size;

				if (k == 0)
					F12[index1].re = -2.0 * sin(x) * cos(y);
				else
					F12[index1].re = 2.0 * cos(x) * sin(y);

				F12[index1].im = 0.0;
			}
		}

		//  get FFT.
		fftwnd_one(planN, F12, NULL);

		//  save them.
		index1 = 0;
		index2 = k;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				F_hat_Re[index2] = F12[index1].re / N2;
				F_hat_Im[index2] = F12[index1].im / N2;

				index1 += 1;
				index2 += 2;
			}
	}

	// Set F to be zero for the locations (KK, i) and (i, KK).
	zeroboundary_hat(F_hat_Re, F_hat_Im, 2);

	//  free the assigned arrays.
	free(F12);
}

void Get_F_hat_t(double *F_hat_Re, double *F_hat_Im, double time,
		int iterations) {
	double mesh_size = PI2 / N;

	int i, j, k, index1, index2;
	double x, y, a, b, p1, p2;

	fftw_complex *F12 = (fftw_complex *) calloc(N2, sizeof(fftw_complex));
	a = 0.05;
	b = 0.05;
	p1 = 20.;
	p2 = 21;
	if (iterations == 0) {
		cout << "forcing (a,b,p1,p2) = (" << a << "," << b << "," << p1 << ","
				<< p2 << ")" << endl;
		cout << "forcing increase (1+sin(2*pi*t/200))" << endl;
	}
	for (k = 0; k < 2; k++) {
		//  the external force.
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				index1 = i * N + j;
				x = i * mesh_size;
				y = j * mesh_size;

				if (k == 0)
					F12[index1].re = -2.0 * sin(x) * cos(y) * (1. + a * sin(
							(PI2 / p1) * time) * sin(x)) * (1 + sin(PI2 * time
									/ 200));
				else
					F12[index1].re = 2.0 * cos(x) * sin(y) * (1. + b * sin((PI2
							/ p2) * time) * sin(y)) * (1
									+ sin(PI2 * time / 200));

				F12[index1].im = 0.0;
			}
		}

		//  get FFT.
		fftwnd_one(planN, F12, NULL);

		//  save them.
		index1 = 0;
		index2 = k;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				F_hat_Re[index2] = F12[index1].re / N2;
				F_hat_Im[index2] = F12[index1].im / N2;

				index1 += 1;
				index2 += 2;
			}
	}

	// Set F to be zero for the locations (KK, i) and (i, KK).
	zeroboundary_hat(F_hat_Re, F_hat_Im, 2);

	//  free the assigned arrays.
	free(F12);
}

void Get_hat(double *S, double *S_hat_Re, double *S_hat_Im, int k_dim) {
	int k;
	int i, j, index1, index3;

	fftw_complex *in_Ri = (fftw_complex*) calloc(N2, sizeof(fftw_complex));

	for (k = 0; k < k_dim; k++) {
		index1 = 0;
		index3 = k;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				in_Ri[index1].re = S[index3];
				in_Ri[index1].im = 0;

				index1 += 1;
				index3 += k_dim;
			}

		//  get FFT.
		fftwnd_one(planN, in_Ri, NULL);

		//  save them.
		index1 = 0;
		index3 = k;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				S_hat_Re[index3] = in_Ri[index1].re / N2; // REMEMBER!!! DIVIDED BY N*N!!!
				S_hat_Im[index3] = in_Ri[index1].im / N2;

				index1 += 1;
				index3 += k_dim;
			}
	} // end for k.

	//  Free the assigned arrays.
	free(in_Ri);
}

void Get_P_U_hat(double *P_hat_Re, double *P_hat_Im, double *U_hat_Re,
		double *U_hat_Im, double *S_hat_Re, double *S_hat_Im, double *F_hat_Re,
		double *F_hat_Im, int Freq[]) {
	int i, j, k, index1, index2, index3;
	double denominator, Re_part, Im_part;

	int ki, kj;
	int ki2, kj2, kij;

	index1 = 0;
	index2 = 0;
	index3 = 0;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			ki = Freq[i];
			kj = Freq[j];

			ki2 = ki * ki;
			kj2 = kj * kj;
			kij = ki * kj;

			if ((ki == 0) && (kj == 0))
				denominator = 1.0;
			else
				denominator = ki2 + kj2;

			Re_part = Beta * (S_hat_Re[index3] * ki2 + 2 * S_hat_Re[index3 + 1]
			                                                        * kij + S_hat_Re[index3 + 2] * kj2) - ki * F_hat_Im[index2]
			                                                                                                            - kj * F_hat_Im[index2 + 1];

			Im_part = Beta * (S_hat_Im[index3] * ki2 + 2 * S_hat_Im[index3 + 1]
			                                                        * kij + S_hat_Im[index3 + 2] * kj2) + ki * F_hat_Re[index2]
			                                                                                                            + kj * F_hat_Re[index2 + 1];

			P_hat_Re[index1] = Re_part / denominator;
			P_hat_Im[index1] = Im_part / denominator;

			for (k = 0; k < 2; k++) // k = 0 -- U1,  k = 1 -- U2.
			{
				if (k == 0)
					kij = ki;
				else
					kij = kj;

				Re_part = P_hat_Im[index1] * kij - Beta * (S_hat_Im[index3 + k]
				                                                    * ki + S_hat_Im[index3 + k + 1] * kj) - F_hat_Re[index2
				                                                                                                     + k];

				Im_part = -P_hat_Re[index1] * kij + Beta
						* (S_hat_Re[index3 + k] * ki + S_hat_Re[index3 + k + 1]
						                                        * kj) - F_hat_Im[index2 + k];

				U_hat_Re[index2 + k] = Re_part / denominator;
				U_hat_Im[index2 + k] = Im_part / denominator;
			}
			index1 += 1;
			index2 += 2;
			index3 += 3;
		}

	// Set "P_hat and U_hat" to be zero for the locations (KK, i) and (i, KK).
	zeroboundary_hat(P_hat_Re, P_hat_Im, 1);

	zeroboundary_hat(U_hat_Re, U_hat_Im, 2);
}

void Get_U_S(double *U, double *U_hat_Re, double *U_hat_Im, int k_dim) {
	int i, j, k;
	int index1, index2;

	fftw_complex *in_Ri = (fftw_complex*) calloc(N2, sizeof(fftw_complex));

	//  Get U or S.
	for (k = 0; k < k_dim; k++) {
		index1 = 0;
		index2 = k;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				in_Ri[index1].re = U_hat_Re[index2];
				in_Ri[index1].im = -U_hat_Im[index2];

				index1 += 1;
				index2 += k_dim;
			}

		//  get inverse FFT .
		fftwnd_one(planN, in_Ri, NULL);

		index1 = 0;
		index2 = k;
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				U[index2] = in_Ri[index1].re;

				index1 += 1;
				index2 += k_dim;
			}

	} // for k.

	//  Free the assigned arrays.
	free(in_Ri);
}

void Get_UgradS_convolution(double *Re_UgradS_hat, double *Im_UgradS_hat,
		double *S_hat_Re, double *S_hat_Im, double *U_hat_Re, double *U_hat_Im,
		int Freq[]) {
	int i, j;
	int kU, kS;
	int index, index1, index2, indexU, indexS;

	int beg, end;

	int shift;
	int ki, kj, kij;

	double *A_re = (double*) calloc(N2, sizeof(double));
	double *A_im = (double*) calloc(N2, sizeof(double));
	double *A = (double*) calloc(N2, sizeof(double));
	double *B_re = (double*) calloc(N2, sizeof(double));
	double *B_im = (double*) calloc(N2, sizeof(double));
	double *B = (double*) calloc(N2, sizeof(double));
	double *C_re = (double*) calloc(N2, sizeof(double));
	double *C_im = (double*) calloc(N2, sizeof(double));
	double *C = (double*) calloc(N2, sizeof(double));

	beg = (int) (floor(N / 3.) + 1);
	end = (int) (N - floor(N / 3.));

	//  get a^{hat}*b^{hat} by using the aliasing technique.
	for (kU = 0; kU < 2; kU++) {
		for (kS = 0; kS < 3; kS++) {
			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++) {
					index1 = i * N + j;

					indexU = 2 * index1 + kU;
					A_re[index1] = U_hat_Re[indexU];
					A_im[index1] = U_hat_Im[indexU];
				}

			Hou_filter(A_re, A_im, 1);

			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++) {
					index1 = i * N + j;

					if (kU == 0)
						kij = i;
					else
						kij = j;

					indexS = 3 * index1 + kS;

					B_re[index1] = -S_hat_Im[indexS] * Freq[kij];
					B_im[index1] = S_hat_Re[indexS] * Freq[kij];
				}

			Hou_filter(B_re, B_im, 1);

			Get_U_S(A, A_re, A_im, 1);

			Get_U_S(B, B_re, B_im, 1);

			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++) {
					index1 = i * N + j;

					C[index1] = A[index1] * B[index1];
				}

			Get_hat(C, C_re, C_im, 1);

			zeroboundary_hat(C_re, C_im, 1);

			//  get the goal convolution sum of the original two arrays.
			shift = 3 * kU + kS;
			indexS = shift;
			index = 0;
			for (i = 0; i < N2; i++) {
				Re_UgradS_hat[indexS] = C_re[index];
				Im_UgradS_hat[indexS] = C_im[index];

				index += 1;
				indexS += 6;
			}

		} // for kS.
	} // for kU.

	//  Free the assigned arrays.
	free(A_re);
	free(A_im);
	free(A);
	free(B_re);
	free(B_im);
	free(B);
	free(C_re);
	free(C_im);
	free(C);

}

void Get_gradUS_convolution(double *Re_gradUS_hat, double *Im_gradUS_hat,
		double *S_hat_Re, double *S_hat_Im, double *U_hat_Re, double *U_hat_Im,
		int Freq[]) {
	int i, j;
	int kU, kS, kxy;
	int index, index1, index2, indexU, indexS, index8;

	int kij;

	int beg, end;

	int shift;

	double *A_re = (double*) calloc(N2, sizeof(double));
	double *A_im = (double*) calloc(N2, sizeof(double));
	double *A = (double*) calloc(N2, sizeof(double));
	double *B_re = (double*) calloc(N2, sizeof(double));
	double *B_im = (double*) calloc(N2, sizeof(double));
	double *B = (double*) calloc(N2, sizeof(double));
	double *C_re = (double*) calloc(N2, sizeof(double));
	double *C_im = (double*) calloc(N2, sizeof(double));
	double *C = (double*) calloc(N2, sizeof(double));

	beg = (int) (floor(N / 3.) + 1);
	end = (int) (N - floor(N / 3.));

	//  get a^{hat}*b^{hat} by using the aliasing technique.
	for (kxy = 0; kxy < 2; kxy++) {
		for (kU = 0; kU < 2; kU++) {
			for (kS = 0; kS < 2; kS++) {
				for (i = 0; i < N; i++)
					for (j = 0; j < N; j++) {
						index1 = i * N + j;
						indexU = index1 * 2 + kU;

						if (kxy == 0)
							kij = i;
						else
							kij = j;

						A_re[index1] = -U_hat_Im[indexU] * Freq[kij];
						A_im[index1] = U_hat_Re[indexU] * Freq[kij];
					}

				Hou_filter(A_re, A_im, 1);

				for (i = 0; i < N; i++)
					for (j = 0; j < N; j++) {
						index1 = i * N + j;

						indexS = 3 * index1 + kxy + kS;

						B_re[index1] = S_hat_Re[indexS];
						B_im[index1] = S_hat_Im[indexS];
					}

				Hou_filter(B_re, B_im, 1);

				Get_U_S(A, A_re, A_im, 1);
				Get_U_S(B, B_re, B_im, 1);

				for (i = 0; i < N; i++)
					for (j = 0; j < N; j++) {
						index1 = i * N + j;

						C[index1] = A[index1] * B[index1];
					}

				Get_hat(C, C_re, C_im, 1);
				zeroboundary_hat(C_re, C_im, 1);

				//  get the goal convolution sum of the original two arrays.
				shift = 4 * kxy + 2 * kU + kS;
				index8 = shift;
				index = 0;
				for (i = 0; i < N2; i++) {
					Re_gradUS_hat[index8] = C_re[index];
					Im_gradUS_hat[index8] = C_im[index];

					index += 1;
					index8 += 8;
				}

			} // for kS.
		} // for kU.
	} // for kxy.

	//  Free the assigned arrays.
	free(A_re);
	free(A_im);
	free(A);
	free(B_re);
	free(B_im);
	free(B);
	free(C_re);
	free(C_im);
	free(C);

}

void Hou_filter(double *Ahat_re, double *Ahat_im, int k_dim) {
	int i, j, k, index1;
	double k1, k2;
	double hou_cutoff;

	for (k = 0; k < k_dim; k++) {
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				index1 = (i * N + j) * k_dim + k;

				k1 = (double) (KK - abs(i - KK));
				k2 = (double) (KK - abs(j - KK));
				hou_cutoff = exp(-36. * pow(sqrt(k1 * k1 + k2 * k2) / (N / 2.),
						36.));
				Ahat_re[index1] = hou_cutoff * Ahat_re[index1];
				Ahat_im[index1] = hou_cutoff * Ahat_im[index1];
			}
	}

}

void Add_Shat(double *mid_Shat_Re, double *mid_Shat_Im, double *S_hat_Re,
		double *S_hat_Im) {
	int i, k;
	int index;

	for (k = 0; k < 3; k++) {
		index = k;
		for (i = 0; i < N2; i++) {
			mid_Shat_Re[index] += S_hat_Re[index];
			mid_Shat_Im[index] += S_hat_Im[index];

			index += 3;
		}
	}
}
void update_S_hat(double *Re_newS_hat, double *Im_newS_hat, double *S_hat_Re,
		double *S_hat_Im, double *U_hat_Re, double *U_hat_Im, int Freq[]) {
	int i, j;

	int index1, index3, index6, index8;
	double slope;

	double *Re_gradUS_hat = (double*) calloc(N * N * 8, sizeof(double));
	double *Im_gradUS_hat = (double*) calloc(N * N * 8, sizeof(double));

	double *Re_UgradS_hat = (double*) calloc(N * N * 6, sizeof(double));
	double *Im_UgradS_hat = (double*) calloc(N * N * 6, sizeof(double));

	Get_gradUS_convolution(Re_gradUS_hat, Im_gradUS_hat, S_hat_Re, S_hat_Im,
			U_hat_Re, U_hat_Im, Freq);

	Get_UgradS_convolution(Re_UgradS_hat, Im_UgradS_hat, S_hat_Re, S_hat_Im,
			U_hat_Re, U_hat_Im, Freq);

	index1 = 0;
	index3 = 0;
	index6 = 0;
	index8 = 0;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {

			//  For S1.
			Re_newS_hat[index3] = -(Re_UgradS_hat[index6]
			                                      + Re_UgradS_hat[index6 + 3]) + 2 * (Re_gradUS_hat[index8]
			                                                                                        + Re_gradUS_hat[index8 + 4]) - reciprocal_Wi
			                                                                                        * S_hat_Re[index3];

			Im_newS_hat[index3] = -(Im_UgradS_hat[index6]
			                                      + Im_UgradS_hat[index6 + 3]) + 2 * (Im_gradUS_hat[index8]
			                                                                                        + Im_gradUS_hat[index8 + 4]) - reciprocal_Wi
			                                                                                        * S_hat_Im[index3];

			//  For S2.
			Re_newS_hat[index3 + 1] = -(Re_UgradS_hat[index6 + 1]
			                                          + Re_UgradS_hat[index6 + 4]) + (Re_gradUS_hat[index8 + 1]
			                                                                                        + Re_gradUS_hat[index8 + 5] + Re_gradUS_hat[index8 + 2]
			                                                                                                                                    + Re_gradUS_hat[index8 + 6]) - reciprocal_Wi
			                                                                                                                                    * S_hat_Re[index3 + 1];

			Im_newS_hat[index3 + 1] = -(Im_UgradS_hat[index6 + 1]
			                                          + Im_UgradS_hat[index6 + 4]) + (Im_gradUS_hat[index8 + 1]
			                                                                                        + Im_gradUS_hat[index8 + 5] + Im_gradUS_hat[index8 + 2]
			                                                                                                                                    + Im_gradUS_hat[index8 + 6]) - reciprocal_Wi
			                                                                                                                                    * S_hat_Im[index3 + 1];

			//  For S3.
			Re_newS_hat[index3 + 2] = -(Re_UgradS_hat[index6 + 2]
			                                          + Re_UgradS_hat[index6 + 5]) + 2 * (Re_gradUS_hat[index8
			                                                                                            + 3] + Re_gradUS_hat[index8 + 7]) - reciprocal_Wi
			                                                                                            * S_hat_Re[index3 + 2];

			Im_newS_hat[index3 + 2] = -(Im_UgradS_hat[index6 + 2]
			                                          + Im_UgradS_hat[index6 + 5]) + 2 * (Im_gradUS_hat[index8
			                                                                                            + 3] + Im_gradUS_hat[index8 + 7]) - reciprocal_Wi
			                                                                                            * S_hat_Im[index3 + 2];

			index1 += 1;
			index3 += 3;
			index6 += 6;
			index8 += 8;
		}
	Re_newS_hat[0] = -(Re_UgradS_hat[0] + Re_UgradS_hat[3]) + 2
			* (Re_gradUS_hat[0] + Re_gradUS_hat[4]) - reciprocal_Wi
			* (S_hat_Re[0] - 1.0);

	Re_newS_hat[2] = -(Re_UgradS_hat[2] + Re_UgradS_hat[5]) + 2
			* (Re_gradUS_hat[3] + Re_gradUS_hat[7]) - reciprocal_Wi
			* (S_hat_Re[2] - 1.0);

	//  Free the assigned arrays.
	free(Re_gradUS_hat);
	free(Im_gradUS_hat);

	free(Re_UgradS_hat);
	free(Im_UgradS_hat);
}

void Save_U_S(double *U, char filename[], int k_dim) {
	int i, j, k, index;
	FILE *fp;

	fp = fopen(filename, "w");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			index = (i * N + j) * k_dim;
			for (k = 0; k < k_dim; k++) {
				fprintf(fp, "%25.16e", U[index + k]);
				//              fprintf(fp, "   ");
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void Save_tr(double *tr, char filename[], int fivept) {
	int i;
	FILE *fp;

	fp = fopen(filename, "w");
	for (i = 0; i < (tr_pts * tr_pts * 2 * fivept); i++) {
		fprintf(fp, "%25.16e", tr[i]);
	}
	fprintf(fp, "\n");

	fclose(fp);
}

void Euler_Step(double *S_hat_Re, double *S_hat_Im, double *RHS_Re,
		double *RHS_Im, double nnu, double ddt, int Freq[]) {
	int i, j, sl, index;
	double coef, ksq;

	for (sl = 0; sl < 3; sl++)
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				index = 3 * (i * N + j) + sl;

				ksq = Freq[i] * Freq[i] + Freq[j] * Freq[j];

				coef = ((NU) * (ksq) * (dt)) / 2;

				S_hat_Re[index] = (1 - coef) / (1 + coef) * S_hat_Re[index]
				                                                     + (1 / (1 + coef)) * dt * RHS_Re[index];

				S_hat_Im[index] = (1 - coef) / (1 + coef) * S_hat_Im[index]
				                                                     + (1 / (1 + coef)) * dt * RHS_Im[index];

			}

}

void ABCN_Step(double *S_hat_Re, double *S_hat_Im, double *RHS_Re,
		double *RHS_Im, double *old_RHS_Re, double *old_RHS_Im, double nnu,
		double ddt, int Freq[]) {
	int i, j, sl, index;
	double coef, ksq;

	for (sl = 0; sl < 3; sl++)
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++) {
				index = 3 * (i * N + j) + sl;

				ksq = Freq[i] * Freq[i] + Freq[j] * Freq[j];

				coef = ((NU) * (ksq) * (dt)) / 2;

				S_hat_Re[index] = (1 - coef) / (1 + coef) * S_hat_Re[index]
				                                                     + (1 / (1 + coef)) * (dt / 2) * (3 * RHS_Re[index]
				                                                                                                 - old_RHS_Re[index]);

				S_hat_Im[index] = (1 - coef) / (1 + coef) * S_hat_Im[index]
				                                                     + (1 / (1 + coef)) * (dt / 2) * (3 * RHS_Im[index]
				                                                                                                 - old_RHS_Im[index]);

			}

}

void Get_Frequency(int Freq[])//, double *ex_Freq1, double *ex_Freq2 )
{
	int i;

	for (i = 0; i < KK + 1; i++)
		Freq[i] = i;
	for (i = KK + 1; i < N; i++)
		Freq[i] = i - N;

}

void zeroboundary_hat(double *A_hat_Re, double *A_hat_Im, int k_dim) {
	//  This function is used to set the values of A_hat to be zero at (KK, i) and (i, KK).
	int i, k;
	int index;

	for (k = 0; k < k_dim; k++) {
		for (i = 0; i < N; i++) {
			index = (KK * N + i) * k_dim + k;

			A_hat_Re[index] = 0;
			A_hat_Im[index] = 0;

			index = (i * N + KK) * k_dim + k;

			A_hat_Re[index] = 0;
			A_hat_Im[index] = 0;
		}
	}
}

void readParameters(int argc, char **argv) {
	char **margv;
	int jjj;

	if (argc != 1) {
		margv = argv;
		margv++;
		while (*margv != NULL) {
			if (strcmp(*margv, "-n") == 0) {
				margv++;
				if (margv == NULL || (jjj = sscanf(*margv, "%d", &N)) != 1) {
					fprintf(stderr,
							"error in parameters: -n size , size is integer\n");
					exit(1);
				}
			}
			else if (strcmp(*margv, "-h") == 0){
				margv++;
				if(margv == NULL || (jjj = sscanf(*margv, "%lf", &h)) != 1){
					fprintf(stderr,"error in parameters: -h configuration space width , value is double\n");
					exit(1);
				}
				offset = h/2;
			}
			else if(strcmp(*margv, "-D") == 0){
				margv++;
				if(margv == NULL || (jjj = sscanf(*margv, "%lf", &D)) != 1){
					fprintf(stderr, "error in parameters: -D value , value is double\n");
					exit(1);
				}
			}
			else if(strcmp(*margv, "-H") == 0){
				margv++;
				if(margv == NULL || (jjj = sscanf(*margv, "%lf", &H)) != 1){
					fprintf(stderr, "error in parameters: -H value , value is double\n");
					exit(1);
				}
			}
			else if(strcmp(*margv, "-Q0") == 0){
				margv++;
				if(margv == NULL || (jjj = sscanf(*margv, "%lf", &Q0)) != 1){
					fprintf(stderr, "error in parameters: -Q0 value , value is double\n");
					exit(1);
				}
			}
			else if (strcmp(*margv, "-trpts") == 0) {
				margv++;
				if (margv == NULL || (jjj = sscanf(*margv, "%d", &tr_pts))
						!= 1) {
					fprintf(stderr,
							"error in parameters: -trpts size , size is integer\n");
					exit(1);
				}
			}
			else if (strcmp(*margv, "-wi") == 0) {
				margv++;
				if (margv == NULL || (jjj = sscanf(*margv, "%lf", &Wi))
						!= 1) {
					fprintf(stderr,
							"error in parameters: -wi value , value is double\n");
					exit(1);
				}
			}
			else if (strcmp(*margv, "-nu") == 0) {
				margv++;
				if (margv == NULL || (jjj = sscanf(*margv, "%lf",
						&NU)) != 1) {
					fprintf(stderr,
							"error in parameters: -NU value , value is double\n");
					exit(1);
				}
			}
			else if (strcmp(*margv, "-dt") == 0) {
				margv++;
				if (margv == NULL || (jjj = sscanf(*margv,
						"%lf", &dt)) != 1) {
					fprintf(stderr,
							"error in parameters: -dt value , value is double\n");
					exit(1);
				}
				dtOn = false;
			}
			else if (strcmp(*margv, "-beta") == 0) {
				margv++;
				if (margv == NULL || (jjj = sscanf(*margv,
						"%lf", &Beta)) != 1) {
					fprintf(stderr,
							"error in parameters: -beta value , value is double\n");
					exit(1);
				}
				betaOn = false;
			}
			else if (strcmp(*margv, "-pi") == 0) {
				pertid = true;

			}
			else if (strcmp(*margv, "-td") == 0) {
				timedep_F = true;
			}
			else if (strcmp(*margv, "-tr") == 0) {
				tracer = true;
			}
			else if (strcmp(*margv, "-restart")
					== 0) {
				start_w_zero = false;
			}
			else if (strcmp(*margv,
					"-restart_tr") == 0) {
				restart_tr = true;
			}
			else if (strcmp(*margv,
					"-mod_off")
					== 0) {
				mod_off = true;
			}
			else if (strcmp(*margv,
					"-rand_tr")
					== 0) {
				rand_tr = true;
			}
			else if (strcmp(
					*margv,
					"-fivept_on")
					== 0) {
				fivept_on
				= true;
			}
			else if (strcmp(
					*margv,
					"-turnfivept_on")
					== 0) {
				turnfivept_on
				= true;
			}
			else if (strcmp(
					*margv,
					"-lastsaved")
					== 0) {
				margv++;
				if (margv
						== NULL
						|| (jjj
								= sscanf(
										*margv,
										"%lf",
										&lastsaved))
										!= 1) {
					fprintf(
							stderr,
							"error in parameters: -lastsaved value , value is double\n");
					exit(
							1);
				}
			}
			else if (strcmp(
					*margv,
					"-finaltime")
					== 0) {
				margv++;
				if (margv
						== NULL
						|| (jjj
								= sscanf(
										*margv,
										"%lf",
										&finalTime))
										!= 1) {
					fprintf(
							stderr,
							"error in parameters: -finaltime value , value is double\n");
					exit(
							1);
				}
			}
			else if (strcmp(
					*margv,
					"-savetime")
					== 0) {
				margv++;
				if (margv
						== NULL
						|| (jjj
								= sscanf(
										*margv,
										"%lf",
										&saveTime))
										!= 1) {
					fprintf(
							stderr,
							"error in parameters: -savetime value , value is double\n");
					exit(
							1);
				}
			}
			else {
				fprintf(
						stderr,
						"unknown argument  %s\n",
						*margv);
				fprintf(
						stderr,
						" arguments are:\n");
				fprintf(
						stderr,
						" -n size (default 512)\n");
				fprintf(
						stderr,
						" -trpts size (default 64)\n");
				fprintf(
						stderr,
						" -wi value (default 0.1)\n");
				fprintf(
						stderr,
						" -nu value (default 0.0)\n");
				fprintf(
						stderr,
						" -dt value (default .02/(pow(2,(log2(N)-6))) )\n");
				fprintf(
						stderr,
						" -beta value (default  1./(2.*Wi)) )\n");
				fprintf(
						stderr,
						" -pi  (sets pertid = true, default false )\n");
				fprintf(
						stderr,
						" -td  (sets timedep_F = true, default false )\n");
				fprintf(
						stderr,
						" -tr  (sets tracer = true, default false )\n");
				fprintf(
						stderr,
						" -restart  (sets start_w_zero = false, default true )\n");
				fprintf(
						stderr,
						" -restart_tr  (sets restart_tr = true, default false )\n");
				fprintf(
						stderr,
						" -mod_off (sets mod_off = true, default false )\n");
				fprintf(
						stderr,
						" -rand_tr (sets rand_tr = true, default false )\n");
				fprintf(
						stderr,
						" -fivept_on (sets fivept_on = true, default false )\n");
				fprintf(
						stderr,
						" -turnfivept_on (sets turnfivept_on = true, default false )\n");
				fprintf(
						stderr,
						" -lastsaved value (default 1.0)\n");
				fprintf(
						stderr,
						" -finaltime value (default 4.0)\n");
				fprintf(
						stderr,
						" -savetime value (default 1.0)\n");

				exit(
						1);
			}

			margv++;
		}
	}

	return;
}

void bcucof(double y[], double y1[], double y2[], double y12[], double d1,
		double d2, double **c) {
	static int wt[16][16] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 3, 0, 0,
			0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0,
			1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0,
			3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0,
			1, 0, 0, 1, -3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, 9, -9, 9, -9, 6, 3,
			-3, -6, 6, -6, -3, 3, 4, 2, 1, 2, -6, 6, -6, 6, -4, -2, 2, 4, -3,
			3, 3, -3, -2, -1, -1, -2, 2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, -6, 6,
			-6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1, 4, -4, 4, -4, 2,
			2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1 };
	int l, k, j, i;
	double xx, d1d2, cl[16], x[16];

	d1d2 = d1 * d2;
	for (i = 1; i <= 4; i++) {
		x[i - 1] = y[i - 1];
		x[i + 3] = y1[i - 1] * d1;
		x[i + 7] = y2[i - 1] * d2;
		x[i + 11] = y12[i - 1] * d1d2;
	}
	for (i = 0; i <= 15; i++) {
		xx = 0.0;
		for (k = 0; k <= 15; k++)
			xx += wt[i][k] * x[k];
		cl[i] = xx;
	}
	l = 0;
	for (i = 1; i <= 4; i++)
		for (j = 1; j <= 4; j++)
			c[i][j] = cl[l++];
}

void bcuint(double y[], double y1[], double y2[], double y12[], double x1l,
		double x1u, double x2l, double x2u, double x1, double x2, double *ansy) {
	void bcucof(double y[], double y1[], double y2[], double y12[], double d1,
			double d2, double **c);
	double **matrix(int nrl, int nrh, int ncl, int nch);
	void free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
	int i, j;
	double t, u, d1, d2, **c;
	double h;

	h = PI2 / N;
	c = matrix(1, 4, 1, 4);
	d1 = h;
	d2 = h;

	bcucof(y, y1, y2, y12, d1, d2, c);
	for (i = 1; i <= 4; i++)
		for (j = 1; j <= 4; j++)
			t = (x1 - x1l) / d1;
	u = (x2 - x2l) / d2;
	//cout<<"before"<<endl;
	*ansy = 0.0;
	for (i = 4; i >= 1; i--) {
		*ansy = t * (*ansy) + ((c[i][4] * u + c[i][3]) * u + c[i][2]) * u
				+ c[i][1];
	}
	// cout<<"bef fm"<< endl;
	free_matrix(c, 1, 4, 1, 4);

}
double **matrix(int nrl, int nrh, int ncl, int nch) {
	int i;
	double **m;

	m = (double **) malloc((unsigned) (nrh - nrl + 1) * sizeof(double*));
	m -= nrl;

	for (i = nrl; i <= nrh; i++) {
		m[i] = (double *) malloc((unsigned) (nch - ncl + 1) * sizeof(double));
		m[i] -= ncl;
	}
	return m;

}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch) {
	int i;

	for (i = nrh; i >= nrl; i--)
		free((char*) (m[i] + ncl));
	//cout<<"i = "<<i<<endl;
	free((char*) (m + nrl));
}

void Euler_tr(double *tr, double *update, int fivept) {
	int i, len, t;

	len = (fivept * 2 * tr_pts * tr_pts);
	for (i = 0; i < len; i++) {
		tr[i] = tr[i] + dt * update[i];
		// cout<<update[i]<<endl;
	}

	if (!mod_off) {
		for (t = 0; t < len; t++) {
			if (tr[t] < 0.)
				tr[t] = tr[t] + PI2;
			if (tr[t] >= PI2)
				tr[t] = tr[t] - PI2;
		}
	}

}

void AB_tr(double *tr, double *update, double *old_update, int fivept) {
	int i, len, t;

	len = (fivept * 2 * tr_pts * tr_pts);
	for (i = 0; i < len; i++) {
		tr[i] = tr[i] + (.5) * dt * (3 * update[i] - old_update[i]);
	}
	if (!mod_off) {
		for (t = 0; t < len; t++) {
			if (tr[t] < 0.)
				tr[t] = tr[t] + PI2;
			if (tr[t] >= PI2)
				tr[t] = tr[t] - PI2;
		}
	}
}

void Get_mixd2U(double *gradU_hat_re, double *gradU_hat_im, double *mixd2U,
		int Freq[]) {
	int i, j, index;
	double *grad2U_hat_re = (double*) calloc(N2 * 8, sizeof(double));
	double *grad2U_hat_im = (double*) calloc(N2 * 8, sizeof(double));
	double *grad2U = (double*) calloc(N2 * 8, sizeof(double));

	Get_grad_hat(gradU_hat_re, gradU_hat_im, grad2U_hat_re, grad2U_hat_im,
			Freq, 4);
	Get_U_S(grad2U, grad2U_hat_re, grad2U_hat_im, 8);
	index = 0;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			mixd2U[index] = grad2U[4 * (i * N + j) + 1];
			mixd2U[index + 1] = grad2U[4 * (i * N + j) + 5];

			index = index + 2;
		}

	free(grad2U_hat_re);
	free(grad2U_hat_im);
	free(grad2U);
}

void Get_grad_hat(double *A_hat_re, double *A_hat_im, double *gradA_hat_re,
		double *gradA_hat_im, int Freq[], int k_dim) {
	int i, j;
	int kxy, kA;
	int index, indexA, indexgA;

	int kij, l;

	l = 0;

	for (kA = 0; kA < k_dim; kA++) {
		for (kxy = 0; kxy < 2; kxy++) {
			for (i = 0; i < N; i++)
				for (j = 0; j < N; j++) {
					if (kxy == 0)
						kij = i;
					else
						kij = j;
					index = (i * N + j);
					indexA = k_dim * index + kA;
					indexgA = 2 * k_dim * index + l;

					gradA_hat_re[indexgA] = -A_hat_im[indexA] * Freq[kij];
					gradA_hat_im[indexgA] = A_hat_re[indexA] * Freq[kij];

				}
			l = l + 1;
		} // for kxy.
	} // for kA.
}

void Get_update(double *tr, double *U, double *U_hat_Re, double *U_hat_Im,
		double *update, int Freq[], int fivept) {

	double *tr_mod = (double*) calloc((tr_pts * tr_pts * 2 * fivept),
			sizeof(double));
	int len, t, kx, ky, kxh, kyh, index;
	double h, xh, yh, x0, y0;
	int i;
	double *ansy = (double*) calloc(1, sizeof(double));
	double x1l, x1u, x2l, x2u, x1, x2;
	double d1, d2, **c;
	double *y = (double*) calloc(4, sizeof(double));
	double *y1 = (double*) calloc(4, sizeof(double));
	double *y2 = (double*) calloc(4, sizeof(double));
	double *y12 = (double*) calloc(4, sizeof(double));
	// double  *U_hat_Re   = (double*) calloc( N2*2, sizeof(double) );
	//double  *U_hat_Im   = (double*) calloc( N2*2, sizeof(double) );
	double *gradU = (double*) calloc(N2 * 4, sizeof(double));
	double *gradU_hat_re = (double*) calloc(N2 * 4, sizeof(double));
	double *gradU_hat_im = (double*) calloc(N2 * 4, sizeof(double));
	double *mixd2U = (double*) calloc(N2 * 2, sizeof(double));

	len = (tr_pts * tr_pts * 2 * fivept);
	h = PI2 / N;
	Get_grad_hat(U_hat_Re, U_hat_Im, gradU_hat_re, gradU_hat_im, Freq, 2);
	Get_U_S(gradU, gradU_hat_re, gradU_hat_im, 4);
	Get_mixd2U(gradU_hat_re, gradU_hat_im, mixd2U, Freq);

	index = 0;

	if (mod_off) {

		for (t = 0; t < len; t++) {
			tr_mod[t] = fmod(tr[t], PI2);
			// fmod misbehaves for negative numbers!  Force it positive again.
			if (tr_mod[t] < 0)
				tr_mod[t] += PI2;
		}
	}
	else {
		for (t = 0; t < len; t++) {
			tr_mod[t] = tr[t];
		}
	}

	for (t = 0; t < len / 2; t++) {

		x0 = h * floor(tr_mod[index] / h + EPS);
		y0 = h * floor(tr_mod[index + 1] / h + EPS);

		//if(x0==PI2) x0=0.0;
		//if(y0==PI2) y0=0.0;
		xh = x0 + h;
		yh = y0 + h;

#ifdef BECCA_DEBUG        
		cerr << "index = " << index << " ";
		cerr << "tr_mod[index] = " << tr_mod[index] << " ";
		cerr << "tr[index] = " << tr[index] << " ";
		cerr << "h = " << h << endl;
#endif
		kx = int(floor(tr_mod[index] / h + EPS));
		ky = int(floor(tr_mod[index + 1] / h + EPS));

		kx = kx % N;
		ky = ky % N;

		if (kx == N - 1)
			kxh = 0;
		else
			kxh = kx + 1;
		if (ky == N - 1)
			kyh = 0;
		else
			kyh = ky + 1;

		//for u
		y[0] = U[UXY(kx,ky)]; //u(x0,y0)
		y1[0] = gradU[4 * (kx * N + ky)]; //ux(x0,y0)
		y2[0] = gradU[4 * (kx * N + ky) + 1]; //uy(x0,y0)
		y12[0] = mixd2U[2 * (kx * N + ky)]; //uxy(x0,y0)

		y[1] = U[2 * (kxh * N + ky)]; //u(x0+h,y0)
		y1[1] = gradU[4 * (kxh * N + ky)]; //etc
		y2[1] = gradU[4 * (kxh * N + ky) + 1];
		y12[1] = mixd2U[2 * (kxh * N + ky)];

		y[2] = U[2 * (kxh * N + kyh)]; //u(x0+h,y0+h)
		y1[2] = gradU[4 * (kxh * N + kyh)];
		y2[2] = gradU[4 * (kxh * N + kyh) + 1];
		y12[2] = mixd2U[2 * (kxh * N + kyh)];

		y[3] = U[2 * (kx * N + kyh)]; //u(x0,y0+h)
		y1[3] = gradU[4 * (kx * N + kyh)];
		y2[3] = gradU[4 * (kx * N + kyh) + 1];
		y12[3] = mixd2U[2 * (kx * N + kyh)];

		bcuint(y, y1, y2, y12, x0, xh, y0, yh, tr_mod[index],
				tr_mod[index + 1], ansy);

		update[index] = *ansy;

		// now for v:
		y[0] = U[2 * (kx * N + ky) + 1]; //v(x0,y0)
		y1[0] = gradU[4 * (kx * N + ky) + 2]; //vx(x0,y0)
		y2[0] = gradU[4 * (kx * N + ky) + 3]; //vy(x0,y0)
		y12[0] = mixd2U[2 * (kx * N + ky) + 1]; //vxy(x0,y0)

		y[1] = U[2 * (kxh * N + ky) + 1]; //v(x0+h,y0)
		y1[1] = gradU[4 * (kxh * N + ky) + 2]; //etc
		y2[1] = gradU[4 * (kxh * N + ky) + 3];
		y12[1] = mixd2U[2 * (kxh * N + ky) + 1];

		y[2] = U[2 * (kxh * N + kyh) + 1]; //v(x0+h,y0+h)
		y1[2] = gradU[4 * (kxh * N + kyh) + 2];
		y2[2] = gradU[4 * (kxh * N + kyh) + 3];
		y12[2] = mixd2U[2 * (kxh * N + kyh) + 1];

		y[3] = U[2 * (kx * N + kyh) + 1]; //v(x0,y0+h)
		y1[3] = gradU[4 * (kx * N + kyh) + 2];
		y2[3] = gradU[4 * (kx * N + kyh) + 3];
		y12[3] = mixd2U[2 * (kx * N + kyh) + 1];

		bcuint(y, y1, y2, y12, x0, xh, y0, yh, tr_mod[index],
				tr_mod[index + 1], ansy);
		update[index + 1] = *ansy;

		index += 2;

	}
	free(tr_mod);
	free(gradU);
	free(gradU_hat_re);
	free(gradU_hat_im);
	free(mixd2U);

}

