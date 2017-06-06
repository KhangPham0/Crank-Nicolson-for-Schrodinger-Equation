#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string.h>
#include <array>
#include <iostream>
#include <complex>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <string>
#include <iomanip> 
#include <time.h>

const double PI{ M_PI };
const double h_bar{ 1.0 };
const double m1{ 1.0 };//{ .5 };
const double m2{ 1.0 };//{ 5. };
const double sig{ .05 };
const int grid_point{ 101 };

//variables for harmonic potential well
const double K1{ 1.0E6 };
const double K2{ 1.0E6 };
const double omega1{ std::sqrt(K1 / m1) };
const double omega2{ std::sqrt(K2 / m2) };

//void Probability(int);

typedef std::complex<double> compx;
const compx I{ 0.,1. };
const compx one{ 1.,0. };
const compx two{ 2.,0. };

const double k1 = 110.;
const double k2 = -110.;
const double dt = 4.2E-06;
const double dx = .005;
const double X01 = .35*grid_point*dx;
const double X02 = .7*grid_point*dx;
const double V_0 = 40000.;
const double al = .062;
double w[3];
double pot[grid_point][grid_point];
const int N1 = grid_point - 2, N2 = grid_point - 2;
const int Nt = 3000;// 40000;
const int spacing = 10; //file spacing

compx value[grid_point][grid_point][2];
compx tmp[grid_point][grid_point];
std::array<double, grid_point> rho_1;
std::array<double, grid_point> rho_2;
std::array<compx, grid_point - 1> alpha;
std::array<compx, grid_point> beta;
//compx alpha[grid_point - 1];
//compx beta[grid_point];
const compx a = -one / (two*m1), b = -one / (two*m2);
const double r = dt / (dx*dx);
const compx A = (I*r*a / two), B = (I*r*b / two), C = I*dt / two;
const compx mid1 = one - two*A;
double Rho[grid_point][grid_point];

std::ofstream file1, file2, file3;

double harm_coor(int index) {
	return (-.5*grid_point*dx) + (index*dx);
}

double coor(int index) {
	return index*dx;
}

void Output(int n) {
	file1.open("x_" + std::to_string(n) + ".dat");
	file2.open("PsiSquared_" + std::to_string(n) + ".dat");
	for (int i = 1; i <= N1; i++) {
		if (rho_1[i] < 1.e-20)rho_1[i] = 1.e-20;
		if (rho_2[i] < 1.e-20)rho_2[i] = 1.e-20;
		file1 << i << "\t" << rho_1[i] << "\t" << rho_2[i] << "\n";
		for (int j = 1; j <= N2; j++) {
			if (Rho[i][j] < 1.e-20) Rho[i][j] = 1.e-20;
			file2 << i << "\t" << j << "\t" << Rho[i][j] << "\n";
		}
	}
	file1.close();
	file2.close();
}




compx R(int dir, int x1, int x2) {
	if (dir == 0) {
		return (one - C*pot[x1][x2])*value[x1][x2][dir] - B*(value[x1][x2 + 1][dir] - two * value[x1][x2][dir] + value[x1][x2 - 1][dir]);
	}
	else {
		return value[x1][x2][dir] - A*(value[x1 + 1][x2][dir] - two * value[x1][x2][dir] + value[x1 - 1][x2][dir]);
	}
}

compx mid2(int x1, int x2) {
	return one - two*B + C*pot[x1][x2];
}

void Probability(int n) {
	double Ptot = 0.;
	int p = 1, k;
	for (int i = 1; i <= N1; i++) {
		k = 1;
		if (p == 3) p = 1;
		for (int j = 1; j <= N2; j++) {
			if (k == 3) k = 1;
			Ptot = Ptot + w[k] * w[p] * abs(value[i][j][0])*abs(value[i][j][0]);
			k += 1;
		}
		p += 1;
	}

	Rho[0][0] = 0.; Rho[N1 + 1][N2 + 1] = 0.;
	for (int i = 1; i <= N1; i++) {
		for (int j = 1; j <= N2; j++) {
			Rho[i][j] = abs(value[i][j][0])*abs(value[i][j][0]) / Ptot;
		}
	}

	rho_1.fill(0.);
	rho_2.fill(0.);
	for (int i = 1; i <= N1; i++) {
		k = 1;
		for (int j = 1; j <= N2; j++) {
			if (k == 3) k = 1;
			rho_1[i] = rho_1[i] + w[k] * Rho[i][j];
			rho_2[i] = rho_2[i] + w[k] * Rho[j][i];
			k += 1;
		}
	}
	if (n % spacing == 0)
		Output(n);
}




//Use this when boudaries are being forced to equal to 0
void SolveSE() {

	alpha.fill(0.);
	beta.fill(0.);
	for (int n = 1; n <= Nt; n++) {
		std::cout << n << "\n";
		//Solving for the x1 direction
		for (int x2 = 1; x2 <= N2; x2++) {

			alpha[1] = A / mid1;

			beta[1] = R(0, 1, x2) / mid1;

			for (int l = 2; l < N2 - 1; l++) {
				alpha[l] = A / (mid1 - A*alpha[l - 1]);

				beta[l] = (R(0, l, x2) - A*beta[l - 1]) / (mid1 - A*alpha[l - 1]);
			}

			beta[N2] = (R(0, N2, x2) - A*beta[N2 - 1]) / (mid1 - A*alpha[N2 - 1]);

			//Backward run
			value[N1][x2][1] = beta[N2];
			for (int l = N1 - 1; l >= 1; l--) {
				value[l][x2][1] = beta[l] - alpha[l] * value[l + 1][x2][1];
			}

		}


		//Solving for the x2 direction
		for (int x1 = 1; x1 < N1; x1++) {

			alpha[1] = B / mid2(x1, 1);

			beta[1] = R(1, x1, 1) / mid2(x1, 1);

			for (int m = 2; m < N1 - 1; m++) {
				alpha[m] = B / (mid2(x1, m) - B*alpha[m - 1]);

				beta[m] = (R(1, x1, m) - B*beta[m - 1]) / (mid2(x1, m) - B*alpha[m - 1]);
			}

			beta[N1] = (R(1, x1, N1) - B*beta[N1 - 1]) / (mid2(x1, N1) - B*alpha[N1 - 1]);

			//Backward run
			value[x1][N2][0] = beta[N1];
			for (int m = N2 - 1; m >= 1; m--) {
				value[x1][m][0] = beta[m] - alpha[m] * value[x1][m + 1][0];
			}
		}

		Probability(n);

	}
}

void potential() {
	file3.open("potential.dat");
	//double con = 0.0; //change this to 1 to turn on the couple harmonic potential, at 0 this is just a harmonic well
	//double x1 = 0., harm_x1 = -0.5*grid_point*dx;
	for (int i = 0; i < grid_point; i++) {
		//x1 += dx;
		//harm_x1 += dx;
		//double x2 = 0., harm_x2 = -0.5*grid_point*dx;
		for (int j = 0; j < grid_point; j++) {
			//x2 += dx;
			//harm_x2 += dx;
			//harmonic well
			pot[i][j] = .5*m1*omega1*omega1*harm_coor(i)*harm_coor(i) +.5*m1*omega2*omega2*harm_coor(j)*harm_coor(j);

			//square well
			//if (abs(i - j)*dx <= al) pot[i][j] = V_0;
			//else pot[i][j] = 0.0;

			file3 << i << "\t" << j << "\t" << pot[i][j] << "\n";
		}
	}
	file3.close();
}


void Initialize() {
	compx x1_psi0, x2_psi0, x1_psi1, x2_psi1, x1_psi2, x2_psi2, psi00, psi10, psi12;
	double cons = 1. / (4.*sig*sig);
	double A1 = m1*omega1 / (PI*h_bar), A2 = m2*omega2 / (PI*h_bar);
	//double x1 = 0., harm_x1 = -.5*grid_point*dx;
	for (int l = 1; l <= N1; l++) {
		//x1 += dx;
		//harm_x1 += dx;
		//double x2 = 0., harm_x2 = -.5*grid_point*dx;
		double expt1 = sqrt(A1*PI)*harm_coor(l);
		x1_psi0 = pow(A1, .25)*exp(-expt1*expt1*0.5);
		x1_psi1 = pow(A1, .25)*sqrt(2.)*expt1*exp(-expt1*expt1*0.5);
		x1_psi2 = pow(A1, .25)*(1. / sqrt(2.))*(2.*expt1*expt1 - 1.)*exp(-expt1*expt1*0.5);
		for (int m = 0; m <= N2; m++) {
			//x2 += dx;
			//harm_x2 += dx;
			double expt2 = sqrt(A2*PI)*harm_coor(m);
			x2_psi0 = pow(A2, .25)*exp(-expt2*expt2*0.5);
			x2_psi1 = pow(A2, .25)*sqrt(2.)*expt2*exp(-expt2*expt2*0.5);
			x2_psi2 = pow(A2, .25)*(1. / sqrt(2.))*(2.*expt2*expt2 - 1.)*exp(-expt2*expt2*0.5);
			psi00 = x1_psi0*x2_psi0;
			psi10 = x1_psi1*x2_psi0;
			psi12 = x1_psi1*x2_psi2;
			value[l][m][0] = psi00+psi10+psi12+x1_psi0*x2_psi2; //harmonic potential
			//value[l][m][0] = exp(I*(k1*coor(l) + k2*coor(m)))*exp(-cons*((coor(l) - X01)*(coor(l) - X01) + (coor(m) - X02)*(coor(m) - X02))); //square well
		}
	}
	for (int j = 0; j <= N1 + 1; j++) {
		value[N1 + 1][j][0] = 0.;
		value[0][j][0] = 0.;
	}
	for (int i = 1; i <= N2; i++) {
		value[i][N2 + 1][0] = 0.;
		value[i][0][0] = 0.;
	}
	Probability(0);
}

int main() {
	time_t timei, timef;
	timei = time(NULL); //set a timer to see how long it takes to execute the program
	w[0] = dx / 3.;
	w[1] = 4.*dx / 3.;
	w[2] = 2.*dx / 3.;

	potential();
	Initialize();
	SolveSE();
	timef = time(NULL);
	std::cout << timef - timei << "\n"; //timer stops
	getchar();
}