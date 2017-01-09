#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <vector>
#include <string.h>
#include <assert.h>
#include <array>
#include <iostream>
#include <complex>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>

#define PI M_PI
#define scale 1.5/.5
#define m1 .5
#define m2 5.
#define L .5//1.5
#define sig .05
#define grid_point 101

typedef std::complex<double> compx;
#define I compx(0.,1.)
#define one compx(1.,0.)
#define two compx(2.,0.)

class wave_function {
public:
	wave_function();
	wave_function(bool);
	double k1 = 110.;
	double k2 = -110.;
	double x_01 = .25*grid_point;
	double x_02 = .6*grid_point;
	double dt = 4.2E-6;
	double dx = L / grid_point;
	double V_0 = 80000.;
	double al = .062;
	double totSteps = 1000 * dt;
	double w[3]; /////////simpson's weight
	std::vector<std::vector<compx>> value;
	std::vector<compx> density_x1;
	std::vector<compx> density_x2;
	void solve_triag();
	double potential(int, int);
	double real_space(int);
	compx sec_space_deriv(char, int, int);
	void rho();
	void normalize();
};

wave_function::wave_function() {
	value.resize(grid_point);
	density_x1.resize(grid_point);
	density_x2.resize(grid_point);
	for (int l = 0; l < grid_point; l++) {
		value[l].resize(grid_point);
		density_x1[l] = 0;
		density_x2[l] = 0;
		for (int m = 0; m < grid_point; m++)
			value[l][m] = compx(0., 0.);
	}
}

wave_function::wave_function(bool init) {
	if (init) {
		w[0] = dx / 3.;
		w[1] = 4.*dx / 3.;
		w[2] = 2.*dx / 3.;
		value.resize(grid_point);
		density_x1.resize(grid_point);
		density_x2.resize(grid_point);
		for (int l = 0; l < grid_point; l++) {
			value[l].resize(grid_point);
			density_x1[l] = 0.;
			density_x2[l] = 0.;
			for (int m = 0; m < grid_point; m++) {
				value[l][m] = compx(0., 0.);
			}
		}
		for (int l = 1; l < grid_point - 1; l++) {
			for (int m = 1; m < grid_point - 1; m++) {
				value[l][m] = exp(I*compx(k1, 0)*compx(real_space(l), 0))*exp(-pow(real_space(l) - real_space(x_01), 2.) / (4.*sig*sig))
					*exp(I*compx(k2, 0)*compx(real_space(m), 0))*exp(-pow(real_space(m) - real_space(x_02), 2.) / (4.*sig*sig));
			}
		}
		normalize();
	}
	else if (!init) {
		value.resize(grid_point - 1);
		for (int l = 0; l < grid_point - 1; l++) {
			value[l].resize(grid_point - 1);
			for (int m = 0; m < grid_point - 1; m++)
				value[l][m] = compx(0., 0.);
		}
	}
}

//Use this when boudaries are being forced to equal to 0
void wave_function::solve_triag() {

	compx a = -one / (two*m1), b = -one / (two*m2);
	double r = dt / (dx*dx);
	compx A = I*r*a / two, B = I*r*b / two, C = I*dt / two, mid = one - two*A;
	std::vector<compx> alpha(grid_point - 1);
	std::vector<compx> beta(grid_point);
	wave_function tmp; //hold the n* values

	for (int a = 0; a < grid_point; a++) {
		for (int b = 0; b < grid_point; b++) {
			tmp.value[a][b] = value[a][b];
		}
	}

	alpha[0] = 0.0, alpha[grid_point - 2] = 0.0;
	beta[0] = 0.0, beta[grid_point - 1] = 0.0;

	//std::cout << "copied" << std::endl;
	//Solving for the x1 direction first
	for (int x2 = 1; x2 < grid_point - 1; x2++) {

		alpha[1] = A / mid;

		beta[1] = ((one - C*potential(1, x2))*value[1][x2] - B*sec_space_deriv('y', 1, x2)) / mid;//((one - D)*value[0][0] - E) / F;

																								  //std::cout << "beginning" << std::endl;
																								  //Forward run
		for (int l = 2; l < grid_point - 2; l++) {
			alpha[l] = A / (mid - A*alpha[l - 1]);

			beta[l] = ((one - C*potential(l, x2)) * value[l][x2] - B*sec_space_deriv('y', l, x2)
				- A*beta[l - 1]) / (mid - A*alpha[l - 1]);
			//std::cout << i << std::endl;
		}
		//std::cout << "first loop" << std::endl;

		beta[grid_point - 2] = ((one - C*potential(grid_point - 2, x2)) * value[grid_point - 2][x2] - B*sec_space_deriv('y', grid_point - 2, x2)
			- A*beta[grid_point - 3]) / (mid - A*alpha[grid_point - 3]);

		//Backward run
		tmp.value[grid_point - 2][x2] = beta[grid_point - 2];
		for (int l = grid_point - 3; l >= 1; l--) { //-2 because -1 is calculated 1 line above
			tmp.value[l][x2] = beta[l] - alpha[l] * tmp.value[l + 1][x2];
		}
		//std::cout << y << "   The end" << std::endl;
	}
	tmp.normalize();
	//std::cout << "x direction done" << std::endl;
	//Solving for the x2 direction
	for (int x1 = 1; x1 < grid_point - 1; x1++) {

		alpha[1] = B / (one - two*B + C*potential(x1, 1));

		beta[1] = (tmp.value[x1][1] - A*tmp.sec_space_deriv('x', x1, 1)) / (one - two*B + C*potential(x1, 1));
		//std::cout << x << "   beginning2" << std::endl;
		//Forward run
		for (int m = 2; m < grid_point - 2; m++) {
			alpha[m] = B / (one - two * B + C*potential(x1, m) - B*alpha[m - 1]);
			//else
			//	alpha[j] = compx(0., 0.);
			beta[m] = (tmp.value[x1][m] - A*tmp.sec_space_deriv('x', x1, m) - B*beta[m - 1]) / (one - two*B + C*potential(x1, m) - B*alpha[m - 1]);
			//std::cout << j << std::endl;
		}
		//std::cout << "first loop2" << std::endl;

		beta[grid_point - 2] = (tmp.value[x1][grid_point - 2] - A*tmp.sec_space_deriv('x', x1, grid_point - 2) - B*beta[grid_point - 3]) /
			(one - two*B + C*potential(x1, grid_point - 2) - B*alpha[grid_point - 3]);

		//Backward run
		value[x1][grid_point - 2] = beta[grid_point - 2];
		for (double m = grid_point - 3; m >= 1; m--) { //-2 because of the same reason above
			value[x1][m] = beta[m] - alpha[m] * value[x1][m + 1];
		}
		//std::cout << x << "   The end2" << std::endl;
	}
	normalize();
}

double wave_function::potential(int x1, int x2) {

	//return .5*m1*real_space(x1)*real_space(x1) + .5*m1*9.*pow(real_space(x2) - real_space(x1), 2.) + .5*m2*real_space(x2)*real_space(x2);
	//square well
	if (al - abs(real_space(x1) - real_space(x2)) > 0.)
		return V_0;
	else
		return 0.0;

	//gaussian
	//return V_0*exp(-pow(abs(real_space(x1)-real_space(x2)),2.)/(2.*al*al));
}

double wave_function::real_space(int index) {
	return scale*(index*dx);
}

compx wave_function::sec_space_deriv(char x1_or_x2, int l, int m) {
	if (x1_or_x2 == 'x') {
		return value[l + 1][m] - two * value[l][m] + value[l - 1][m];
	}
	else if (x1_or_x2 == 'y') {
		return value[l][m + 1] - two * value[l][m] + value[l][m - 1];
	}
}

void wave_function::rho() {
	for (int k = 0; k < grid_point; k++) {
		density_x1[k] = 0.;
		density_x2[k] = 0.;
	}

	int k;
	for (int i = 0; i < grid_point; i++) {
		k = 1;
		for (int j = 0; j < grid_point; j++) {
			if (k == 3) k = 1;
			density_x1[i] += w[k]*pow(abs(value[i][j]), 2.);
			density_x2[i] += w[k]*pow(abs(value[j][i]), 2.);
			k += 1;
		}
	}
}
void wave_function::normalize() {
	compx sum = 0;
	for (int i1 = 0; i1 < grid_point; i1++) {
		for (int i2 = 0; i2 < grid_point; i2++) {
			sum += pow(abs(value[i1][i2]), 2.0);
		}
	}
	sum *= dx * dx;
	compx amplitude = one / pow(sum, .5);
	for (int i1 = 0; i1 < grid_point; i1++) {
		for (int i2 = 0; i2 < grid_point; i2++) {
			value[i1][i2] *= amplitude;
		}
	}
}

int main() {
	double check = 0.0;
	wave_function v(true);
	std::ofstream file1; //density of both particles
	std::ofstream file2; //density individually
	std::ofstream file5; //potential
	file5.open("potential.dat");
	for (int i = 0; i < grid_point; i++)
		for (int j = 0; j < grid_point; j++)
			file5 << v.real_space(i) << "\t" << v.real_space(j) << "\t" << v.potential(i, j) << std::endl;
	int index = 0;
	for (double k = v.dt; k <= v.totSteps; k += v.dt) {
		std::cout << index << std::endl;
		v.rho();
		if (index % 10 == 0) {
			file1.open("data_" + std::to_string(index) + ".dat");
			file1 << "Time   " << k - v.dt << std::endl
				<< "x" << "\t" << "y" << "\t" << "imag" << "\t" << "real" << "\t" << "abs" << std::endl;
			file2.open("x_" + std::to_string(index) + ".dat");
			file2 << "Time   " << k - v.dt << std::endl
				<< "Coordinate" << "\t" << "x1" << "\t" << "x2" << std::endl;
			//file3.open("x2_" + to_string_with_precision(k - v.dt, 6) + ".dat");
			for (int i = 0; i < grid_point; i++) {
				for (int j = 0; j < grid_point; j++) {
					file1 << v.real_space(i) << "\t"
						<< v.real_space(j) << "\t"
						<< imag(v.value[i][j]) << "\t"
						<< real(v.value[i][j]) << "\t"
						<< abs(v.value[i][j]) << std::endl;
				}
			}
			for (int k = 0; k < grid_point; k++) {
				file2 << v.real_space(k) << "\t"
					<< abs(v.density_x1[k]) << "\t"
					<< abs(v.density_x2[k]) << std::endl;
			}
			file1.close();
			file2.close();

		}
		v.solve_triag();
		index++;
	}
	getchar();
}