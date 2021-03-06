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
#include <math.h>
#include <iomanip> //this is for setprecision

#define PI M_PI
#define h_bar 1.0
#define M 1.0 
#define L .5
#define sig .05
#define grid_point 100
#define K1 1.0E6
#define K2 1.0E6
#define omega1 std::sqrt(K1/M)
#define omega2 std::sqrt(K2/M)


typedef std::complex<double> compx;
#define I compx(0.,1.)
#define one compx(1.,0.)
#define two compx(2.,0.)

class wave_function {
public:
	wave_function();
	wave_function(bool);
	double dt = 4.2E-06;
	double dx = L / grid_point;
	std::vector<std::vector<compx>> value;
	void solve_triag();
	double potential(int, int);
	double real_space(int);
	double density_x1[grid_point];
	double density_x2[grid_point];
	compx sec_space_deriv(char, int, int);
	void normalize();
	void rho();
};

wave_function::wave_function() {
	value.resize(grid_point);
	for (int l = 0; l < grid_point; l++) {
		value[l].resize(grid_point);
		for (int m = 0; m < grid_point; m++)
			value[l][m] = compx(0., 0.);
	}
}

wave_function::wave_function(bool init) {
	if (init) {
		compx psi00, psi10, psi12;
		value.resize(grid_point);
		for (int l = 0; l < grid_point; l++) {
			value[l].resize(grid_point);
			for (int m = 0; m < grid_point; m++) {
				/*value[l][m] = exp(I*compx(k1, 0)*compx(real_space(l), 0))*exp(-pow(real_space(l) - real_space(x_01), 2.) / (4.*sig*sig))
				*exp(I*compx(k2, 0)*compx(real_space(m), 0))*exp(-pow(real_space(m) - real_space(x_02), 2.) / (4.*sig*sig));*/
				psi00 = pow(M*omega1 / (PI*h_bar), .25)*exp(-M*omega1*pow(real_space(l), 2.) / (2.*h_bar))*pow(M*omega2 / (PI*h_bar), .25)*exp(-M*omega2*pow(real_space(m), 2.) / (2.*h_bar));
				//psi10 = pow(M*omega1 / (PI*h_bar), .25)*pow(2., .5)*pow(M*omega1 / (h_bar), .5)*real_space(l)*exp(-M*omega1*pow(real_space(l), 2.) / (2.*h_bar))*pow(M*omega2 / (PI*h_bar), .25)*exp(-M*omega2*pow(real_space(m), 2.) / (2.*h_bar));
				value[l][m] = psi00;// +psi10;
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

//Use this when boudaries are periodic
void wave_function::solve_triag() {

	compx a = -h_bar / (two*M), b = -h_bar / (two*M);
	double r = dt / (dx*dx);
	compx A = (I*r*a / two), B = (I*r*b / two), C = I*dt / two;
	compx mid = one - two*A;
	compx x_N, q_N, tvalue;
	std::vector<compx> alpha(grid_point - 2);
	std::vector<compx> beta1(grid_point - 1);//for x_1
	std::vector<compx> beta2(grid_point - 1);//for x_2
	wave_function tmp; //hold the n* values
	wave_function x_1, x_2; //hold the solutions of the tridiagonal "condensed" system with grid_point - 1 unknown

							//Solving for the x1 direction first
	for (int x2 = 0; x2 < grid_point; x2++) {
		//std::cout << " x2  " << x2 << std::endl;
		alpha[0] = A / mid;

		compx cur_pot = potential(0, x2);
		compx der2 = sec_space_deriv('y', 0, x2);

		beta1[0] = ((one - C*cur_pot)*value[0][x2] - B*der2) / mid;
		beta2[0] = -A / mid;

		//std::cout << "beginning" << std::endl;
		//Forward run
		for (int l = 1; l < grid_point - 2; l++) {
			alpha[l] = A / (mid - A*alpha[l - 1]);

			beta1[l] = ((one - C*potential(l, x2))*value[l][x2] - B*sec_space_deriv('y', l, x2)
				- A*beta1[l - 1]) / (mid - A*alpha[l - 1]);
			beta2[l] = -A*beta2[l - 1] / (mid - A*alpha[l - 1]);
			//std::cout << i << std::endl;
		}
		//std::cout << "potential  " << potential(1, x2) << "  here   " << beta1[1] << std::endl;
		beta1[grid_point - 2] = ((one - C*potential(grid_point - 2, x2))*value[grid_point - 2][x2] - B*sec_space_deriv('y', grid_point - 2, x2)
			- A*beta1[grid_point - 3]) / (mid - A*alpha[grid_point - 3]);
		beta2[grid_point - 2] = (-A - A*beta2[grid_point - 3]) / (mid - A*alpha[grid_point - 3]);
		//std::cout << "first loop" << std::endl;

		//Backward run
		x_1.value[grid_point - 2][x2] = beta1[grid_point - 2];
		x_2.value[grid_point - 2][x2] = beta2[grid_point - 2];
		for (int l = grid_point - 3; l >= 0; l--) {
			x_1.value[l][x2] = beta1[l] - alpha[l] * x_1.value[l + 1][x2];
			x_2.value[l][x2] = beta2[l] - alpha[l] * x_2.value[l + 1][x2];
		}

		compx q_N = (one - C*potential(grid_point - 1, x2))*value[grid_point - 1][x2] - B*sec_space_deriv('y', grid_point - 1, x2);
		compx x_N = (q_N - A*x_1.value[0][x2] - A*x_1.value[grid_point - 2][x2])
			/ (mid + A*x_2.value[0][x2] + A*x_2.value[grid_point - 2][x2]);

		//Solve for tmp values (aka n*)
		tmp.value[grid_point - 1][x2] = x_N;
		for (int l = 0; l < grid_point - 1; l++) {
			tmp.value[l][x2] = x_1.value[l][x2] + x_2.value[l][x2] * x_N;
			tvalue = tmp.value[l][x2];
			//std::cout << l << "  " << x2 << std::endl;
		}
		//std::cout << tmp.value[0][x2] << std::endl;
		//std::cout << x2 << "   The end" << std::endl;
	}

	tmp.normalize();

	//std::cout << "x direction done" << std::endl;
	//Solving for the x2 direction
	for (int x1 = 0; x1 < grid_point; x1++) {

		alpha[0] = B / (one - two*B + C*potential(x1, 0));

		beta1[0] = (tmp.value[x1][0] - A*tmp.sec_space_deriv('x', x1, 0)) / (one - two*B + C*potential(x1, 0));
		beta2[0] = -B / (one - two*B + C*potential(x1, 0));
		//std::cout << x1 << "   beginning2" << std::endl;
		//Forward run
		for (int m = 1; m < grid_point - 2; m++) {
			alpha[m] = B / (one - two*B + C*potential(x1, m) - B*alpha[m - 1]);

			beta1[m] = (tmp.value[x1][m] - A*tmp.sec_space_deriv('x', x1, m) - B*beta1[m - 1])
				/ (one - two*B + C*potential(x1, m) - B*alpha[m - 1]);
			beta2[m] = -B*beta2[m - 1] / (one - two*B + C*potential(x1, m) - B*alpha[m - 1]);
		}
		beta1[grid_point - 2] = (tmp.value[x1][grid_point - 2] - A*tmp.sec_space_deriv('x', x1, grid_point - 2) - B*beta1[grid_point - 3])
			/ (one - two*B + C*potential(x1, grid_point - 2) - B*alpha[grid_point - 3]);
		beta2[grid_point - 2] = (-B - B*beta2[grid_point - 3]) / (one - two*B + C*potential(x1, grid_point - 2) - B*alpha[grid_point - 3]);
		//std::cout << "first loop2" << std::endl;

		//Backward run
		x_1.value[x1][grid_point - 2] = beta1[grid_point - 2];
		x_2.value[x1][grid_point - 2] = beta2[grid_point - 2];
		for (double m = grid_point - 3; m >= 0; m--) { //-2 because of the same reason above
			x_1.value[x1][m] = beta1[m] - alpha[m] * x_1.value[x1][m + 1];
			x_2.value[x1][m] = beta2[m] - alpha[m] * x_2.value[x1][m + 1];
		}

		q_N = tmp.value[x1][grid_point - 1] - A*tmp.sec_space_deriv('x', x1, grid_point - 1);
		x_N = (q_N - B*x_1.value[x1][0] - B*x_1.value[x1][grid_point - 2])
			/ ((one - two*B + C*potential(x1, grid_point - 1)) + B*x_2.value[x1][0] + B*x_2.value[x1][grid_point - 2]);

		value[x1][grid_point - 1] = x_N;
		for (int m = 0; m < grid_point - 1; m++) {
			value[x1][m] = x_1.value[x1][m] + x_2.value[x1][m] * x_N;
			//std::cout << x1 << "  " << m << std::endl;
		}

	}

	normalize();
}

double wave_function::potential(int x1, int x2) {

	return .5*M*omega1*omega1*(real_space(x1))*(real_space(x1)) + .5*M*omega1*omega1*(real_space(x2))*(real_space(x2)) + .5*M*omega2*omega2*pow(real_space(x2)-real_space(x1),2.);
	//square well
	/*if (al - abs(real_space(x1) - real_space(x2)) > 0.)
	return V_0;
	else
	return 0.0;*/

	//gaussian
	//return V_0*exp(-pow(abs(real_space(x1)-real_space(x2)),2.)/(2.*al*al));
}

double wave_function::real_space(int index) {
	return (-L / 2.0) + (index*dx);
}


//x stands for x1 or particle 1 and y stands for x2 or particle 2
compx wave_function::sec_space_deriv(char x1_or_x2, int l, int m) {
	if (x1_or_x2 == 'x') {
		if (l == 0)
			return value[l + 1][m] - two * value[l][m] + value[grid_point - 1][m];
		else if (l == grid_point - 1)
			return value[0][m] - two * value[l][m] + value[l - 1][m];
		else
			return value[l + 1][m] - two * value[l][m] + value[l - 1][m];
	}
	else if (x1_or_x2 == 'y') {
		if (m == 0)
			return value[l][m + 1] - two * value[l][m] + value[l][grid_point - 1];
		else if (m == grid_point - 1)
			return value[l][0] - two * value[l][m] + value[l][m - 1];
		else
			return value[l][m + 1] - two * value[l][m] + value[l][m - 1];
	}
}

//Copied from stack overflow to deal with string for now 
template <typename NEW>
std::string to_string_with_precision(const NEW a_value, const int n = 6)
{
	std::ostringstream out;
	out << std::fixed << std::setprecision(n) << a_value;
	return out.str();
}


//normalization for the function above rho
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

void wave_function::rho() {
	for (int k = 0; k < grid_point; k++) {
		density_x1[k] = 0.;
		density_x2[k] = 0.;
	}

	for (int i = 0; i < grid_point; i++) {
		for (int j = 0; j < grid_point; j++) {
			density_x1[i] += pow(abs(value[i][j]), 2.)*dx;
			density_x2[i] += pow(abs(value[j][i]), 2.)*dx;
		}
	}
}

int main() {
	double check = 0.0;
	wave_function v(true);
	std::ofstream file1; //density of both particles
	std::ofstream file2;
	std::ofstream file5; //potential
	file5.open("potential.dat");
	for (int i = 0; i < grid_point; i++)
		for (int j = 0; j < grid_point; j++)
			file5 << v.real_space(i) << "\t" << v.real_space(j) << "\t" << v.potential(i, j) << std::endl;
	int index = 0;
	for (double k = v.dt; k <= 2100 * v.dt; k += v.dt) {
		std::cout << index << std::endl;
		v.rho();
		if (index % 10 == 0) {
			file1.open("data_" + to_string_with_precision(index, 0) +".dat");
			file1 << "Time   " << k - v.dt << std::endl
				<< "x" << "\t" << "y" << "\t" << "imag" << "\t" << "real" << "\t" << "abs" << std::endl;
			file2.open("x_" + std::to_string(index) + ".dat");
			file2 << "Time   " << k - v.dt << std::endl
				<< "Coordinate" << "\t" << "x1" << "\t" << "x2" << std::endl;
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
			file2.close();
			file1.close();

		}
		v.solve_triag();
		index++;
	}
	getchar();
}