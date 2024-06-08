/*

This program uses Runge Kutta 4 for model the orbit of Jupiter and the Earth
around the sun. The initial goal of this exercise was to model how the orbits
are affected by the mass of the sun. You can see this for yourself by changing
mass_Sun in main(). mass_Sun = 1 means we are using a modeled sun that is 100%
the mass of the actual sun.

INPUT: None

OUTPUTS 3 files:
	Earth_periodicity.txt: 	a 2d timeseries of the position of the Earth relative to
													some starting position. Formatted in two columns
													'time xposition'
	orbit_Earth.txt: 		a 2d plot of the Earth's orbit around the sun.
													Formatted in two columns 'xposition yposition'
	orbit_Jupiter.txt: 	a 2d plot of Jupiters's orbit around the sun.
													Formatted in two columns 'xposition yposition'


*/


#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <stdlib.h>

using namespace std;


// For one orbit (not taking into account forces from both sun and secondary planet)

// void rk4(vector<double> &system, vector<double> &dydx, const double x, const double h,
// 	vector<double> &yout, void derivs(double, vector<double> &, vector<double> &, vector<double> &, vector<double> &))
// {
// 	int n=system.size();
// 	vector<double> dym(n),dyt(n),yt(n);
// 	double h_2=h*0.5;
// 	double h_6 = h/6.0;
// 	//double xh = x + h_2;
// 	for (int i=0;i<n;i++) yt[i]=system[i]+h_2*dydx[i];

// 	derivs(x+h_2,yt,dyt); //derivs(xh,yt,dyt);
// 	for (int i=0;i<n;i++) yt[i]=system[i]+h_2*dyt[i];
// 	derivs(x+h_2,yt,dym); //derivs(xh,yt,dym);
// 	for (int i=0;i<n;i++) {
// 		yt[i]=system[i]+h*dym[i];
// 		dym[i] += dyt[i];
// 	}
// 	derivs(x+h,yt,dyt); //derivs(x+h,yt,dyt);
// 	for (int i=0;i<n;i++)
// 		yout[i]=system[i]+h_6*(dydx[i]+dyt[i]+2.0*dym[i]);

// }

// void derivsPlanet(double x, vector<double>& system, vector<double>& dydx) {
// 	double value,gm;
// 	gm = 1.32167e8;
// 	//y[0] = x;
// 	value = gm*pow(system[0]*system[0] + system[2]*system[2], -1.5);
// 	//value=1;
// 	dydx[0] = system[1];
// 	dydx[1] = -system[0] * value;
// 	dydx[2] = system[3];
// 	dydx[3] = -system[2] * value;
// }


void rk4(vector<double> &Earth, vector<double> &Jupiter, vector<double> &Earth_dydx,
	vector<double> &Jupiter_dydx, const double x, const double h,
	vector<double> &Earth_UPDATED, vector<double> &Jupiter_UPDATED,
	void derivs(double, vector<double> &, vector<double> &, vector<double> &, vector<double> &, float, float), 
	const float Earth_scale, const float Jupiter_scale) {

	int n = Earth.size();
	vector<double> Earth_dym(n),Earth_dyt(n),Earth_yt(n);
	vector<double> Jupiter_dym(n),Jupiter_dyt(n),Jupiter_yt(n);
	double h_2=h*0.5;
	double h_6 = h/6.0;

	for (int i = 0; i < n; i++) {
		Earth_yt[i]=Earth[i]+h_2*Earth_dydx[i];
		Jupiter_yt[i]=Jupiter[i]+h_2*Jupiter_dydx[i];
	}

	derivs(x+h_2,Earth_yt,Jupiter_yt,Earth_dyt,Jupiter_dyt,Earth_scale,Jupiter_scale);

	for (int i = 0; i < n; i++) {
		Earth_yt[i]=Earth[i]+h_2*Earth_dyt[i];
		Jupiter_yt[i]=Jupiter[i]+h_2*Jupiter_dyt[i];
	}

	derivs(x+h_2,Earth_yt,Jupiter_yt,Earth_dym,Jupiter_dym,Earth_scale,Jupiter_scale);

	for (int i = 0; i < n; i++) {
		Earth_yt[i]=Earth[i]+h*Earth_dym[i];
		Earth_dym[i] += Earth_dyt[i];
		Jupiter_yt[i]=Jupiter[i]+h*Jupiter_dym[i];
		Jupiter_dym[i] += Jupiter_dyt[i];
	}

	derivs(x+h,Earth_yt,Jupiter_yt,Earth_dyt,Jupiter_dyt,Earth_scale,Jupiter_scale); //derivs(x+h,yt,dyt);

	for (int i=0;i<n;i++) {
		Earth_UPDATED[i]=Earth[i]+h_6*(Earth_dydx[i]+Earth_dyt[i]+2.0*Earth_dym[i]);
		Jupiter_UPDATED[i]=Jupiter[i]+h_6*(Jupiter_dydx[i]+Jupiter_dyt[i]+2.0*Jupiter_dym[i]);
	}

}

void derivsPlanet(double t, vector<double>& Earth, vector<double>& Jupiter, 
	vector<double>& Earth_dydx, vector<double>& Jupiter_dydx, const float Earth_scale, const float Jupiter_scale) {

	double G = 39.483;// AU^(3/2) / (D(M)^1/2) where D is solar day and M is solar mass

	// masses are written in terms of fractions of Sun's mass

	double mass_Sun = 1;
	// double mass_Sun = 2;
	double mass_J = (1/1048.) * Jupiter_scale;
	// cout << "Jupiter mass " << mass_J << endl;

	//double mass_E = 1/332.;
	double mass_E = (1/332950.) * Earth_scale;
	// cout << "Earth mass " << mass_E << endl;



	double sun_acceleration_1 = G*mass_Sun * pow(Earth[0]*Earth[0] + Earth[2]*Earth[2], -1.5);
	double sun_acceleration_2 = G*mass_Sun * pow(Jupiter[0]*Jupiter[0] + Jupiter[2]*Jupiter[2], -1.5);

	double planetary_acceleration_1, planetary_acceleration_2;

	// acceleration of Earth because of Jupiter
	planetary_acceleration_1 = G*mass_J*pow(pow(Earth[0]*Earth[0]-Jupiter[0]*Jupiter[0],2) + pow(Earth[2]*Earth[2]-Jupiter[2]*Jupiter[2],2),-1.5);
	// acceleration of Jupiter because of Earth
	planetary_acceleration_2 = G*mass_E*pow(pow(Earth[0]*Earth[0]-Jupiter[0]*Jupiter[0],2) + pow(Earth[2]*Earth[2]-Jupiter[2]*Jupiter[2],2),-1.5);

	double value_1, value_2;
	value_1 = sun_acceleration_1 + planetary_acceleration_1; // on Earth by the Sun and Jupiter
	value_2 = sun_acceleration_2 + planetary_acceleration_2; // on Jupiter by the Sun and Earth

	//cout << Earth[0] << " " << Earth[2] << " " << value_1 << endl;

	Earth_dydx[0] = Earth[1];
	Earth_dydx[1] = -Earth[0] * value_1;
	Earth_dydx[2] = Earth[3];
	Earth_dydx[3] = -Earth[2] * value_1;

	Jupiter_dydx[0] = Jupiter[1];
	Jupiter_dydx[1] = -Jupiter[0] * value_2;
	Jupiter_dydx[2] = Jupiter[3];
	Jupiter_dydx[3] = -Jupiter[2] * value_2;
}

int main(int argc, char *argv[]) {
	float Earth_scale = 1.;  // default scale of the Earth
	float Jupiter_scale = 1.;  // default scale of the Earth

	cout << argc << endl;
	if (argc > 1) {  // a new scale for Earth has been input
		Earth_scale = atof(argv[1]);
	} if (argc > 2) {  // new scale for Jupiter has been provided
		Jupiter_scale = atof(argv[2]);
	}

	// double AU = 149597870.691; // km	(always a good value to have on hand)
	vector<double> Earth;	// Earth = [x0,vx0,y0,vy0]
	// starting positions used for Runge Kutta 4
	Earth.push_back(0.983);
	Earth.push_back(0);
	Earth.push_back(0);
	Earth.push_back(6.392);

	vector<double> Jupiter;	// Jupiter = [x0,vx0,y0,vy0]
	Jupiter.push_back(4.95);
	Jupiter.push_back(0);
	Jupiter.push_back(0);
	Jupiter.push_back(2.894);

	vector<double> Earth_dydx(4,0.0);
	vector<double> Jupiter_dydx(4,0.0);

	derivsPlanet(0,Earth,Jupiter,Earth_dydx,Jupiter_dydx,Earth_scale,Jupiter_scale);

	double h = 0.001;  // time step size for RK4 'h'
	vector<double> Earth_UPDATED(Earth);
	vector<double> Jupiter_UPDATED(Jupiter);

	ofstream file_a("orbit_Earth.txt");
	ofstream file_b("orbit_Jupiter.txt");

	double t = 0;

	int Earth_orbits = 0;
	int Jupiter_orbits = 0;
	ofstream period("Earth_periodicity.txt");
	for (double i = 0; i < 5000000; i+=100) {
		t = i*h;
		rk4(Earth,Jupiter,Earth_dydx,Jupiter_dydx,t,h,Earth_UPDATED,Jupiter_UPDATED,derivsPlanet,Earth_scale,Jupiter_scale);
		Earth = Earth_UPDATED;
		Jupiter = Jupiter_UPDATED;
		file_a << Earth[0] << " " << Earth[2] << endl;
		file_b << Jupiter[0] << " " << Jupiter[2] << endl;
		period << t << " " << Earth[0] << endl;
		if (abs(Earth[2]) < 0.01) Earth_orbits++;
		if (abs(Jupiter[2]) < 0.005) Jupiter_orbits++;

		derivsPlanet(t,Earth,Jupiter,Earth_dydx,Jupiter_dydx,Earth_scale,Jupiter_scale);

	}
	file_a.close();
	file_b.close();
	period.close();

	cout << Earth_orbits/2. << endl;
	cout << Jupiter_orbits/2. << endl;

  

}
