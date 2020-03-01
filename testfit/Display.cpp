#include "Display.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <shapefil.h>
#include <fstream>
#include <string>
#include <iostream>

using namespace std;
using namespace Eigen;

Display::Display()
{
}


Display::~Display()
{
}

void Display::generateShapeByParaFile(const string &parafile, const string &outputFile)
{
	ifstream ifs(parafile);
	vector<Vector3d> output_points;
	vector<double> X0, Y0, S, Mu, Ka, Ps, Ks, Yt, Z0;

	string str;
	while (getline(ifs, str))
	{
		vector<string> sVec;
		boost::split(sVec, str, boost::is_any_of(", "), boost::token_compress_on);

		X0.push_back(stod(sVec[0]));
		Y0.push_back(stod(sVec[1]));
		S.push_back(stod(sVec[2]));
		Mu.push_back(stod(sVec[3]));
		Ka.push_back(stod(sVec[4]));
		Ps.push_back(stod(sVec[5]));
		Ks.push_back(stod(sVec[6]));
		Yt.push_back(stod(sVec[7]));
		Z0.push_back(stod(sVec[8]));
	}

	for (int i = 0; i < X0.size(); ++i)
	{
		// 遍历每一段
		int segment_point_number = int(S[i] / 0.001);
		double mu_ = Mu[i];
		double ka_ = Ka[i];
		double ps_ = Ps[i];
		double x0_ = X0[i];
		double y0_ = Y0[i];
		double z0_ = Z0[i];
		double ks_ = Ks[i];
		double yt_ = Yt[i];

		double accu_x = 0.0; double accu_y = 0.0;
		for (int p = 0; p < segment_point_number; ++p)
		{
			double s_ = 0.001 * p;
			double cos_ = cos(mu_ + s_*ka_ + 0.5*ps_*s_*s_);
			double sin_ = sin(mu_ + s_*ka_ + 0.5*ps_*s_*s_);

			accu_x += cos_ * 0.001; accu_y += sin_ * 0.001;

			double z_plus = s_ * ks_ + 0.5*s_*s_*yt_;

			Vector3d point_;
			point_(0) = x0_ + accu_x; point_(1) = y0_ + accu_y; point_(2) = z0_ + z_plus;
			output_points.push_back(point_);
		}
	}

	int nSHPType = SHPT_POLYGONZ;
	const int point_number = output_points.size();

	SHPHandle	hSHPHandle;
	SHPObject	*psShape;
	double		*x;
	double		*y;
	double		*z;
	double		*m;
	int			i;

	x = (double *)malloc(output_points.size() * sizeof(double));
	y = (double *)malloc(output_points.size() * sizeof(double));
	z = (double *)malloc(output_points.size() * sizeof(double));
	m = (double *)malloc(output_points.size() * sizeof(double));
	string path_t = outputFile;
	hSHPHandle = SHPCreate(path_t.c_str(), nSHPType);

	for (i = 0; i < point_number; i++)
	{
		x[i] = output_points[i](0);
		y[i] = output_points[i](1);
		z[i] = output_points[i](2);
		m[i] = 100.0;
	}
	psShape = SHPCreateObject(nSHPType, -1, 0, NULL, NULL,
		point_number, x, y, z, m);
	SHPWriteObject(hSHPHandle, -1, psShape);
	SHPDestroyObject(psShape);
	SHPClose(hSHPHandle);

	free(x); free(y); free(z); free(m);

	std::cout << "input para file: " << parafile << endl;
	std::cout << "output shape file: " << outputFile << endl;
	std::cout << "________DONE________";
}

void Display::generateShapeByGreek(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
	const vector<Vector3d> &start_points, const Vector3d &mean_coor)
{
	vector<Vector3d> output_points;
	for (int i = 0; i < greek.size(); ++i)
	{
		// 遍历每一段
		int segment_point_number = int(u[i].back() / 0.001);
		double mu_ = greek[i](0);
		double ka_ = greek[i](1);
		double ps_ = greek[i](2);
		double x0_ = start_points[i](0);
		double y0_ = start_points[i](1);
		double z0_ = start_points[i](2);
		double ks_ = greek_v[i](0);
		double yt_ = greek_v[i](1);

		double accu_x = 0.0; double accu_y = 0.0;
		for (int p = 0; p < segment_point_number; ++p)
		{
			double s_ = 0.001 * p;
			double cos_ = cos(mu_ + s_*ka_ + 0.5*ps_*s_*s_);
			double sin_ = sin(mu_ + s_*ka_ + 0.5*ps_*s_*s_);

			accu_x += cos_ * 0.001; accu_y += sin_ * 0.001;

			double z_plus = s_ * ks_ + 0.5*s_*s_*yt_;

			Vector3d point_;
			point_(0) = x0_ + accu_x + mean_coor(0); point_(1) = y0_ + accu_y + mean_coor(1); point_(2) = z0_ + z_plus + mean_coor(2);
			output_points.push_back(point_);
		}
	}

	int nSHPType = SHPT_POLYGONZ;
	const int point_number = output_points.size();

	SHPHandle	hSHPHandle;
	SHPObject	*psShape;
	double		*x;
	double		*y;
	double		*z;
	double		*m;
	int			i;

	x = (double *)malloc(output_points.size() * sizeof(double));
	y = (double *)malloc(output_points.size() * sizeof(double));
	z = (double *)malloc(output_points.size() * sizeof(double));
	m = (double *)malloc(output_points.size() * sizeof(double));
	string path_t = outputFile;
	hSHPHandle = SHPCreate(path_t.c_str(), nSHPType);

	for (i = 0; i < point_number; i++)
	{
		x[i] = output_points[i](0);
		y[i] = output_points[i](1);
		z[i] = output_points[i](2);
		m[i] = 100.0;
	}
	psShape = SHPCreateObject(nSHPType, -1, 0, NULL, NULL,
		point_number, x, y, z, m);
	SHPWriteObject(hSHPHandle, -1, psShape);
	SHPDestroyObject(psShape);
	SHPClose(hSHPHandle);

	free(x); free(y); free(z); free(m);
}

void Display::ZgenerateShapeByGreek(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
	const vector<Vector3d> &start_points, const Vector3d &mean_coor,
	const vector<double> &z_start, const vector<double> &base_u, const vector<double> &base_v)
{
	vector<Vector3d> output_points;

	int v_seglabel = 0;

	for (int i = 0; i < greek.size(); ++i)
	{
		// 遍历每一段
		int segment_point_number = int(u[i].back() / 0.001 + 0.5) + 1;
		double mu_ = greek[i](0);
		double ka_ = greek[i](1);
		double ps_ = greek[i](2);
		double x0_ = start_points[i](0);
		double y0_ = start_points[i](1);

		double accu_x = 0.0; double accu_y = 0.0;
		for (int p = 0; p < segment_point_number; ++p)
		{
			double s_ = 0.001 * p;
			double cos_ = cos(mu_ + s_*ka_ + 0.5*ps_*s_*s_);
			double sin_ = sin(mu_ + s_*ka_ + 0.5*ps_*s_*s_);

			accu_x += cos_ * 0.001; accu_y += sin_ * 0.001;

			// 确定所属的z分段
			if (v_seglabel + 1 < z_start.size() && s_ + base_u[i] >= base_v[v_seglabel + 1])
				++v_seglabel;

			double ks_ = greek_v[v_seglabel](0);
			double yt_ = greek_v[v_seglabel](1);
			
			double Vs_ = s_ + base_u[i] - base_v[v_seglabel];
			double z_plus = Vs_ * ks_ + 0.5*Vs_*Vs_*yt_;

			Vector3d point_;
			point_(0) = x0_ + accu_x + mean_coor(0); point_(1) = y0_ + accu_y + mean_coor(1); point_(2) = z_plus + z_start[v_seglabel] + mean_coor(2);
			output_points.push_back(point_);
		}
	}

	int nSHPType = SHPT_POLYGONZ;
	const int point_number = output_points.size();

	SHPHandle	hSHPHandle;
	SHPObject	*psShape;
	double		*x;
	double		*y;
	double		*z;
	double		*m;
	int			i;

	x = (double *)malloc(output_points.size() * sizeof(double));
	y = (double *)malloc(output_points.size() * sizeof(double));
	z = (double *)malloc(output_points.size() * sizeof(double));
	m = (double *)malloc(output_points.size() * sizeof(double));
	string path_t = outputFile;
	hSHPHandle = SHPCreate(path_t.c_str(), nSHPType);

	for (i = 0; i < point_number; i++)
	{
		x[i] = output_points[i](0);
		y[i] = output_points[i](1);
		z[i] = output_points[i](2);
		m[i] = 100.0;
	}
	psShape = SHPCreateObject(nSHPType, -1, 0, NULL, NULL,
		point_number, x, y, z, m);
	SHPWriteObject(hSHPHandle, -1, psShape);
	SHPDestroyObject(psShape);
	SHPClose(hSHPHandle);

	free(x); free(y); free(z); free(m);
}