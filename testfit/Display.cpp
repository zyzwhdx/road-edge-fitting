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

void Display::ZgenerateShapeByZPara(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
	const vector<Vector3d> &start_points, const Vector3d &mean_coor,
	const vector<ZPara> &zpara, const vector<vector<Vector3d> > &data)
{
	ofstream ofs("curvature.txt");
	ofstream ofs_dir("direction.txt");
	ofstream ofs_knot("knots.txt");

	double step_ = 0.01;
	vector<Vector3d> output_points;

	int v_seglabel;
	int knot_label;

	for (int i = 0; i < greek.size(); ++i)
	{
		// 遍历每一段
		ZPara para = zpara[i];
		v_seglabel = 0; knot_label = 0;
		int segment_point_number = int(u[i].back() / step_ + 0.5) + 1;
		double mu_ = greek[i](0);
		double ka_ = greek[i](1);
		double ps_ = greek[i](2);
		double x0_ = start_points[i](0);
		double y0_ = start_points[i](1);

		double accu_x = 0.0; double accu_y = 0.0;
		for (int p = 0; p < segment_point_number; ++p)
		{
			double s_ = step_ * p;
			double cos_ = cos(mu_ + s_*ka_ + 0.5*ps_*s_*s_);
			double sin_ = sin(mu_ + s_*ka_ + 0.5*ps_*s_*s_);

			accu_x += cos_ * step_; accu_y += sin_ * step_;

			// 确定所属的z分段
			if (v_seglabel + 1 < para.u_start.size() && s_ > para.u_start[v_seglabel + 1])
				++v_seglabel;

			double ks_ = para.para[v_seglabel](1);
			double yt_ = para.para[v_seglabel](2);
			double Vs_ = s_ - para.u_start[v_seglabel];

			double z_plus = Vs_ * ks_ + 0.5*Vs_*Vs_*yt_;

			Vector3d point_;
			point_(0) = x0_ + accu_x + mean_coor(0); point_(1) = y0_ + accu_y + mean_coor(1); point_(2) = z_plus + para.para[v_seglabel](0) + mean_coor(2);
			output_points.push_back(point_);

			// 输出曲率
			if (p % 10 == 0)
			{
				double curvature = s_*ps_ + ka_;
				ofs << to_string(point_(0)) << "," << to_string(point_(1)) << "," << to_string(point_(2)) << "," << to_string(abs(curvature * 1000000)) << endl;

				double direction = mu_ + s_*ka_ + 0.5*ps_*s_*s_;
				ofs_dir << to_string(point_(0)) << "," << to_string(point_(1)) << "," << to_string(point_(2)) << "," << to_string(abs(direction * 1000)) << endl;
			}
		}
		for (int j = 0; j < u[i].size(); j++)
		{
			double up = u[i][j];
			double xx = CosIntegral(greek[i], up);
			double yy = SinIntegral(greek[i], up);
			Vector3d point_;
			point_(0) = x0_ + xx + mean_coor(0); point_(1) = y0_ + yy + mean_coor(1); point_(2) = 0.0;
			Vector3d diff = point_ - data[i][j] - mean_coor;
			double distance2d = sqrt(diff(0)*diff(0) + diff(1)*diff(1));

			point_(2) = findNearestZ(output_points, point_);

			ofs_knot << to_string(point_(0)) << "," << to_string(point_(1)) << "," << to_string(point_(2)) << "," << to_string(distance2d) << endl;
		}
	}

	ofs.close(); ofs_dir.close(); ofs_knot.close();

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

double Display::CosIntegral(const Vector3d &para, const double &up)
{
	double g_iterator = up;
	double step_ = 0.01; //ATTENTION

	double x_predicted = cos(para[0])*0.5*step_;

	for (long i = 0; i*step_ < g_iterator; ++i)
		x_predicted += 2 * 0.5 * step_*cos(para[0] + para[1] * i*step_ + 0.5*para[2] * i*i*step_*step_);

	double f_step = g_iterator - floor(100 * g_iterator) / 100.0;
	x_predicted += 0.5 * f_step*cos(para[0] + para[1] * g_iterator + 0.5*para[2] * g_iterator*g_iterator);

	return x_predicted;
}

double Display::SinIntegral(const Vector3d &para, const double &up)
{
	double g_iterator = up;
	double step_ = 0.01; //ATTENTION

	double x_predicted = sin(para[0])*0.5*step_;

	for (long i = 0; i*step_ < g_iterator; ++i)
		x_predicted += 2 * 0.5 * step_*sin(para[0] + para[1] * i*step_ + 0.5*para[2] * i*i*step_*step_);

	double f_step = g_iterator - floor(100 * g_iterator) / 100.0;
	x_predicted += 0.5 * f_step*sin(para[0] + para[1] * g_iterator + 0.5*para[2] * g_iterator*g_iterator);

	return x_predicted;
}

double Display::findNearestZ(const vector<Vector3d> &points, const Vector3d &ref)
{
	double rst = points[0](2);
	Vector3d dif = points[0] - ref;
	double minDis = dif(0)*dif(0) + dif(1)*dif(1);

	for (int i = 0; i < points.size(); i++)
	{
		dif = points[i] - ref;
		double dis = dif(0)*dif(0) + dif(1)*dif(1);

		if (dis < minDis)
		{
			minDis = dis;
			rst = points[i](2);
		}
	}
	return rst;
}