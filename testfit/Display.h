#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Solvepara.h"
using std::string;
using std::vector;
using Eigen::Vector3d;
using Eigen::Vector2d;

class Display
{
public:
	// 根据matlab输出的para文件输出shape文件
	void generateShapeByParaFile(const string &parafile, const string &outputFile);

	// 根据计算得到的参数输出shp文件
	void generateShapeByGreek(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<Vector3d> &start_points, const Vector3d &mean_coor);

	// 根据计算得到的参数输出shp文件，加入垂直方向，通过u的转换
	void ZgenerateShapeByGreek(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<Vector3d> &start_points, const Vector3d &mean_coor,
		const vector<double> &z_start, const vector<double> &base_u, const vector<double> &base_v);

	void ZgenerateShapeByZPara(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<Vector3d> &start_points, const Vector3d &mean_coor,
		const vector<ZPara> &zpara);

	Display();
	~Display();
};

