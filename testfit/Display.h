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
	// ����matlab�����para�ļ����shape�ļ�
	void generateShapeByParaFile(const string &parafile, const string &outputFile);

	// ���ݼ���õ��Ĳ������shp�ļ�
	void generateShapeByGreek(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<Vector3d> &start_points, const Vector3d &mean_coor);

	// ���ݼ���õ��Ĳ������shp�ļ������봹ֱ����ͨ��u��ת��
	void ZgenerateShapeByGreek(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<Vector3d> &start_points, const Vector3d &mean_coor,
		const vector<double> &z_start, const vector<double> &base_u, const vector<double> &base_v);

	void ZgenerateShapeByZPara(const string &outputFile, const vector<vector<double> > &u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<Vector3d> &start_points, const Vector3d &mean_coor,
		const vector<ZPara> &zpara);

	Display();
	~Display();
};

