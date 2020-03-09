#include <fstream>
#include <string>
#include <iostream>
#include "Display.h"
#include "Solvepara.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <ctime>
#include <vector>
using namespace std;
using namespace Eigen;

void main()
{
	Solvepara solvepara;

	vector<Vector3d> total_points;
	int total_point_number;
	// 从文件中读取点坐标
	Vector3d mean_coor = solvepara.readPointsFromTxt("./midmark.txt", total_points, total_point_number);

	vector<vector<Vector3d> > points_segmented;
	vector<vector<double> > u_init, u_medi, u_rst;
	vector<Vector3d> start_points;
	// 分段
	solvepara.segmentPoints(200.0, total_points, points_segmented, u_init, start_points);
	int segment_number = points_segmented.size();

	// 曲线参数迭代初值
	vector<Vector3d> greek_init(segment_number), greek_medi, greek_rst;
	vector<Vector2d> greek_v_init(segment_number), greek_v_medi, greek_v_rst;
	for (int i = 0; i < segment_number; ++i)
	{
		greek_init[i](0) = 1.90; greek_init[i](1) = 0.0; greek_init[i](2) = 0.0;
		greek_v_init[i](0) = 0.0; greek_v_init[i](1) = 0.0;
	}

	clock_t start, end;
	start = clock();

	u_medi = u_init; greek_medi = greek_init; greek_v_medi = greek_v_init;

	for (int major_iter = 0; major_iter < 2; major_iter++)
	{
		solvepara.panel_.major_iter_num = major_iter + 1; // 告诉优化器这是第几次循环
		solvepara.solveParaByFixedU(u_medi, greek_medi, greek_v_medi, points_segmented, greek_rst, greek_v_rst);
		greek_medi = greek_rst;
		greek_v_medi = greek_v_rst;
		solvepara.solveUByFixedPara(u_medi, greek_medi, greek_v_medi, points_segmented, u_rst);
		u_medi = u_rst;
	}

	end = clock();
	double endtime = (double)(end - start) / CLOCKS_PER_SEC;
	cout << endl << "Total time:" << endtime << "s" << endl;		//s为单位

	// 高程方向优化
	//vector<double> global_u = solvepara.translateUtoGlobal(u_medi);
	//vector<vector<double> > Vu_init; vector<vector<Vector3d> > Vpoints_segmented; vector<Vector3d> Vstart_points;
	//solvepara.segmentPointsByKnownU(30.0, total_points, global_u, Vpoints_segmented, Vu_init, Vstart_points);
	//vector<double> Zstart = solvepara.ZsolveParaByFixedU(Vu_init, Vpoints_segmented, greek_v_rst);
	//vector<double> base_u = solvepara.ccltGlobalBaseU(u_medi);  //水平线形u转换
	//vector<double> base_v = solvepara.ccltGlobalBaseU(Vu_init); //垂直线形转换

	Display display;
	const string output_shape = "curve_vec.shp";
	//display.ZgenerateShapeByGreek(output_shape, u_medi, greek_medi, greek_v_rst, start_points, mean_coor, Zstart, base_u, base_v);
	vector<ZPara> zpara = solvepara.optimizeHeight(points_segmented, u_medi);
	display.ZgenerateShapeByZPara(output_shape, u_medi, greek_medi, greek_v_rst, start_points, mean_coor, zpara);


	ofstream ofs("u.txt");
	for (int i = 0; i < u_medi.size(); i++)
	{
		for (int j = 0; j < u_medi[i].size(); j++)
		{
			ofs << i << "," << j << "," << u_medi[i][j] << endl;
		}
	}

	int pp;
	std::cin >> pp;
}
