#include "Solvepara.h"
#include <ceres/ceres.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <ctime>
#include <fstream>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
using ceres::AutoDiffCostFunction;
using ceres::CauchyLoss;
using ceres::CostFunction;
using ceres::LossFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

Solvepara::Solvepara()
{
}


Solvepara::~Solvepara()
{
}

Vector3d Solvepara::readPointsFromTxt(const string &filename, vector<Vector3d> &points, int &point_number_total)
{
	points.clear();
	point_number_total = 0;

	ifstream ifs(filename);
	string str;
	vector<string> vec;

	vector<double> xs, ys, zs;

	while (getline(ifs, str))
	{
		vector<string> sVec;
		boost::split(sVec, str, boost::is_any_of(", "), boost::token_compress_on);

		Vector3d cp;
		cp(0) = atof(sVec[0].c_str()); xs.push_back(cp(0));
		cp(1) = atof(sVec[1].c_str()); ys.push_back(cp(1));
		cp(2) = atof(sVec[2].c_str()); zs.push_back(cp(2));

		points.push_back(cp); ++point_number_total;
	}
	ifs.close();

	//求平均值并重心化
	Vector3d mean_coor;
	double sumx = std::accumulate(std::begin(xs), std::end(xs), 0.0); mean_coor(0) = sumx / point_number_total;
	double sumy = std::accumulate(std::begin(ys), std::end(ys), 0.0); mean_coor(1) = sumy / point_number_total;
	double sumz = std::accumulate(std::begin(zs), std::end(zs), 0.0); mean_coor(2) = sumz / point_number_total;

	for (size_t i = 0; i < point_number_total; i++)
	{
		points[i](0) -= mean_coor(0);
		points[i](1) -= mean_coor(1);
		points[i](2) -= mean_coor(2);
	}
	//排序
	sort(points.begin(), points.end(), compPointByY);
	return mean_coor;
}

void Solvepara::segmentPoints(const float &threshold, const vector<Vector3d> &points, vector<vector<Vector3d> > &points_segmented, 
	vector<vector<double> > &u_init, vector<Vector3d> &start_points)
{
	int segment_label = 0;
	vector<double> distance_accumu(points.size());
	distance_accumu[0] = 0.0; //从每段的初始点到i点的距离和

	vector<Vector3d> points_seg;
	vector<double> u_seg;
	start_points.push_back(points.front());
	
	for (size_t i = 0; i < points.size() - 1; i++)
	{
		// i到i+1点的距离
		double dist2nextP = (points[i](0) - points[i + 1](0))*(points[i](0) - points[i + 1](0))
			+ (points[i](1) - points[i + 1](1))*(points[i](1) - points[i + 1](1));
		dist2nextP = sqrt(dist2nextP);
		distance_accumu[i + 1] = dist2nextP + distance_accumu[i];

		u_seg.push_back(distance_accumu[i]);
		points_seg.push_back(points[i]);

		// 以第i点为连接点,最后一段至少有3个点
		if (distance_accumu[i + 1] > threshold && i + 3 < points.size())
		{
			u_init.push_back(u_seg);
			points_segmented.push_back(points_seg);
			u_seg.clear();
			points_seg.clear();
			u_seg.push_back(0.0);
			points_seg.push_back(points[i]);
			start_points.push_back(points[i]);
			distance_accumu[i + 1] -= distance_accumu[i];
			distance_accumu[i] = 0.0;
		}
		else if (distance_accumu[i + 1] > threshold && i + 3 >= points.size())
		{
			for (int p = i + 1; p < points.size()-1; ++p)
			{
				dist2nextP = (points[i](0) - points[i + 1](0))*(points[i](0) - points[i + 1](0))
					+ (points[i](1) - points[i + 1](1))*(points[i](1) - points[i + 1](1));
				dist2nextP = sqrt(dist2nextP);
				distance_accumu[i + 1] = dist2nextP + distance_accumu[i];

				u_seg.push_back(distance_accumu[i]);
				points_seg.push_back(points[i]);
			}
			u_seg.push_back(distance_accumu.back());
			points_seg.push_back(points.back());
			u_init.push_back(u_seg);
			points_segmented.push_back(points_seg);

			break;
		}

		if (i == points.size() - 2)
		{
			u_seg.push_back(distance_accumu[i+1]);
			points_seg.push_back(points[i+1]);
			u_init.push_back(u_seg);
			points_segmented.push_back(points_seg);
		}
	}
}

void Solvepara::segmentPointsByKnownU(const float &threshold, const vector<Vector3d> &points, const vector<double> &u, vector<vector<Vector3d> > &points_segmented,
	vector<vector<double> > &u_init, vector<Vector3d> &start_points)
{
	int segment_label = 0;
	vector<Vector3d> points_seg;
	vector<double> u_seg;
	start_points.push_back(points.front());

	double dist_accu = 0.0;
	int low = 0;
	int hig = 0;

	for (int i = 1; i < points.size() - 1; ++i)
	{
		int a = floor(u[i] / threshold);
		int b = floor(u[i+1] / threshold);
		if ( a != b && i+3 < points.size())
		{
			hig = i;
			for (int j = low; j <= hig; ++j)
			{
				u_seg.push_back(u[j] - u[low]);
				points_seg.push_back(points[j]);
			}
			start_points.push_back(points_seg.front());
			u_init.push_back(u_seg); u_seg.clear();
			points_segmented.push_back(points_seg); points_seg.clear();
			low = i;
		}
		else if (a != b && i + 3 >= points.size() )
		{
			hig = points.size();
			for (int j = low; j <= hig; ++j)
			{
				u_seg.push_back(u[j] - u[low]);
				points_seg.push_back(points[j]);
			}
			start_points.push_back(points_seg.front());
			u_init.push_back(u_seg); u_seg.clear();
			points_segmented.push_back(points_seg); points_seg.clear();
			low = i;
			return;
		}
	}

	for (int j = low; j <= points.size(); ++j)
	{
		u_seg.push_back(u[j] - u[low]);
		points_seg.push_back(points[j]);
	}
	start_points.push_back(points_seg.front());
	u_init.push_back(u_seg);
	points_segmented.push_back(points_seg);
}

void Solvepara::solveUByFixedPara(const vector<vector<double> > &input_u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
	const vector<vector<Vector3d> > &data_points, vector<vector<double> > &output_u)
{
	int segment_number = input_u.size();
	output_u.resize(segment_number);
	
	for (int i = 0; i < segment_number; ++i)
	{
		int segment_point_number = input_u[i].size();
		output_u[i].resize(segment_point_number);
		output_u[i][0] = 0.0;
		// 里程数u*1000即为对应的数组下标ref
		vector<Vector3d> coordinates = ccltIntegralVector(data_points[i][0], greek[i], greek_v[i], input_u[i].back());
		for (int j = 1; j < segment_point_number; ++j)
		{
			long pos = ccltNearestPointIndex(data_points[i][j], coordinates);
			output_u[i][j] = pos * 0.001;
		}
		cout << "---------UPDATE U: segment " << i + 1 << " / " << segment_number << "------------" << endl << endl;
	}
	cout << endl;
}

void Solvepara::solveParaByFixedU(const vector<vector<double> > &input_u, const vector<Vector3d> &greek_input, const vector<Vector2d> &greek_v_input,
	const vector<vector<Vector3d> > &data_points, vector<Vector3d> &output_greek, vector<Vector2d> &output_greek_v)
{
	int segment_number = input_u.size(); //有多少段
	output_greek.resize(segment_number); output_greek_v.resize(segment_number);

	// 遍历每一段
	for (int i = 0; i < segment_number; i++)
	{
		int seg_point_number = data_points[i].size(); //一段里有多少点
		// 第一次大循环直接采用init值作为初值
		double greek_i[3] = { greek_input[i](0), greek_input[i](1), greek_input[i](2) };
		double greek_solve[3] = { greek_input[i](0), greek_input[i](1), greek_input[i](2) };

		// 每段solve一次
		Problem problem;

		// Configure the loss function.
		LossFunction* loss = NULL;

		// Add the residuals.
		double xx, yy;
		double x00 = data_points[i][0](0);
		double y00 = data_points[i][0](1);

		// 对每个观测加一个cost
		for (int j = 0; j < seg_point_number; ++j)
		{
			double g_iterator = input_u[i][j];
			xx = data_points[i][j](0);
			yy = data_points[i][j](1);
			CostFunction *cost =
				new AutoDiffCostFunction<DistanceFromCurveCost, 2, 3>(
				new DistanceFromCurveCost(xx, yy, g_iterator, x00, y00));
			problem.AddResidualBlock(cost, loss, greek_solve);
		}

		problem.SetParameterLowerBound(greek_solve, 0, 0.00);	problem.SetParameterUpperBound(greek_solve, 0, 2 * M_PI);
		problem.SetParameterLowerBound(greek_solve, 1, -0.05);	problem.SetParameterUpperBound(greek_solve, 1, 0.05);
		problem.SetParameterLowerBound(greek_solve, 2, -0.0025);problem.SetParameterUpperBound(greek_solve, 2, 0.0025);

		Solver::Options options;
		options.max_num_iterations = 1000;
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = true;
		Solver::Summary summary;

		std::cout << "---------- segment:  " << i + 1 << " / " << segment_number << "  ------------" << endl;
		Solve(options, &problem, &summary);

		std::cout << summary.BriefReport() << "\n";
		std::cout << to_string(greek_i[0]) << "--->" << to_string(greek_solve[0]) << endl;
		std::cout << to_string(greek_i[1]) << "--->" << to_string(greek_solve[1]) << endl;
		std::cout << to_string(greek_i[2]) << "--->" << to_string(greek_solve[2]) << endl;
		std::cout << "---------- segment:  " << i + 1 << " / " << segment_number << "  ------------" << endl;

		std::cout << endl << endl;

		output_greek[i](0) = greek_solve[0];
		output_greek[i](1) = greek_solve[1];
		output_greek[i](2) = greek_solve[2];
	}
}

vector<double> Solvepara::translateUtoGlobal(const vector<vector<double> > &u)
{
	if (u.size() < 2)
		return u[0];

	vector<double > rst;
	for (int i = 0; i < u[0].size(); i++)
		rst.push_back(u[0][i]);

	for (int i = 1; i < u.size(); ++i)
	{
		double offset = rst.back();
		for (int j = 1; j < u[i].size(); ++j)
		{
			double u_tmp = u[i][j] + offset;
			rst.push_back(u_tmp);
		}
	}

	return rst;
}


vector<Vector3d> Solvepara::ccltIntegralVector(const Vector3d &start_point, const Vector3d &greek, const Vector2d &greek_v, const double &length)
{
	double step = 0.001;
	double gene_length = floor((length + 20.0) * 100 + 0.5) / 100.0;

	vector<double> ref; //离散的求积分节点
	for (long i = 0; i*step < gene_length; ++i)
		ref.push_back(i*step);
	ref.push_back(gene_length);
	long vector_size = ref.size();

	vector<Vector2d> triangle_val(vector_size);
	vector<Vector3d> coordi_point(vector_size);
	for (long i = 0; i < vector_size; ++i)
	{
		triangle_val[i](0) = cos(greek[0] + greek[1] * ref[i] + 0.5*greek[2] * ref[i] * ref[i]);
		triangle_val[i](1) = sin(greek[0] + greek[1] * ref[i] + 0.5*greek[2] * ref[i] * ref[i]);
		if (i > 0)
		{
			coordi_point[i](0) = coordi_point[i - 1](0) + 0.5*step*(triangle_val[i - 1](0) + triangle_val[i](0));
			coordi_point[i](1) = coordi_point[i - 1](1) + 0.5*step*(triangle_val[i - 1](1) + triangle_val[i](1));
			coordi_point[i](2) = coordi_point[0](2) + ref[i] * greek_v(0) + 0.5*greek_v(1)*ref[i] * ref[i];
		}
		else
		{
			coordi_point[i](0) = start_point(0);
			coordi_point[i](1) = start_point(1);
			coordi_point[i](2) = start_point(2);
		}
	}
	return coordi_point;
}

long Solvepara::ccltNearestPointIndex(const Vector3d &ref_point, const vector<Vector3d> &curve)
{
	long pos = 0;
	double least_distance = (ref_point(0) - curve[0](0))*(ref_point(0) - curve[0](0))
		+ (ref_point(1) - curve[0](1))*(ref_point(1) - curve[0](1));

	for (long i = 1; i < curve.size(); ++i)
	{
		double distance = (ref_point(0) - curve[i](0))*(ref_point(0) - curve[i](0))
			+ (ref_point(1) - curve[i](1))*(ref_point(1) - curve[i](1));
		if (distance <= least_distance)
		{
			least_distance = distance;
			pos = i;
		}
	}
	return pos;
}

bool compPointByY(const Vector3d &a, const Vector3d &b)
{
	if (a(1) < b(1))
		return true;
	else
		return false;
}

vector<double> Solvepara::ZsolveParaByFixedU(const vector<vector<double> > &input_u, const vector<vector<Vector3d> > &data_points, vector<Vector2d> &output_greek_v)
{
	// 每一段起始z
	vector<double> start_z;

	int segment_number = input_u.size();
	output_greek_v.clear();
	output_greek_v.resize(segment_number);

	for (int i = 0; i < segment_number; i++)
	{
		int seg_point_number = input_u[i].size();

		double greek_i[2] = { 0.0, 0.0 };
		double greek_solve[2] = { 0.0, 0.0 };

		// 每段solve一次
		Problem problem;

		// Configure the loss function.
		LossFunction* loss = NULL;

		// Add the residuals.
		double z00 = data_points[i][0](2);
		start_z.push_back(z00);

		// 对每个观测加一个cost
		for (int j = 0; j < seg_point_number; ++j)
		{
			double zz = data_points[i][j](2);
			double u = input_u[i][j];
			CostFunction *cost =
				new AutoDiffCostFunction<ZDistanceFromCurveCost, 1, 2>(
				new ZDistanceFromCurveCost(zz, u, z00));
			problem.AddResidualBlock(cost, loss, greek_solve);
		}
		problem.SetParameterLowerBound(greek_solve, 0, -0.2);	problem.SetParameterUpperBound(greek_solve, 0, 0.2);
		problem.SetParameterLowerBound(greek_solve, 1, -0.01);	problem.SetParameterUpperBound(greek_solve, 1, 0.01);

		Solver::Options options;
		options.max_num_iterations = 1000;
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = false;
		Solver::Summary summary;

		std::cout << "------V------ segment:  " << i + 1 << " / " << segment_number << "  ------V-------" << endl;
		Solve(options, &problem, &summary);

		std::cout << summary.BriefReport() << "\n";
		std::cout << to_string(greek_i[0]) << "--->" << to_string(greek_solve[0]) << endl;
		std::cout << to_string(greek_i[1]) << "--->" << to_string(greek_solve[1]) << endl;
		std::cout << "------V------ segment:  " << i + 1 << " / " << segment_number << "  ------V-------" << endl;

		std::cout << endl << endl;

		output_greek_v[i](0) = greek_solve[0];
		output_greek_v[i](1) = greek_solve[1];
	}
	return start_z;
}

vector<double> Solvepara::ccltGlobalBaseU(const vector<vector<double> > &u)
{
	if (u.size() < 2)
		return u[0];

	vector<double > rst;
	rst.push_back(0.0);

	for (int i = 1; i < u.size(); ++i)
	{
		rst.push_back(rst[i - 1] + u[i - 1].back());
	}

	return rst;
}
