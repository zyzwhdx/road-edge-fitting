#include "Solvepara.h"
#include <ceres/ceres.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <ctime>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <ceres/loss_function.h>
#include <opencv2/core/types_c.h>
#include <opencv2/core/core_c.h>
#include <cmath>
#include <math.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
using namespace std;
using namespace pcl;
using Eigen::Vector3d;
using ceres::AutoDiffCostFunction;
using ceres::CauchyLoss;
using ceres::CostFunction;
using ceres::LossFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

Solvepara::Solvepara()
{
	this->panel_.major_iter_num = 0;
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
	this->mean_point = mean_coor;
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

	for (int i = 0; i < start_points.size() - 1; ++i)
	{
		int p_num = points_segmented[i].size();
		Vector3d sum = points_segmented[i][p_num - 1] + points_segmented[i + 1][1];
		Vector3d mean = sum / 2.0; //连接点

		// 这一段加在段尾
		Vector3d fin_p = points_segmented[i].back(); // 这一段的最后一点
		Vector3d dif = mean - fin_p;
		double distance = sqrt(dif(0)*dif(0) + dif(1)*dif(1) + dif(2)*dif(2)) + u_init[i].back();
		points_segmented[i].push_back(mean);
		u_init[i].push_back(distance);

		// 下一段加在段头
		Vector3d fst_p = points_segmented[i + 1][0];
		Vector3d sec_p = points_segmented[i + 1][1];
		double distance_ori = u_init[i + 1][1];
		points_segmented[i + 1][0] = mean;
		dif = mean - sec_p;
		distance = sqrt(dif(0)*dif(0) + dif(1)*dif(1) + dif(2)*dif(2));
		for (int j = 1; j < u_init[i + 1].size(); j++)
		{
			u_init[i + 1][j] = u_init[i + 1][j] - distance_ori + distance;
		}
		start_points[i + 1] = mean;
	}

	ofstream ofs("segment_points.txt");
	for (int i = 0; i < points_segmented.size(); i++)
	{
		for (int j = 0; j < points_segmented[i].size(); j++)
		{
			ofs << to_string(points_segmented[i][j](0) + mean_point(0)) << "," << to_string(points_segmented[i][j](1) + mean_point(1))
				<< "," << to_string(points_segmented[i][j](2) + mean_point(2)) << ","
				<< i << "," << u_init[i][j] << endl;
		}
	}
	ofs.close();
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


/************************************************************************/
/*   ATTENTION STEP SIZE                                                */
/************************************************************************/
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
			output_u[i][j] = pos * 0.01; //ATTENTION
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

	double curvature_cache;

	// 遍历每一段
	for (int i = 0; i < segment_number; i++)
	{
		int seg_point_number = data_points[i].size(); //一段里有多少点
		// 第一次大循环通过计算得到初值
		double greek_i[3] = { greek_input[i](0), greek_input[i](1), greek_input[i](2) };
		double greek_solve[3] = { greek_input[i](0), greek_input[i](1), greek_input[i](2) };
		if (panel_.major_iter_num < 2)
		{
			greek_i[0] = ccltBeginEndDir2D(data_points, i)(0); greek_i[1] = 0.0; greek_i[2] = 0.0;
			greek_solve[0] = greek_i[0]; greek_solve[1] = 0.0; greek_solve[2] = 0.0;
		}

		// 每段solve一次
		Problem problem;

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
			problem.AddResidualBlock(cost, new CauchyLoss(0.1), greek_solve);
		}

		//连接点约束，加权，仅当迭代次数大于1的时候加
		if (this->panel_.major_iter_num > 1)
		{
			double g_iterator = input_u[i].back();
			xx = data_points[i].back()(0);
			yy = data_points[i].back()(1);

			CostFunction *cost =
				new AutoDiffCostFunction<WeightedDistanceFromCurveCost_0, 2, 3>(
				new WeightedDistanceFromCurveCost_0(xx, yy, g_iterator, x00, y00, seg_point_number));
			problem.AddResidualBlock(cost, NULL, greek_solve);
		}

		// 连接方向约束，加权，仅当迭代次数大于1的时候加
		// 选参考点的方式有所不同
		Vector2d refdir;
		if (this->panel_.major_iter_num > 1)
		{
			Vector2d dir = ccltBeginEndDir2D(data_points, i); //开始和结尾的方向约束
			refdir = dir;
			double length = input_u[i].back();

			CostFunction *cost =
				new AutoDiffCostFunction<DirectionConnection2D, 2, 3>(
				new DirectionConnection2D(dir(0), dir(1), length));
			problem.AddResidualBlock(cost, NULL, greek_solve);
		}

		// 曲率约束
		if (this->panel_.major_iter_num > 1 && i > 0)
		{
			CostFunction *cost =
				new AutoDiffCostFunction<CurvatureCost, 1, 3>(
				new CurvatureCost(curvature_cache));
			problem.AddResidualBlock(cost, NULL, greek_solve);
		}

		// 待优化参数范围
		problem.SetParameterLowerBound(greek_solve, 0, 0.00);	problem.SetParameterUpperBound(greek_solve, 0, 2 * M_PI);
		problem.SetParameterLowerBound(greek_solve, 1, -0.05); problem.SetParameterUpperBound(greek_solve, 1, 0.05);
		problem.SetParameterLowerBound(greek_solve, 2, -0.0025); problem.SetParameterUpperBound(greek_solve, 2, 0.0025);

		Solver::Options options;
		options.max_num_iterations = 20;
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = true;
		Solver::Summary summary;

		std::cout << "---------- segment:  " << i + 1 << " / " << segment_number << "  ------------" << endl;
		std::cout << "major iteration: " << panel_.major_iter_num << endl;
		Solve(options, &problem, &summary);

		// 参考方向和计算方向对比
		std::cout << "ref_dir: " << to_string(refdir(0)) << "  rst_dir: " << to_string(greek_solve[0]) << endl;
		std::cout << "ref_dir: " << to_string(refdir(1)) << "  rst_dir: " << to_string(greek_solve[0] + greek_solve[1] * input_u[i].back()
			+ 0.5*greek_solve[2] * input_u[i].back() *input_u[i].back()) << endl;

		std::cout << summary.BriefReport() << "\n";
		std::cout << to_string(greek_i[0]) << "--->" << to_string(greek_solve[0]) << endl;
		std::cout << to_string(greek_i[1]) << "--->" << to_string(greek_solve[1]) << endl;
		std::cout << to_string(greek_i[2]) << "--->" << to_string(greek_solve[2]) << endl;
		std::cout << "length: " << to_string(input_u[i].back()) << endl;
		std::cout << "---------- segment:  " << i + 1 << " / " << segment_number << "  ------------" << endl;

		std::cout << endl << endl;

		output_greek[i](0) = greek_solve[0];
		output_greek[i](1) = greek_solve[1];
		output_greek[i](2) = greek_solve[2];

		curvature_cache = greek_solve[1] + input_u[i].back()*greek_solve[2];
	}
}

Vector2d Solvepara::ccltBeginEndDir2D(const vector<vector<Vector3d> > &data_points, const int &i)
{
	Vector2d rst;
	int segment_number = data_points.size();
	vector<Vector3d> dirP;
	if (i == 0) //第一段
	{
		for (int m = 0; m < 10; m++)
			dirP.push_back(data_points[i][m]);
		double dir = PCAccltDirectionConstraint2D(dirP);
		rst(0) = dir;
		cout << "PCA dir: " << to_string(dir) << endl;

		dirP.clear();
		for (int m = 0; m < 6; m++)
			dirP.push_back(data_points[i][data_points[i].size() - 6 + m]);
		for (int m = 0; m < 5; m++)
			dirP.push_back(data_points[i + 1][m + 1]);
		dir = PCAccltDirectionConstraint2D(dirP);
		rst(1) = dir;
		cout << "PCA dir: " << to_string(dir) << endl;
	}
	else if (i == segment_number - 1) //最后一段
	{
		if (data_points[i].size() < 12) // 如果最后一段点数太少
		{
			for (int p = 0; p < data_points[i].size(); p++)
				dirP.push_back(data_points[i][p]);
			rst(0) = PCAccltDirectionConstraint2D(dirP);
			rst(1) = rst(0);
			return rst;
			cout << "PCA dir: " << to_string(rst(0)) << endl;
			cout << "PCA dir: " << to_string(rst(0)) << endl;
		}

		for (int m = 0; m < 6; m++)
			dirP.push_back(data_points[i - 1][data_points[i - 1].size() - 6 + m]);
		for (int m = 0; m < 5; m++)
			dirP.push_back(data_points[i][m + 1]);
		double dir = PCAccltDirectionConstraint2D(dirP);
		rst(0) = dir;
		cout << "PCA dir: " << to_string(dir) << endl;

		dirP.clear();
		for (int m = 0; m < 10; m++)
			dirP.push_back(data_points[i][data_points[i].size() - 10 + m]);
		dir = PCAccltDirectionConstraint2D(dirP);
		rst(1) = dir;
		cout << "PCA dir: " << to_string(dir) << endl;
	}
	else
	{
		for (int m = 0; m < 6; m++)
			dirP.push_back(data_points[i - 1][data_points[i - 1].size() - 6 + m]);
		for (int m = 0; m < 5; m++)
			dirP.push_back(data_points[i][m + 1]);
		double dir = PCAccltDirectionConstraint2D(dirP);
		rst(0) = dir;
		cout << "PCA dir: " << to_string(dir) << endl;

		dirP.clear();
		for (int m = 0; m < 6; m++)
			dirP.push_back(data_points[i][data_points[i].size() - 6 + m]);
		for (int m = 0; m < 5; m++)
			dirP.push_back(data_points[i + 1][m + 1]);
		dir = PCAccltDirectionConstraint2D(dirP);
		rst(1) = dir;
		cout << "PCA dir: " << to_string(dir) << endl;
	}
	return rst;
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
	double step = 0.01; //ATTENTION
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

vector<ZPara> Solvepara::optimizeHeight(const vector<vector<Vector3d> > &points, const vector<vector<double> > &u)
{
	/************************************************************************/
	/*                       分段                                           */
	/************************************************************************/
	ofstream ofs("hpara.txt");
	vector<ZPara> para;
	for (int i = 0; i < points.size(); ++i)
	{
		vector<Vector3d> pointsseg = points[i]; //一大段里所有的点
		vector<double > useg = u[i]; //一大段里所有的u
		unsigned int pointnum = pointsseg.size(); //一大段里的点数
		ZPara seg_para; //一大段的参数

		vector<vector<Vector2d> > UHpair; //每一大段下分小段
		vector<Vector2d> pair_seg; //每一小段的pair
		pair_seg.push_back(Vector2d(0.0, pointsseg.front()(2)));
		for (int j = 1; j < pointnum; ++j)
		{
			pair_seg.push_back(Vector2d(useg[j], pointsseg[j](2)));
			if (pair_seg.size() > 19)
			{
				Vector2d cache_h = pair_seg.back();
				UHpair.push_back(pair_seg);
				pair_seg.clear();
				pair_seg.push_back(cache_h);
			}
		}

		//最后一小段的情况
		if (pair_seg.size() > 4)
		{
			UHpair.push_back(pair_seg);
		}
		else if (pair_seg.size() == 1)
		{		
		}
		else
		{
			for (int j = 1; j < pair_seg.size(); ++j)
				UHpair.back().push_back(pair_seg[j]);
		}

		// output
		for (int m = 0; m < UHpair.size(); ++m)
		{
			for (int j = 0; j < UHpair[m].size(); ++j)
			{
				ofs << i << "," << m << "," << j << "," << UHpair[m][j](0) << "," << UHpair[m][j](1) << endl;
			}
		}
		ofs << endl;

		// 对于每一小段有一个初始高程
		for (int m = 0; m < UHpair.size(); ++m)
		{
			// 每段solve一次
			Problem problem;
			double greek[2] = { 0.0, 0.0 };
			double u00 = UHpair[m][0](0);
			seg_para.u_start.push_back(u00);
			double z00 = UHpair[m][0](1);
			Vector3d greek_v(z00, 0.0, 0.0);
			for (int n = 0; n < UHpair[m].size(); ++n)
			{
				double zz = UHpair[m][n](1);
				double uu = UHpair[m][n](0) - u00;
				CostFunction *cost =
					new AutoDiffCostFunction<WeightedZDistanceFromCurveCost_0, 1, 2>(
					new WeightedZDistanceFromCurveCost_0(zz, uu, z00, 1.0));
				problem.AddResidualBlock(cost, NULL, greek);

				if (n == UHpair[m].size()-1)
				{
					CostFunction *cost =
						new AutoDiffCostFunction<WeightedZDistanceFromCurveCost_0, 1, 2>(
						new WeightedZDistanceFromCurveCost_0(zz, uu, z00, 10.0));
					problem.AddResidualBlock(cost, NULL, greek);
				}
			}
			problem.SetParameterLowerBound(greek, 0, -0.2);  problem.SetParameterUpperBound(greek, 0, 0.2);
			problem.SetParameterLowerBound(greek, 1, -0.01); problem.SetParameterUpperBound(greek, 1, 0.01);
			Solver::Options options;
			options.max_num_iterations = 20;
			options.linear_solver_type = ceres::DENSE_QR;
			options.minimizer_progress_to_stdout = false;

			Solver::Summary summary;
			Solve(options, &problem, &summary);

			std::cout << summary.BriefReport() << "\n";

			greek_v(1) = greek[0]; greek_v(2) = greek[1];
			
			seg_para.para.push_back(greek_v);
		}

		para.push_back(seg_para);
	}
	ofs.close();
	return para;
}

//求直线方向2D
double Solvepara::ccltDirectionConstraint2D(const vector<Vector3d> &points)
{
	Problem problem;
	double para[2] = { 1.0, 0.0 };

	for (int i = 0; i < points.size(); ++i)
	{
		double xx = points[i](0);
		double yy = points[i](1);

		CostFunction *cost =
			new AutoDiffCostFunction<DirectionCost2D, 1, 2>(
			new DirectionCost2D(xx, yy));
		problem.AddResidualBlock(cost, new CauchyLoss(0.1), para);
	}

	Solver::Options options;
	options.max_num_iterations = 1000;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = false;

	Solver::Summary summary;
	Solve(options, &problem, &summary);

	if (para[0] < 0)
		para[0] += M_PI;

	return para[0];
}

double Solvepara::PCAccltDirectionConstraint2D(const vector<Vector3d> &points)
{
	CvMat *pData = cvCreateMat(points.size(), 3, CV_32FC1);
	CvMat *pMean = cvCreateMat(1, 3, CV_32FC1);
	CvMat *pEigVals = cvCreateMat(1, 3, CV_32FC1);
	CvMat *pEigVecs = cvCreateMat(3, 3, CV_32FC1);

	for (int k = 0; k < points.size(); ++k) 
	{
		cvmSet(pData, k, 0, points[k](0));
		cvmSet(pData, k, 1, points[k](1));
		cvmSet(pData, k, 2, 0.0);
	}
	cvCalcPCA(pData, pMean, pEigVals, pEigVecs, CV_PCA_DATA_AS_ROW);

	double xx = cvmGet(pEigVecs, 0, 0);
	double yy = cvmGet(pEigVecs, 0, 1);
	double zz = cvmGet(pEigVecs, 0, 2);

	double rst = atan(yy / xx);
	if (rst < 0)
	{
		rst += M_PI;
	}

	return rst;
}

// ransac提取圆弧
void Solvepara::ransacExtractCircles(const vector<Vector3d> &points)
{
	ofstream ofs("ransac.txt");
	// 拷贝点云
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_f(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_circle(new pcl::PointCloud<pcl::PointXYZ>);

	for (int i = 0; i < points.size(); i++)
	{
		PointXYZ p;
		p.x = points[i](0); p.y = points[i](1); p.z = 0.0;
		cloud->points.push_back(p);
	}
	cloud->width = cloud->points.size();
	cloud->height = 1;
	cloud->points.resize(cloud->width * cloud->height);

	//ransac
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
	pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_line(new pcl::PointCloud<pcl::PointXYZ>());

	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_CIRCLE2D);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setMaxIterations(1000);
	seg.setDistanceThreshold(0.1);

	int label = 0;
	while (cloud->points.size() > 10)
	{
		// Segment the largest planar component from the remaining cloud
		seg.setInputCloud(cloud);
		seg.segment(*inliers, *coefficients);
		if (inliers->indices.size() == 0)
		{
			std::cout << "Could not estimate a planar model for the given dataset." << std::endl;
			break; // 如果剩下的点中没有一个圆弧,则剩下的点都是outliers
		}

		// Extract the planar inliers from the input cloud
		pcl::ExtractIndices<pcl::PointXYZ> extract;
		extract.setInputCloud(cloud);
		extract.setIndices(inliers);
		extract.setNegative(false);

		// Get the points associated with the planar surface
		extract.filter(*cloud_circle);

		// Remove the planar inliers, extract the rest
		extract.setNegative(true);
		extract.filter(*cloud_f);
		*cloud = *cloud_f;

		for (int i = 0; i < cloud_circle->points.size(); i++)
		{
			PointXYZ pxy = cloud_circle->points[i];
			ofs << to_string(pxy.x + mean_point(0)) << "," << to_string(pxy.y + mean_point(1)) << "," << 0.0 << "," << label << endl;
		}
		label++;
	}
	ofs.close();
	cout << "ransac circle extraction complete" << endl;
}