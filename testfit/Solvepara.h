#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <ceres/ceres.h>
#include <glog/logging.h>
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using namespace std;
using namespace Eigen;


//extern double g_iterator; // 即此时的u值，用于判断迭代终止
//extern Vector3d start_point_solve_greek;

struct Panel
{
	unsigned int major_iter_num; //现在是第几个大循环
};

struct ZPara
{
	vector<double> u_start; //每一大段中的每个小段的起点，第一个是0.0
	vector<Vector3d> para;  //每一小段的参数
};

class Solvepara
{
public:
	Solvepara();
	~Solvepara();

	// 控制面板，传参用
	Panel panel_;

	Vector3d mean_point;
	
	// 读取所有点，平移，根据y坐标排序。 输出平均坐标
	Vector3d readPointsFromTxt(const string &filename, vector<Vector3d> &points, int &point_number_total);

	// 分段
	void segmentPoints(const float &threshold, //分段长度阈值
		const vector<Vector3d> &points, //输入边界点
		vector<vector<Vector3d> > &points_segmented, //输出边界段
		vector<vector<double> > &u_init, //边界点对应的初始u
		vector<Vector3d> &start_points); //每一段的起始点

	void segmentPointsByKnownU(const float &threshold, const vector<Vector3d> &points, const vector<double> &u, vector<vector<Vector3d> > &points_segmented,
		vector<vector<double> > &u_init, vector<Vector3d> &start_points);

	// 计算para，优化
	void solveParaByFixedU(const vector<vector<double> > &input_u, const vector<Vector3d> &greek_input, const vector<Vector2d> &greek_v_input,
		const vector<vector<Vector3d> > &data_points, vector<Vector3d> &output_greek, vector<Vector2d> &output_greek_v);
	
	// 高程方向
	vector<double> ZsolveParaByFixedU(const vector<vector<double> > &input_u, const vector<vector<Vector3d> > &data_points, vector<Vector2d> &output_greek_v);

	// 计算U
	void solveUByFixedPara(const vector<vector<double> > &input_u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<vector<Vector3d> > &data_points, vector<vector<double> > &output_u);

	// 把优化后分段的u转换为全局的u
	vector<double> translateUtoGlobal(const vector<vector<double> > &u);

	// 每一段高程优化
	vector<ZPara> optimizeHeight(const vector<vector<Vector3d> > &points, const vector<vector<double> > &u);

	//
	vector<double> ccltGlobalBaseU(const vector<vector<double> > &u);

private:
	vector<Vector3d> ccltIntegralVector(const Vector3d &start_point, const Vector3d &greek, const Vector2d &greek_v, const double &length);
	long ccltNearestPointIndex(const Vector3d &ref_point, const vector<Vector3d> &curve);
};

// 没有权重的控制点
class DistanceFromCurveCost {
public:
	DistanceFromCurveCost(double xx, double yy, double u, double x0, double y0) : xx_(xx), yy_(yy), u_(u), x0_(x0), y0_(y0){}
	template <typename T> bool operator()(
		const T* const greek,		// 待优化的参数greek
		T* residual) const {		
		// 积分计算
		double g_iterator = u_;
		double start_point_solve_greek_x = x0_;
		double start_point_solve_greek_y = y0_;

		T x_predicted = ceres::cos(greek[0])*0.5*0.001 + start_point_solve_greek_x;
		T y_predicted = ceres::sin(greek[0])*0.5*0.001 + start_point_solve_greek_y;

		for (long i = 0; i*0.001 < g_iterator; ++i)
		{
			x_predicted += 2 * 0.5 * 0.001*ceres::cos(greek[0] + greek[1] * T(i)*0.001 + 0.5*greek[2] * T(i)*T(i)*0.001*0.001);
			y_predicted += 2 * 0.5 * 0.001*ceres::sin(greek[0] + greek[1] * T(i)*0.001 + 0.5*greek[2] * T(i)*T(i)*0.001*0.001);
		}

		double f_step = g_iterator - floor(1000 * g_iterator) / 1000.0;
		x_predicted += 0.5 * T(f_step)*ceres::cos(greek[0] + greek[1] * T(g_iterator) + 0.5*greek[2] * T(g_iterator)*T(g_iterator));
		y_predicted += 0.5 * T(f_step)*ceres::sin(greek[0] + greek[1] * T(g_iterator) + 0.5*greek[2] * T(g_iterator)*T(g_iterator));

		residual[0] = T(xx_) - x_predicted;
		residual[1] = T(yy_) - y_predicted;
		return true;
	}

private:
	// u是积分终止里程
	double xx_, yy_, u_, x0_, y0_;
};

// 用于控制连接点连接约束
// 有权重的控制点
class WeightedDistanceFromCurveCost_0 {
public:
	WeightedDistanceFromCurveCost_0(double xx, double yy, double u, double x0, double y0, double w) : xx_(xx), yy_(yy), u_(u), x0_(x0), y0_(y0), w_(w){}
	template <typename T> bool operator()(
		const T* const greek,		// 待优化的参数greek
		T* residual) const {
		// 积分计算
		double g_iterator = u_;
		double start_point_solve_greek_x = x0_;
		double start_point_solve_greek_y = y0_;

		T x_predicted = ceres::cos(greek[0])*0.5*0.001 + start_point_solve_greek_x;
		T y_predicted = ceres::sin(greek[0])*0.5*0.001 + start_point_solve_greek_y;

		for (long i = 0; i*0.001 < g_iterator; ++i)
		{
			x_predicted += 2 * 0.5 * 0.001*ceres::cos(greek[0] + greek[1] * T(i)*0.001 + 0.5*greek[2] * T(i)*T(i)*0.001*0.001);
			y_predicted += 2 * 0.5 * 0.001*ceres::sin(greek[0] + greek[1] * T(i)*0.001 + 0.5*greek[2] * T(i)*T(i)*0.001*0.001);
		}

		double f_step = g_iterator - floor(1000 * g_iterator) / 1000.0;
		x_predicted += 0.5 * T(f_step)*ceres::cos(greek[0] + greek[1] * T(g_iterator) + 0.5*greek[2] * T(g_iterator)*T(g_iterator));
		y_predicted += 0.5 * T(f_step)*ceres::sin(greek[0] + greek[1] * T(g_iterator) + 0.5*greek[2] * T(g_iterator)*T(g_iterator));

		residual[0] = (T(xx_) - x_predicted)*ceres::sqrt(T(w_));
		residual[1] = (T(yy_) - y_predicted)*ceres::sqrt(T(w_));
		return true;
	}

private:
	// u是积分终止里程
	// w是权重
	double xx_, yy_, u_, x0_, y0_, w_;
};

class ZDistanceFromCurveCost {
public:
	ZDistanceFromCurveCost(double zz, double u, double z0) : zz_(zz), u_(u), z0_(z0){}
	template <typename T> bool operator()(
		const T* const greek,		// 待优化的参数greek
		T* residual) const {

		T z_predicted = T(z0_) + T(u_)*greek[0] + 0.5*T(u_)*T(u_)*greek[1];

		residual[0] = T(zz_) - z_predicted;
		return true;
	}

private:
	// u是积分终止里程
	double zz_, u_, z0_;
};

class WeightedZDistanceFromCurveCost_0 {
public:
	WeightedZDistanceFromCurveCost_0(double zz, double u, double z0, double w) : zz_(zz), u_(u), z0_(z0), w_(w){}
	template <typename T> bool operator()(
		const T* const greek,		// 待优化的参数greek
		T* residual) const {

		T z_predicted = T(z0_) + T(u_)*greek[0] + 0.5*T(u_)*T(u_)*greek[1];

		residual[0] = (T(zz_) - z_predicted)*ceres::sqrt(T(w_));
		return true;
	}

private:
	// u是积分终止里程
	double zz_, u_, z0_, w_;
};

bool compPointByY(const Vector3d &a, const Vector3d &b);
