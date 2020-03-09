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

	// 计算每一段结束或者末尾的时候的direction指向，返回角度值，范围0-2*pi
	double ccltDirectionConstraint2D(const vector<Vector3d> &points);

	// 和上面一个函数目的相同，用PCA
	double PCAccltDirectionConstraint2D(const vector<Vector3d> &points);

	Vector2d ccltBeginEndDir2D(const vector<vector<Vector3d> > &data_points, const int &i);
};

/************************************************************************/
/*      STEP SIZE                  !!!!                                 */
/************************************************************************/
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
		double step_ = 0.01; //ATTENTION

		T x_predicted = ceres::cos(greek[0])*0.5*step_ + start_point_solve_greek_x;
		T y_predicted = ceres::sin(greek[0])*0.5*step_ + start_point_solve_greek_y;

		for (long i = 0; i*step_ < g_iterator; ++i)
		{
			x_predicted += 2 * 0.5 * step_*ceres::cos(greek[0] + greek[1] * T(i)*step_ + 0.5*greek[2] * T(i)*T(i)*step_*step_);
			y_predicted += 2 * 0.5 * step_*ceres::sin(greek[0] + greek[1] * T(i)*step_ + 0.5*greek[2] * T(i)*T(i)*step_*step_);
		}

		double f_step = g_iterator - floor(100 * g_iterator) / 100.0;
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
		double step_ = 0.01; //ATTENTION

		T x_predicted = ceres::cos(greek[0])*0.5*step_ + start_point_solve_greek_x;
		T y_predicted = ceres::sin(greek[0])*0.5*step_ + start_point_solve_greek_y;

		for (long i = 0; i*step_ < g_iterator; ++i)
		{
			x_predicted += 2 * 0.5 * step_*ceres::cos(greek[0] + greek[1] * T(i)*step_ + 0.5*greek[2] * T(i)*T(i)*step_*step_);
			y_predicted += 2 * 0.5 * step_*ceres::sin(greek[0] + greek[1] * T(i)*step_ + 0.5*greek[2] * T(i)*T(i)*step_*step_);
		}

		double f_step = g_iterator - floor(100 * g_iterator) / 100.0;
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

// 直线拟合
class DirectionCost2D {
public:
	DirectionCost2D(double xx, double yy) : xx_(xx), yy_(yy){}
	template <typename T> bool operator()(
		const T* const para,		// 待优化的参数,直线参数2dimension
		T* residual) const {

		T y_predicted = ceres::tan(para[0])*T(xx_) + para[1];

		residual[0] = T(yy_) - y_predicted;
		return true;
	}

private:
	// u是积分终止里程
	double xx_, yy_;
};

// 两段间的方向连接
class DirectionConnection2D {
public:
	DirectionConnection2D(double dir1, double dir2, double length) : dir1_(dir1), dir2_(dir2), length_(length){}

	template <typename T> bool operator()(
		const T* const para,		// 待优化的参数，缓和曲线的参数3dimension
		T* residual) const {

		T dir_begin = para[0];
		T dir_end = para[0] + para[1] * T(length_) + para[2] * T(length_)*T(length_)*0.5;

		residual[0] = (T(dir1_) - dir_begin)*ceres::sqrt(100);
		residual[1] = (T(dir2_) - dir_end)*ceres::sqrt(100);

		return true;
	}

private:
	// u是积分终止里程
	double dir1_, dir2_, length_;
};

bool compPointByY(const Vector3d &a, const Vector3d &b);
