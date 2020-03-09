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


//extern double g_iterator; // ����ʱ��uֵ�������жϵ�����ֹ
//extern Vector3d start_point_solve_greek;

struct Panel
{
	unsigned int major_iter_num; //�����ǵڼ�����ѭ��
};

struct ZPara
{
	vector<double> u_start; //ÿһ����е�ÿ��С�ε���㣬��һ����0.0
	vector<Vector3d> para;  //ÿһС�εĲ���
};

class Solvepara
{
public:
	Solvepara();
	~Solvepara();

	// ������壬������
	Panel panel_;

	Vector3d mean_point;
	
	// ��ȡ���е㣬ƽ�ƣ�����y�������� ���ƽ������
	Vector3d readPointsFromTxt(const string &filename, vector<Vector3d> &points, int &point_number_total);

	// �ֶ�
	void segmentPoints(const float &threshold, //�ֶγ�����ֵ
		const vector<Vector3d> &points, //����߽��
		vector<vector<Vector3d> > &points_segmented, //����߽��
		vector<vector<double> > &u_init, //�߽���Ӧ�ĳ�ʼu
		vector<Vector3d> &start_points); //ÿһ�ε���ʼ��

	void segmentPointsByKnownU(const float &threshold, const vector<Vector3d> &points, const vector<double> &u, vector<vector<Vector3d> > &points_segmented,
		vector<vector<double> > &u_init, vector<Vector3d> &start_points);

	// ����para���Ż�
	void solveParaByFixedU(const vector<vector<double> > &input_u, const vector<Vector3d> &greek_input, const vector<Vector2d> &greek_v_input,
		const vector<vector<Vector3d> > &data_points, vector<Vector3d> &output_greek, vector<Vector2d> &output_greek_v);
	
	// �̷߳���
	vector<double> ZsolveParaByFixedU(const vector<vector<double> > &input_u, const vector<vector<Vector3d> > &data_points, vector<Vector2d> &output_greek_v);

	// ����U
	void solveUByFixedPara(const vector<vector<double> > &input_u, const vector<Vector3d> &greek, const vector<Vector2d> &greek_v,
		const vector<vector<Vector3d> > &data_points, vector<vector<double> > &output_u);

	// ���Ż���ֶε�uת��Ϊȫ�ֵ�u
	vector<double> translateUtoGlobal(const vector<vector<double> > &u);

	// ÿһ�θ߳��Ż�
	vector<ZPara> optimizeHeight(const vector<vector<Vector3d> > &points, const vector<vector<double> > &u);

	//
	vector<double> ccltGlobalBaseU(const vector<vector<double> > &u);

private:
	vector<Vector3d> ccltIntegralVector(const Vector3d &start_point, const Vector3d &greek, const Vector2d &greek_v, const double &length);
	long ccltNearestPointIndex(const Vector3d &ref_point, const vector<Vector3d> &curve);

	// ����ÿһ�ν�������ĩβ��ʱ���directionָ�򣬷��ؽǶ�ֵ����Χ0-2*pi
	double ccltDirectionConstraint2D(const vector<Vector3d> &points);

	// ������һ������Ŀ����ͬ����PCA
	double PCAccltDirectionConstraint2D(const vector<Vector3d> &points);

	Vector2d ccltBeginEndDir2D(const vector<vector<Vector3d> > &data_points, const int &i);
};

/************************************************************************/
/*      STEP SIZE                  !!!!                                 */
/************************************************************************/
// û��Ȩ�صĿ��Ƶ�
class DistanceFromCurveCost {
public:
	DistanceFromCurveCost(double xx, double yy, double u, double x0, double y0) : xx_(xx), yy_(yy), u_(u), x0_(x0), y0_(y0){}
	template <typename T> bool operator()(
		const T* const greek,		// ���Ż��Ĳ���greek
		T* residual) const {		
		// ���ּ���
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
	// u�ǻ�����ֹ���
	double xx_, yy_, u_, x0_, y0_;
};

// ���ڿ������ӵ�����Լ��
// ��Ȩ�صĿ��Ƶ�
class WeightedDistanceFromCurveCost_0 {
public:
	WeightedDistanceFromCurveCost_0(double xx, double yy, double u, double x0, double y0, double w) : xx_(xx), yy_(yy), u_(u), x0_(x0), y0_(y0), w_(w){}
	template <typename T> bool operator()(
		const T* const greek,		// ���Ż��Ĳ���greek
		T* residual) const {
		// ���ּ���
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
	// u�ǻ�����ֹ���
	// w��Ȩ��
	double xx_, yy_, u_, x0_, y0_, w_;
};

class ZDistanceFromCurveCost {
public:
	ZDistanceFromCurveCost(double zz, double u, double z0) : zz_(zz), u_(u), z0_(z0){}
	template <typename T> bool operator()(
		const T* const greek,		// ���Ż��Ĳ���greek
		T* residual) const {

		T z_predicted = T(z0_) + T(u_)*greek[0] + 0.5*T(u_)*T(u_)*greek[1];

		residual[0] = T(zz_) - z_predicted;
		return true;
	}

private:
	// u�ǻ�����ֹ���
	double zz_, u_, z0_;
};

class WeightedZDistanceFromCurveCost_0 {
public:
	WeightedZDistanceFromCurveCost_0(double zz, double u, double z0, double w) : zz_(zz), u_(u), z0_(z0), w_(w){}
	template <typename T> bool operator()(
		const T* const greek,		// ���Ż��Ĳ���greek
		T* residual) const {

		T z_predicted = T(z0_) + T(u_)*greek[0] + 0.5*T(u_)*T(u_)*greek[1];

		residual[0] = (T(zz_) - z_predicted)*ceres::sqrt(T(w_));
		return true;
	}

private:
	// u�ǻ�����ֹ���
	double zz_, u_, z0_, w_;
};

// ֱ�����
class DirectionCost2D {
public:
	DirectionCost2D(double xx, double yy) : xx_(xx), yy_(yy){}
	template <typename T> bool operator()(
		const T* const para,		// ���Ż��Ĳ���,ֱ�߲���2dimension
		T* residual) const {

		T y_predicted = ceres::tan(para[0])*T(xx_) + para[1];

		residual[0] = T(yy_) - y_predicted;
		return true;
	}

private:
	// u�ǻ�����ֹ���
	double xx_, yy_;
};

// ���μ�ķ�������
class DirectionConnection2D {
public:
	DirectionConnection2D(double dir1, double dir2, double length) : dir1_(dir1), dir2_(dir2), length_(length){}

	template <typename T> bool operator()(
		const T* const para,		// ���Ż��Ĳ������������ߵĲ���3dimension
		T* residual) const {

		T dir_begin = para[0];
		T dir_end = para[0] + para[1] * T(length_) + para[2] * T(length_)*T(length_)*0.5;

		residual[0] = (T(dir1_) - dir_begin)*ceres::sqrt(100);
		residual[1] = (T(dir2_) - dir_end)*ceres::sqrt(100);

		return true;
	}

private:
	// u�ǻ�����ֹ���
	double dir1_, dir2_, length_;
};

bool compPointByY(const Vector3d &a, const Vector3d &b);
