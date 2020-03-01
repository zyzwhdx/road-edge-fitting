#define GLOG_NO_ABBREVIATED_SEVERITIES
#include "point.h"
#include <fstream>
#include <string>
#include <vector>
#include "dataIo.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <numeric>
#include <opencv2/core/core.hpp>
#include <ceres/ceres.h>
#include <pcl/io/pcd_io.h>
#include <string>
#include <pcl/kdtree/kdtree_flann.h>
using namespace std;
using namespace pcl;

void readVertices(string path, vector<vector<zPoint> > &vertices, vector<vector<zSegment> > &segments, double &ox, double &oy)
{
	vertices.clear();
	segments.clear();
	ifstream ifs(path);
	string str;
	vector<string> vec;

	//第一行
	getline(ifs, str);
	boost::split(vec, str, boost::is_any_of(", "), boost::token_compress_on);
	ox = atof(vec[0].c_str());
	oy = atof(vec[1].c_str());

	//第二行
	getline(ifs, str);
	int num = atoi(str.c_str());
	vertices.resize(num);
	segments.resize(num);

	//后面的若干行
	while (getline(ifs, str))
	{
		vector<string> sVec;
		boost::split(sVec, str, boost::is_any_of(", "), boost::token_compress_on);

		int cLabel = atoi(sVec[3].c_str());
		zPoint pt;
		pt.x = atof(sVec[0].c_str());
		pt.y = atof(sVec[1].c_str());
		vertices[cLabel].push_back(pt);
	}

	for (int i = 0; i < vertices.size(); i++)
	{
		segments[i].resize(vertices[i].size());
		for (int j = 0; j < segments[i].size() - 1; j++)
		{
			segments[i][j].src = vertices[i][j];
			segments[i][j].dst = vertices[i][j + 1];
		}
		segments[i][segments[i].size() - 1].src = vertices[i][segments[i].size() - 1];
		segments[i][segments[i].size() - 1].dst = vertices[i][0];
	}
}

void readCheckPoints(string path, vector<zPoint> &checkpts, double &ox, double &oy)
{
	checkpts.clear();
	ifstream ifs(path);
	string str;
	
	while (getline(ifs, str))
	{
		vector<string> sVec;
		boost::split(sVec, str, boost::is_any_of(", "), boost::token_compress_on);

		zPoint pt;
		pt.x = atof(sVec[0].c_str()) - ox;
		pt.y = atof(sVec[1].c_str()) - oy;
		checkpts.push_back(pt);
	}
}

double ccltDistancePt2Sg(zPoint p, zPoint m, zPoint n)
{
	//pm点乘mn
	double pmDmn = (m.x - p.x)*(n.x - m.x) + (m.y - p.y)*(n.y - m.y);
	//pn点乘mn
	double pnDmn = (n.x - p.x)*(n.x - m.x) + (n.y - p.y)*(n.y - m.y);

	//pm边长
	double pm = sqrt((m.x - p.x)*(m.x - p.x) + (m.y - p.y)*(m.y - p.y));
	//pn边长
	double pn = sqrt((n.x - p.x)*(n.x - p.x) + (n.y - p.y)*(n.y - p.y));
	//mn边长
	double mn = sqrt((n.x - m.x)*(n.x - m.x) + (n.y - m.y)*(n.y - m.y));


	double small_edge = pm;
	if (pn < pm)
	{
		small_edge = pn;
	}

	if (pmDmn * pnDmn >= 1e-6)
	{
		return small_edge;
	}
	else
	{
		double x1 = m.x;
		double y1 = m.y;
		double x2 = n.x;
		double y2 = n.y;
		double x3 = p.x;
		double y3 = p.y;
		double s = fabs(x1*y2 + x2*y3 + x3*y1 - x1*y3 - x2*y1 - x3*y2);
		return fabs(s / mn);
	}
}

void ccltCheckDistances(vector<zPoint> checks, vector<vector<zSegment> > segments, vector<vector<double> > &distances)
{
	distances.clear();
	distances.resize(checks.size());
	for (int i = 0; i < checks.size(); i++)
	{
		zPoint ckPoint = checks[i];
		for (int m = 0; m < segments.size();m++)
		{
			for (int j = 0; j < segments[m].size(); j++)
			{
				double distance = ccltDistancePt2Sg(ckPoint, segments[m][j].src, segments[m][j].dst);
				distances[i].push_back(distance);
			}
		}
		sort(distances[i].begin(), distances[i].end());
	}
}

void preRead(string path, eString &pts, eString &pts2, double &ox, double &oy)
{
	ifstream ifs(path);
	string str;

	vector<double> xx, yy;
	while (getline(ifs, str))
	{
		vector<string> sVec;
		boost::split(sVec, str, boost::is_any_of(", "), boost::token_compress_on);

		ePoint pt;
		pt.x = atof(sVec[0].c_str()); xx.push_back(pt.x);
		pt.y = atof(sVec[1].c_str()); yy.push_back(pt.y);
		pt.z = atof(sVec[2].c_str());
		pts.push_back(pt);
	}

	ox = medianNum(xx);
	oy = medianNum(yy);

	for (size_t i = 0; i < pts.size(); i++)
	{
		pts[i].id = i;
		pts[i].x -= ox;
		pts[i].y -= oy;
	}

	eString all = pts;
	eString right;
	eString seeds;

	eEdge segments;

	while (all.size() > 0)
	{
		right.push_back(all.front());
		seeds.push_back(all.front());
		all.erase(all.begin());

		while (seeds.size() > 0)
		{
			for (int j = 0; j < all.size(); j++)
			{
				if (ccltPtDistance2D(seeds.front(), all[j]) < 8)
				{
					right.push_back(all[j]);
					seeds.push_back(all[j]);
					all.erase(all.begin() + j);
				}
			}
			seeds.erase(seeds.begin());
		}
		segments.push_back(right);
		right.clear();
	}

	sort(segments.begin(), segments.end(), compSize);

	pts = segments[0];
	pts2 = segments[1];

	ofstream ii("right.txt");
	for (int i = 0; i < pts.size(); i++)
	{
		ii << to_string(pts[i].x) << "," << to_string(pts[i].y) << endl;
	}

	ofstream iii("left.txt");
	for (int i = 0; i < pts2.size(); i++)
	{
		iii << to_string(pts2[i].x) << "," << to_string(pts2[i].y) << endl;
	}
	ii.close();
	iii.close();
}

void findCenter(eString &left, eString &right, eString &center)
{
	if (right.size() > left.size())
	{
		eString tmp = right;
		right = left;
		left = tmp;
	}

	for (int i = 0; i < left.size(); i++)
	{
		vector<eDistance > distances;
		for (int j = 0; j < right.size(); j++)
		{
			eDistance dis;
			dis.xx = right[j].x;
			dis.yy = right[j].y;
			dis.zz = right[j].z;
			dis.dis = ccltPtDistance2D(left[i], right[j]);
			distances.push_back(dis);
		}
		sort(distances.begin(), distances.end(), compdis);
		
		eDistance nearest = distances.front();
		ePoint ct;
		ct.x = 0.5*(nearest.xx + left[i].x);
		ct.y = 0.5*(nearest.yy + left[i].y);
		ct.z = 0.5*(nearest.zz + left[i].z);

		center.push_back(ct);
	}
	//ofstream ii("center.txt");
	//for (int i = 0; i < center.size(); i++)
	//{
	//	ii << to_string(center[i].x) << "," << to_string(center[i].y) << "," << to_string(center[i].z) << endl;
	//}
	//ii.close();
}


double ccltPtDistance2D(ePoint a, ePoint b)
{
	double dx = a.x - b.x;
	double dy = a.y - b.y;
	double d2 = dx*dx + dy*dy;
	return sqrt(d2);
}


void readEdges(string path, eString &pts, double &ox, double &oy, extent_dir &edir)
{
	ifstream ifs(path);
	string str;

	vector<double> xx, yy;
	while (getline(ifs, str))
	{
		vector<string> sVec;
		boost::split(sVec, str, boost::is_any_of(", "), boost::token_compress_on);

		ePoint pt;
		pt.x = atof(sVec[0].c_str()); xx.push_back(pt.x);
		pt.y = atof(sVec[1].c_str()); yy.push_back(pt.y);
		pt.z = atof(sVec[2].c_str());
		pts.push_back(pt);
	}

	ox = medianNum(xx);
	oy = medianNum(yy);

	double range_x = fabs(xx.back() - xx.front());
	double range_y = fabs(yy.back() - yy.front());

	if (range_y > range_x)
	{
		//按y值从小到大排序
		sort(pts.begin(), pts.end(), compy);
		edir = e_y;
	}
	else
	{
		//按x值从小到大排序
		sort(pts.begin(), pts.end(), compx);
		edir = e_x;
	}

	for (size_t i = 0; i < pts.size(); i++)
	{
		pts[i].id = i;
		pts[i].x -= ox;
		pts[i].y -= oy;
	}
}

void readEdges(eString &center, extent_dir &edir)
{
	vector<double > xx, yy;
	for (int i = 0; i < center.size(); i++)
	{
		xx.push_back(center[i].x);
		yy.push_back(center[i].y);
	}
	sort(xx.begin(), xx.end());
	sort(yy.begin(), yy.end());

	double range_x = fabs(xx.back() - xx.front());
	double range_y = fabs(yy.back() - yy.front());

	if (range_y > range_x)
	{
		//按y值从小到大排序
		sort(center.begin(), center.end(), compy);
		edir = e_y;
	}
	else
	{
		//按x值从小到大排序
		sort(center.begin(), center.end(), compx);
		edir = e_x;
	}

	for (size_t i = 0; i < center.size(); i++)
	{
		center[i].id = i;
	}
}

bool compy(ePoint &a, ePoint &b)
{
	if (a.y < b.y)
	{
		return true;
	}
	return false;
}

bool compx(ePoint &a, ePoint &b)
{
	if (a.x < b.x)
	{
		return true;
	}
	return false;
}

bool compid(ePoint &a, ePoint &b)
{
	if (a.id < b.id)
	{
		return true;
	}
	return false;
}

bool compdis(eDistance &a, eDistance &b)
{
	if (a.dis < b.dis)
	{
		return true;
	}
	return false;
}

bool compSize(eString &a, eString &b)
{
	if (a.size()<b.size())
	{
		return true;
	}
	return false;
}

void iterativeBreakLines(eString &oPts, eEdge &whole_edge, int th)
{
	if (oPts.size() < 30)
	{
		whole_edge.push_back(oPts);
		return;
	}
	int pos;

	if (ccltMaxOffset(oPts, pos) < th)
	{
		whole_edge.push_back(oPts);
		return;
	}

	eString a;
	eString b;
	for (size_t i = 0; i <= pos - oPts.front().id; i++)
	{
		a.push_back(oPts[i]);
	}
	for (size_t i = pos - oPts.front().id + 1; i < oPts.size(); i++)
	{
		b.push_back(oPts[i]);
	}
	iterativeBreakLines(a, whole_edge, th);
	iterativeBreakLines(b, whole_edge, th);
}


double ccltMaxOffset(eString &pts, int &pos_maxoft)
{
	if (pts.size() < 4)
		return 0.0;

	sort(pts.begin(), pts.end(), compid);

	zPoint source;
	zPoint target;

	source.x = pts[0].x;
	source.y = pts[0].y;
	target.x = pts[pts.size() - 1].x;
	target.y = pts[pts.size() - 1].y;

	double max_ost = 0.0;
	int max_id = 0;

	for (size_t i = 1; i < pts.size() - 1; i++)
	{
		zPoint cPt;
		cPt.x = pts[i].x;
		cPt.y = pts[i].y;
		double distance = ccltDistancePt2Sg(cPt, source, target);
		if (distance > max_ost)
		{
			max_ost = distance;
			max_id = pts[i].id;
		}
	}
	pos_maxoft = max_id;
	return max_ost;
}

int zrand()
{
	srand(int(time(0)));
	return rand() % 99 + 1;
}

double medianNum(vector<double> &a)
{
	int N = a.size();
	sort(a.begin(), a.end());
	if (N % 2 != 0)
	{
		return a.at(N / 2);
	}
	else
	{
		return (a.at(N / 2) + a.at(N / 2 - 1)) / 2.0;
	}
}

void ccltNormalPara(eString &a, normalPara &para, const extent_dir &edir)
{
	vector<double> xx;
	if (edir == e_y)
	{
		for (int i = 0; i < a.size(); i++)
		{
			xx.push_back(a.at(i).y);
		}
	}
	else
	{
		for (int i = 0; i < a.size(); i++)
		{
			xx.push_back(a.at(i).x);
		}
	}


	double sum = accumulate(begin(xx), end(xx), 0.0);
	double mean = sum / xx.size();

	double accum = 0.0;
	for_each(begin(xx), end(xx), [&](const double d)
	{
		accum += (d - mean)*(d - mean);
	});

	double stv = sqrt(accum / (xx.size() - 1));

	para.ave = mean;
	para.stv = stv;
}

struct function_residual
{
	function_residual(double x_, double y_) : x(x_), y(y_){}

	template<typename T>
	bool operator()(const T* const abcde, T* residual) const
	{
		residual[0] = T(y) - (abcde[0] * T(x)*T(x)*T(x)*T(x) +
			abcde[1] * T(x)*T(x)*T(x) +
			abcde[2] * T(x)*T(x) +
			abcde[3] * T(x) +
			abcde[4]);//这里记得要将x和y进行T类型转换！！！
		return true;//这里记得return!!!!!!!!
	}

	const double x, y;
};

void solveFit(const vector<double> &xx, const vector<double> &yy, vector<double> &result)
{
	double abcde[5] = { 0, 0, 0, 0, 0 };
	ceres::Problem problem;//构建问题

	long N = xx.size();
	for (int j = 0; j < N; j++)
	{
		//这里也可以直接 auto cost_function = ...。如果用下方形式，一定记得是CostFunction*指针
		ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<function_residual, 1, 5>(new function_residual(xx[j], yy[j]));
		//这里CostFunction和AutoDiffCostFunction和function_residual都是类，类的传递用指针。
		//其中AutoDiffCostFunction还是模板类，要传递函数误差的类型，误差维度，估计参数维度

		//这里调用AddResidualBlock()方法，将代价函数，核函数，待估计参数传入
		problem.AddResidualBlock(cost_function, nullptr, abcde);
	}

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	//此函数定义：CERES_EXPORT void Solve(const Solver::Options& options, Problem* problem, Solver::Summary* summary);
	ceres::Solve(options, &problem, &summary);//注意这里是solve，不是solver了。options为引用，&problem,和　&summary为指针

	cout << summary.BriefReport() << endl;
	cout << "abcde=";
	for (auto item : abcde)
		cout << item << " ";

	for (int i = 0; i < 5;i++)
	{
		result.push_back(abcde[i]);
	}
}

void fitString(const eString &pts, polyPara &rst, polyPara &zrst, const normalPara &para, const extent_dir &edir)
{
	vector<double > xx;
	vector<double > yy;
	vector<double > zz;
	if (edir == e_y)
	{
		for (long i = 0; i < pts.size(); i++)
		{
			yy.push_back(pts[i].x);
			xx.push_back((pts[i].y - para.ave) / para.stv);
			zz.push_back(pts[i].z);
		}
	}
	else
	{
		for (long i = 0; i < pts.size(); i++)
		{
			xx.push_back((pts[i].x - para.ave) / para.stv);
			yy.push_back(pts[i].y);
			zz.push_back(pts[i].z);
		}
	}

	vector<double> result;
	vector<double> zresult;
	solveFit(xx, yy, result);
	solveFit(xx, zz, zresult);
	rst.a = result[0];
	rst.b = result[1];
	rst.c = result[2];
	rst.d = result[3];
	rst.e = result[4];

	zrst.a = zresult[0];
	zrst.b = zresult[1];
	zrst.c = zresult[2];
	zrst.d = zresult[3];
	zrst.e = zresult[4];
}

void fitEdge(const eEdge &edge, const double &ox, const double &oy, const extent_dir &dir, PointCloud<PointXYZI> &Vector)
{
	PointCloud<PointXYZI> oc;
	PointXYZI poi;
	
	vector<ePoint > links;
	int lbuf = 3;

	for (int i = 0; i < edge.size();i++)
	{
		eString str = edge[i];
		normalPara nPara;
		polyPara pPara;
		polyPara zPara;

		ccltNormalPara(str, nPara, dir);

		fitString(str, pPara, zPara, nPara, dir);

		if (dir == e_y)
		{
			double xx = edge[i].front().y;
			double end_x = edge[i].back().y;
			double mean = nPara.ave;
			double stv = nPara.stv;

			double a = pPara.a;
			double b = pPara.b;
			double c = pPara.c;
			double d = pPara.d;
			double e = pPara.e;

			double m = zPara.a;
			double n = zPara.b;
			double o = zPara.c;
			double p = zPara.d;
			double q = zPara.e;

			if (i == 0)
			{
				end_x -= lbuf;
				double x = (end_x - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP;
				adP.x = y;
				adP.y = end_x;
				adP.z = z;

				links.push_back(adP);
			}
			else if (i == edge.size() - 1)
			{
				xx += lbuf;
				double x = (xx - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP;

				adP.x = y;
				adP.y = xx;
				adP.z = z;

				links.push_back(adP);
			}
			else
			{
				end_x -= lbuf;
				xx += lbuf;

				double x = (end_x - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP2;
				adP2.x = y;
				adP2.y = end_x;
				adP2.z = z;

				x = (xx - mean) / stv;
				y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP1;

				adP1.x = y;
				adP1.y = xx;
				adP1.z = z;

				links.push_back(adP1);
				links.push_back(adP2);
			}

			while (xx <= end_x)
			{
				double x = (xx - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				poi.x = y;
				poi.y = xx;
				poi.z = z;
				poi.intensity = 100;
				xx += 0.01;
				oc.push_back(poi);
			}
		}		
		else
		{
			double xx = edge[i].front().x;
			double end_x = edge[i].back().x;
			double mean = nPara.ave;
			double stv = nPara.stv;

			double a = pPara.a;
			double b = pPara.b;
			double c = pPara.c;
			double d = pPara.d;
			double e = pPara.e;

			double m = zPara.a;
			double n = zPara.b;
			double o = zPara.c;
			double p = zPara.d;
			double q = zPara.e;

			if (i == 0)
			{
				end_x -= lbuf;
				double x = (end_x - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP;
				adP.x = end_x;
				adP.y = y;
				adP.z = z;

				links.push_back(adP);
			}
			else if (i == edge.size() - 1)
			{
				xx += lbuf;
				double x = (xx - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP;

				adP.x = xx;
				adP.y = y;
				adP.z = z;

				links.push_back(adP);
			}
			else
			{
				end_x -= lbuf;
				xx += lbuf;

				double x = (end_x - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP2;
				adP2.x = end_x;
				adP2.y = y;
				adP2.z = z;

				x = (xx - mean) / stv;
				y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				ePoint adP1;

				adP1.x = xx;
				adP1.y = y;
				adP1.z = z;

				links.push_back(adP1);
				links.push_back(adP2);
			}

			while (xx <= end_x)
			{
				double x = (xx - mean) / stv;
				double y = a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
				double z = m*x*x*x*x + n*x*x*x + o*x*x + p*x + q;

				poi.x = xx;
				poi.y = y;
				poi.z = z;
				poi.intensity = 100;
				xx += 0.01;
				oc.push_back(poi);
			}
		}
	}
	for (int i = 0; i < edge.size() - 1; i++)
	{
		addPoint(oc, links[i * 2], links[i * 2 + 1]);
	}

	oc.width = oc.points.size();
	oc.height = 1;
	oc.is_dense = false;

	Vector = oc;

	dataIo ioo;
	ioo.writePointCloudIntoLasFile("vector.las", oc, ox, oy);
}

void addPoint(PointCloud<PointXYZI> &cloud, ePoint src, ePoint dst)
{
	double dx = fabs(src.x - dst.x);
	double dy = fabs(src.y - dst.y);
	if (dx > dy)
	{
		int iteration = dx * 100;

		for (int i = 1; i < iteration; i++)
		{
			PointXYZI pt;
			pt.x = src.x + 0.01*i;
			pt.y = src.y + 1.0*i*(dst.y - src.y) / iteration;
			pt.z = src.z + 1.0*i*(dst.z - src.z) / iteration;
			pt.intensity = 100;

			cloud.push_back(pt);
		}
	}
	else
	{
		int iteration = dy * 100;

		for (int i = 1; i < iteration; i++)
		{
			PointXYZI pt;
			pt.y = src.y + 0.01*i;
			pt.x = src.x + 1.0*i*(dst.x - src.x) / iteration;
			pt.z = src.z + 1.0*i*(dst.z - src.z) / iteration;
			pt.intensity = 100;

			cloud.push_back(pt);
		}
	}
}

void ccltLine(ePoint outP, const PointCloud<PointXYZI> &Vector, line2D &Line, line2D &cLine)
{
	eString pts;

	PointCloud<PointXY>::Ptr cloud(new PointCloud<PointXY>);
	cloud->width = Vector.points.size();
	cloud->height = 1;
	cloud->points.resize(cloud->width);

	for (size_t i = 0; i < cloud->points.size(); i++)
	{
		cloud->points[i].x = Vector.points[i].x;
		cloud->points[i].y = Vector.points[i].y;
	}

	KdTreeFLANN<PointXY> kdtree;
	kdtree.setInputCloud(cloud);

	PointXY sPt;
	sPt.x = outP.x;
	sPt.y = outP.y;

	int K = 1;
	vector<int> idx(K);
	vector<float> distance(K);

	kdtree.nearestKSearch(sPt, K, idx, distance);

	cout << endl;
	cout << sPt.x << "," << sPt.y << endl;
	cout << Vector.points[idx[0]].x << "," << Vector.points[idx[0]].y << "," << Vector.points[idx[0]].z << endl;

	ePoint cPt;
	cPt.x = Vector.points[idx[0]].x;
	cPt.y = Vector.points[idx[0]].y;
	cPt.z = Vector.points[idx[0]].z;

	sPt.x = Vector.points[idx[0]].x;
	sPt.y = Vector.points[idx[0]].y;
	K = 2;
	vector<int> idx2(K);
	vector<float> distance2(K);

	kdtree.nearestKSearch(sPt, K, idx2, distance2);

	//ePoint cPt;
	//cPt.x = Vector.points[idx2[0]].x;
	//cPt.y = Vector.points[idx2[0]].y;
	//cPt.z = Vector.points[idx2[0]].z;
	
	pts.push_back(cPt);
	cPt.x = Vector.points[idx2[1]].x;
	cPt.y = Vector.points[idx2[1]].y;
	cPt.z = Vector.points[idx2[1]].z;
	pts.push_back(cPt);

	if (fabs(pts[0].x - pts[1].x) < 1e-8)
	{
		Line.a = 1.0;
		Line.b = 0.0;
		Line.c = -pts[0].x;

		cLine.a = 0.0;
		cLine.b = 1.0;
		cLine.c = -pts[0].y;
	}
	else if (fabs(pts[0].y - pts[1].y) < 1e-8)
	{
		cLine.a = 1.0;
		cLine.b = 0.0;
		cLine.c = -pts[0].x;

		Line.a = 0.0;
		Line.b = 1.0;
		Line.c = -pts[0].y;
	}
	else
	{
		Line.a = (pts[1].y - pts[0].y) / (pts[1].x - pts[0].x);
		Line.b = -1;
		Line.c = -Line.a*pts[1].x - Line.b*pts[1].y;

		cLine.a = 1.0 / Line.a;
		cLine.b = -1;
		cLine.c = -cLine.a*pts[0].x - cLine.b*pts[0].y;
	}

	cout << "a = " << Line.a << endl;
	cout << "b = " << Line.b << endl;
	cout << "c = " << Line.c << endl;

	cout << pts[0].x << "," << pts[0].y << endl;
	cout << pts[1].x << "," << pts[1].y << endl;
}
