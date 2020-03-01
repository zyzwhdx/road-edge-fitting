#ifndef _POINT_H_
#define _POINT_H_

#include <vector>
#include <string>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
using namespace pcl;
using namespace std;

struct zPoint
{
	double x;
	double y;
};

struct zSegment
{
	zPoint src;
	zPoint dst;
};

struct ePoint
{
	int id;
	double x;
	double y;
	double z;
	ePoint()
	{
		x = y = z = 0.0;
		id = 0;
	}
};

struct eDistance 
{
	double xx;
	double yy;
	double zz;
	double dis;
};

struct normalPara
{
	double ave;
	double stv;

	normalPara()
	{
		ave = 0.0;
		stv = 1.0;
	}
};

struct polyPara
{
	double a;
	double b;
	double c;
	double d;
	double e;
	polyPara()
	{
		a = b = c = d = e = 0.0;
	}
};

struct line2D 
{
	double a;
	double b;
	double c;
	line2D()
	{
		a = 1;
		b = 1;
		c = 0;
	}
};

enum extent_dir
{
	e_x, e_y
};

//一串边缘点定义一个string
typedef vector<ePoint> eString;

//一串string定义一个edge
typedef vector<eString> eEdge;

void readVertices(string path, vector<vector<zPoint> > &vertices, vector<vector<zSegment> > &segments, double &ox, double &oy);

void readCheckPoints(string path, vector<zPoint> &checkpts, double &ox, double &oy);

double ccltDistancePt2Sg(zPoint p, zPoint m, zPoint n);

double ccltPtDistance2D(ePoint a, ePoint b);

void ccltCheckDistances(vector<zPoint> checks, vector<vector<zSegment> > segments, vector<vector<double> > &distances);

double ccltMaxOffset(eString &pts, int &pos_maxoft);

void preRead(string path, eString &pts, eString &pts2, double &ox, double &oy);

void findCenter(eString &left, eString &right, eString &center);

void readEdges(string path, eString &pts, double &ox, double &oy, extent_dir &edir);

void readEdges(eString &center, extent_dir &edir);

void iterativeBreakLines(eString &oPts, eEdge &whole_edge, int th);

bool compy(ePoint &a, ePoint &b);

bool compx(ePoint &a, ePoint &b);

bool compid(ePoint &a, ePoint &b);

bool compdis(eDistance &a, eDistance &b);

bool compSize(eString &a, eString &b);

double medianNum(vector<double> &a);

int zrand();

void ccltNormalPara(eString &a, normalPara &para, const extent_dir &edir);

void fitString(const eString &pts, polyPara &rst, polyPara &zrst, const normalPara &para, const extent_dir &edir);

void fitEdge(const eEdge &edge, const double &ox, const double &oy, const extent_dir &dir, PointCloud<PointXYZI> &Vector);

void solveFit(const vector<double> &xx, const vector<double> &yy, vector<double> &result);

void addPoint(PointCloud<PointXYZI> &cloud, ePoint src, ePoint dst);

//Line是切线，cLine是法线
void ccltLine(ePoint outP, const PointCloud<PointXYZI> &Vector, line2D &Line, line2D &cLine);

void findCrossSection(const PointCloud<PointXYZI> &cloud, const line2D &line, const line2D &cLine, const int &width, const int &thick);

#endif