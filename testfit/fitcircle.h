#ifndef _FITCIRCLE_H_
#define _FITCIRCLE_H_
#include <vector>
#include <complex>
#include "point.h"
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

using namespace pcl;
using namespace std;

typedef complex<double> f_POINT;

bool FitCircle(const std::vector<f_POINT> &points, double &cent_x, double &cent_y, double &radius);

void ccltCircle(const vector<vector<zPoint> > &vertices, PointCloud<PointXYZI> &curvatures);

#endif