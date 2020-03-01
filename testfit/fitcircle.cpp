#include "fitcircle.h"
#include <vector>
#include <complex>
#include "point.h"
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

using namespace std;
using namespace pcl;

bool FitCircle(const std::vector<f_POINT> &points, double &cent_x, double &cent_y, double &radius)
{
	cent_x = 0.0f;
	cent_y = 0.0f;
	radius = 0.0f;
	if (points.size() < 3)
	{
		return false;
	}

	double sum_x = 0.0f, sum_y = 0.0f;
	double sum_x2 = 0.0f, sum_y2 = 0.0f;
	double sum_x3 = 0.0f, sum_y3 = 0.0f;
	double sum_xy = 0.0f, sum_x1y2 = 0.0f, sum_x2y1 = 0.0f;

	int N = points.size();
	for (int i = 0; i < N; i++)
	{
		double x = points[i].real();
		double y = points[i].imag();
		double x2 = x * x;
		double y2 = y * y;
		sum_x += x;
		sum_y += y;
		sum_x2 += x2;
		sum_y2 += y2;
		sum_x3 += x2 * x;
		sum_y3 += y2 * y;
		sum_xy += x * y;
		sum_x1y2 += x * y2;
		sum_x2y1 += x2 * y;
	}

	double C, D, E, G, H;
	double a, b, c;

	C = N * sum_x2 - sum_x * sum_x;
	D = N * sum_xy - sum_x * sum_y;
	E = N * sum_x3 + N * sum_x1y2 - (sum_x2 + sum_y2) * sum_x;
	G = N * sum_y2 - sum_y * sum_y;
	H = N * sum_x2y1 + N * sum_y3 - (sum_x2 + sum_y2) * sum_y;
	a = (H * D - E * G) / (C * G - D * D);
	b = (H * C - E * D) / (D * D - G * C);
	c = -(a * sum_x + b * sum_y + sum_x2 + sum_y2) / N;

	cent_x = a / (-2);
	cent_y = b / (-2);

	radius = sqrt(a * a + b * b - 4 * c) / 2;
	return true;
}

void ccltCircle(const vector<vector<zPoint> > &vertices, PointCloud<PointXYZI> &curvatures)
{
	int num = 100;
	for (size_t i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].size() < 100)
		{
			continue;
		}
		for (size_t j = 0; j < vertices[i].size(); j++)
		{
			vector<f_POINT> points;
			points.clear();
			points.push_back(f_POINT(vertices[i][j].x, vertices[i][j].y));
			for (size_t m = 1; m <= num; m++)
			{
				if (j + m >= vertices[i].size())
				{
					f_POINT a(vertices[i][j + m - vertices[i].size()].x, vertices[i][j + m - vertices[i].size()].y);
					f_POINT b(vertices[i][j - m].x, vertices[i][j - m].y);
					points.push_back(a);
					points.push_back(b);
				}
				else if (j - m < 0)
				{
					f_POINT a(vertices[i][j + m].x, vertices[i][j + m].y);
					f_POINT b(vertices[i][j - m + vertices[i].size()].x, vertices[i][j - m + vertices[i].size()].y);
					points.push_back(a);
					points.push_back(b);
				}
				else
				{
					f_POINT a(vertices[i][j + m].x, vertices[i][j + m].y);
					f_POINT b(vertices[i][j - m].x, vertices[i][j - m].y);
					points.push_back(a);
					points.push_back(b);
				}
			}
			double a, b, r;
			FitCircle(points, a, b, r);

			PointXYZI p;
			p.x = vertices[i][j].x;
			p.y = vertices[i][j].y;
			p.z = 0;
			p.intensity = 1000000 / r;

			curvatures.points.push_back(p);
		}
	}
	curvatures.width = curvatures.points.size();
	curvatures.height = 1;
	curvatures.is_dense = false;
}