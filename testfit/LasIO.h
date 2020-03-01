#pragma once
#ifndef _LASIO_H_
#define _LASIO_H_
#include <string>
#include <vector>
#include <liblas\liblas.hpp>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

struct MinBound
{
	double minX;
	double minY;
	double minZ;
};
struct MaxBound
{
	double maxX;
	double maxY;
	double maxZ;
};

/*xyzi点结构*/
struct Point3DI
{
	double z_X;
	double z_Y;
	double z_Z;
	short z_Intensity;

	Point3DI(double dx=0.0,double dy=0.0,double dz=0.0, short inten=0)
	{
		z_X = dx;
		z_Y = dy;
		z_Z = dz;
		z_Intensity = inten;
	}
};

struct Point3DRGBI : Point3DI
{
	short z_R;
	short z_G;
	short z_B;

	Point3DRGBI()
		:Point3DI()
	{
		z_R = 0;
		z_G = 0;
		z_B = 0;
	}

	Point3DRGBI(double dx, double dy, double dz, short inten, short r = 0, short g = 0, short b = 0)
		:Point3DI(dx, dy, dz, inten)
	{
		z_R = r;
		z_G = g;
		z_B = b;
	}
};

class LasIO
{
	typedef pcl::PointCloud<pcl::PointXYZI> pclXYZI;
public:
	//构造函数读取并显示header信息
	LasIO(std::string pack = "datapath.txt");
	//当sig不等于0时 向数据文件目录下写las_pcd.pcd文件
	pclXYZI transToPcd(int sig = 0, std::string outname = "las_pcd.pcd");
	//void writePCD2Las(pclXYZI cloud, std::string outname = "pcd2las.las");
	
	long getPointCount() const { return this->m_pointCount; }
	std::vector<Point3DI> getPointsVec() const { return this->m_pointsVec; }
	MinBound getMinBound() const { return this->m_minbound; }
	MaxBound getMaxBound() const { return this->m_maxbound; }
	std::string getDataDir() const { return this->m_dataDir; }
	void writePCD2Las(pclXYZI cloud, std::string outname);
	~LasIO();

private:
	std::string m_dataDir;
	//点数
	long m_pointCount;
	//最小xyz坐标
	MinBound m_minbound;
	MaxBound m_maxbound;
	//点坐标和强度
	std::vector<Point3DI> m_pointsVec;
	pclXYZI pclPoints;
	liblas::Header m_header;
};
#endif

