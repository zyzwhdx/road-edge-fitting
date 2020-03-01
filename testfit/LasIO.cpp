#include "LasIO.h"
#include <liblas\liblas.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <boost/filesystem.hpp>
#include <windows.h>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
typedef pcl::PointCloud<pcl::PointXYZI> pclXYZI;

LasIO::LasIO(std::string pack)
{
	try
	{
		ifstream fileMenu;
		fileMenu.open(pack, 0);
		string filePathAndName;
		//暂时只能读一个数据文件
		if (fileMenu.is_open()) {
			while (getline(fileMenu, filePathAndName))
			{
				cout << "Data file is: " << filePathAndName << endl;
				m_dataDir = boost::filesystem::path(filePathAndName).parent_path().string();
			}
		}
		else {
			cout << "Cannot find datapath.txt in current directory." << endl;
		}

		/************************************************************************/

		std::ifstream ifs;
		ifs.open(filePathAndName.c_str(), std::ios::in | std::ios::binary);
		if (!ifs)
		{
			cout << "null" << endl;
		}
		liblas::ReaderFactory f;

		//若文件名有错 此处会有异常
		liblas::Reader reader = f.CreateWithStream(ifs);
		liblas::Header const& header = reader.GetHeader();
		m_header = header;

		m_minbound.minX = header.GetMinX();
		m_minbound.minY = header.GetMinY();
		m_minbound.minZ = header.GetMinZ();
		m_maxbound.maxX = header.GetMaxX();
		m_maxbound.maxY = header.GetMaxY();
		m_maxbound.maxZ = header.GetMaxZ();

		m_pointCount = header.GetPointRecordsCount();

		std::cout << std::setiosflags(std::ios::fixed);
		std::cout << std::setprecision(2);

		//std::cout << "X scale factor              : " << header.GetScaleX() << std::endl;
		//std::cout << "Y scale factor              : " << header.GetScaleY() << std::endl;
		//std::cout << "Z scale factor              : " << header.GetScaleZ() << std::endl;
		//std::cout << "X offset                    : " << header.GetOffsetX() << std::endl;
		//std::cout << "Y offset                    : " << header.GetOffsetY() << std::endl;
		//std::cout << "Z offset                    : " << header.GetOffsetZ() << std::endl;
		//std::cout << "Max X                       : " << header.GetMaxX() << std::endl;
		//std::cout << "Max Y                       : " << header.GetMaxY() << std::endl;
		//std::cout << "Max Z                       : " << header.GetMaxZ() << std::endl;
		//std::cout << "Min X                       : " << header.GetMinX() << std::endl;
		//std::cout << "Min Y                       : " << header.GetMinY() << std::endl;
		//std::cout << "Min Z                       : " << header.GetMinZ() << std::endl;

		cout << "Points count: " << m_pointCount << endl;
		cout << "MinX: " << m_minbound.minX << endl;
		cout << "MinY: " << m_minbound.minY << endl;
		cout << "MinZ: " << m_minbound.minZ << endl;
		/************************************************************************/
		m_pointsVec.clear();
		reader.Seek(0);
		//在此处根据需要加百分比进度条 太慢了
		//long i = 0;
		while (reader.ReadNextPoint())
		{
			const liblas::Point &curPt = reader.GetPoint();
			Point3DI tmpPt(curPt.GetX(), curPt.GetY(), curPt.GetZ(), curPt.GetIntensity());
			m_pointsVec.push_back(tmpPt);
			//i++;
			//printf("%.2lf%%\r", i * 100.0 / m_pointCount);
		}
		cout << "Successfully read " << m_pointsVec.size() << " points into vector.\n";
		ifs.close();
	}
	catch (std::out_of_range &exc)
	{
		std::cout << exc.what() << std::endl;
	}
}

/* 当sig！=0时 向数据目录下写一个pcd文件                                                                     */
pclXYZI LasIO::transToPcd(int sig, std::string outname)
{
	pclPoints.width = m_pointCount;
	pclPoints.height = 1;
	pclPoints.is_dense = false;
	pclPoints.resize(pclPoints.width * pclPoints.height);
	for (long i=0;i<pclPoints.width;i++)
	{
		pclPoints.points[i].x = m_pointsVec[i].z_X - m_minbound.minX;
		pclPoints.points[i].y = m_pointsVec[i].z_Y - m_minbound.minY;
		pclPoints.points[i].z = m_pointsVec[i].z_Z;
		pclPoints.points[i].intensity = m_pointsVec[i].z_Intensity;
	}
	//因为太慢 不用写pcd的方法，如果有需要则用CC
	if (sig != 0)
	{
		pcl::io::savePCDFileASCII(m_dataDir + outname, pclPoints);
		std::cerr << "Saved " << pclPoints.points.size() << " data points to las_pcd.pcd." << std::endl;
	}
	return pclPoints;
}

//void LasIO::writePCD2Las(pclXYZI cloud, std::string outname)
//{
//	std::string filePathName = m_dataDir + outname;
//
//	std::vector<Point3DRGBI> Pts;
//	Pts.clear();
//	for (size_t i = 0; i < cloud.points.size(); i++)
//	{
//		Point3DRGBI tPnt(cloud.points[i].x + m_minbound.minX, cloud.points[i].y + m_minbound.minY, cloud.points[i].z, cloud.points[i].intensity, 180);
//		Pts.push_back(tPnt);
//	}
//
//	std::ofstream Ofs;
//	Ofs.open(filePathName.c_str(), std::ios::out | std::ios::binary);
//	liblas::WriterFactory f;
//
//	//liblas::Header Header;
//	//Header.SetVersionMajor(1);
//	//Header.SetVersionMinor(2);
//	//Header.SetDataFormatId(liblas::ePointFormat2);
//	//Header.SetPointRecordsCount(Pts.size());
//	//Header.SetScale(0.001, 0.001, 0.001);
//	//Header.SetOffset(0, 0, 0);
//	liblas::Writer Writer = f.CreateWithStream(Ofs, m_header);
//
//	for (int i = 0; i < Pts.size(); i++)
//	{
//		liblas::Point TempPnt(&m_header);
//		liblas::Color Color(Pts[i].z_R, Pts[i].z_G, Pts[i].z_B);
//		TempPnt.SetCoordinates(Pts[i].z_X, Pts[i].z_Y, Pts[i].z_Z);
//		TempPnt.SetIntensity(Pts[i].z_Intensity);
//		TempPnt.SetColor(Color);
//		Writer.WritePoint(TempPnt);
//	}
//	Ofs.close();
//	cout << "successfully write " << Pts.size() << " points to " << filePathName << endl;
//}


LasIO::~LasIO()
{
}
