#include "dataIo.h"

#include <iostream>


#include <pcl/io/pcd_io.h>
#include <fstream>


#include <liblas/liblas.hpp>
#include <liblas/version.hpp>
#include <liblas/point.hpp>
#include <pcl/io/ply_io.h>


using namespace  std;


bool dataIo::readPointCloudFromPcdFileA(const std::string &fileName, pointCloudXYZI &pointCloud)
{
	if (pcl::io::loadPCDFile<pcl::PointXYZI>(fileName, pointCloud) == -1) //* load the file
	{
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return false;
	}

	return true;
}

bool dataIo::readPointCloudFromLasFile(const std::string &fileName, pointCloudXYZI &pointCloud)
{
	if (!fileName.substr(fileName.rfind('.')).compare(".las"))
	{
		std::ifstream ifs;
		ifs.open(fileName, std::ios::in | std::ios::binary);
		if (ifs.bad())
		{
			cout << "未发现匹配项" << endl;
			return false;
		}
		liblas::ReaderFactory f;
		liblas::Reader reader = f.CreateWithStream(ifs);

		const liblas::Header& header = reader.GetHeader();
		liblas::PointFormatName format = header.GetDataFormatId();//获得点的格式类型,数据类型0和1没有颜色值，其他数据类型有颜色值;
		bool has_color = false;
		if (liblas::ePointFormat0 == format || liblas::ePointFormat1 == format)
			has_color = false;
		else
			has_color = true;

		/*while循环中遍历所有的点;*/
		while (reader.ReadNextPoint())
		{
			const liblas::Point& p = reader.GetPoint();
			pcl::PointXYZI  las_point;
			/*将重心化后的坐标和强度值赋值给PCL中的点;*/
			las_point.x = p.GetX();
			las_point.y = p.GetY();
			las_point.z = p.GetZ();
			las_point.intensity = p.GetIntensity();
			pointCloud.push_back(las_point);

		}
	}

	return true;
}

bool dataIo::writePointCloudIntoLasFile(const std::string &fileName, pointCloudXYZI &pointCloud)
{
	pointCloudBound bound;
	getCloudBound(pointCloud, bound);

	std::ofstream ofs;
	ofs.open(fileName, std::ios::out | std::ios::binary);
	if (ofs.is_open())
	{
		liblas::Header header;
		header.SetDataFormatId(liblas::ePointFormat2);
		header.SetVersionMajor(1);
		header.SetVersionMinor(2);
		header.SetMin(bound.minx , bound.miny, bound.minz);
		header.SetMax(bound.maxx, bound.maxy, bound.maxz);
		header.SetOffset((bound.minx + bound.maxx) / 2.0, (bound.miny + bound.maxy) / 2.0, (bound.minz + bound.maxz) / 2.0);
		header.SetScale(0.000001, 0.000001, 0.0001);
		header.SetPointRecordsCount(pointCloud.points.size());

		liblas::Writer writer(ofs, header);
		liblas::Point pt(&header);

		for (int i = 0; i < pointCloud.points.size(); i++)
		{
			pt.SetCoordinates(double(pointCloud.points[i].x),double(pointCloud.points[i].y),double(pointCloud.points[i].z));
			pt.SetIntensity(pointCloud.points[i].intensity);
			writer.WritePoint(pt);
		}
		ofs.flush();
		ofs.close();

		return true;
	}

	return false;
}

bool dataIo::writePointCloudIntoLasFile(const std::string &fileName, pointCloudXYZI &pointCloud, double ox, double oy)
{
	pointCloudBound bound;
	getCloudBound(pointCloud, bound);

	std::ofstream ofs;
	ofs.open(fileName, std::ios::out | std::ios::binary);
	if (ofs.is_open())
	{
		liblas::Header header;
		header.SetDataFormatId(liblas::ePointFormat2);
		header.SetVersionMajor(1);
		header.SetVersionMinor(2);
		header.SetMin(bound.minx + ox, bound.miny + oy, bound.minz);
		header.SetMax(bound.maxx + ox, bound.maxy + oy, bound.maxz);
		header.SetOffset(ox + (bound.minx + bound.maxx) / 2.0, oy + (bound.miny + bound.maxy) / 2.0, (bound.minz + bound.maxz) / 2.0);
		header.SetScale(0.001, 0.001, 0.0001);
		header.SetPointRecordsCount(pointCloud.points.size());

		liblas::Writer writer(ofs, header);
		liblas::Point pt(&header);

		for (int i = 0; i < pointCloud.points.size(); i++)
		{
			pt.SetCoordinates(double(pointCloud.points[i].x) + ox, double(pointCloud.points[i].y) + oy, double(pointCloud.points[i].z));
			pt.SetIntensity(pointCloud.points[i].intensity);
			writer.WritePoint(pt);
		}
		ofs.flush();
		ofs.close();

		return true;
	}

	return false;
}


bool dataIo::readPointCloudFromPlyFileA(const std::string &fileName, pointCloudXYZI &pointCloud)
{
	if (pcl::io::loadPLYFile(fileName, pointCloud) == -1) //* load the file
	{
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return false;
	}
	return true;
}

bool dataIo::readPointCloudFromPlyFileA(const std::string &fileName, pointCloudXYZ &pointCloud)
{
	if (pcl::io::loadPLYFile(fileName, pointCloud) == -1) //* load the file
	{
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return false;
	}
	return true;
}

bool dataIo::readPointCloudFromPlyFileA(const std::string &fileName, pcl::PolygonMesh &mesh)
{
	if (pcl::io::loadPLYFile(fileName, mesh) == -1) //* load the file
	{
		PCL_ERROR("Couldn't read file test_pcd.pcd \n");
		return false;
	}
	return true;
}
