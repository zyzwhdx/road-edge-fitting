#ifndef DATAIO
#define DATAIO


#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include <string>
#include <pcl/PolygonMesh.h>


typedef  pcl::PointCloud<pcl::PointXYZI>::Ptr      pointCloudXYZIPtr;
typedef  pcl::PointCloud<pcl::PointXYZI>           pointCloudXYZI;


typedef  pcl::PointCloud<pcl::PointXYZ>::Ptr      pointCloudXYZPtr;
typedef  pcl::PointCloud<pcl::PointXYZ>           pointCloudXYZ;


typedef  pcl::PointCloud<pcl::PointXY>::Ptr      pointCloudXYPtr;
typedef  pcl::PointCloud<pcl::PointXY>           pointCloudXY;

class dataIo
{
public:

	struct pointCloudBound
	{
		double minx;
		double maxx;
		double miny;
		double maxy;
		double minz;
		double maxz;
		pointCloudBound()
		{
			minx = maxx = miny = maxy = minz = maxz = 0.0;
		}
	};

	bool readPointCloudFromPcdFileA(const std::string &fileName,pointCloudXYZI &pointCloud);

	bool readPointCloudFromLasFile(const std::string &fileName, pointCloudXYZI &pointCloud);

	bool readPointCloudFromPlyFileA(const std::string &fileName, pointCloudXYZ &pointCloud);
	bool readPointCloudFromPlyFileA(const std::string &fileName, pointCloudXYZI &pointCloud);
	bool readPointCloudFromPlyFileA(const std::string &fileName, pcl::PolygonMesh &mesh);

	
	bool writePointCloudIntoLasFile(const std::string &fileName, pointCloudXYZI &pointCloud);

	bool writePointCloudIntoLasFile(const std::string &fileName, pointCloudXYZI &pointCloud, double ox, double oy);

	void getCloudBound(pointCloudXYZI & cloud, pointCloudBound & bound)
	{
		double min_x = 0.0;
		double min_y = 0.0;
		double min_z = 0.0;
		double max_x = 0.0;
		double max_y = 0.0;
		double max_z = 0.0;

		for (size_t i = 0; i<cloud.size(); ++i)
		{
			//»ñÈ¡±ß½ç
			if (min_x>cloud.points[i].x)
				min_x = cloud.points[i].x;
			if (min_y > cloud.points[i].y)
				min_y = cloud.points[i].y;
			if (min_z > cloud.points[i].z)
				min_z = cloud.points[i].z;
			if (max_x < cloud.points[i].x)
				max_x = cloud.points[i].x;
			if (max_y < cloud.points[i].y)
				max_y = cloud.points[i].y;
			if (max_z < cloud.points[i].z)
				max_z = cloud.points[i].z;
		}
		bound.minx = min_x;
		bound.maxx = max_x;
		bound.miny = min_y;
		bound.maxy = max_y;
		bound.minz = min_z;
		bound.maxz = max_z;
	}

protected:
	
private:
	
};
#endif