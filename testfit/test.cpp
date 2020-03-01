/************************************************************************/
/* Created by Joey Zhou on December 19th, 2018                          */
/* VS2015 liblas OpenCV310 PCL1.8.1                                     */
/************************************************************************/

#include <liblas\liblas.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include "LasIO.h"
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
//#include "Grid.h"
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <vector>
#include "dataIo.h"
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>

using namespace pcl;
using namespace std;

/*r邻域强度中值滤波*/
int midFilter(PointCloud<PointXYZI>::ConstPtr cloud, PointCloud<PointXYZI> &outCloud, int flag = 0, string outPathName = "d:/mid_filtered.pcd", float r = 0.5);
/*r邻域强度均值滤波*/
int aveFilter(PointCloud<PointXYZI>::ConstPtr cloud, PointCloud<PointXYZI> &outCloud, int flag = 0, string outPathName = "d:/ave_filtered.pcd", float r = 0.5);
//void pcd2las(PointCloud<PointXYZI> cloud, std::string outname);
/*计算k邻域高程标准差，输出为强度量*/
int stdevHeight(PointCloud<PointXYZI>::ConstPtr cloud, PointCloud<PointXYZI> &outCloud, int flag = 0, string outPathName = "stdev_height.las", int k = 10);

//intensity of put pointcloud need to be a label but true intensity
//return percentage of plane inlier points
float ransacPlnDet(PointCloud<PointXYZ>::ConstPtr inCloud, PointCloud<PointXYZ>::Ptr outCloud,PointIndices::Ptr inliers, float threshold = 0.03);


int mainlt(int argc, char *argv[])
{
	LasIO las;
	pcl::PointCloud<pcl::PointXYZI> PCcloud = las.transToPcd(0);
	//采用深层复制，所以从此PCcloud作废，以下均使用智能指针
	MinBound bound = las.getMinBound();
	pcl::PointCloud<pcl::PointXYZI>::ConstPtr cloud = PCcloud.makeShared();
	//PCcloud实体和cloud指针均不变，作为数据源

	/************************************************************************/
	/*                       中值滤波去噪                                   */
	/************************************************************************/

	////创建原始点云副本
	//PointCloud<PointXYZI> midValFiltered = PCcloud;

	////滤波
	//midFilter(cloud, midValFiltered);

	/************************************************************************/
	/*                       均值滤波去噪                                   */
	/************************************************************************/
	//PointCloud<PointXYZI> aveValFiltered = PCcloud;

	//aveFilter(cloud, aveValFiltered);

	/************************************************************************/
	/*                       高程标准差                                     */
	/************************************************************************/

	PointCloud<PointXYZI> stdZ = PCcloud;

	stdevHeight(cloud, stdZ, 1);

	system("pause");
	return 0;
}

//中值滤波，r邻域
int midFilter(PointCloud<PointXYZI>::ConstPtr cloud, PointCloud<PointXYZI> &outCloud, int flag, string outPathName, float r)
{
	KdTreeFLANN<pcl::PointXYZI> kdtree;

	kdtree.setInputCloud(cloud);
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;

	float radius = r;

	for (size_t p = 0; p < cloud->points.size(); p++)
	{
		//对每个点取r邻域，有大于1个邻近点则中值滤波
		PointXYZI searchPoint = cloud->points[p];
		if (kdtree.radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 1)
		{
			vector<float> iWindow;
			float mid = 0;
			iWindow.clear();
			for (size_t i = 0; i < pointIdxRadiusSearch.size(); ++i)
			{
				iWindow.push_back(cloud->points[pointIdxRadiusSearch[i]].intensity);
			}

			sort(iWindow.begin(), iWindow.end());

			if (iWindow.size() % 2 == 0)//偶数个取平均值
			{
				mid = (iWindow.at(iWindow.size() / 2) + iWindow.at(iWindow.size() / 2 - 1)) / 2;
			}
			else
			{
				mid = iWindow.at(iWindow.size() / 2);
			}
			//cout << abs(searchPoint.intensity - mid) << endl;
			outCloud.points[p].intensity = mid;
		}
		else
		{
			outCloud.points[p].intensity = 0;
		}
	}
	cout << "Successfully commit middle value filter " << endl;
	if (flag != 0)
	{
		pcl::PCDWriter writer;
		writer.write(outPathName, outCloud);
		cout << "Successfully write filtered cloud to " << outPathName << endl;
	}
	return 0;
}


//均值滤波，r邻域
int aveFilter(PointCloud<PointXYZI>::ConstPtr cloud, PointCloud<PointXYZI> &outCloud, int flag, string outPathName, float r)
{
	KdTreeFLANN<pcl::PointXYZI> kdtree;

	kdtree.setInputCloud(cloud);
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;

	float radius = r;

	for (size_t p = 0; p < cloud->points.size(); p++)
	{
		//对每个点取r邻域，有大于1个邻近点则均值滤波
		PointXYZI searchPoint = cloud->points[p];
		if (kdtree.radiusSearch(searchPoint, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 1)
		{
			vector<float> iWindow;
			float mid = 0;
			iWindow.clear();
			for (size_t i = 0; i < pointIdxRadiusSearch.size(); ++i)
			{
				iWindow.push_back(cloud->points[pointIdxRadiusSearch[i]].intensity);
			}

			float sum = accumulate(std::begin(iWindow), std::end(iWindow), 0.0);
			float mean = sum / iWindow.size();

			outCloud.points[p].intensity = mean;
		}
		else
		{
			outCloud.points[p].intensity = 0;
		}
	}
	cout << "Successfully commit mean value filter " << endl;
	if (flag != 0)
	{
		pcl::PCDWriter writer;
		writer.write(outPathName, outCloud);
		cout << "Successfully write filtered cloud to " << outPathName << endl;
	}
	return 0;
}

/********************************************************************************/
/* 推测因为版本问题，liblas和boost，所以无法使用writeLas功能……，用txt代替吧   */
/********************************************************************************/
//void pcd2las(PointCloud<PointXYZI> cloud, std::string outname)
//{
//	std::string filePathName = outname;
//
//	std::vector<Point3DRGBI> Pts;
//	Pts.clear();
//	for (size_t i = 0; i < cloud.points.size(); i++)
//	{
//		Point3DRGBI tPnt(cloud.points[i].x, cloud.points[i].y, cloud.points[i].z, cloud.points[i].intensity, 180);
//		Pts.push_back(tPnt);
//	}
//
//	std::ofstream Ofs(filePathName, std::ios::out | std::ios::binary);
//
//	liblas::Header Header;
//	Header.SetVersionMajor(1);
//	Header.SetVersionMinor(2);
//	Header.SetDataFormatId(liblas::ePointFormat2);
//	Header.SetPointRecordsCount(Pts.size());
//	Header.SetScale(0.001, 0.001, 0.001);
//	Header.SetOffset(0, 0, 0);
//
//	liblas::Writer Writer(Ofs, Header);
//	liblas::Point TempPnt(&Header);
//	for (int i = 0; i < Pts.size(); i++)
//	{
//		liblas::Color Color(Pts[i].z_R, Pts[i].z_G, Pts[i].z_B);
//		TempPnt.SetCoordinates(Pts[i].z_X, Pts[i].z_Y, Pts[i].z_Z);
//		TempPnt.SetIntensity(Pts[i].z_Intensity);
//		TempPnt.SetColor(Color);
//		Writer.WritePoint(TempPnt);
//	}
//	Ofs.flush();
//	Ofs.close();
//	cout << "successfully write " << Pts.size() << " points to " << filePathName << endl;
//}

int stdevHeight(PointCloud<PointXYZI>::ConstPtr cloud, PointCloud<PointXYZI> &outCloud, int flag, string outPathName, int k)
{
	KdTreeFLANN<pcl::PointXYZI> kdtree;

	kdtree.setInputCloud(cloud);
	int kNeib = k;
	std::vector<int> pointIdxKNNSearch(kNeib);
	std::vector<float> pointKNNSquaredDistance(kNeib);

	for (size_t p = 0; p < cloud->points.size(); p++)
	{
		PointXYZI searchPoint = cloud->points[p];
		if (kdtree.nearestKSearch(searchPoint, kNeib, pointIdxKNNSearch, pointKNNSquaredDistance) > 1)
		{
			vector<float> sta;
			sta.clear();
			for (int i = 0; i<pointIdxKNNSearch.size();i++)
			{
				sta.push_back(cloud->points[pointIdxKNNSearch[i]].z);
			}
			float sum = accumulate(begin(sta), end(sta), 0.0);
			float mean = sum / sta.size(); //均值

			double accum = 0.0;
			for_each(begin(sta), end(sta), [&](const double d) {
				accum += (d - mean)*(d - mean);
			});
			double stdev = sqrt(accum / (sta.size() - 1)); //标准差

			outCloud.points[p].intensity = stdev * 1000;
		}
		else
		{
			outCloud.points[p].intensity = 0;
		}
	}
	cout << "Successfully calculated height standard deviation " << endl;

	if (flag != 0)
	{
		dataIo dio;
		dio.writePointCloudIntoLasFile(outPathName, outCloud);
		cout << "Successfully write calculation result to intensity cloud: " << outPathName << endl;
	}
	return 0;
}

float ransacPlnDet(PointCloud<PointXYZ>::ConstPtr inCloud, PointCloud<PointXYZ>::Ptr outCloud, PointIndices::Ptr inliers, float threshold)
{
	//这个函数下一步计划把第三个参数也换成点云输出去，输出indice现在看非常恶心……
	//plane parameters
	ModelCoefficients::Ptr coef(new ModelCoefficients);
	//PointIndices::Ptr inliers(new PointIndices);
	SACSegmentation<PointXYZ> seg;

	inliers->indices.clear();

	seg.setOptimizeCoefficients(true);
	seg.setModelType(SACMODEL_PLANE);
	seg.setMethodType(SAC_RANSAC);
	seg.setDistanceThreshold(threshold);
	seg.setInputCloud(inCloud);

	seg.segment(*inliers, *coef);

	ExtractIndices<PointXYZ> extract;
	extract.setInputCloud(inCloud);
	extract.setIndices(inliers);
	extract.setNegative(true);
	extract.filter(*outCloud);
	
	float percentage = inliers->indices.size() * 1.0 / inCloud->points.size();
	cout << "Successfully extracted a plane: " << percentage * 100 << "%" << endl;
	return percentage;
}