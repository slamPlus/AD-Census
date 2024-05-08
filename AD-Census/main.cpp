/* -*-c++-*- AD-Census - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/AD-Census
* Describe	: main
*/
#include <iostream>
#include <string>
#include <fstream>
#include "ADCensusStereo.h"
#include "adcensus_util.h"
#include <chrono>
using namespace std::chrono;

// opencv library
#include <opencv2/opencv.hpp>
#ifdef _DEBUG
#pragma comment(lib,"opencv_world451d.lib")
#else
#pragma comment(lib,"opencv_world451.lib")
#endif

/*保存视差图*/
void SaveDisparityMap(const float32* disp_map, const sint32& width, const sint32& height, const std::string& path);
/*保存视差点云*/
void SaveDisparityCloud(const uint8* img_bytes, const float32* disp_map,
	const sint32& width, const sint32& height, const std::string& path, bool cpt_err = true);
/*匹配精度评估*/
void AssessDisparityMap(const std::string txt_path, const float32* disp_map, const sint32& width, const sint32& height, bool show_map = true);

/**
* \brief
* \param argv 3
* \param argc argc[1]:左影像路径 argc[2]: 右影像路径 argc[3]: 最小视差[可选，默认0] argc[4]: 最大视差[可选，默认64]
* \param eg. ..\Data\cone\im2.png ..\Data\cone\im6.png 0 64
* \param eg. ..\Data\Cloth3\view1.png ..\Data\Cloth3\view5.png 0 128
* \return
*/
int main(int argv, char** argc)
{
	if (argv < 3) {
		std::cout << "参数过少，请至少指定左右影像路径！" << std::endl;
		return -1;
	}

	printf("Image Loading...");
	//···············································································//
	// 读取影像
	/*std::string path_left = "E:/3d_reconstruction/rectifyImageL.jpg";
	std::string path_right = "E:/3d_reconstruction/rectifyImageR.jpg";
	std::string path_dispTxt = "E:/3d_reconstruction/f5/disparityMap_6.txt";*/
	//std::string path_left = "E:/3d_reconstruction/20221020/20221020153407_569_L/rectifyImageL_.png";
	//std::string path_right = "E:/3d_reconstruction/20221020/20221020153407_569_L/rectifyImageR_.png";

	std::ifstream infile;
	std::string s;
	vector<std::string> paths;
	infile.open("文件名.txt");


	while (getline(infile, s))
	{
		std::istringstream sin(s); //将整行字符串line读入到字符串流istringstream中
		//vector<std::string> fields; //声明一个字符串向量
		std::string field;
		while (getline(sin, field, '\n')) //将字符串流sin中的字符读入到field字符串中，以逗号为分隔符
		{
			paths.push_back(field); //将刚刚读取的字符串添加到向量fields中
		}
	}

	for (int i = 0; i < 18; i++) {
		std::cout << paths[i] << std::endl << std::endl;
		//std::string path_left = "E:/3d_reconstruction/20221020/"+ paths[i] + "/rectifyImageL__ROI.png";
		//std::string path_right = "E:/3d_reconstruction/20221020/" + paths[i] + "/rectifyImageR__ROI.png";
		std::string path_left = "E:/3d_reconstruction/rectifyImageL.jpg";
		std::string path_right = "E:/3d_reconstruction/rectifyImageR.jpg";
		//std::string path_left = "E:/3d_reconstruction/data/230806/watermelon1/l/images_vis/11.png";
		//std::string path_right = "E:/3d_reconstruction/data/230806/watermelon1/r/images_vis/11_.png";

		cv::Mat img_left = cv::imread(path_left, cv::IMREAD_COLOR);
		cv::Mat img_right = cv::imread(path_right, cv::IMREAD_COLOR);

		if (img_left.data == nullptr || img_right.data == nullptr) {
			std::cout << "读取影像失败！" << std::endl;
			return -1;
		}
		if (img_left.rows != img_right.rows || img_left.cols != img_right.cols) {
			std::cout << "左右影像尺寸不一致！" << std::endl;
			return -1;
		}


		//···············································································//
		
		// 图像预处理
		//cv::imshow("预处理前", img_left);

		cv::Mat tl, tr;
		adcensus_util::equalizeHist_color(img_left, tl, false);
		adcensus_util::equalizeHist_color(img_right, tr, false);

		cv::bilateralFilter(tl, img_left, 9, 50, 25);
		cv::bilateralFilter(tr, img_right, 9, 50, 25);

		cv::imshow("预处理后", img_left);
		cv::imshow("预处理后r", img_right);
		cv::waitKey(0);
		cv::destroyAllWindows();

		const sint32 width = static_cast<uint32>(img_left.cols);		// 直方图均衡化
		const sint32 height = static_cast<uint32>(img_right.rows);

		// 左右影像的彩色数据
		auto bytes_left = new uint8[width * height * 3];
		auto bytes_right = new uint8[width * height * 3];
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				bytes_left[i * 3 * width + 3 * j] = img_left.at<cv::Vec3b>(i, j)[0];
				bytes_left[i * 3 * width + 3 * j + 1] = img_left.at<cv::Vec3b>(i, j)[1];
				bytes_left[i * 3 * width + 3 * j + 2] = img_left.at<cv::Vec3b>(i, j)[2];
				bytes_right[i * 3 * width + 3 * j] = img_right.at<cv::Vec3b>(i, j)[0];
				bytes_right[i * 3 * width + 3 * j + 1] = img_right.at<cv::Vec3b>(i, j)[1];
				bytes_right[i * 3 * width + 3 * j + 2] = img_right.at<cv::Vec3b>(i, j)[2];
			}
		}
		printf("Done!\n");

		// AD-Census匹配参数设计
		ADCensusOption ad_option;

		// 计算 视差范围 min-max disparity
		/*std::string str_dist = paths[i].substr(paths[i].find("_") + 1, 3);
		std::stringstream geek(str_dist);
		float32 f_dist;
		geek >> f_dist;
		ad_option.min_disparity = std::max((700. - float(f_dist)), 0.);
		ad_option.max_disparity = ad_option.min_disparity + 256;*/
		printf("w = %d, h = %d, d = [%d, %d]\n\n", width, height, ad_option.min_disparity, ad_option.max_disparity);

		// 定义AD-Census匹配类实例
		ADCensusStereo ad_census;

		printf("AD-Census Initializing...\n");
		auto start = steady_clock::now();
		//···············································································//
		// 初始化
		if (!ad_census.Initialize(width, height, ad_option)) {
			std::cout << "AD-Census初始化失败！" << std::endl;
			return -2;
		}
		auto end = steady_clock::now();
		auto tt = duration_cast<milliseconds>(end - start);
		printf("AD-Census Initializing Done! Timing :	%lf s\n\n", tt.count() / 1000.0);

		printf("AD-Census Matching...\n");
		// disparity数组保存子像素的视差结果
		auto disparity = new float32[uint32(width * height)]();

		start = steady_clock::now();
		//···············································································//
		// 匹配
		if (!ad_census.Match(bytes_left, bytes_right, disparity)) {
			std::cout << "AD-Census匹配失败！" << std::endl;
			return -2;
		}
		end = steady_clock::now();
		tt = duration_cast<milliseconds>(end - start);
		printf("\nAD-Census Matching...Done! Timing :	%lf s\n", tt.count() / 1000.0);

		//···············································································//

		// 保存视差图
		SaveDisparityMap(disparity, width, height, path_left);
		// 误匹配像素评价
		std::string path_dispTxt = "E:/3d_reconstruction/f5/disparityMap_6.txt";
		//AssessDisparityMap(path_dispTxt, disparity, width, height, true);
		// 保存点云
		SaveDisparityCloud(bytes_left, disparity, width, height, path_left, false);

		//cv::destroyAllWindows();

		//···············································································//
		// 释放内存
		delete[] disparity;
		disparity = nullptr;
		delete[] bytes_left;
		bytes_left = nullptr;
		delete[] bytes_right;
		bytes_right = nullptr;
	}

	system("pause");
	return 0;
}


void SaveDisparityMap(const float32* disp_map, const sint32& width, const sint32& height, const std::string& path)
{
	// 保存视差图
	const cv::Mat disp_mat = cv::Mat(height, width, CV_8UC1);
	float32 min_disp = float32(width), max_disp = -float32(width);
	for (sint32 i = 0; i < height; i++) {
		for (sint32 j = 0; j < width; j++) {
			const float32 disp = abs(disp_map[i * width + j]);
			if (disp != Invalid_Float) {
				min_disp = std::min(min_disp, disp);
				max_disp = std::max(max_disp, disp);
			}
		}
	}
	for (sint32 i = 0; i < height; i++) {
		for (sint32 j = 0; j < width; j++) {
			const float32 disp = abs(disp_map[i * width + j]);
			if (disp == Invalid_Float) {
				disp_mat.data[i * width + j] = 0;
			}
			else {
				disp_mat.data[i * width + j] = static_cast<uchar>((disp - min_disp) / (max_disp - min_disp) * 255);
			}
		}
	}

	cv::imwrite(path + "-d.png", disp_mat);
	cv::Mat disp_color;
	applyColorMap(disp_mat, disp_color, cv::COLORMAP_MAGMA);
	cv::imwrite(path + "-c.png", disp_color);
}


void SaveDisparityCloud(const uint8* img_bytes, const float32* disp_map,
	const sint32& width, const sint32& height, const std::string& path, bool cpt_err)
{
	// 保存视差点云(x,y,disp,r,g,b)
	FILE* fp_disp_cloud = nullptr;
	fopen_s(&fp_disp_cloud, (path + "-cloud.txt").c_str(), "w");
	if (fp_disp_cloud) {
		for (sint32 i = 0; i < height; i++) {
			for (sint32 j = 0; j < width; j++) {
				const float32 disp = abs(disp_map[i * width + j]);
				if (disp == Invalid_Float) {
					fprintf_s(fp_disp_cloud, "%f %f %f %d %d %d\n", float32(j), float32(i),
						0, img_bytes[i * width * 3 + 3 * j + 2], img_bytes[i * width * 3 + 3 * j + 1], img_bytes[i * width * 3 + 3 * j]);
					continue;
				}
				fprintf_s(fp_disp_cloud, "%f %f %f %d %d %d\n", float32(j), float32(i),
					disp, img_bytes[i * width * 3 + 3 * j + 2], img_bytes[i * width * 3 + 3 * j + 1], img_bytes[i * width * 3 + 3 * j]);
			}
		}
		fclose(fp_disp_cloud);
	}

	float32 min_disp = float32(width), max_disp = -float32(width);
	float32 max_gt_disp = -float32(width);
	for (sint32 i = 0; i < height; i++) {
		for (sint32 j = 0; j < width; j++) {
			const float32 disp = abs(disp_map[i * width + j]);
			if (disp != Invalid_Float) {
				min_disp = std::min(min_disp, disp);
				max_disp = std::max(max_disp, disp);
			}
		}
	}

	float32 min_invd = min_disp + 0.2 * (max_disp - min_disp);
	float32 max_invd = min_disp + 0.9 * (max_disp - min_disp);

	// 计算点云误差
	if (cpt_err) {
		std::cout << "检查下相机参数 和 GroundTruth路径" << std::endl;
		const std::string txt_path = "E:\\3d_reconstruction\\f5\\heartDepthMap_6.txt";
		float32 u0 = 165.964371;
		float32 v0 = 154.498138;
		float32 fx = 391.656525;
		float32 fy = 426.835144;
		float32 baseline = 5;
		float32 doffs = 190.896454 - 165.964371;

		uint32 i = 0, j = 0;
		float32 err_sum = 0.0;
		uint32 err_cnt = 0;
		float32 delta = 25;

		std::string strline;
		std::ifstream fin(txt_path);
		while (getline(fin, strline))
		{
			std::istringstream sin(strline); //将整行字符串line读入到字符串流istringstream中
			std::string xyd_str;
			vector<std::string> xyd;
			while (getline(sin, xyd_str, ' ')) //将字符串流sin中的字符读入到field字符串中
			{
				xyd.push_back(xyd_str);
			}

			// ground_truth : x, y, z
			float32 x_;
			float32 y_;
			float32 z_;
			std::stringstream s0(xyd[0]);	s0 >> x_;
			std::stringstream s1(xyd[1]);	s1 >> y_;
			std::stringstream s2(xyd[2]);	s2 >> z_;

			float32 disp = abs(disp_map[i * width + j]); 
			if (disp == Invalid_Float || disp >= max_invd || disp <= min_invd) {
				// 更新 i, j
				j++;
				if (j == width) {
					i++;
					j = 0;
				}
				continue;
			}

			float32 z = fx * baseline / (disp + doffs); // Zc = baseline * f / (d + doffs)
			float32 x = (j - u0) * z / fx; // Xc向右，Yc向下为正
			float32 y = (i - v0) * z / fy;

			// 计算距离
			float32 dist_ij = sqrt((x - x_) * (x - x_) + (y - y_) * (y - y_) + (z - z_) * (z - z_));
			if (dist_ij > delta) {
				err_cnt++;
			}
			err_sum = err_sum + dist_ij;

			// 更新 i, j
			j++;
			if (j == width) {
				i++;
				j = 0;
			}

		}

		float32 err_mean = err_sum / (width * height);
		float32 err_rate = 1.0 * err_cnt / (width * height);
		std::cout << "平均误差：" << err_mean << "mm" << std::endl;
		std::cout << "错误率：" << err_rate << std::endl;

	}

}

void AssessDisparityMap(const std::string txt_path, const float32* disp_map, 
	const sint32& width, const sint32& height, bool show_map)
{
	// 读取标准视差图
	cv::Mat disp_GT = cv::Mat::zeros(height, width, CV_32FC1);
	cv::Mat disp_GT_norm = cv::Mat::zeros(height, width, CV_32FC1);
	std::string strline;
	int i = 0, j = 0;
	std::ifstream fin(txt_path);
	while (getline(fin, strline))
	{
		if (strline == "") {
			break;
		}
		std::stringstream geek(strline);
		float32 tt;
		geek >> tt;
		disp_GT.at<float>(i, j) = float(tt);
		
		j++;
		if (j == width) {
			i++;
			j = 0;
		}
	}
	
	cv::Mat disp_mat_8u = cv::Mat(height, width, CV_8UC1);		// uchar类型，归一化后的图像
	cv::Mat disp_mat_32f = cv::Mat(height, width, CV_32FC1);	// 单精度浮点，归一化后的图像便于显示
	float32 min_disp = float32(width), max_disp = -float32(width);
	float32 max_gt_disp = -float32(width);
	for (sint32 i = 0; i < height; i++) {
		for (sint32 j = 0; j < width; j++) {
			const float32 disp = abs(disp_map[i * width + j]);
			const float32 disp_gt = abs(disp_GT.at<float>(i, j));
			if (disp != Invalid_Float) {
				min_disp = std::min(min_disp, disp);
				max_disp = std::max(max_disp, disp);
				max_gt_disp = std::max(max_gt_disp, disp_gt);
			}
		}
	}
	cv::normalize(disp_GT, disp_GT_norm, 0, 1, cv::NormTypes::NORM_MINMAX);

	float32 min_invd = min_disp + 0.2 * (max_disp - min_disp);
	float32 max_invd = min_disp + 0.9 * (max_disp - min_disp);
	for (sint32 i = 0; i < height; i++) {
		for (sint32 j = 0; j < width; j++) {
			const float32 disp = abs(disp_map[i * width + j]);
			if (disp == Invalid_Float || disp >= max_invd || disp <= min_invd) {
				disp_mat_8u.data[i * width + j] = 0;
				disp_mat_32f.data[i * width + j] = 0.;
			}
			else {
				disp_mat_8u.data[i * width + j] = static_cast<uchar>(255 * (disp - min_disp) / (max_disp - min_disp));
				disp_mat_32f.at<float>(i, j) = float(1.0 *  (disp - min_invd) / (max_invd - min_invd));
			}
		}
	}

	// 统计误匹配点
	float32 err_sum = 0.0;
	float32 delta = 0.1;	// 视差值偏差阈值
	int cnt = 0;
	int Num = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const float32 disp = float32(disp_mat_32f.at<float>(i, j));
			const float32 d0 = float32(disp_GT_norm.at<float>(i, j));
			if (d0 == 0 || disp >= width || disp <= -width) {
				continue;
			}
			Num++;
			err_sum = err_sum + (disp - d0) * (disp - d0);
			if (abs(disp - d0) > delta) {
				cnt++;
			}
		}
	}
	// 输出评估结果
	float32 err_mean = sqrt(1.0 * err_sum / Num);
	float32 Bn = 1.0 * cnt / Num;
	std::cout << std::endl << std::endl;
	std::cout << "平均匹配误差RMS Error：" << err_mean << ", "<< max_gt_disp * err_mean << std::endl;
	std::cout << "平均误匹配率：" << Bn << std::endl;

	// 显示视差图
	if (show_map) {
		cv::Mat disp_color;
		applyColorMap(disp_mat_8u, disp_color, cv::COLORMAP_JET);

		cv::imshow("disp_mat", disp_mat_32f);
		cv::imshow("disp_mat-color", disp_color);
		cv::imshow("ground truth", disp_GT_norm);
		cv::imwrite("E:/3d_reconstruction/disp_GT.png", disp_GT);
		cv::waitKey(0);
	}
}
