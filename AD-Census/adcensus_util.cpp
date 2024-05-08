/* -*-c++-*- AD-Census - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/AD-Census
* Describe	: implement of adcensus_util
*/

#include "adcensus_util.h"
#include <cassert>

void adcensus_util::census_transform_color2(const uint8* source,
	vector<std::string>& census, const sint32& width, const sint32& height,
	const sint32 window_size_w, const sint32& window_size_h, bool sparse_window)
{
	if (source == nullptr || census.empty() || width <= window_size_w || height <= window_size_h) {
		return;
	}

	sint32 ws_w2 = sint32((window_size_w - 1) / 2);
	sint32 ws_h2 = sint32((window_size_h - 1) / 2);
	
	sint32 Np = ((window_size_w + 1) / 2) * ((window_size_h + 1) / 2);	// 默认使用稀疏窗口
	sint32 step = 2;
	if (!sparse_window) {
		step = 1;
		Np = window_size_h * window_size_w;
	}
	
	// 逐像素计算census值
	for (sint32 i = ws_w2; i < height - ws_w2; i++) {
		for (sint32 j = ws_h2; j < width - ws_h2; j++) {

			// 中心像素值
			const uint8 b_center = source[i * width * 3 + j * 3];
			const uint8 g_center = source[i * width * 3 + j * 3 + 1];
			const uint8 r_center = source[i * width * 3 + j * 3 + 2];

			// 遍历窗口内邻域像素
			double* dis_mc = new double[Np];
			int Np_i = 0;
			for (sint32 r = -ws_w2; r <= ws_w2; r = r + step) {
				for (sint32 c = -ws_h2; c <= ws_h2; c = c + step) {
					// 计算p、q两点间颜色的曼哈顿距离
					const uint8 bb = source[(i + r) * width * 3 + (j + c) * 3];
					const uint8 gg = source[(i + r) * width * 3 + (j + c) * 3 + 1];
					const uint8 rr = source[(i + r) * width * 3 + (j + c) * 3 + 2];
					dis_mc[Np_i] = double(abs(bb - b_center) + abs(gg - g_center) + abs(rr - r_center));
					Np_i++;
				}
			}
			// 计算 d_mean, d_std
			double dis_mean = 0., dis_std = 0., dis_sum = 0.;
			for (int ind = 0; ind < Np; ind++) {
				dis_sum = dis_sum + dis_mc[ind];
			}
			dis_mean = 1.0 * dis_sum / Np;
			dis_sum = 0.;
			for (int ind = 0; ind < Np; ind++) {
				dis_sum = dis_sum + (dis_mc[ind] - dis_mean) * (dis_mc[ind] - dis_mean);
			}
			dis_std = sqrt(dis_sum / Np);

			// 逐一比较像素值与中心像素值的的大小，计算p点的census值
			std::string cen_val = "";
			for (int ind = 0; ind < Np; ind++) {
				double dis_m = abs(dis_mc[ind] - dis_mean);
				if (dis_m < dis_std)	cen_val = cen_val + "1";
				else cen_val = cen_val + "0";
			}

			// 中心像素的census值
			census[i * width + j] = cen_val;
		}
	}
}

void adcensus_util::census_transform_color(const uint8* source, 
	vector<uint64>& census, const sint32& width, const sint32& height,
	const sint32 window_size_w, const sint32& window_size_h)
{
	sint32 ws_w2 = sint32((window_size_w - 1) / 2);
	sint32 ws_h2 = sint32((window_size_h - 1) / 2);
	sint32 Np = window_size_w * window_size_h;
	if (Np > 64) {
		printf("减小窗口大小为 5 * 5 \n");
		ws_w2 = 5;
		ws_h2 = 5;
		Np = 25;
	}
	if (source == nullptr || census.empty() || width <= window_size_w || height <= window_size_h) {
		return;
	}

	// 逐像素计算census值
	for (sint32 i = ws_w2; i < height - ws_w2; i++) {
		for (sint32 j = ws_h2; j < width - ws_h2; j++) {

			// 中心像素值
			const uint8 b_center = source[i * width * 3 + j * 3];
			const uint8 g_center = source[i * width * 3 + j * 3 + 1];
			const uint8 r_center = source[i * width * 3 + j * 3 + 2];

			// 遍历窗口内邻域像素
			uint64* dis_mc = new uint64[Np];
			int Np_i = 0;
			for (sint32 r = -ws_w2; r <= ws_w2; r++) {
				for (sint32 c = -ws_h2; c <= ws_h2; c++) {
					// 计算p、q两点间颜色的曼哈顿距离
					const uint8 bb = source[(i + r) * width * 3 + (j + c) * 3];
					const uint8 gg = source[(i + r) * width * 3 + (j + c) * 3 + 1];
					const uint8 rr = source[(i + r) * width * 3 + (j + c) * 3 + 2];
					dis_mc[Np_i] = uint64(abs(bb - b_center) + abs(gg - g_center) + abs(rr - r_center));
					Np_i++;
				}
			}
			// 计算 d_mean, d_std
			double dis_mean = 0., dis_std = 0., dis_sum = 0.;
			for (int ind = 0; ind < Np; ind++) {
				dis_sum = dis_sum + dis_mc[ind];
			}
			dis_mean = 1.0 * dis_sum / Np;
			dis_sum = 0.;
			for (int ind = 0; ind < Np; ind++) {
				dis_sum = dis_sum + (dis_mc[ind] - dis_mean) * (dis_mc[ind] - dis_mean);
			}
			dis_std = sqrt(dis_sum / Np);

			// 逐一比较像素值与中心像素值的的大小，计算p点的census值
			uint64 cen_val = 0u;
			for (int ind = 0; ind < Np; ind++) {
				double dis_m = dis_mc[ind] - dis_mean;
				if (dis_m < dis_std)	cen_val = (cen_val << 1) + 1;
				else cen_val = cen_val << 1;
			}

			// 中心像素的census值
			census[i * width + j] = cen_val;
		}
	}
}

uint8 adcensus_util::HammingDist(const std::string& x, const std::string& y)
{
	if ((x.size() != y.size()) || (x.size() <= 0) || (y.size() <= 0))
	{
		return 1.0f;
	}
	sint32 len = x.size();
	uint64 dist = 0;

	for (sint32 i = 0; i < len; i++) {
		dist += (x[i] != y[i]) ? 1 : 0;
	}
	
	//sint32 n = len / 64;	// 切割成n个64位整数来做
	//for (sint32 i = 0; i < n; i++) {
	//	std::string s1 = x.substr(i * 64, (i + 1) * 64 - 1);
	//	std::string s2 = y.substr(i * 64, (i + 1) * 64 - 1);
	//	uint64 temp1 = strtol(s1.data(), NULL, 2);
	//	uint64 temp2 = strtol(s2.data(), NULL, 2);
	//	dist = dist + Hamming64(temp1, temp2);
	//}
	//std::string s1 = x.substr(n * 64, x.length());
	//std::string s2 = y.substr(n * 64, y.length());
	//uint64 temp1 = strtol(s1.data(), NULL, 2);
	//uint64 temp2 = strtol(s2.data(), NULL, 2);
	//dist = dist + Hamming64(temp1, temp2);

	return static_cast<uint8>(dist);
}

void adcensus_util::census_transform2(const uint8* source, vector<std::string>& census, 
	const sint32& width, const sint32& height)
{
	if (source == nullptr || census.empty() || width <= 9 || height <= 7) {
		return;
	}

	// 逐像素计算census值
	for (sint32 i = 4; i < height - 4; i++) {
		for (sint32 j = 3; j < width - 3; j++) {

			// 中心像素值
			const uint8 gray_center = source[i * width + j];

			// 遍历大小为9x7的窗口内邻域像素，逐一比较像素值与中心像素值的的大小，计算census值
			std::string census_val = "";
			for (sint32 r = -4; r <= 4; r++) {
				for (sint32 c = -3; c <= 3; c++) {
					const uint8 gray = source[(i + r) * width + j + c];
					
					if (gray < gray_center) {
						census_val += "1";
					}
					else {
						census_val += "0";
					}
				}
			}

			// 中心像素的census值
			census[i * width + j] = census_val;
		}
	}
}

void adcensus_util::census_transform(const uint8* source, vector<uint64>& census,
	const sint32& width, const sint32& height)
{
	if (source == nullptr || census.empty() || width <= 9 || height <= 7) {
		return;
	}

	// 逐像素计算census值
	for (sint32 i = 4; i < height - 4; i++) {
		for (sint32 j = 3; j < width - 3; j++) {

			// 中心像素值
			const uint8 gray_center = source[i * width + j];

			// 遍历大小为9x7的窗口内邻域像素，逐一比较像素值与中心像素值的的大小，计算census值
			uint64 census_val = 0u;

			for (sint32 r = -4; r <= 4; r++) {
				for (sint32 c = -3; c <= 3; c++) {
					const uint8 gray = source[(i + r) * width + j + c];

					if (gray < gray_center) {
						census_val = (census_val << 1) + 1;
					}
					else {
						census_val = census_val << 1;
					}
				}
			}

			// 中心像素的census值
			census[i * width + j] = census_val;
		}
	}
}

uint8 adcensus_util::Hamming64(const uint64& x, const uint64& y)
{
	uint64 dist = 0, val = x ^ y;

	// Count the number of set bits
	while (val) {
		++dist;
		val &= val - 1;
	}

	return static_cast<uint8>(dist);
}

void adcensus_util::MedianFilter(const float32* in, float32* out, const sint32& width, const sint32& height, const sint32 wnd_size)
{
	const sint32 radius = wnd_size / 2;
	const sint32 size = wnd_size * wnd_size;
	
	std::vector<float32> wnd_data;
	wnd_data.reserve(size);

	for (sint32 y = 0; y < height; y++) {
		for (sint32 x = 0; x < width; x++) {
			wnd_data.clear();
			for (sint32 r = -radius; r <= radius; r++) {
				for (sint32 c = -radius; c <= radius; c++) {
					const sint32 row = y + r;
					const sint32 col = x + c;
					if (row >= 0 && row < height && col >= 0 && col < width) {
						wnd_data.push_back(in[row * width + col]);
					}
				}
			}
			std::sort(wnd_data.begin(), wnd_data.end());
			if (!wnd_data.empty()) {
				out[y * width + x] = wnd_data[wnd_data.size() / 2];
			}
		}
	}
}

void adcensus_util::equalizeHist_color(cv::Mat src, cv::Mat& dst, bool flag=true)
{
	// 分割通道
	vector<cv::Mat>channels;
	split(src, channels);
	if (flag)	merge(channels, dst);

	cv::Mat blue, green, red;
	blue = channels.at(0);
	green = channels.at(1);
	red = channels.at(2);
	// 分别对BGR通道做直方图均衡化
	equalizeHist(blue, blue);
	equalizeHist(green, green);
	equalizeHist(red, red);
	// 合并通道
	if (flag == false)
	merge(channels, dst);

}
