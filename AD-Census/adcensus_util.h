/* -*-c++-*- AD-Census - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/AD-Census
* Describe	: header of adcensus_util
*/

#pragma once
#include <algorithm>
#include <iostream>
#include <string.h>
#include "adcensus_types.h"

#include <opencv2/opencv.hpp>

namespace adcensus_util
{
	/**
	* \brief census变换
	* \param source	输入，影像数据
	* \param census	输出，census值数组，预分配空间
	* \param width	输入，影像宽
	* \param height	输入，影像高
	* \param window_size_w	输入，窗口宽
	* \param window_size_h	输入，窗口高
	*/
	void census_transform_color2(const uint8* source, vector<std::string>& census, 
		const sint32& width, const sint32& height, 
		const sint32 window_size_w, const sint32& window_size_h, bool sparse_window = true);

	void census_transform_color(const uint8* source, vector<uint64>& census, const sint32& width, const sint32& height, const sint32 window_size_w, const sint32& window_size_h);

	uint8 HammingDist(const std::string& x, const std::string& y);

	/// <summary>
	/// 原始census变换方案，只针对灰色图像，窗口固定为 7*9
	/// </summary>
	/// <param name="source"></param>
	/// <param name="census"></param>
	/// <param name="width"></param>
	/// <param name="height"></param>
	void census_transform2(const uint8* source, vector<std::string>& census, const sint32& width, const sint32& height);
	void census_transform(const uint8* source, vector<uint64>& census, const sint32& width, const sint32& height);

	// Hamming距离
	uint8 Hamming64(const uint64& x, const uint64& y);

	/**
	* \brief 中值滤波
	* \param in				输入，源数据
	* \param out			输出，目标数据
	* \param width			输入，宽度
	* \param height			输入，高度
	* \param wnd_size		输入，窗口宽度
	*/
	void MedianFilter(const float32* in, float32* out, const sint32& width, const sint32& height, const sint32 wnd_size);


	void equalizeHist_color(cv::Mat src, cv::Mat& dst, bool flag);

}