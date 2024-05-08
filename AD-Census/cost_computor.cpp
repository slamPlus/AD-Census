/* -*-c++-*- AD-Census - Copyright (C) 2020.
* Author	: Yingsong Li(Ethan Li) <ethan.li.whu@gmail.com>
* https://github.com/ethan-li-coding/AD-Census
* Describe	: implement of class CostComputor
*/

#include "cost_computor.h"
#include "adcensus_util.h"
#include <core.hpp>
#include <iostream>

CostComputor::CostComputor() :
	width_(0), height_(0), img_left_(nullptr), img_right_(nullptr), census_method_(0),
	lambda_ad_(0), lambda_census_(0), min_disparity_(0), max_disparity_(0),
	census_window_width_(3), census_window_height_(3), Ta_(0.1), Th_(0.3), lambda_(0.5),
	is_initialized_(false) { }

CostComputor::~CostComputor()
{
	
}

bool CostComputor::Initialize(
	const sint32& width, const sint32& height, 
	const sint32& min_disparity, const sint32& max_disparity)
{
	width_ = width;
	height_ = height;
	min_disparity_ = min_disparity;
	max_disparity_ = max_disparity;

	const sint32 img_size = width_ * height_;
	const sint32 disp_range = max_disparity_ - min_disparity_;
	if (img_size <= 0 || disp_range <= 0) {
		is_initialized_ = false;
		return false;
	}

	// 灰度数据（左右影像）
	gray_left_.resize(img_size);
	gray_right_.resize(img_size);
	// census数据（左右影像）
	census_left_.resize(img_size);
	census_right_.resize(img_size);
	// 初始代价数据
	cost_init_.resize(img_size * disp_range);

	is_initialized_ = !gray_left_.empty() && !gray_right_.empty() && !census_left_.empty() && !census_right_.empty() && !cost_init_.empty();
	return is_initialized_;
}

void CostComputor::SetData(const uint8* img_left, const uint8* img_right)
{
	img_left_ = img_left;
	img_right_ = img_right;
}

void CostComputor::SetParams(
	const uint8 census_method, const sint32& lambda_ad, const sint32& lambda_census,
	const sint32& census_window_width, const sint32& census_window_height,
	const float32& Ta, const float32& Th, const float32& lambda)
{
	census_method_ = census_method;
	lambda_ad_ = lambda_ad;
	lambda_census_ = lambda_census;
	Ta_ = Ta;
	Th_ = Th;
	lambda_ = lambda;	
	census_window_height_ = census_window_height;
	census_window_width_ = census_window_width;

	std::cout << std::endl;
	if (census_method_ == 0) {
		std::cout << "AD_Census：lambda_ad = " << lambda_ad_ << ", lambda_census = " << lambda_census_ << std::endl;
		std::cout << "census window size: w * h = " << " 9 * 7 " << std::endl;
	}
	else if (census_method_ == 1) {
		std::cout << "Adaptive Method：Ta = " << Ta_ << ", Th = " << Th_ << ", lambda = " << lambda_ << std::endl;
		std::cout << "census window size: w * h = " << census_window_width_ << " * " << census_window_height_ << std::endl;
	}

}

void CostComputor::ComputeGray()
{
	// 彩色转灰度
	for (sint32 n = 0; n < 2; n++) {
		const auto color = (n == 0) ? img_left_ : img_right_;
		auto& gray = (n == 0) ? gray_left_ : gray_right_;
		for (sint32 y = 0; y < height_; y++) {
			for (sint32 x = 0; x < width_; x++) {
				const auto b = color[y * width_ * 3 + 3 * x];
				const auto g = color[y * width_ * 3 + 3 * x + 1];
				const auto r = color[y * width_ * 3 + 3 * x + 2];
				gray[y * width_ + x] = uint8(r * 0.299 + g * 0.587 + b * 0.114);
			}
		}
	}
}

void CostComputor::CensusTransform()
{
	// 左右影像census变换
	if (census_method_ == 0) {
		adcensus_util::census_transform(&gray_left_[0], census_left_, width_, height_);
		adcensus_util::census_transform(&gray_right_[0], census_right_, width_, height_);
	}
	
	else if (census_method_ == 1) {
		adcensus_util::census_transform_color(img_left_, census_left_, width_, height_, census_window_width_, census_window_height_);
		adcensus_util::census_transform_color(img_right_, census_right_, width_, height_, census_window_width_, census_window_height_);
	}
	
}

void CostComputor::ComputeCost()
{
	//float32 cost_census = static_cast<float32>(adcensus_util::HammingDist("1001", "1000"));
	const sint32 disp_range = max_disparity_ - min_disparity_;
	// 预设参数
	const auto lambda_ad = lambda_ad_;
	const auto lambda_census = lambda_census_;

	// 计算代价
	for (sint32 y = 0; y < height_; y++) {
		for (sint32 x = 0; x < width_; x++) {
			const auto bl = img_left_[y * width_ * 3 + 3 * x];
			const auto gl = img_left_[y * width_ * 3 + 3 * x + 1];
			const auto rl = img_left_[y * width_ * 3 + 3 * x + 2];
			const uint64 census_val_l = census_left_[y * width_ + x];
			// 逐视差计算代价值
			for (sint32 d = min_disparity_; d < max_disparity_; d++) {
				auto& cost = cost_init_[y * width_ * disp_range + x * disp_range + (d - min_disparity_)];
				const sint32 xr = x - d;

				// 超出右图图像范围
				if (xr < 0 || xr >= width_) {
					cost = 1.0f;
					continue;
				}

				const auto br = img_right_[y * width_ * 3 + 3 * xr];
				const auto gr = img_right_[y * width_ * 3 + 3 * xr + 1];
				const auto rr = img_right_[y * width_ * 3 + 3 * xr + 2];

				// 背景区域匹配代价设为1
				if ((bl == 0 && gl == 0 && rl == 0) || (br == 0 && gr == 0 && rr == 0)) {
					cost = 1.0f;
					continue;
				}

				// ad代价
				const float32 cost_ad = (abs(bl - br) + abs(gl - gr) + abs(rl - rr)) / 3.0f / 255.0f;

				// census代价
				const uint64 census_val_r = census_right_[y * width_ + xr];
				float32 cost_census;
				// ad-census代价
				if (census_method_ == 0) {
					cost_census = static_cast<float32>(adcensus_util::Hamming64(census_val_l, census_val_r));
					//if (census_val_l.length() != 0)	cost_census = cost_census / (census_val_l.length());
					cost = 1 - exp(-cost_ad / lambda_ad) + 1 - exp(-cost_census / lambda_census);
				}
				else if (census_method_ == 1) {
					cost_census = static_cast<float32>(adcensus_util::Hamming64(census_val_l, census_val_r));
					//if (census_val_l.length() != 0)	cost_census = cost_census / (census_val_l.length());
					cost = lambda_ * std::min(Ta_, cost_ad) + (1 - lambda_) * std::min(Th_, cost_census);
				}
			}
		}
	}
}

void CostComputor::Compute()
{
	if(!is_initialized_) {
		return;
	}

	// 计算灰度图
	if (census_method_ == 0) {
		ComputeGray();
	}
	
	// census变换
	CensusTransform();

	// 代价计算
	ComputeCost();
}

float32* CostComputor::get_cost_ptr()
{
	if (!cost_init_.empty()) {
		return &cost_init_[0];
	}
	else {
		return nullptr;
	}
}
