#include "cost_aggregator_2.h"
#include <iostream>

Aggregator::Aggregator() : width_(0), height_(0), img_left_(nullptr), img_right_(nullptr),
cost_init_(nullptr),
min_disparity_(0), max_disparity_(0), window_size_w_(7), window_size_h_(7), 
c1_(10), c2_(5), d1_(10), d2_(15), lambda_c_(0.05), lambda_d_(0.3), spars_window_(true),
is_initialized_(false) { }

Aggregator::~Aggregator()
{

}

bool Aggregator::Initialize(const sint32& width, const sint32& height, const sint32& min_disparity, const sint32& max_disparity)
{
	width_ = width;
	height_ = height;
	min_disparity_ = min_disparity;
	max_disparity_ = max_disparity;

	const sint32 img_size = width_ * height_;
	const sint32 disp_range = max_disparity_ - min_disparity_;
	if (img_size <= 0 || disp_range <= 0) {
		is_initialized_ = false;
		return is_initialized_;
	}

	// 为聚合代价数组分配内存
	cost_aggr_.clear();
	cost_aggr_.resize(img_size * disp_range);

	is_initialized_ = !cost_aggr_.empty();

	return is_initialized_;
}

void Aggregator::SetData(const uint8* img_left, const uint8* img_right, const float32* cost_init)
{
	img_left_ = img_left;
	img_right_ = img_right;
	cost_init_ = cost_init;
}

void Aggregator::SetParams(const uint32& window_size_w, const uint32& window_size_h,
	const float32& c1, const float32& c2, const float32& d1, const float32& d2,
	const float32& lambda_c, const float32& lambda_d, const bool& spars_window)
{
	spars_window_ = spars_window;
	window_size_w_ = window_size_w;
	window_size_h_ = window_size_h;
	c1_ = c1;
	c2_ = c2;
	d1_ = d1;
	d2_ = d2;
	lambda_c_ = lambda_c;
	lambda_d_ = lambda_d;

	std::cout << std::endl;
	std::cout << "Testing Adaptive Aggregation of Weights... " << std::endl;
	std::cout << "Use sparse window: " << spars_window_ << std::endl;
	std::cout << "window_size: w * h = " << window_size_w_ << " * " << window_size_h_ << std::endl;
	std::cout << "c1 = " << c1_ << ", c2 = " << c2_ << ", d1 = " << d1_ << ", d2 = " << d2_ << std::endl;
	std::cout << "lambda_c = " << lambda_c_ << ", lambda_d = " << lambda_d_ << std::endl;
}


void Aggregator::Aggregate()
{
	uint8 step = 1;
	uint32 Np = window_size_h_ * window_size_w_;
	if (spars_window_) {
		step = 2;
		Np = ((window_size_h_ + 1) / 2) * ((window_size_w_ + 1) / 2);
	}
	sint32 disp_range = max_disparity_ - min_disparity_;
	sint32 win_w2 = (window_size_w_ - 1) / 2;
	sint32 win_h2 = (window_size_h_ - 1) / 2;

	//初始化权重数据
	wight_l_.resize(width_ * height_ * Np);
	wight_r_.resize(width_ * height_ * Np);

	// 遍历窗口内的像素，记录权重值
	for (sint32 y = 0; y < height_; y++) {
		for (sint32 x = 0; x < width_; x++) {

			sint32 Np_i = 0;
			for (sint32 ind_y = -win_h2; ind_y <= win_h2; ind_y += step) {
				for (sint32 ind_x = -win_w2; ind_x <= win_w2; ind_x += step) {
					if (y + ind_y < 0 || y + ind_y >= height_ || x + ind_x < 0 || x + ind_x >= width_) {
						wight_l_[y * width_ * Np + x * Np + Np_i] = 0.f;
						wight_r_[y * width_ * Np + x * Np + Np_i] = 0.f;
						Np_i++; 
						continue;
					}
					float32 dis_cil = 0;
					float32 dis_cir = 0;
					for (sint32 c = 0; c < 3; c++) {
						float32 I_pl = img_left_[y * width_ * 3, x * 3, c];
						float32 I_ql = img_left_[(y + ind_y) * width_ * 3, (x + ind_x) * 3, c];
						float32 I_pr = img_right_[y * width_ * 3, x * 3, c];
						float32 I_qr = img_right_[(y + ind_y) * width_ * 3, (x + ind_x) * 3, c];
						dis_cil = dis_cil + pow(I_pl - I_ql, 2);
						dis_cir = dis_cir + pow(I_pr - I_qr, 2);
					}
					float32 dis_cl = sqrt(dis_cil);
					float32 dis_cr = sqrt(dis_cir);
					float32 dis_dl = sqrt(pow(ind_x, 2) + pow(ind_y, 2));
					float32 dis_dr = sqrt(pow(ind_x, 2) + pow(ind_y, 2));
					// 计算p、q两点的权重值
					float32 w_pql = ComputeW(dis_cl, dis_dl);
					float32 w_pqr = ComputeW(dis_cr, dis_dr);
					wight_l_[y * width_ * Np + x * Np + Np_i] = w_pql;
					wight_r_[y * width_ * Np + x * Np + Np_i] = w_pqr;
					Np_i++;
				}
			}
		}
	}

	// 遍历每个像素，求聚合代价
	for (sint32 y = 0; y < height_; y++) {
		for (sint32 x = 0; x < width_; x++) {

			for (sint32 d = min_disparity_; d < max_disparity_; d++) {
				sint32 xr = x - d;
				/*if (min_disparity_ != 0) {
					d = d - min_disparity_;
				}*/
				if (xr < 0 || xr >= width_) {
					cost_aggr_[y * width_ * disp_range + x * disp_range + d] = Large_Float;
					continue;
				}
				//// 背景像素聚合代价无穷大
				//if ((img_left_[y * width_ * 3 + x * 3] == 0 && img_left_[y * width_ * 3 + x * 3 + 1] == 0
				//	&& img_left_[y * width_ * 3 + x * 3 + 2] == 0) || (img_right_[y * width_ * 3 + xr * 3] == 0 &&
				//		img_right_[y * width_ * 3 + xr * 3 + 1] == 0 && img_right_[y * width_ * 3 + xr * 3 + 2] == 0)) {
				//	cost_aggr_[y * width_ * disp_range + x * disp_range + d] = Large_Float;
				//	continue;
				//}
				float32 E1 = 0.f;
				float32 E2 = 0.f;
				sint32 Np_i = 0;
				for (sint32 ind_y = -win_h2; ind_y <= win_h2; ind_y += step) {
					for (sint32 ind_x = -win_w2; ind_x <= win_w2; ind_x += step) {
						sint32 yq = y + ind_y;
						sint32 xq = x + ind_x;
						if (yq < 0 || yq >= height_ || xq < 0 || xq >= width_) {
							Np_i++;
							continue;
						}
						float32 cqd = cost_init_[yq * width_ * disp_range + xq * disp_range + d];
						float32 wl = wight_l_[y * width_ * Np + x * Np + Np_i];
						float32 wr = wight_r_[y * width_ * Np + xr * Np + Np_i];
						Np_i++;

						E1 = E1 + wl * wr * cqd;
						E2 = E2 + wl * wr;
					}
				}

				// E(i, j, d) = E(pl, d) = E1/E2
				if (E2 != 0) {
					cost_aggr_[y * width_ * disp_range + x * disp_range + d] = E1 / E2;
				}
				// E(i, j, d) = C(i, j, d)
				else {
					cost_aggr_[y * width_ * disp_range + x * disp_range + d] = 
						cost_init_[y * width_ * disp_range + x * disp_range + d];
				}
			}
		}
	}

}


float32* Aggregator::get_cost_ptr()
{
	if (!cost_aggr_.empty()) {
		return &cost_aggr_[0];
	}
	else {
		return nullptr;
	}
}

float32* Aggregator::get_wightL_ptr()
{
	if (!wight_l_.empty()) {
		return &wight_l_[0];
	}
	else {
		return nullptr;
	}
}

float32* Aggregator::get_wightR_ptr()
{
	if (!wight_r_.empty()) {
		return &wight_r_[0];
	}
	else {
		return nullptr;
	}
}

float32 Aggregator::ComputeW(float32 Dc, float32 Dd)
{
	if (Dc <= c1_ && Dd <= d1_) {
		return 1.0f;
	}
	else if (Dc <= c2_ && Dd > d1_ && Dd <= d2_) {
		return 1.0f;
	}
	else {
		float32 w = exp(-Dc / lambda_c_ - Dd / lambda_d_);
		return w;
	}
}
