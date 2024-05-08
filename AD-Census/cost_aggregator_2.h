#ifndef AD_CENSUS_CROSS_AGGREGATOR_2_H_
#define AD_CENSUS_CROSS_AGGREGATOR_2_H_

#include "adcensus_types.h"
#include <algorithm>

/**
 * \brief 基于自适应权重的匹配代价聚合
 */
class Aggregator {
public:
	Aggregator();
	~Aggregator();

	/**
	 * \brief 初始化代价聚合器
	 * \param width		影像宽
	 * \param height	影像高
	 * \return true:初始化成功
	 */
	bool Initialize(const sint32& width, const sint32& height, const sint32& min_disparity, const sint32& max_disparity);

	/**
	 * \brief 设置代价聚合器的数据
	 * \param img_left		// 左影像数据，三通道
	 * \param img_right		// 右影像数据，三通道
	 * \param cost_init		// 初始代价数组
	 */
	void SetData(const uint8* img_left, const uint8* img_right, const float32* cost_init);

	/**
	 * \brief 设置代价聚合器的参数
	 * \param 
	 */
	void Aggregator::SetParams(const uint32& window_size_w, const uint32& window_size_h,
		const float32& c1, const float32& c2, const float32& d1, const float32& d2,
		const float32& lambda_c, const float32& lambda_d, const bool& spars_window = true);

	/** \brief 聚合 */
	void Aggregate();

	/** \brief 获取聚合代价数组指针 */
	float32* get_cost_ptr();

	/** \brief 获取权重数据数组指针 */
	float32* get_wightL_ptr();
	float32* get_wightR_ptr();

private:
	/** \brief 图像尺寸 */
	sint32	width_;
	sint32	height_;

	/** \brief 影像数据 */
	const uint8* img_left_;
	const uint8* img_right_;

	/** \brief 初始代价数组指针 */
	const float32* cost_init_;
	/** \brief 聚合代价数组 */
	vector<float32> cost_aggr_;

	/** \brief 权重数据 */
	vector<float32> wight_l_;
	vector<float32> wight_r_;


	sint32  min_disparity_;			// 最小视差
	sint32	max_disparity_;			// 最大视差

	/** \brief 是否成功初始化标志	*/
	bool is_initialized_;

	/** \brief 是否使用稀疏窗	*/
	bool spars_window_;

	/** \brief 窗口尺寸	*/
	uint32 window_size_w_; 
	uint32 window_size_h_;
	
	/** \brief 设定的颜色差值的阈值	*/
	float32 c1_;
	float32 c2_;

	/** \brief 设定的距离的阈值	*/
	float32 d1_;
	float32 d2_;

	/** \brief 调节颜色差值与距离对权重的影响	*/
	float32 lambda_c_; 
	float32 lambda_d_;

	float32 ComputeW(float32 Dc, float32 Dd);
};
#endif