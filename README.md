# AD-Census

Improved the AD-Census algorithm, ① adaptive weight cost aggregation algorithm, ② implemented pre-processing and simplified original post-processing filtering algorithm.\
实现了改进的 AD-Census 算法，① 基于自适应权重的代价聚合算法，② 实现了预处理和简化的后处理滤波算法

![图片1](https://s2.loli.net/2024/05/17/vgyYkBfP3j76WlU.png)



## Start

**Implementation Details:** This project requires OpenCV (451 version has been tested feasible).

**Hyper-parameter：** You can set hyper-parameters in the ***’adcensus_types.h‘*** file. The specific meanings of these parameters can be understood by reading the references [1], [2].


## Reference 

[1] Mei X , Sun X , Zhou M , et al. <b>On building an accurate stereo matching system on graphics hardware</b>[C]// IEEE International Conference on Computer Vision Workshops. IEEE, 2012.

[2] 马波涛. 基于双目立体视觉的心脏软组织三维重构技术研究[D].电子科技大 学,2017.

This work is based on repository: https://github.com/ethan-li-coding/AD-Census, we appreciate this work.
