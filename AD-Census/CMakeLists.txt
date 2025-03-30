cmake_minimum_required(VERSION 3.10)
project(AD-Census-Stereo)

# 设置C++标准
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找OpenCV依赖
find_package(OpenCV REQUIRED)

# 添加可执行文件
add_executable(AD-Census-Stereo     main.cpp
    ADCensusStereo.cpp
    cross_aggregator.cpp
    cost_computor.cpp
    scanline_optimizer.cpp
    multistep_refiner.cpp
    adcensus_util.cpp
)

# 包含头文件目录
target_include_directories(AD-Census-Stereo PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OpenCV_INCLUDE_DIRS}
)

# 链接OpenCV库
target_link_libraries(AD-Census-Stereo
    ${OpenCV_LIBS}
)
