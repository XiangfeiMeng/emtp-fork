cmake_minimum_required(VERSION 3.16)
project(EMTP_Simple CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -march=native")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# 1) Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# 2) SuiteSparse (你原来的写法保留)
find_library(SUITESPARSE_KLU klu)
find_library(SUITESPARSE_AMD amd)
find_library(SUITESPARSE_BTF btf)
find_library(SUITESPARSE_COLAMD colamd)
find_library(SUITESPARSE_CONFIG suitesparseconfig)

# 3) 让 Eigen 使用外部 BLAS（LAPACKE 不是必须，先不强制）
add_compile_definitions(EIGEN_USE_BLAS)

# 4) 用 FindBLAS / FindLAPACK（Ubuntu 22.04 安装 libopenblas-dev 后就能找到）
find_package(BLAS REQUIRED)        # 会解析到 OpenBLAS
find_package(LAPACK)               # 可选；若要让 Eigen 用 LAPACK 再加 EIGEN_USE_LAPACKE

# 如果你真的需要 LAPACKE（C 接口），可以在装了 liblapacke-dev 后启用：
# add_compile_definitions(EIGEN_USE_LAPACKE)
# find_package(LAPACKE)  # 注意：老版本 CMake 可能没有 FindLAPACKE 模块

# 5) 头文件
include_directories(src)

# 6) 可执行文件
add_executable(emtp_sim
    main.cpp
    Control.cpp
    Grid.cpp
    Simulation.cpp
    Curve.cpp
    Devices.cpp
    converter.cpp
    KluSolverEigen.cpp
    transformer.cpp
    Generator.cpp  
)

# 7) 链接库（注意这里的 target 是 emtp_sim，不是 emtp）
target_link_libraries(emtp_sim PRIVATE
    Eigen3::Eigen
    ${SUITESPARSE_KLU}
    ${SUITESPARSE_AMD}
    ${SUITESPARSE_BTF}
    ${SUITESPARSE_COLAMD}
    ${SUITESPARSE_CONFIG}
)

# 链接 BLAS（优先使用现代的导入目标）
if (TARGET BLAS::BLAS)
  target_link_libraries(emtp_sim PRIVATE BLAS::BLAS)
else()
  target_link_libraries(emtp_sim PRIVATE ${BLAS_LIBRARIES})
endif()

# （可选）链接 LAPACK
if (LAPACK_FOUND)
  if (TARGET LAPACK::LAPACK)
    target_link_libraries(emtp_sim PRIVATE LAPACK::LAPACK)
  else()
    target_link_libraries(emtp_sim PRIVATE ${LAPACK_LIBRARIES})
  endif()
endif()

# 8) 输出目录（可选）
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
