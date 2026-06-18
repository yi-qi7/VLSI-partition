# ============================================================
#  TopoPart Makefile — ICCAD 2021 多 FPGA 拓扑驱动划分框架
#
#  用法 (在 VS Developer Command Prompt 中)：
#    mingw32-make              # 纯 CPU 版本
#    mingw32-make CUDA=1       # GPU 加速版本
#    mingw32-make clean        # 清理
# ============================================================

# ---- 构建目录 ----
BUILD_DIR = ../build
TARGET    = $(BUILD_DIR)/topopart.exe

# ---- TopoPart C++ 源文件 ----
SRC = topopart_main.cpp      \
      topopart_types.cpp     \
      topopart_utils.cpp     \
      topopart_candidate.cpp \
      topopart_coarsen.cpp   \
      topopart_initial.cpp   \
      topopart_refine.cpp    \
      topopart_validate.cpp

OBJ   = $(SRC:.cpp=.obj)
OBJS  = $(addprefix $(BUILD_DIR)/, $(OBJ))

# ---- TopoPart CUDA 源文件（可选，设置 CUDA=1 启用） ----
CU_SRC  = topopart_cuda.cu
CU_OBJ  = $(CU_SRC:.cu=.obj)
CU_OBJS = $(addprefix $(BUILD_DIR)/, $(CU_OBJ))

# ---- 编译器 (MSVC cl.exe) ----
CXX      = cl
CXXFLAGS = /nologo /utf-8 /W3 /std:c++17 /Zi /EHsc /bigobj /O2

# ---- CUDA 配置（CUDA=1 时生效） ----
ifdef CUDA
  CUDA_PATH  = C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.3
  NVCC       = "$(CUDA_PATH)/bin/nvcc"
  NVCCFLAGS  = -arch=sm_120 -std=c++17 -O2
  NVCCFLAGS += -Xcompiler /utf-8,/W3,/bigobj,/EHsc
  CUDA_INC   = /I"$(CUDA_PATH)/include"
  NVCC_INC   = -I"$(CUDA_PATH)/include"
  LDFLAGS    = -L"$(CUDA_PATH)/lib/x64" -lcudart

  CXXFLAGS  += /DTOPOPART_USE_CUDA $(CUDA_INC)
  ALL_OBJS   = $(OBJS) $(CU_OBJS)
  LINK_CMD   = $(NVCC) $(NVCCFLAGS) -o $@ $(ALL_OBJS) $(LDFLAGS)
else
  ALL_OBJS   = $(OBJS)
  LINK_CMD   = link /nologo /DEBUG /out:$@ $(ALL_OBJS)
endif

# ============================================================
#  目标
# ============================================================
.PHONY: all clean

all: $(TARGET)

# ---- 链接 ----
$(TARGET): $(ALL_OBJS)
	$(LINK_CMD)

# ---- C++ 编译 (MSVC) ----
$(BUILD_DIR)/topopart_%.obj: topopart_%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) /c /Fo$@ $<

# ---- CUDA 编译 (nvcc) ----
ifdef CUDA
$(BUILD_DIR)/topopart_cuda.obj: topopart_cuda.cu | $(BUILD_DIR)
	$(NVCC) $(NVCCFLAGS) $(NVCC_INC) -c $< -o $@
endif

# ---- 构建目录 ----
$(BUILD_DIR):
	@if not exist "$(BUILD_DIR)" mkdir "$(BUILD_DIR)"

# ---- 清理 ----
clean:
	@if exist "$(BUILD_DIR)" rmdir /s /q "$(BUILD_DIR)"