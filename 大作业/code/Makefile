# ============================================================
#  多 FPGA 拓扑驱动划分框架 Makefile
#  用法：
#    mingw32-make              # 编译
#    mingw32-make clean        # 清理
# ============================================================

# ---- 构建目录 ----
BUILD_DIR = ../build
TARGET    = $(BUILD_DIR)/main.exe

# ---- C++ 源文件列表 ----
SRC = main.cpp      \
      types.cpp     \
      utils.cpp     \
      candidate.cpp \
      coarsen.cpp   \
      initial.cpp   \
      refine.cpp    \
      validate.cpp

OBJ   = $(SRC:.cpp=.o)
OBJS  = $(addprefix $(BUILD_DIR)/, $(OBJ))

# ---- 编译器配置 (MinGW g++) ----
CXX      = g++
CXXFLAGS = -Wall -Wextra -std=c++17 -g -O2

ALL_OBJS   = $(OBJS)
LINK_CMD   = $(CXX) -o $@ $(ALL_OBJS)

# ============================================================
#  构建目标
# ============================================================
.PHONY: all clean

all: $(TARGET)

# ---- 链接 ----
$(TARGET): $(ALL_OBJS)
	$(LINK_CMD)

# ---- C++ 编译 ----
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# ---- 构建目录创建 ----
$(BUILD_DIR):
	@if not exist "$(BUILD_DIR)" mkdir "$(BUILD_DIR)"

# ---- 清理 ----
clean:
	@if exist "$(BUILD_DIR)" rmdir /s /q "$(BUILD_DIR)"
