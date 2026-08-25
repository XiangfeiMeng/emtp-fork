#!/bin/bash
# EMTP 仿真自动化执行脚本（Linux/Mac/WSL 适用）

echo -e "\033[36m===== EMTP 仿真自动化脚本开始执行 =====\033[0m"

# 定义错误处理函数
error_exit() {
    echo -e "\033[31m❌ 脚本执行失败: $1\033[0m"
    exit 1
}

# 1. 切换到上级目录（当前是 build/bin，上级是 build）
echo -e "\n\033[32m[步骤1/5] 切换到上级目录（build 目录）...\033[0m"
cd .. || error_exit "切换到 build 目录失败，请检查当前路径：$(pwd)"
echo -e "✅ 当前目录: $(pwd)"

# 2. 执行 ninja 构建
echo -e "\n\033[32m[步骤2/5] 执行 ninja 构建...\033[0m"
if ! command -v ninja &> /dev/null; then
    error_exit "未找到 ninja 命令，请先安装并添加到环境变量"
fi
ninja || error_exit "ninja 执行失败，退出码: $?"
echo -e "✅ ninja 构建完成"

# 3. 切换到 bin 目录（回到最初的执行目录）
echo -e "\n\033[32m[步骤3/5] 切换到 bin 目录...\033[0m"
cd bin || error_exit "切换到 bin 目录失败，请检查路径：$(pwd)"
echo -e "✅ 当前目录: $(pwd)"

# 4. 执行 emtp_sim 程序
echo -e "\n\033[32m[步骤4/5] 执行 emtp_sim 仿真程序...\033[0m"
if [ ! -f "./emtp_sim" ]; then
    error_exit "未找到 emtp_sim 可执行文件（当前目录：$(pwd)）"
fi
./emtp_sim || error_exit "emtp_sim 执行失败，退出码: $?"
echo -e "✅ emtp_sim 仿真完成"

# 5. 执行 Python 绘图脚本（curve.py 需在 bin 目录，若在上级目录则改为 python ../curve.py）
echo -e "\n\033[32m[步骤5/5] 执行 curve.py 绘图脚本...\033[0m"
if ! command -v python &> /dev/null; then
    if ! command -v python3 &> /dev/null; then
        error_exit "未找到 Python 环境，请先安装并添加到环境变量"
    fi
    PY_CMD="python3"
else
    PY_CMD="python"
fi

if [ ! -f "./curve.py" ]; then
    error_exit "未找到 curve.py 脚本（当前目录：$(pwd)）"
fi

$PY_CMD curve.py || error_exit "curve.py 执行失败，退出码: $?"
echo -e "✅ curve.py 绘图完成"

echo -e "\n\033[32m===== 所有步骤执行完成！=====\033[0m"
