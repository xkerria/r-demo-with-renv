# 加载 ieugwasr 包
library(ieugwasr)

# 获取 OpenGWAS API 状态
api_status()

# 获取 OpenGWAS 用户信息
user()

# 获取研究 ieu-a-22 的列表
gwasinfo(id="ieu-a-22")

# 获取研究 ieu-a-22 的高频数据，并赋值给变量 exposure_schizophrenia
exposure_schizophrenia <- tophits(id="ieu-a-22", clump=T, r2=0.001, kb=10000)

# head 函数用于输出数据集的前面若干行数据，默认为 6 行
head(exposure_schizophrenia)

# 输出数据集的前 10 行数据
head(exposure_schizophrenia, 10)

# 输出数据集的全部数据
print(exposure_schizophrenia)

