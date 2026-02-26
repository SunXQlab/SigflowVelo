import scanpy as sc
import pandas as pd
import magic
import os
import time

# --- 接口：输入文件 ---
# 脚本 3 的输出文件
input_file = "result\script_3_result\Mono-Macro_count_file.csv"

# --- 接口：输出文件 ---
# 脚本 4 需要的输入文件
output_file = "result/script_3_result/Mono-Macro_scimpute_count.csv"

print(f"--- 步骤 1/5: 正在加载原始 counts: {input_file} ---")
if not os.path.exists(input_file):
    print(f"错误: 输入文件未找到: {input_file}")
    print("请确保脚本 3 已成功运行，并正确生成了此文件。")
else:
    # index_col=0 将第一列（细胞名）设为索引
    df_count = pd.read_csv(input_file, index_col=0)

    # df_count 是 (细胞 x 基因)
    print(f"加载数据成功，形状 (Cells x Genes): {df_count.shape}")

    # 创建 AnnData 对象
    adata = sc.AnnData(df_count)

    # --- 步骤 2/5: 运行 MAGIC (这是一个基于扩散的快速插补) ---
    # 1. 初始化 MAGIC operator
    print("--- 步骤 2/5: 正在初始化 MAGIC... ---")
    # n_jobs=-1 告诉 magic 使用所有可用的 CPU 核心
    magic_op = magic.MAGIC(n_jobs=-1)

    # 2. 拟合-转换
    print("--- 步骤 3/5: 正在运行 MAGIC 插补... ---")
    print("(这将充分利用您的 CPU，请在任务管理器中查看)")
    start_time = time.time()

    # MAGIC 会自动处理 n_jobs，通常会使用所有可用的核心
    adata_magic = magic_op.fit_transform(adata.X)

    end_time = time.time()
    print(f"插补完成。耗时: {end_time - start_time:.2f} 秒。")

    # 3. 提取结果
    # 它的形状是 (细胞 x 基因)
    print(f"插补矩阵形状: {adata_magic.shape}")

    df_imputed = pd.DataFrame(
        adata_magic, 
        index=df_count.index, 
        columns=df_count.columns
    )


    # 确保 'result' 文件夹存在
    os.makedirs("result", exist_ok=True)

    # --- 步骤 5/5: 保存为 CSV ---
    print(f"--- 步骤 5/5: 正在保存插补后的矩阵到: {output_file} ---")
    # index=True 保存行索引 (基因名)
    df_imputed.to_csv(output_file, index=True)

    print("--- Python 插补流程全部完成。---")