import os
import pandas as pd

base_dir = "./cell2location_map"

results = []

for subdir in os.listdir(base_dir):
    sample_dir = os.path.join(base_dir, subdir)
    if not os.path.isdir(sample_dir):
        continue
    
    for file in os.listdir(sample_dir):
        if file.endswith("_spatialDM_global.csv"):
            csv_path = os.path.join(sample_dir, file)
            
            df = pd.read_csv(csv_path, sep="\t|,", engine="python")
            
            subset = df[df.iloc[:, 0] == "SPP1_CD44"]
            
            if not subset.empty:
                subset = subset.copy()
                subset["file"] = file
                results.append(subset)
                
if results:
    merged_df = pd.concat(results, ignore_index=True)
    print(f"✅ 共提取到 {len(merged_df)} 条 SPP1_CD44 记录")
    display(merged_df)
    merged_df.to_csv(os.path.join(base_dir, "SPP1_CD44_summary.csv"), index=False)
    print("💾 已保存为:", os.path.join(base_dir, "SPP1_CD44_summary.csv"))
else:
    print("⚠️ 没有找到包含 SPP1_CD44 的行")

merged_df = merged_df.rename(columns={"Unnamed: 0": "geneInter"})

merged_df["file_clean"] = merged_df["file"].str.replace(r"^GSM\d+_", "", regex=True)

split_cols = merged_df["file_clean"].str.split("_", expand=True)

merged_df["sample"] = split_cols[0]
merged_df["ligand"] = split_cols[1]
merged_df["ligand_ct"] = split_cols[2]
merged_df["receptor"] = split_cols[3]
merged_df["receptor_ct"] = split_cols[4]

merged_df = merged_df[["geneInter", "global_I", "pval", "sample", "ligand", "ligand_ct", "receptor", "receptor_ct", "file"]]

merged_df.to_csv("./cell2location_map/SPP1_CD44_summary_parsed.csv", index=False)
print("✅ 已保存解析后的结果文件：SPP1_CD44_summary_parsed.csv")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import FormatStrFormatter

df = pd.read_csv("./cell2location_map/SPP1_CD44_summary_parsed.csv")

df["ligand_ct"] = df["ligand_ct"].str.replace("Dutcal", "Ductal")
df["receptor_ct"] = df["receptor_ct"].str.replace("Dutcal", "Ductal")
df["pair"] = df["ligand_ct"] + " → " + df["receptor_ct"]

df["log_pval"] = -np.log10(df["pval"].replace(0, 1e-10))

pair_palette = {
    "Macrophage → Ductal": "#9ecae1",
    "Macrophage → T cells": "#fdae6b",
    "T cells → Ductal": "#a1d99b",
    "T cells → Macrophage": "#fcaeae"
}

sns.set(style="white", context="talk")

plt.figure(figsize=(20, 6))

bubble = sns.scatterplot(
    data=df,
    x="sample",
    y="global_I",
    hue="pair",
    palette=pair_palette,
    size="log_pval",
    sizes=(50, 300),
    alpha=0.7,
    edgecolor=None,
    linewidth=0
)

bubble.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

bubble.set_xticklabels(sorted(df["sample"].unique()), rotation=45, ha="right", fontsize=17)
bubble.tick_params(axis='both', which='major', labelsize=18, width=1.2)

plt.xlabel("Sample", fontsize=18, labelpad=10)
plt.ylabel("Global Moran's I", fontsize=18, labelpad=10)
plt.title("SPP1–CD44 Interaction Strength Across Samples",
          fontsize=20, pad=25)

handles, labels = bubble.get_legend_handles_labels()

pair_handles = handles[1:1+len(pair_palette)]
pair_labels = labels[1:1+len(pair_palette)]

size_handles = handles[1+len(pair_palette):]
size_labels = [lab for lab in labels[2+len(pair_palette):]]
size_labels = [f"{float(lab):.2f}" for lab in size_labels]

plt.legend().remove()

plt.subplots_adjust(right=0.75)
plt.tight_layout(rect=[0, 0, 0.75, 1])

fig = plt.gcf()

legend1 = fig.legend(
    pair_handles, pair_labels,
    title="celltype pairs",
    title_fontsize=18,
    fontsize=15,
    bbox_to_anchor=(0.95, 1),
    loc="upper right",
    frameon=True,
    borderaxespad=0.
)

legend2 = fig.legend(
    size_handles, size_labels,
    title="-log10(p_val)",
    title_fontsize=18,
    fontsize=15,
    bbox_to_anchor=(0.95, 0.75),
    loc="upper right",
    frameon=True,
    borderaxespad=0.
)

plt.savefig("SPP1_CD44_interaction_bubble_Style2_light_splitLegend_bigFont_1decimal.pdf",
            dpi=300, bbox_inches="tight", pad_inches=0.1)
plt.show()