# 3、Methscan

```bash
cd /share/LCZX_Data/data/allcools/25110891_IR01_Met
tar -xzf allcools.tar.gz
```

```bash
cd /share/LCZX_Data/data/allcools

cat > convert_to_cov_v2.py << 'EOF'
import os
import gzip
import multiprocessing as mp
from glob import glob
import time

def process_single_file(allc_path):
    # 1. 更加严谨的路径替换，确保生成 .cov.gz
    base_dir = os.path.dirname(allc_path)
    file_name = os.path.basename(allc_path)
    new_name = file_name.replace('_allc.gz', '.cov.gz').replace('.allc.gz', '.cov.gz')
    cov_path = os.path.join(base_dir, new_name)

    if os.path.exists(cov_path):
        return True

    try:
        # 设置 compresslevel=1 提高并行写入速度
        with gzip.open(allc_path, 'rt') as fin, \
             gzip.open(cov_path, 'wt', compresslevel=1) as fout:
            for line in fin:
                if line.startswith('chr\t') or line.startswith('chrom'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                strand = parts[2]
                context = parts[3]
                mc = int(parts[4])
                cov = int(parts[5])

                # 只保留 CG 上下文
                if not context.startswith('CG'):
                    continue

                # 🌟 正负链合并逻辑 (CpG dyad collapse)
                if strand == '+':
                    start = pos
                elif strand == '-':
                    start = pos - 1
                else:
                    continue

                unmc = cov - mc
                meth_perc = (mc / cov * 100.0) if cov > 0 else 0.0

                # 标准 Bismark Cov 格式：chr, start, end, perc, mc, unmc
                out_line = f"{chrom}\t{start}\t{start}\t{meth_perc:.6f}\t{mc}\t{unmc}\n"
                fout.write(out_line)
        return True
    except Exception as e:
        print(f"Error processing {allc_path}: {e}")
        return False

if __name__ == '__main__':
    input_dir = "methscan_input"
    if not os.path.exists(input_dir):
        print(f"❌ 错误：找不到目录 {input_dir}")
        exit(1)

    all_files = glob(f"{input_dir}/*allc.gz")
    print(f"🧬 发现 {len(all_files)} 个 ALLC 文件，准备开始转换...")

    cores = max(1, mp.cpu_count() - 2)
    print(f"🚀 使用 {cores} 个核心并行加速...")

    start_time = time.time()
    with mp.Pool(cores) as pool:
        pool.map(process_single_file, all_files)

    print(f"✅ 转换完成！耗时: {(time.time() - start_time)/60:.2f} 分钟。")
EOF

# 挂在后台全速运行
nohup python convert_to_cov_v2.py > format_convert_v2.log 2>&1 &

tail -f format_convert_v2.log
ls -1 methscan_input/*.cov.gz | wc -l
```

```bash
mkdir -p compact_data
nohup methscan prepare /share/LCZX_Data/25110891_IR01_Met/allcools/methscan_input/*.cov.gz ./compact_data > prepare_final.log 2>&1 &
tail -f prepare_final.log
```

```bash
# 确保你在分析目录下
cd /share/LCZX_Data/25110891_IR01_Met

# 后台运行 profile 分析
# 修正了开头的 r，补全了 .csv，并添加了日志输出
nohup methscan profile \
    --strand-column 6 \
    /share/LCZX_Data/ref/human_hg38_TSS.bed \
    ./compact_data \
    TSS_profile.csv > profile_analysis.log 2>&1 &

tail -f profile_analysis.log
```

```bash
# 1. 确保进入你的工作目录
cd /share/LCZX_Data/25110891_IR01_Met

# 2. 创建或清空过滤后的文件夹
mkdir -p filtered_data

# 3. 后台运行过滤任务
# --min-sites 60000: 保证数据量
# --min-meth 20 --max-meth 85: 覆盖人类正常的甲基化区间
nohup methscan filter \
    --min-sites 60000 \
    --min-meth 20 \
    --max-meth 85 \
    ./compact_data \
    ./filtered_data > cell_filter_v2.log 2>&1 &

cat cell_filter_v2.log
```

```bash
# 1. 删除之前那个路径混乱的目录（防止里面有嵌套文件夹）
rm -rf ./smoothed_data

# 2. 重新复制一份干净的 filtered_data 过来
cp -r ./filtered_data ./smoothed_data

# 3. 再次确认一下文件是不是就在根目录下（不应该有嵌套）
ls -l ./smoothed_data/column_header.txt

# 4. 重新启动 Smooth 任务
nohup methscan smooth ./smoothed_data > smoothing.log 2>&1 &

tail -f smoothing.log
```

```bash
# 1. 创建 VMR 输出目录
mkdir -p scan_results

# 2. 运行 scan (寻找高变区域)
# 使用 4 个线程并行计算（这个步骤支持多线程！）
# 输入必须是平滑后的目录 ./smoothed_data
nohup methscan scan \
    --threads 4 \
    ./smoothed_data \
    ./scan_results/VMRs.bed > scan.log 2>&1 &

tail -f scan.log
```

```bash
mkdir -p ./VMR_matrix

# 3. 开启 32 线程全力运行
nohup methscan matrix \
    --threads 32 \
    ./scan_results/VMRs.bed \
    ./smoothed_data \
    ./VMR_matrix > matrix_32threads.log 2>&1 &

tail -f  matrix_32threads.log

# 检查 CSV 的行数（应该等于 VMR 数量 + 1 行表头）
zcat ./VMR_matrix/methylation_fractions.csv.gz | wc -l

# 查看文件
ls -lth ./VMR_matrix/
du -sh ./VMR_matrix/*
```

下游

```bash
library(tidyverse)
library(irlba)

meth_mtx <- read.csv("share/LCZX_Data/25110891_RI01_Met/VMR_matrix/mean_shrunken_residuals.csv.gz", row.names=1) %>% as.matrix()

# PCA that iteratively imputes missing values
prcomp_iterative <- function(x, n=10, n_iter=50, min_gain=0.001, ...) {
  mse <- rep(NA, n_iter)
  na_loc <- is.na(x)
  x[na_loc] = 0  # zero is our first guess

  for (i in 1:n_iter) {
    prev_imp <- x[na_loc]  # what we imputed in the previous round
    # PCA on the imputed matrix
    pr <- prcomp_irlba(x, center = F, scale. = F, n = n, ...)
    # impute missing values with PCA
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    # compare our new imputed values to the ones from the previous round
    mse[i] = mean((prev_imp - new_imp) ^ 2)
    # if the values didn't change a lot, terminate the iteration
    gain <- mse[i] / max(mse, na.rm = T)
    if (gain < min_gain) {
      message(paste(c("\n\nTerminated after ", i, " iterations.")))
      break
    }
  }
  pr$mse_iter <- mse[1:i]
  pr
}


pca <- meth_mtx %>%
  scale(center = T, scale = F) %>%
  prcomp_iterative(n = 5)  # increase this value to e.g. 15 for real data sets

pca_tbl <- as_tibble(pca$x) %>%
  add_column(cell=rownames(meth_mtx))

p <- pca_tbl %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.6, size = 1) +  # 建议增加透明度，防止点太多重叠
  coord_fixed() +
  labs(title="PCA based on VMR methylation")
ggsave("PCA_VMR_methylation.png", plot = p, width = 8, height = 7, dpi = 300)

library(uwot)  # R package for UMAP

umap_obj <- uwot::umap(pca$x, min_dist=0.05, n_neighbors=5, seed=2, ret_nn=T)
umap_tbl <- umap_obj$embedding %>%
  magrittr::set_colnames(c("UMAP1", "UMAP2")) %>%
  as_tibble() %>%
  add_column(cell=rownames(meth_mtx))

umap_tbl %>%
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point() +
  coord_fixed() +
  labs(title="UMAP based on VMR methylation")

library(igraph)  # R package for graph manipulation, also implements the Leiden algorithm

# get the edges of the neighbor graph from the UMAP object
neighbor_graph_edges <-
  tibble(from = rep(1:nrow(umap_obj$nn$euclidean$idx), times=ncol(umap_obj$nn$euclidean$idx)),
         to = as.vector(umap_obj$nn$euclidean$idx),
         weight = as.vector(umap_obj$nn$euclidean$dist)) %>%
  filter(from != to) %>%
  mutate(from = rownames(meth_mtx)[from],
         to = rownames(meth_mtx)[to])

# run Leiden clustering
clust_obj <- neighbor_graph_edges %>%
  igraph::graph_from_data_frame(directed=F) %>%
  igraph::cluster_leiden(resolution_parameter = .5)  # adjust the resolution parameter to your needs

# put the clustering results into a data frame (tibble) for plotting
clust_tbl <- tibble(
  leiden_cluster = as.character(clust_obj$membership),
  cell = clust_obj$names
) %>%
  full_join(umap_tbl, by="cell")

clust_tbl %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = leiden_cluster)) +
  geom_point() +
  coord_fixed()

/share/LCZX_Data/25110891_IR01_Met/VMR_matrix/mean_shrunken_residuals.csv.gz
```