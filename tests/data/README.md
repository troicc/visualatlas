# 测试数据

## uspto31k_steps.csv

从 mech-USPTO-31K 中随机抽取的 40 个基元步，每条反应抽一步。各列：

| 列 | 内容 |
| --- | --- |
| `step` | 这一步的 MechSMILES |
| `expected_product` | 期望产物：不带原子映射的规范 SMILES，保留不参与反应的立体中心 |
| `chrimp_product` | ChRIMP 自己应用箭头得到的产物 |
| `stereo_dropped_by_chrimp` | ChRIMP 的产物是否把立体标记弄丢了 |

**来源：**
- Chen, Babazade, Kim, Han, Jung, *A large-scale reaction dataset of mechanistic pathways of organic reactions*, Sci. Data 11, 863 (2024)，DOI 10.1038/s41597-024-03709-y。
- 数据发布在 figshare（10.6084/m9.figshare.24797220.v2），许可证为 CC BY 4.0。
- MechSMILES 格式的转换取自 ChRIMP（github.com/schwallergroup/ChRIMP，MIT）的 `data/preprocessed/mech_uspto_31k.csv`。

**期望产物怎么核实的：** 每一行都满足下面三条。
1. `chrimp_product` 与数据集中同一反应下一步的反应物逐字相同（比较时用规范 SMILES，去掉原子映射）。
2. 把 `expected_product` 的立体标记去掉后，与 `chrimp_product` 相同。
3. 这一步没有任何一根箭头碰到立体中心，所以原有的立体标记应原样保留。

**ChRIMP 丢了立体标记：** 40 行中有 7 行，`stereo_dropped_by_chrimp = 1`。
- 原始的 mech-USPTO-31K 只存整条反应和箭头列表，立体标记都在：6155/31364 条反应带立体中心。
- 中间体是 ChRIMP 转换时用自己的箭头应用生成的，生成时立体标记被丢掉了。

**许可证：** PMechDB 和 RMechDB 的许可证是 CC BY-NC-ND，不允许发布转换后的版本，所以这两个库的数据一律不放进本仓库。
