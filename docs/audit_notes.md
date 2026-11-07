# 审计笔记

本页只记行号、错误类型和判断，不复制数据行：RMechDB 和 PMechDB 的许可证是 CC BY-NC-ND。要看原文，用本地下载的 `rmechdb_data/all.csv` 对照行号（文件没有表头，第 1 行就是数据）。

## RMechDB `all.csv`（2026-11-05 下载，2026-11-07 审计）

**命令：** `python scripts/audit_rmechdb.py <all.csv> out.csv`

**结果：** 共 5,426 行。
- 解析成功 5,421 行，其中 5,407 行通过，14 行不一致；
- 5 行箭头编码解析不了。

### 按步骤类型

| 类型 | 标出 / 总数 |
| --- | --- |
| abstraction（夺氢） | 3 / 1,696 |
| recombine（偶联） | 1 / 1,401 |
| homolyze（均裂） | 0 / 1,122 |
| addition（加成） | 8 / 589 |
| retroaddition（逆加成） | 1 / 295 |
| resonance（共振） | 0 / 199 |
| ha resonance | 1 / 119 |

### 逐行判断

「已核」表示逐根箭头手算过；「初判」表示看过、方向明确，但没有逐根手算。

| 行 | 类型 | 核对器报的问题 | 判断 | 程度 |
| --- | --- | --- | --- | --- |
| 279 | addition | half_bond, next_mismatch | 少一根箭头，留下 15–14 半键 | 初判 |
| 393 | addition | next_mismatch | 写出的产物里，自由基和 H 的位置与箭头推出的不符 | 初判 |
| 394 | addition | no_electrons, half_bond | 箭头从炔碳 20 出发，它没有非键电子；应从键 20,21 出发 | 已核 |
| 400 | addition | stereo_mismatch | 没被箭头碰到的环并联碳，在产物里丢了立体标记 | 已核 |
| 404 | addition | half_bond, next_mismatch | 箭头编号指向了错误的原子，产生 21–5、22–5 两个半键 | 初判 |
| 476 | abstraction | no_electrons, half_bond | 鱼钩从 H 原子出发；应从键 S–H（9,7）出发 | 已核 |
| 483 | addition | next_mismatch | 产物与箭头不符，原因未细查 | 初判 |
| 526 | addition | half_bond | 第二根箭头的终点应是原子 2，写成了键 1,5 | 已核 |
| 540 | — | 解析失败 | 编码 `1=2` 格式错误 | 已核 |
| 541 | retroaddition | half_bond | 少一根把电子留在 C10 上的鱼钩 | 已核 |
| 1859 | addition | no_electrons, half_bond | 箭头从 =CH₂ 碳出发，它没有非键电子 | 初判 |
| 2046 | ha resonance | over_limit | 单电子移动的方向写反了：写出的产物是 N⁺/O⁻，箭头却是 O → N | 已核 |
| 2160 | abstraction | half_bond | 终点应是原子 2，写成了键 2,3 | 已核 |
| 2516 | recombine | over_limit | 第二根鱼钩用了 O9，应是 O8 | 已核 |
| 2928、2929 | — | 解析失败 | 编码 `21-22` 应为 `21,22` | 已核 |
| 3970、4852 | — | 解析失败 | 编号被截断：`20,21-2` 应为 `20,21-21` | 已核 |
| 4325 | abstraction | over_limit | 写出的产物里 C11–O10 键断了，但箭头里没有断这根键 | 初判 |

### 开发中修掉的误报（都是核对器的错，不是数据的错）

| 修正 | 修正后通过数 |
| --- | --- |
| 起点：读入时没有标出单电子，`[O]` 被当成 OH | 1,642 → 4,254 |
| 编码分号后的空格把前一根箭头切丢；鱼钩起点规则过严；单电子约定两边不一致 | 4,254 → 5,406 |
| 清除立体时误伤了没被碰到的双键 | 5,406 → 5,407 |

## ChRIMP 转换版 mech-USPTO-31K

抽查 40 步，有 7 步的立体中心在 ChRIMP 生成的下一步中消失，详见 `tests/data/README.md`。全量统计待做（第 3 周）。
