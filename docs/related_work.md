# 相关工作笔记（留给以后写该领域综述用）

2026-11-07 起记。每条标明核实状态：
- **已核实**：读过原文、摘要或源码；
- **未核实**：只见过检索摘要或二手描述，引用前必须再查。

## 与本项目最接近的三项

### Phan, *Lewis-labeled graphs: curly arrows and fishhooks as executable electron transfers*

- **出处：** arXiv 2607.26088，2026-07-27。**已核实**（读了摘要）。
- **做了什么：**
  - 在带孤对、单电子、σ/π 标签的分子图上，把弯箭头和鱼钩当作可执行的电子转移。
  - 在 39,732 个参考端点上复现了产物。
  - 严格重放时，160 个人工审过的基元步全部接受，1,120 个人为制造的错误全部拒绝。
- **对应的软件：** SynKit 1.6.x（JCIM 2025，10.1021/acs.jcim.5c02123）。PyPI 页面写明能核对用户给出的电子流机理。正文**未核实**。
- **摘要没有提到的：** 金属、立体化学、画图。
- **与本项目的关系：** 核对器的核心思路相同（箭头即电子转移，重放后与写出的结构比对）。实验三应把它当作对照组。
- **作者评价（2026-11-07）：** 文中的图画得很丑。

### MolDraw 的 MolScript（Scidart Academy）

- **出处：** 博客 https://www.moldraw.com/blog/molscript-reaction-mechanism-language.html ，2026-10-17。**已核实**（读了博客原文）。
- **原文宣称：**
  - 弯箭头自动摆放，避开键、标签、电荷、电子点和其他箭头；几根箭头指向同一个目标时会散开；摆不下时照样画出来，并报一条警告。
  - 每一步重放电子，算出的结构必须与写出的下一格一致；检查电荷和原子守恒、八隅律、鱼钩是否单电子、共振。
  - 路线图里有有机金属的氧化态和金属–配体键，支持画催化循环。
  - 打包成 AI 技能，配了 41 条请求的测评集，参考脚本得 100%。
- **没有的：**
  - 不开源，开源列为「下一步」；
  - 没有同行评审；
  - 没有和其他画法对比，也没有违规率这类数字；
  - 没提立体化学。
- **作者评价（2026-11-07）：** 还不够好看，尤其是自由基的箭头。

### MechSMILES / ChRIMP（Neukomm, Jončev, Schwaller）

- **出处：** arXiv 2512.05722；代码 github.com/schwallergroup/ChRIMP（MIT）。**已核实**（读了源码并本地运行过）。
- **核对：**
  - `MoleculeSet.check_step_legality` 判断两个状态之间是否存在一种合法的单步箭头组合；
  - 另有原子与电荷守恒的预检；
  - `MechSmiles.check_validity()` 是占位函数。
- **电子数上限：**
  - 第二周期 8；P、S、Cl 12；第三周期其余元素 18；Br、I 按 8；
  - 第 4–5 周期（含过渡金属）18；
  - 只支持到 Ba，La 以后直接报错「未实现」。
- **自由基：** 代码里有单电子进攻和均裂，论文里列为以后的工作。
- **画图：** 用 RDKit 画每一格，弯箭头是贝塞尔曲线，不避让。实际跑出来的问题：箭头穿过苯环的键，Na⁺ 与 H–H 叠在一起，Br 标签压在键上，多步图里试剂漂在结构上方。
- **数据问题（本项目发现）：** 生成中间体时会丢掉立体标记。抽查 40 步中有 7 步如此，见 `tests/data/README.md`。

## 机理数据集

| 数据集 | 出处 | 规模 | 许可证 | 核实 |
| --- | --- | --- | --- | --- |
| mech-USPTO-31K | Chen et al., Sci. Data 11, 863 (2024)，10.1038/s41597-024-03709-y | 31,364 条反应，只存整条反应和箭头列表，不存中间体 | figshare 10.6084/m9.figshare.24797220.v2，CC BY 4.0 | 已核实 |
| PMechDB | 10.1021/acs.jcim.3c01810 | — | CC BY-NC-ND 4.0 | 已核实许可证 |
| RMechDB | Tavakoli et al., JCIM 2023，10.1021/acs.jcim.2c01359 | `all.csv` 中 5,425 步；均裂 1,122、偶联 1,401、夺氢 1,697、加成 592、逆加成 295、共振 199、夺氢共振 119 | CC BY-NC-ND | 已核实（本地数据，2026-11-05 下载） |
| FlowER | Nature 2025，10.1038/s41586-025-09426-9 | — | — | 未核实 |
| oMeBench | arXiv 2510.07731 | — | — | 未核实 |

## 其他相关（都未核实，来自检索）

| 工作 | 出处 | 内容 |
| --- | --- | --- |
| Andersen et al., *Representing catalytic mechanisms with rule composition* | arXiv 2201.04515 (2022) | 把箭头推动写成图重写规则，用于 M-CSA |
| Flamm, Müller, Stadler, *Every atom-atom map can be explained by electron pushing diagrams* | arXiv 2311.13492 | 电子推动的理论 |
| Chen & Baldi, *No Electron Left Behind*（Reaction Explorer） | JCIM 2009，10.1021/ci900157k | 基于规则的电子流步骤，用于教学 |
| Kayala et al., *Learning to Predict Chemical Reactions* | JCIM 2011，10.1021/ci200207y | 逐步更新键级和电荷 |
| Bran, Neukomm, … Schwaller | arXiv 2503.08537 | 用 LLM 评判每一步电子推动是否合理（不是形式核对） |
| FukuyamaBench | arXiv 2607.12771 (2026) | LLM 机理基准 |
| MAELLE | arXiv 2608.27429 (2026) | 对电子占据做流匹配 |
| *The Mechanisms App*（Winter et al.） | ChemRxiv 2020 | 教学用，即时反馈箭头的起止点 |
| ACE Organic | 作业系统 | 给机理打分，标出不合理的步骤 |
| SMiCRM | arXiv 2407.18338 | 从图片识别弯箭头机理（453 张） |
| RxnScribe、MolScribe | JCIM 2023 | 从图片识别反应和分子 |
| EzMechanism | Nat. Methods 2023，10.1038/s41592-023-02006-7 | 从 M-CSA 提取机理规则 |

## 写综述时要回答的问题

- 机理的表示格式：MechSMILES、RMechDB 箭头编码、Lewis-labeled graphs、MolScript，各能表示什么、不能表示什么。
- 核对：形式核对（电子账本）与「合理性」判断（LLM、动力学）的边界在哪里。
- 画图：目前所有工具在自由基箭头、拥挤的多片段场景、多步拼版上的表现。
