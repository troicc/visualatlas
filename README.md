# visualatlas

把箭头推动机理（MechSMILES 或箭头表）自动画成化学家能快速核对的图，并自动查出其中形式上不成立的步骤。

## 两项能力

| 能力 | 输入 | 输出 | 衡量方式 |
| --- | --- | --- | --- |
| 排版 | 一条机理（单步或多步） | SVG：结构式、弯箭头、多格拼版，没有几何违规 | 用一个与引擎无关的检测器在最终 SVG 上数违规 |
| 核对 | 同一条机理 | 每一步的判定：电子记账、原子与电荷守恒、下一格是否与箭头推出的结构一致、立体中心 | 公开数据集上的审计结果，错误注入后的召回率与误报率 |

覆盖范围：极性、自由基（鱼钩箭头）、立体中心做完整核对；过渡金属只核对原子与电荷守恒，外加价电子数上限，不做完整的 18 电子判定。

## 安装与试用

```bash
pip install -e ".[dev]"
visualatlas-demo demo.svg      # 手写的一步：HO⁻ 加成到丙酮羰基，检查后输出 SVG
pytest -q
```

需要一个与 Arial 字宽一致的字体：macOS 和 Windows 自带的 Arial 即可；Linux 上装 Liberation Sans（`apt install fonts-liberation`）。也可以用环境变量 `VISUALATLAS_FONT` 指定 .ttf 文件。

## 目录

```
src/visualatlas/
  fonts.py            字宽（标签排版用）
  draw/mol.py         绘图模型：原子坐标与标签、键的画法、渲染、几何检查、二维松弛
  draw/geomcheck.py   通用几何检查：键长、键角、环内角、楔形模板、共线
  draw/arrows.py      弯箭头：圆弧，在硬检查下穷举落点与拱高
  draw/canvas.py      画布：放分子、孤对、单电子、弯箭头、反应箭头；汇总全部违规
  stereo.py           从画出的楔形反推立体中心的 CIP 构型
  mechsmiles.py       MechSMILES 解析：带原子映射的分子图和箭头列表
  moves.py            统一的电子移动表示（1 个或 2 个电子）
  rmechdb.py          RMechDB 读取
  check/              核对器：电子账本
scripts/              审计脚本
docs/                 精读指南、审计笔记、相关工作笔记
baselines/            对照组脚本（ChRIMP 自带可视化）
tests/data/           测试数据及其来源和许可证说明
  demo.py             示例
```

## 现状

- **已有（第 1 周）**：绘图层，从作者此前的讲义引擎移植而来。移植时与原引擎在同一组场景上逐字比对，SVG 输出和检查结果完全一致。
- **已有**：
  - MechSMILES 与 RMechDB 的读取；
  - 核对器：电子账本，同时处理极性箭头和鱼钩箭头，并比对立体中心。【标注】这部分由 AI 编写，作者待精读，见 [`docs/STUDY_checker.md`](docs/STUDY_checker.md)。
  - RMechDB 全量审计：5,426 行中 19 行被标出，逐条人工看过，都是数据本身的错误，见 [`docs/audit_notes.md`](docs/audit_notes.md)。
- **还没有**：自动排版。
- 方向和节奏见 [`PLAN.md`](PLAN.md)，设计取舍随做随记在 [`DECISIONS.md`](DECISIONS.md)。

在 2026-12-20 公开之前，仓库保持私有。

## 与相关工作的关系

- [MechSMILES / ChRIMP](https://github.com/schwallergroup/ChRIMP)（MIT）：输入格式和对照组都取自这里。
  - ChRIMP 能做的检查：逐原子电子计数，以及 `MoleculeSet.check_step_legality`，后者判断两个状态之间是否存在一种合法的单步箭头组合。`MechSmiles.check_validity()` 仍是占位函数。
  - 本项目核对的是另一件事：画出来的这几根箭头，推出来的是不是画出来的下一格；不是的话，错在哪一根箭头、哪个原子。
  - 覆盖范围：ChRIMP 的代码对自由基和金属（到第 5 周期）已有初步支持，只是论文里没写。本项目在这一点上的差别是逐步核对立体中心，以及画法。ChRIMP 生成中间体时会丢掉立体标记，见 `tests/data/README.md`。
- 其他相关工作（Lewis-labeled graphs、MolScript 等）记在 [`docs/related_work.md`](docs/related_work.md)。
- 排版和弯箭头几何沿用作者此前为讲义写的绘图规则（ACS 1996 参数），只移植通用部分，不含讲义素材。
