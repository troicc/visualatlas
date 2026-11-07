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
  demo.py             示例
```

## 现状

- **已有（第 1 周）**：绘图层，从作者此前的讲义引擎移植而来。移植时与原引擎在同一组场景上逐字比对，SVG 输出和检查结果完全一致。
- **还没有**：化学模型（RDKit 分子图）、核对器、MechSMILES 导入、自动排版。
- 方向和节奏见 [`PLAN.md`](PLAN.md)，设计取舍随做随记在 [`DECISIONS.md`](DECISIONS.md)。

在 2026-12-20 公开之前，仓库保持私有。

## 与相关工作的关系

- [MechSMILES / ChRIMP](https://github.com/schwallergroup/ChRIMP)（MIT）：输入格式和对照组都取自这里。ChRIMP 的 `check_validity()` 目前是占位函数；它的环境只覆盖极性、闭壳层机理。本项目补上的是：机理的可视化与核对，以及对自由基、金属、立体中心的覆盖。
- 排版和弯箭头几何沿用作者此前为讲义写的绘图规则（ACS 1996 参数），只移植通用部分，不含讲义素材。
