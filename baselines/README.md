# 对照组

## ChRIMP 自带可视化（`chrimp_render.py`）

ChRIMP 装在它自己的环境里，不和本项目的环境混用。它的依赖写在 `requirements.txt` 里，没有写进包的元数据，所以 `pip install -e .` 不会把依赖一起装上。画图只需要下面几个：

```bash
git clone https://github.com/schwallergroup/ChRIMP.git && cd ChRIMP
uv venv --python 3.11 && uv pip install -e . rdkit cairosvg svgutils colorama pandas matplotlib pillow
```

在 macOS 上运行要注意两件事：
- cairosvg 需要系统的 cairo 库，要用环境变量指到 Homebrew 的库目录。
- ChRIMP 的可视化会调用 `plt.show()`。不设 `MPLBACKEND=Agg` 的话，它会弹出窗口并一直卡住。

```bash
MPLBACKEND=Agg DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib \
  ChRIMP/.venv/bin/python baselines/chrimp_render.py ChRIMP out/chrimp 0 3 5
```

**输出：**
- `rxn<编号>.svg`：整条机理画成一张图，每格是 RDKit 单独画好的图，横着排成一行。
- `rxn<编号>_step<k>.svg`：每一步单独一张图。
