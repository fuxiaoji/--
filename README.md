# 高温作业专用服装设计 — 复现

本项目基于附件数据复现 2018 A 题三问的数值结果，核心求解器为 Crank–Nicolson 差分法。

## 文件说明

- `model.py`: 底层 C-N 求解器 + 数据读取 + 简易写 Excel 工具
- `repro_notebook.ipynb`: 按题目分块的 Jupyter 复现
- `run_repro.py`: 最小可运行示例（短时步验证）
- `requirements.txt`: 依赖列表

## 快速开始

建议使用虚拟环境：

```zsh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

运行最小示例：

```zsh
python run_repro.py
```

打开 Notebook：

```zsh
jupyter lab
```

在 `repro_notebook.ipynb` 中依次执行单元，即可得到问题一到三的复现结果，并输出 `problem1.xlsx`。
