# 附录：源代码

本附录包含竞赛论文的全部核心代码，按模块分为三部分：

| 部分 | 文件 | 功能 |
|------|------|------|
| A | `model.py` | 底层 PDE 求解器（Crank-Nicolson 差分 + Thomas 算法） |
| B | `repro_notebook.ipynb` | 原论文三问题复现（参数反演、单变量搜索、二维搜索） |
| C | `enhanced_solution.ipynb` | 增强模型（Cobb-Douglas 效用函数 + 二分法 + SLSQP） |

运行环境：Python 3.11+，依赖 `numpy`, `scipy`, `matplotlib`。

---

## A. 底层 PDE 求解器 — `model.py`

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import math
import numpy as np


@dataclass(frozen=True)
class Layer:
    name: str
    rho: float
    c: float
    k: float
    thickness_mm: float

    @property
    def alpha(self) -> float:
        return self.k / (self.c * self.rho)


@dataclass(frozen=True)
class AppendixData:
    materials: Dict[str, Dict[str, float]]
    thickness_range_mm: Dict[str, Tuple[float, float]]
    measured_time_s: np.ndarray
    measured_temp_c: np.ndarray


def _read_shared_strings(z: ZipFile) -> List[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    strings = []
    for si in root.findall(".//a:si", ns):
        texts = [t.text for t in si.findall(".//a:t", ns) if t.text]
        strings.append("".join(texts))
    return strings


def _iter_sheet_rows(z: ZipFile, sheet_name: str) -> List[List[str | None]]:
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    ns = {
        "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    relmap = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
    shared = _read_shared_strings(z)

    sheet_path = None
    for s in wb.findall(".//a:sheets/a:sheet", ns):
        if s.attrib["name"] == sheet_name:
            rel_id = s.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
            sheet_path = "xl/" + relmap[rel_id]
            break
    if sheet_path is None:
        raise ValueError(f"Sheet {sheet_name} not found")

    root = ET.fromstring(z.read(sheet_path))
    rows = []
    for row in root.findall(".//a:sheetData/a:row", ns):
        vals: List[str | None] = []
        for c in row.findall("a:c", ns):
            t = c.attrib.get("t")
            v = c.find("a:v", ns)
            value: str | None = None if v is None else v.text
            if t == "s" and value is not None:
                value = shared[int(value)]
            vals.append(value)
        rows.append(vals)
    return rows


def load_appendix_xlsx(path: str | Path) -> AppendixData:
    """从竞赛附件 Excel 加载材料参数与实测温度数据。"""
    path = Path(path)
    with ZipFile(path) as z:
        sheet1 = _iter_sheet_rows(z, "附件1")
        sheet2 = _iter_sheet_rows(z, "附件2")

    materials: Dict[str, Dict[str, float]] = {}
    thickness_range: Dict[str, Tuple[float, float]] = {}
    for row in sheet1:
        if not row or row[0] in (None, "分层") or row[1] is None:
            continue
        layer = str(row[0]).strip()
        if layer.endswith("层"):
            name = layer.replace("层", "")
        else:
            name = layer
        try:
            rho = float(row[1])
            c = float(row[2])
            k = float(row[3])
        except (TypeError, ValueError):
            continue
        thickness = str(row[4])
        if "-" in thickness:
            low, high = thickness.split("-")
            thickness_range[name] = (float(low), float(high))
        else:
            thickness_range[name] = (float(thickness), float(thickness))
        materials[name] = {"rho": rho, "c": c, "k": k}

    times: List[float] = []
    temps: List[float] = []
    for row in sheet2:
        if not row or row[0] in (None, "时间 (s)"):
            continue
        try:
            times.append(float(row[0]))
            temps.append(float(row[1]))
        except (TypeError, ValueError):
            continue

    return AppendixData(
        materials=materials,
        thickness_range_mm=thickness_range,
        measured_time_s=np.asarray(times, dtype=float),
        measured_temp_c=np.asarray(temps, dtype=float),
    )


def _build_layers(materials, thickness_mm):
    return [
        Layer("I",  materials["I"]["rho"],  materials["I"]["c"],  materials["I"]["k"],  thickness_mm["I"]),
        Layer("II", materials["II"]["rho"], materials["II"]["c"], materials["II"]["k"], thickness_mm["II"]),
        Layer("III",materials["III"]["rho"],materials["III"]["c"],materials["III"]["k"],thickness_mm["III"]),
        Layer("IV", materials["IV"]["rho"], materials["IV"]["c"], materials["IV"]["k"], thickness_mm["IV"]),
    ]


def _grid_from_layers(layers, dx_mm):
    """划分空间网格，返回各节点的热扩散率 alpha、热导率 k 和层间交界面索引。"""
    dx = dx_mm / 1000.0
    counts = [max(1, int(round(layer.thickness_mm / dx_mm))) for layer in layers]
    total = sum(counts)
    alpha = np.zeros(total, dtype=float)
    k = np.zeros(total, dtype=float)
    interfaces = []
    idx = 0
    for layer, count in zip(layers, counts):
        alpha[idx: idx + count] = layer.alpha
        k[idx: idx + count] = layer.k
        idx += count
        interfaces.append(idx)
    interfaces = interfaces[:-1]  # 去掉最后一个（最外侧不是交界面）
    return alpha, k, interfaces


def _build_matrices(alpha, k, interfaces, dx, dt, h_out, h_in, t_env, t_skin):
    """
    构建 Crank-Nicolson 离散化的三对角矩阵 A 和右端项。
    
    返回 (a, b, c, ba, bb, bc, d)，其中：
      A·T^{n+1} = rhs,  rhs = B·T^n + d
      A: 对角 b, 次对角 a/c;  B: 对角 bb, 次对角 ba/bc.
    
    内部节点：标准 C-N 格式
       -r T_{i-1}^{n+1} + (1+2r) T_i^{n+1} - r T_{i+1}^{n+1}
     =  r T_{i-1}^n     + (1-2r) T_i^n     + r T_{i+1}^n
    其中 r = alpha * dt / (2 * dx^2)
    
    交界面：热流连续条件
       -k_L T_{i-1} + (k_L+k_R) T_i - k_R T_{i+1} = 0
    边界：对流换热条件
       外边界: (k0/dx + h_out) T_0 - (k0/dx) T_1 = h_out * T_env
       内边界: -(kN/dx) T_{N-1} + (kN/dx + h_in) T_N = h_in * T_skin
    """
    n = len(alpha)
    a = np.zeros(n - 1, dtype=float)
    b = np.zeros(n, dtype=float)
    c = np.zeros(n - 1, dtype=float)
    ba = np.zeros(n - 1, dtype=float)
    bb = np.zeros(n, dtype=float)
    bc = np.zeros(n - 1, dtype=float)
    d = np.zeros(n, dtype=float)

    r = alpha * dt / (2.0 * dx * dx)

    for i in range(1, n - 1):
        if i in interfaces:
            k_left = k[i - 1]
            k_right = k[i]
            a[i - 1] = -k_left
            b[i] = k_left + k_right
            c[i] = -k_right
            bb[i] = 0.0
            d[i] = 0.0
        else:
            a[i - 1] = -r[i]
            b[i] = 1 + 2 * r[i]
            c[i] = -r[i]
            ba[i - 1] = r[i]
            bb[i] = 1 - 2 * r[i]
            bc[i] = r[i]

    # Outer boundary (environment)
    k0 = k[0]
    b[0] = k0 / dx + h_out
    c[0] = -k0 / dx
    d[0] = h_out * t_env

    # Inner boundary (skin)
    kN = k[-1]
    a[-1] = -kN / dx
    b[-1] = kN / dx + h_in
    d[-1] = h_in * t_skin

    return a, b, c, ba, bb, bc, d


def _thomas_prepare(a, b, c):
    """Thomas 算法预处理：LU 分解（只做一次）。"""
    n = len(b)
    c_prime = np.zeros(n - 1, dtype=float)
    denom = np.zeros(n, dtype=float)

    denom[0] = b[0]
    c_prime[0] = c[0] / denom[0]
    for i in range(1, n - 1):
        denom[i] = b[i] - a[i - 1] * c_prime[i - 1]
        c_prime[i] = c[i] / denom[i]
    denom[-1] = b[-1] - a[-1] * c_prime[-1]
    return c_prime, denom


def _thomas_solve(a, c_prime, denom, d):
    """Thomas 算法每步回代：O(n) 求解三对角方程组。"""
    n = len(d)
    d_prime = np.zeros(n, dtype=float)
    d_prime[0] = d[0] / denom[0]
    for i in range(1, n):
        d_prime[i] = (d[i] - a[i - 1] * d_prime[i - 1]) / denom[i]
    x = np.zeros(n, dtype=float)
    x[-1] = d_prime[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return x


def solve_heat_cn(
    materials, thickness_mm, t_env, h_out, h_in,
    t_skin=37.0, total_time_s=5400, dx_mm=0.1, dt_s=1.0, store_all=True,
):
    """
    一维四层热传导 Crank-Nicolson 求解器。
    
    参数
    ----
    materials : dict       — 各层材料参数 {rho, c, k}
    thickness_mm : dict    — 各层厚度 (mm)
    t_env : float          — 环境温度 (°C)
    h_out : float          — 外边界对流换热系数 (W/(m²·K))
    h_in : float           — 内边界对流换热系数 (W/(m²·K))
    t_skin : float         — 皮肤内侧温度 (°C)，默认 37
    total_time_s : int     — 总模拟时间 (s)
    dx_mm : float          — 空间步长 (mm)
    dt_s : float           — 时间步长 (s)
    store_all : bool       — True=返回全空间温度场, False=仅皮肤侧
    
    返回
    ----
    times : np.ndarray     — 时间序列
    temps : np.ndarray     — 温度矩阵 (steps+1, n) 或 (steps+1, 1)
    """
    layers = _build_layers(materials, thickness_mm)
    alpha, k, interfaces = _grid_from_layers(layers, dx_mm)
    dx = dx_mm / 1000.0
    n = len(alpha)

    a, b, c, ba, bb, bc, d = _build_matrices(
        alpha=alpha, k=k, interfaces=interfaces, dx=dx, dt=dt_s,
        h_out=h_out, h_in=h_in, t_env=t_env, t_skin=t_skin,
    )

    c_prime, denom = _thomas_prepare(a, b, c)

    steps = int(total_time_s / dt_s)
    times = np.arange(0, steps + 1) * dt_s
    if store_all:
        temps = np.zeros((steps + 1, n), dtype=float)
    else:
        temps = np.zeros((steps + 1, 1), dtype=float)

    temp_prev = np.full(n, t_skin, dtype=float)
    if store_all:
        temps[0, :] = temp_prev
    else:
        temps[0, 0] = temp_prev[-1]

    for step in range(1, steps + 1):
        rhs = bb * temp_prev
        rhs[1:] += ba * temp_prev[:-1]
        rhs[:-1] += bc * temp_prev[1:]
        rhs += d
        temp_next = _thomas_solve(a, c_prime, denom, rhs)
        temp_prev = temp_next
        if store_all:
            temps[step, :] = temp_next
        else:
            temps[step, 0] = temp_next[-1]

    return times, temps


def steady_state_h_in(materials, thickness_mm, t_env, t_skin, t_surface, h_out):
    """
    稳态热平衡求内边界换热系数 h_in。
    
    利用已知稳态皮肤温度 T_surface，通过串联热阻计算热流 q，
    再由 q = h_in * (T_surface - T_skin) 反推 h_in。
    
    这一步将问题一的二维搜索 (h_out, h_in) 降为一维 h_out。
    """
    resistance = 1.0 / h_out
    for layer_name in ["I", "II", "III", "IV"]:
        k = materials[layer_name]["k"]
        d_m = thickness_mm[layer_name] / 1000.0
        resistance += d_m / k
    q = (t_env - t_surface) / resistance
    return q / (t_surface - t_skin)


def rmse(a, b):
    """均方根误差。"""
    return float(math.sqrt(np.mean((a - b) ** 2)))


def write_simple_xlsx(path, header, data):
    """导出不含第三方库的 .xlsx 文件（直接写 OOXML）。"""
    path = Path(path)
    from io import BytesIO

    def col_letter(idx):
        letters = ""
        while idx >= 0:
            letters = chr(idx % 26 + 65) + letters
            idx = idx // 26 - 1
        return letters

    def cell_ref(r, c):
        return f"{col_letter(c)}{r}"

    rows_xml = []
    r = 1
    header_cells = []
    for c, text in enumerate(header):
        header_cells.append(
            f'<c r="{cell_ref(r,c)}" t="inlineStr"><is><t>{text}</t></is></c>'
        )
    rows_xml.append(f'<row r="{r}">{"".join(header_cells)}</row>')

    for i in range(data.shape[0]):
        r += 1
        cells = []
        for c in range(data.shape[1]):
            cells.append(f'<c r="{cell_ref(r,c)}"><v>{data[i,c]}</v></c>')
        rows_xml.append(f'<row r="{r}">{"".join(cells)}</row>')

    sheet_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<sheetData>{"".join(rows_xml)}</sheetData>'
        '</worksheet>'
    )

    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets>'
        '</workbook>'
    )

    workbook_rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        '</Relationships>'
    )

    root_rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        '</Relationships>'
    )

    content_types = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '</Types>'
    )

    buffer = BytesIO()
    with ZipFile(buffer, "w") as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", workbook_xml)
        z.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        z.writestr("xl/worksheets/sheet1.xml", sheet_xml)

    path.write_bytes(buffer.getvalue())
```

---

## B. 原论文复现 — `repro_notebook.ipynb`

```python
# ===== Cell 1: 导入 =====
from pathlib import Path
import numpy as np
from model import (
    load_appendix_xlsx, solve_heat_cn, steady_state_h_in, rmse, write_simple_xlsx,
)

ROOT = Path('.').resolve()
APPENDIX = ROOT / 'CUMCM-2018-Problem-A-Chinese-Appendix.xlsx'

# ===== Cell 2: 加载数据 =====
data = load_appendix_xlsx(APPENDIX)
data.materials, data.thickness_range_mm, data.measured_time_s[:5], data.measured_temp_c[-5:]

# ===== Cell 3: 速度优化开关 =====
FAST_MODE = True
DX_MM = 0.2 if FAST_MODE else 0.1
DT_S = 2.0 if FAST_MODE else 1.0
H_OUT_STEP = 0.5 if FAST_MODE else 0.2
P2_COARSE_STEP = 2.0 if FAST_MODE else 1.0
P2_FINE_STEP = 0.2 if FAST_MODE else 0.1
P3_D2_STEP = 0.2 if FAST_MODE else 0.1
P3_D4_STEP = 0.2 if FAST_MODE else 0.1

# ===== Cell 5: 问题一 — 参数反演 =====
materials = data.materials
thickness = {'I': 0.6, 'II': 6.0, 'III': 3.6, 'IV': 5.0}
t_env = 75.0
t_skin = 37.0
t_surface = float(data.measured_temp_c[-1])

def estimate_h_in(h_out):
    return steady_state_h_in(materials, thickness, t_env, t_skin, t_surface, h_out)

def fit_h_out(h_out_values):
    best = {'h_out': None, 'h_in': None, 'rmse': float('inf')}
    for h_out in h_out_values:
        h_in = estimate_h_in(h_out)
        times, temps = solve_heat_cn(
            materials, thickness, t_env, h_out, h_in, t_skin,
            total_time_s=int(data.measured_time_s[-1]),
            dx_mm=DX_MM, dt_s=DT_S, store_all=False,
        )
        sim_skin = temps[:, 0]
        target = np.interp(times, data.measured_time_s, data.measured_temp_c)
        score = rmse(sim_skin, target)
        if score < best['rmse']:
            best = {'h_out': h_out, 'h_in': h_in, 'rmse': score}
    return best

h_out_grid = np.arange(100.0, 120.01, H_OUT_STEP)
best = fit_h_out(h_out_grid)

# ===== Cell 6: 导出全温度场 =====
h_out = float(best['h_out'])
h_in = float(best['h_in'])
times, temps = solve_heat_cn(
    materials, thickness, t_env, h_out, h_in, t_skin,
    total_time_s=int(data.measured_time_s[-1]),
    dx_mm=DX_MM, dt_s=DT_S, store_all=True,
)

dx_mm = DX_MM
counts = [
    int(round(thickness['I'] / dx_mm)),
    int(round(thickness['II'] / dx_mm)),
    int(round(thickness['III'] / dx_mm)),
    int(round(thickness['IV'] / dx_mm)),
]
n = sum(counts)
x_mm = np.arange(n) * dx_mm
header = ['time_s'] + [f'x_{x:.2f}_mm' for x in x_mm]
write_simple_xlsx(ROOT / 'problem1.xlsx', header, np.column_stack([times, temps]))

# ===== Cell 8: 问题二 — 单变量定步长搜索 =====
t_env = 65.0
thickness_p2 = {'I': 0.6, 'II': 6.0, 'III': 3.6, 'IV': 5.5}
total_time_s = 3600

# 关键：使用问题一反演得到的固定换热系数，不再重新计算 h_in
h_out = float(best['h_out'])
h_in = float(best['h_in'])

def check_constraints(d2_mm):
    thickness_p2['II'] = d2_mm
    _, temps = solve_heat_cn(
        materials, thickness_p2, t_env, h_out, h_in, t_skin,
        total_time_s=total_time_s, dx_mm=DX_MM, dt_s=DT_S, store_all=False,
    )
    skin = temps[:, 0]
    max_temp = float(np.max(skin))
    over_44 = float(np.sum(skin > 44.0))
    over_time = over_44 * DT_S
    ok = (max_temp <= 47.0) and (over_time <= 300.0)
    return {'d2': d2_mm, 'ok': ok, 'max': max_temp, 'over_44_s': over_time}

# 粗搜
coarse = np.arange(10.0, 25.01, P2_COARSE_STEP)
coarse_results = [check_constraints(d2) for d2 in coarse]
feasible_coarse = [r for r in coarse_results if r['ok']]
print('feasible coarse count:', len(feasible_coarse))

# 精搜
fine = np.arange(17.0, 19.01, P2_FINE_STEP)
fine_results = [check_constraints(d2) for d2 in fine]
feasible = [r for r in fine_results if r['ok']]
if feasible:
    best_p2 = min(feasible, key=lambda r: r['d2'])
    print('best_p2:', best_p2)

# ===== Cell 11: 问题三 — 二维区域搜索 =====
t_env = 80.0
total_time_s = 1800

h_out = float(best['h_out'])
h_in = float(best['h_in'])

def settle_time(skin, dt=DT_S, tol=0.05, window=300):
    target = float(skin[-1])
    for i in range(len(skin) - window):
        segment = skin[i:i + window]
        if np.all(np.abs(segment - target) <= tol):
            return i * dt
    return float(len(skin) * dt)

def check_pair(d2_mm, d4_mm):
    thickness_p3 = {'I': 0.6, 'II': d2_mm, 'III': 3.6, 'IV': d4_mm}
    _, temps = solve_heat_cn(
        materials, thickness_p3, t_env, h_out, h_in, t_skin,
        total_time_s=total_time_s, dx_mm=DX_MM, dt_s=DT_S, store_all=False,
    )
    skin = temps[:, 0]
    max_temp = float(np.max(skin))
    over_44 = float(np.sum(skin > 44.0))
    over_time = over_44 * DT_S
    ok = (max_temp <= 47.0) and (over_time <= 300.0)
    return {
        'd2': d2_mm, 'd4': d4_mm, 'ok': ok,
        'max': max_temp, 'over_44_s': over_time, 'settle_s': settle_time(skin),
    }

d2_grid = np.arange(18.0, 21.01, P3_D2_STEP)
d4_grid = np.arange(0.6, 6.41, P3_D4_STEP)
feasible = []
for d2 in d2_grid:
    for d4 in d4_grid:
        r = check_pair(d2, d4)
        if r['ok']:
            feasible.append(r)

best_p3 = min(feasible, key=lambda r: (r['d2'], r['settle_s']))
```

---

## C. 增强模型 — `enhanced_solution.ipynb`

```python
# ===== Cell 1: 导入 =====
import sys; sys.path.insert(0, '..')
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from model import load_appendix_xlsx, solve_heat_cn, steady_state_h_in, rmse

plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

ROOT = Path('..').resolve()
APPENDIX = ROOT / 'CUMCM-2018-Problem-A-Chinese-Appendix.xlsx'
DX_MM = 0.2
DT_S = 2.0

# ===== Cell 2: 问题一反演 =====
data = load_appendix_xlsx(APPENDIX)
materials = data.materials
thickness_p1 = {'I': 0.6, 'II': 6.0, 'III': 3.6, 'IV': 5.0}
t_skin = 37.0
t_surface = float(data.measured_temp_c[-1])

best_h = {'h_out': None, 'h_in': None, 'rmse': float('inf')}
for h_out in np.arange(100.0, 120.01, 0.5):
    h_in = steady_state_h_in(materials, thickness_p1, 75.0, t_skin, t_surface, h_out)
    times_p1, temps_p1 = solve_heat_cn(
        materials, thickness_p1, 75.0, h_out, h_in, t_skin,
        total_time_s=5400, dx_mm=DX_MM, dt_s=DT_S, store_all=False,
    )
    target = np.interp(times_p1, data.measured_time_s, data.measured_temp_c)
    score = rmse(temps_p1[:, 0], target)
    if score < best_h['rmse']:
        best_h = {'h_out': h_out, 'h_in': h_in, 'rmse': score}

H_OUT = float(best_h['h_out'])
H_IN = float(best_h['h_in'])
print(f'Problem 1: h_I={H_OUT:.1f}, h_IV={H_IN:.4f}, RMSE={best_h["rmse"]:.4f}')

# ===== Cell 4: 问题二 PDE 封装 =====
T_ENV_P2 = 65.0
TOTAL_TIME_P2 = 3600
D2_MIN, D2_MAX = 0.6, 25.0
K2 = materials['II']['k']
gamma_min = D2_MIN / K2
gamma_max = D2_MAX / K2

def run_skin_temp(d2, t_env=T_ENV_P2, total_time_s=TOTAL_TIME_P2):
    thickness = {'I': 0.6, 'II': d2, 'III': 3.6, 'IV': 5.5}
    _, temps = solve_heat_cn(
        materials, thickness, t_env, H_OUT, H_IN, t_skin,
        total_time_s=total_time_s, dx_mm=DX_MM, dt_s=DT_S, store_all=False,
    )
    return temps[:, 0]

def max_temp(d2):
    return float(np.max(run_skin_temp(d2)))

def temp_at_time(d2, target_time_s):
    skin = run_skin_temp(d2)
    idx = min(int(target_time_s / DT_S), len(skin) - 1)
    return float(skin[idx])

# ===== Cell 5: 二分法求临界厚度 =====
def bisection_critical(f_target, target_val, a, b, tol=0.01, max_iter=50):
    """二分法求解 f(d) = target_val。f 关于 d 单调递减。"""
    fa = f_target(a)
    fb = f_target(b)
    if not (fb <= target_val <= fa):
        return None  # 约束恒满足或恒不满足
    for _ in range(max_iter):
        m = (a + b) / 2.0
        fm = f_target(m)
        if abs(fm - target_val) < tol:
            return m
        if fm > target_val: a = m
        else: b = m
    return (a + b) / 2.0

D1 = bisection_critical(max_temp, 47.0, D2_MIN, D2_MAX)
D2 = bisection_critical(lambda d: temp_at_time(d, 55*60), 44.0, D2_MIN, D2_MAX)
d2_critical = max(D1 if D1 is not None else D2_MIN,
                  D2 if D2 is not None else D2_MIN)
print(f'd2_critical = {d2_critical:.2f} mm')

# ===== Cell 6: Cobb-Douglas 效用函数 & 优化 =====
ALPHA = 0.5

def comfort(d2):
    return 1.0 - np.sqrt((d2 - D2_MIN) / (D2_MAX - D2_MIN))

def insulation(d2):
    gamma = d2 / K2
    return (gamma - gamma_min) / (gamma_max - gamma_min)

def utility(d2):
    c = comfort(d2)
    r = insulation(d2)
    return (c ** ALPHA) * (r ** (1.0 - ALPHA))

# 安全域内最大化 U
d2_search = np.arange(d2_critical, D2_MAX + 0.01, 0.01)
u_values = utility(d2_search)
best_idx = np.argmax(u_values)
d2_opt = d2_search[best_idx]
u_opt = u_values[best_idx]
C_opt = comfort(d2_opt)
R_opt = insulation(d2_opt)

print(f'd2* = {d2_opt:.2f} mm, C = {C_opt:.4f}, R = {R_opt:.4f}, U = {u_opt:.4f}')

# ===== Cell 8: 问题二可视化 (四象限) =====
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
d2_range = np.linspace(D2_MIN, D2_MAX, 500)
d2_feasible = np.linspace(d2_critical, D2_MAX, 200)

# C(d2)
ax = axes[0, 0]
ax.plot(d2_range, comfort(d2_range), 'b-', linewidth=1.5, alpha=0.3,
        label='Infeasible')
ax.plot(d2_feasible, comfort(d2_feasible), 'b-', linewidth=2.5,
        label='Feasible')
ax.axvline(d2_critical, color='red', linestyle='--', linewidth=2,
           label=f'Critical d2 = {d2_critical:.1f} mm')
ax.set_xlabel('d2 (mm)'); ax.set_ylabel('C')
ax.set_title('Comfort Index C(d2)')
ax.legend(); ax.set_ylim(0, 1.05); ax.grid(True, alpha=0.3)

# R(d2)
ax = axes[0, 1]
ax.plot(d2_range, insulation(d2_range), 'r-', linewidth=1.5, alpha=0.3)
ax.plot(d2_feasible, insulation(d2_feasible), 'r-', linewidth=2.5)
ax.axvline(d2_critical, color='red', linestyle='--', linewidth=2)
ax.set_xlabel('d2 (mm)'); ax.set_ylabel('R')
ax.set_title('Insulation Index R(d2)')
ax.set_ylim(0, 1.05); ax.grid(True, alpha=0.3)

# U(d2) — 区分全局最大与约束最优
ax = axes[1, 0]
ax.plot(d2_range, utility(d2_range), 'purple', linewidth=1.5, alpha=0.3)
ax.plot(d2_feasible, utility(d2_feasible), 'purple', linewidth=2.5)
ax.axvline(d2_critical, color='red', linestyle='--', linewidth=2)
# 全局最大（不可行）
d2_global = d2_range[np.argmax(utility(d2_range))]
u_global = np.max(utility(d2_range))
ax.plot(d2_global, u_global, 'o', color='gray', markersize=12, alpha=0.6,
        label=f'Global max (infeasible): d2={d2_global:.1f}mm')
# 约束最优
ax.plot(d2_opt, u_opt, 'go', markersize=12,
        label=f'Constrained max: d2={d2_opt:.1f}mm, U={u_opt:.3f}')
ax.axvspan(0, d2_critical, alpha=0.08, color='red')
ax.set_xlabel('d2 (mm)'); ax.set_ylabel('U')
ax.set_title(f'Utility U = C^{ALPHA} * R^(1-{ALPHA})')
ax.legend(fontsize=8); ax.set_ylim(0, 0.45); ax.grid(True, alpha=0.3)

# C-R Trade-off
ax = axes[1, 1]
c_feas = comfort(d2_feasible)
r_feas = insulation(d2_feasible)
ax.scatter(comfort(d2_full), insulation(d2_full),
           c=utility(d2_full), cmap='plasma', s=3, alpha=0.15)
ax.scatter(c_feas, r_feas, c=utility(d2_feasible),
           cmap='plasma', s=20, alpha=0.9)
ax.plot(C_opt, R_opt, 'go', markersize=14,
        label=f'Optimal (d2={d2_opt:.1f})')
ax.plot(comfort(d2_global), insulation(d2_global), 'o',
        color='gray', markersize=12, alpha=0.5,
        label=f'Global (infeasible, d2={d2_global:.1f})')
ax.set_xlabel('Comfort C'); ax.set_ylabel('Insulation R')
ax.set_title('C-R Trade-off (highlight = feasible)')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.suptitle('Enhanced Model: Problem 2 Analysis', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('problem2_analysis.png', dpi=150, bbox_inches='tight')

# ===== Cell 9: α 敏感性分析 =====
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
alphas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(alphas)))
optimal_d2 = []
for alpha, color in zip(alphas, colors):
    u_vals = np.array([
        (comfort(d)**alpha) * (insulation(d)**(1-alpha))
        for d in d2_range
    ])
    d2_best = d2_range[np.argmax(u_vals)]
    optimal_d2.append((alpha, d2_best))
    axes[0].plot(d2_range, u_vals, color=color, linewidth=2,
                 label=f'alpha={alpha}')
axes[0].set_xlabel('d2 (mm)'); axes[0].set_ylabel('U')
axes[0].set_title('Utility for Different alpha')
axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.3)

axes[1].plot(d2_range, comfort(d2_range), 'b-', linewidth=2, label='C(d2)')
axes[1].plot(d2_range, insulation(d2_range), 'r-', linewidth=2, label='R(d2)')
axes[1].axvline(d2_critical, color='gray', linestyle='--')
axes[1].set_xlabel('d2 (mm)'); axes[1].set_ylabel('Index')
axes[1].set_title('C(d2) and R(d2)')
axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.3)

alphas_list = [a for a, _ in optimal_d2]
d2_list = [d for _, d in optimal_d2]
axes[2].bar(range(len(alphas_list)), d2_list,
            tick_label=[f'{a:.1f}' for a in alphas_list], color=colors)
axes[2].axhline(d2_critical, color='red', linestyle='--', alpha=0.5)
axes[2].set_xlabel('alpha'); axes[2].set_ylabel('Optimal d2 (mm)')
axes[2].set_title('Optimal d2 vs alpha')

plt.suptitle('Sensitivity Analysis: alpha', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('alpha_sensitivity.png', dpi=150, bbox_inches='tight')

# ===== Cell 11: 问题三参数 =====
T_ENV_P3 = 80.0
TOTAL_TIME_P3 = 1800
D4_MIN, D4_MAX = 0.6, 6.4
K4 = materials['IV']['k']

gamma2_min = D2_MIN / K2; gamma2_max = D2_MAX / K2
gamma4_min = D4_MIN / K4; gamma4_max = D4_MAX / K4
gamma_sum_min = gamma2_min + gamma4_min
gamma_sum_max = gamma2_max + gamma4_max
d_sum_min = D2_MIN + D4_MIN
d_sum_max = D2_MAX + D4_MAX

def run_skin_temp_p3(d2, d4):
    thickness = {'I': 0.6, 'II': d2, 'III': 3.6, 'IV': d4}
    _, temps = solve_heat_cn(
        materials, thickness, T_ENV_P3, H_OUT, H_IN, t_skin,
        total_time_s=TOTAL_TIME_P3, dx_mm=DX_MM, dt_s=DT_S, store_all=False,
    )
    return temps[:, 0]

def comfort_p3(d2, d4):
    return 1.0 - np.sqrt((d2 + d4 - d_sum_min) / (d_sum_max - d_sum_min))

def insulation_p3(d2, d4):
    gamma_sum = d2 / K2 + d4 / K4
    return (gamma_sum - gamma_sum_min) / (gamma_sum_max - gamma_sum_min)

def utility_p3(d2, d4, alpha=ALPHA):
    c = comfort_p3(d2, d4)
    r = insulation_p3(d2, d4)
    if c <= 0 or r <= 0:
        return 0.0
    return (c ** alpha) * (r ** (1.0 - alpha))

# ===== Cell 12: 网格遍历可行域 =====
d2_grid = np.arange(18.0, 21.01, 0.2)
d4_grid = np.arange(0.6, 6.41, 0.2)
feasible_p3 = []
for d2 in d2_grid:
    for d4 in d4_grid:
        skin = run_skin_temp_p3(d2, d4)
        max_t = float(np.max(skin))
        over_44 = float(np.sum(skin > 44.0)) * DT_S
        if max_t <= 47.0 and over_44 <= 300.0:
            feasible_p3.append({
                'd2': d2, 'd4': d4, 'max_T': max_t, 'over_44_s': over_44,
                'C': comfort_p3(d2, d4), 'R': insulation_p3(d2, d4),
                'U': utility_p3(d2, d4),
            })

best_by_u = max(feasible_p3, key=lambda r: r['U'])
print(f'Max U: d2={best_by_u["d2"]:.1f}, d4={best_by_u["d4"]:.1f}')
print(f'  C={best_by_u["C"]:.4f}, R={best_by_u["R"]:.4f}, U={best_by_u["U"]:.4f}')

# ===== Cell 13: SLSQP 精化 =====
def objective(x):
    return -utility_p3(x[0], x[1])

def constraint_max_temp(x):
    skin = run_skin_temp_p3(x[0], x[1])
    return 47.0 - float(np.max(skin))

def constraint_over_44(x):
    skin = run_skin_temp_p3(x[0], x[1])
    return 300.0 - float(np.sum(skin > 44.0)) * DT_S

constraints = [
    {'type': 'ineq', 'fun': constraint_max_temp},
    {'type': 'ineq', 'fun': constraint_over_44},
]
bounds = [(D2_MIN, D2_MAX), (D4_MIN, D4_MAX)]

result = minimize(
    objective, x0=np.array([best_by_u['d2'], best_by_u['d4']]),
    method='SLSQP', bounds=bounds, constraints=constraints,
    options={'ftol': 1e-6, 'maxiter': 50},
)
if result.success:
    d2_slsqp, d4_slsqp = result.x
    print(f'SLSQP: d2={d2_slsqp:.2f}, d4={d4_slsqp:.2f}')

# ===== Cell 15: 问题三可视化 =====
fig, axes = plt.subplots(1, 2, figsize=(15, 6))

U_matrix = np.zeros((len(d2_grid), len(d4_grid)))
for i, d2 in enumerate(d2_grid):
    for j, d4 in enumerate(d4_grid):
        U_matrix[i, j] = utility_p3(d2, d4)

D4_mesh, D2_mesh = np.meshgrid(d4_grid, d2_grid)

ax = axes[0]
contour = ax.contourf(D2_mesh, D4_mesh, U_matrix, levels=20,
                       cmap='viridis', alpha=0.9)
ax.scatter([r['d2'] for r in feasible_p3],
           [r['d4'] for r in feasible_p3],
           c='red', s=8, alpha=0.6, label='Feasible')
ax.plot(best_by_u['d2'], best_by_u['d4'], 'r*', markersize=18,
        label=f'Max U ({best_by_u["d2"]:.1f}, {best_by_u["d4"]:.1f})')
ax.set_xlabel('d2 (mm)'); ax.set_ylabel('d4 (mm)')
ax.set_title('Utility U(d2, d4) with Feasible Region')
ax.legend(fontsize=9); plt.colorbar(contour, ax=ax, label='U')

ax = axes[1]
c_vals = np.array([r['C'] for r in feasible_p3])
r_vals = np.array([r['R'] for r in feasible_p3])
u_vals = np.array([r['U'] for r in feasible_p3])
sc = ax.scatter(c_vals, r_vals, c=u_vals, cmap='plasma',
                s=40, alpha=0.8, edgecolors='none')
ax.plot(best_by_u['C'], best_by_u['R'], 'r*', markersize=18,
        label=f'Max U')
ax.set_xlabel('Comfort C'); ax.set_ylabel('Insulation R')
ax.set_title('C-R Trade-off (Pareto Frontier)')
ax.legend(fontsize=9); plt.colorbar(sc, ax=ax, label='U')

plt.suptitle('Enhanced Model: Problem 3 Analysis', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('problem3_analysis.png', dpi=150, bbox_inches='tight')

# ===== Cell 16: 模型对比柱状图 =====
fig, ax = plt.subplots(figsize=(10, 6))
models = ['Original (min d2)', 'Enhanced (max U)']
p2_vals = [17.8, d2_opt]
p3_d2_vals = [19.3, d2_slsqp]
p3_d4_vals = [6.4, d4_slsqp]
x = np.arange(len(models)); width = 0.25
ax.bar(x - width, p2_vals, width, label='P2: d2 (mm)', color='steelblue')
ax.bar(x, p3_d2_vals, width, label='P3: d2 (mm)', color='coral')
ax.bar(x + width, p3_d4_vals, width, label='P3: d4 (mm)', color='seagreen')
# ... (标注数值、图例)
ax.set_xticks(x); ax.set_xticklabels(models)
ax.set_ylabel('Thickness (mm)')
ax.set_title('Model Comparison: Original vs Enhanced')
ax.legend(); ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig('model_comparison.png', dpi=150, bbox_inches='tight')

# ===== Cell 17: 温度演化验证 =====
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
# Problem 2
for d2 in np.arange(2, 26, 4):
    skin = run_skin_temp(d2)
    axes[0].plot(np.arange(len(skin))*DT_S/60, skin, alpha=0.7,
                 label=f'd2={d2:.0f}')
axes[0].axhline(47, color='red', linestyle='--', label='47 C')
axes[0].axhline(44, color='orange', linestyle='--', label='44 C')
axes[0].axvline(55, color='orange', linestyle=':', label='t=55min')
axes[0].set_title('Problem 2: T(t) for various d2')
axes[0].legend(fontsize=7, ncol=2); axes[0].grid(True, alpha=0.3)
# Problem 3
for d2 in np.arange(19, 23, 1):
    skin = run_skin_temp_p3(d2, 6.4)
    axes[1].plot(np.arange(len(skin))*DT_S/60, skin, alpha=0.7,
                 label=f'd2={d2:.1f}, d4=6.4')
axes[1].axhline(47, color='red', linestyle='--')
axes[1].axhline(44, color='orange', linestyle='--')
axes[1].axvline(25, color='orange', linestyle=':', label='t=25min')
axes[1].set_title('Problem 3: T(t) for various d2, d4=6.4')
axes[1].legend(fontsize=7); axes[1].grid(True, alpha=0.3)
plt.suptitle('Temperature Evolution Verification', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('temperature_evolution.png', dpi=150, bbox_inches='tight')
```

---

## 依赖环境

```
numpy>=1.26
scipy>=1.13
matplotlib>=3.9
```

所有代码纯 Python 实现，不依赖 MATLAB 或商业求解器。PDE 求解器用自写 Thomas 算法（追赶法），不依赖 `scipy.sparse.linalg`。

---

## 代码量统计

| 文件 | 行数 | 主要内容 |
|------|------|----------|
| `model.py` | 397 | PDE 求解器、数据加载、稳态方程、xlsx 导出 |
| `repro_notebook.ipynb` | ~120 (code) | 三问题网格搜索 |
| `enhanced_solution.ipynb` | ~200 (code) | 效用函数、二分法、SLSQP、5 张可视化图 |
| **合计** | **~720** | |
