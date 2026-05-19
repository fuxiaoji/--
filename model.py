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


def _build_layers(materials: Dict[str, Dict[str, float]], thickness_mm: Dict[str, float]) -> List[Layer]:
    return [
        Layer("I", materials["I"]["rho"], materials["I"]["c"], materials["I"]["k"], thickness_mm["I"]),
        Layer("II", materials["II"]["rho"], materials["II"]["c"], materials["II"]["k"], thickness_mm["II"]),
        Layer("III", materials["III"]["rho"], materials["III"]["c"], materials["III"]["k"], thickness_mm["III"]),
        Layer("IV", materials["IV"]["rho"], materials["IV"]["c"], materials["IV"]["k"], thickness_mm["IV"]),
    ]


def _grid_from_layers(layers: List[Layer], dx_mm: float) -> Tuple[np.ndarray, np.ndarray, List[int]]:
    dx = dx_mm / 1000.0
    counts = [max(1, int(round(layer.thickness_mm / dx_mm))) for layer in layers]
    total = sum(counts)
    alpha = np.zeros(total, dtype=float)
    k = np.zeros(total, dtype=float)
    interfaces = []
    idx = 0
    for layer, count in zip(layers, counts):
        alpha[idx : idx + count] = layer.alpha
        k[idx : idx + count] = layer.k
        idx += count
        interfaces.append(idx)
    interfaces = interfaces[:-1]
    return alpha, k, interfaces


def _build_matrices(
    alpha: np.ndarray,
    k: np.ndarray,
    interfaces: List[int],
    dx: float,
    dt: float,
    h_out: float,
    h_in: float,
    t_env: float,
    t_skin: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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


def _thomas_prepare(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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


def _thomas_solve(a: np.ndarray, c_prime: np.ndarray, denom: np.ndarray, d: np.ndarray) -> np.ndarray:
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
    materials: Dict[str, Dict[str, float]],
    thickness_mm: Dict[str, float],
    t_env: float,
    h_out: float,
    h_in: float,
    t_skin: float = 37.0,
    total_time_s: int = 5400,
    dx_mm: float = 0.1,
    dt_s: float = 1.0,
    store_all: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    layers = _build_layers(materials, thickness_mm)
    alpha, k, interfaces = _grid_from_layers(layers, dx_mm)
    dx = dx_mm / 1000.0
    n = len(alpha)

    a, b, c, ba, bb, bc, d = _build_matrices(
        alpha=alpha,
        k=k,
        interfaces=interfaces,
        dx=dx,
        dt=dt_s,
        h_out=h_out,
        h_in=h_in,
        t_env=t_env,
        t_skin=t_skin,
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


def steady_state_h_in(
    materials: Dict[str, Dict[str, float]],
    thickness_mm: Dict[str, float],
    t_env: float,
    t_skin: float,
    t_surface: float,
    h_out: float,
) -> float:
    resistance = 1.0 / h_out
    for layer_name in ["I", "II", "III", "IV"]:
        k = materials[layer_name]["k"]
        d_m = thickness_mm[layer_name] / 1000.0
        resistance += d_m / k

    q = (t_env - t_surface) / resistance
    return q / (t_surface - t_skin)


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(math.sqrt(np.mean((a - b) ** 2)))


def write_simple_xlsx(path: str | Path, header: List[str], data: np.ndarray) -> None:
    path = Path(path)
    from io import BytesIO

    def col_letter(idx: int) -> str:
        letters = ""
        while idx >= 0:
            letters = chr(idx % 26 + 65) + letters
            idx = idx // 26 - 1
        return letters

    def cell_ref(r: int, c: int) -> str:
        return f"{col_letter(c)}{r}"

    rows_xml = []
    r = 1
    header_cells = []
    for c, text in enumerate(header):
        header_cells.append(
            f"<c r=\"{cell_ref(r, c)}\" t=\"inlineStr\"><is><t>{text}</t></is></c>"
        )
    rows_xml.append(f"<row r=\"{r}\">{''.join(header_cells)}</row>")

    for i in range(data.shape[0]):
        r += 1
        cells = []
        for c in range(data.shape[1]):
            cells.append(f"<c r=\"{cell_ref(r, c)}\"><v>{data[i, c]}</v></c>")
        rows_xml.append(f"<row r=\"{r}\">{''.join(cells)}</row>")

    sheet_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
        f"<sheetData>{''.join(rows_xml)}</sheetData>"
        "</worksheet>"
    )

    workbook_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">"
        "<sheets><sheet name=\"Sheet1\" sheetId=\"1\" r:id=\"rId1\"/></sheets>"
        "</workbook>"
    )

    workbook_rels = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" "
        "Target=\"worksheets/sheet1.xml\"/>"
        "</Relationships>"
    )

    root_rels = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" "
        "Target=\"xl/workbook.xml\"/>"
        "</Relationships>"
    )

    content_types = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
        "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
        "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
        "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
        "<Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
        "</Types>"
    )

    buffer = BytesIO()
    with ZipFile(buffer, "w") as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", workbook_xml)
        z.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        z.writestr("xl/worksheets/sheet1.xml", sheet_xml)

    path.write_bytes(buffer.getvalue())
