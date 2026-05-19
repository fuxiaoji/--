from pathlib import Path
import numpy as np

from model import load_appendix_xlsx, solve_heat_cn

ROOT = Path('.').resolve()
APPENDIX = ROOT / 'CUMCM-2018-Problem-A-Chinese-Appendix.xlsx'


def main() -> None:
    data = load_appendix_xlsx(APPENDIX)
    materials = data.materials
    thickness = {
        'I': 0.6,
        'II': 6.0,
        'III': 3.6,
        'IV': 5.0,
    }
    times, temps = solve_heat_cn(
        materials=materials,
        thickness_mm=thickness,
        t_env=75.0,
        h_out=115.0,
        h_in=8.0,
        t_skin=37.0,
        total_time_s=120,
        dx_mm=0.2,
        dt_s=1.0,
        store_all=False,
    )
    print('steps:', len(times))
    print('last skin temp:', float(temps[-1, 0]))


if __name__ == '__main__':
    main()
