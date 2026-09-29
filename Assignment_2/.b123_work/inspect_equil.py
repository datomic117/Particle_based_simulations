import csv
import math
from pathlib import Path
import statistics as st
import sys

p = Path(sys.argv[1])
lines = p.read_text(encoding='utf-8').splitlines()
# A running simulation can leave its final row partially flushed.
rows = [r for r in csv.DictReader(lines) if None not in r.values()]
rows = [r for r in rows if all(r.values())]
complete_ps = (int(rows[-1]['step']) // 5000) * 5
print('complete time (ps):', complete_ps)
for lo in range(0, complete_ps, 10):
    b = [r for r in rows if lo*1000 < int(r['step']) <= (lo+10)*1000]
    if len(b) == 100:
        print(f'{lo}-{lo+10} ps: U={st.mean(float(r["Epot_internal"]) for r in b):.2f}, '
              f'T3={st.mean(float(r["T_CH3_K"]) for r in b):.3f}, '
              f'T2={st.mean(float(r["T_CH2_K"]) for r in b):.3f}')
if complete_ps < 50:
    raise SystemExit()
print('Final complete 50 ps window; ten 5 ps blocks; 95% slope CI uses t(8)=2.306:')
for field in ['Epot_internal', 'T_CH3_K', 'T_CH2_K']:
    blocks = []
    for lo in range(complete_ps-50, complete_ps, 5):
        b = [float(r[field]) for r in rows if lo*1000 < int(r['step']) <= (lo+5)*1000]
        assert len(b) == 50
        blocks.append(st.mean(b))
    t = [2.5+5*i for i in range(10)]
    ssx = sum((x-st.mean(t))**2 for x in t)
    slope = sum((x-st.mean(t))*(y-st.mean(blocks)) for x,y in zip(t,blocks))/ssx
    residual = [y-st.mean(blocks)-slope*(x-st.mean(t)) for x,y in zip(t,blocks)]
    se = math.sqrt(sum(y*y for y in residual)/8/ssx)
    print(field, 'mean', st.mean(blocks), 'block SE', st.stdev(blocks)/math.sqrt(10),
          'slope', slope, '95% half-width', 2.306*se, 'stationary slope', abs(slope)<=2.306*se)
