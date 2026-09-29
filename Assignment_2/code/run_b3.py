"""Build/run one B3 protocol, preserving parameters, hashes and wall time.

Run from any directory: python code/run_b3.py PROTOCOL STAGE [DURATION_PS] [FD_DELTA] [VIRIAL_H]
The selectors are compile-time options in setparameters.c, not simulation input files.
Existing outputs are protected against overwriting.
"""
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

code = Path(__file__).resolve().parent
root = code.parent
protocol, stage = map(int, sys.argv[1:3])
duration = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0
delta = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-6
virial_h = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-6
assert protocol in range(4) and stage > 0 and duration > 0 and delta > 0
name = (f'b3_final_equil_{stage:02d}' if protocol == 0 else
        f'b3_final_fd_{delta:g}' if protocol == 1 else
        'b3_final_nve_1fs' if protocol == 2 else 'b3_final_nve_05fs')
data = root / 'data'
output = data / (name + ('.txt' if protocol == 1 else '.csv'))
if output.exists():
    raise SystemExit(f'Refusing to overwrite {output}')
restart = data / ('b3_equil_long_restart.dat' if protocol == 0 and stage == 1 else
                  f'b3_final_equil_{stage-1 if protocol == 0 else stage:02d}_restart.dat')
if not restart.is_file():
    raise SystemExit(f'Missing input: {restart}')
build = root / '.b123_work'
build.mkdir(exist_ok=True)
exe = build / (name + ('.exe' if sys.platform == 'win32' else ''))
sources = sorted(code.glob('*.c'))
command = ['gcc', '-O3', '-Wall', '-Wextra', f'-DB3_PROTOCOL={protocol}',
           f'-DB3_EQ_STAGE={stage}', f'-DB3_DURATION_PS={duration:g}', f'-DB3_FD_DELTA={delta:g}',
           f'-DB3_VIRIAL_H={virial_h:g}',
           *[p.name for p in sources], '-o', str(exe), '-lm']
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
metadata = {'protocol': protocol, 'stage': stage, 'duration_ps': duration, 'fd_delta': delta, 'virial_h': virial_h,
            'machine': platform.platform(), 'compiler': subprocess.check_output(['gcc', '--version'], text=True).splitlines()[0],
            'compile_argv': command, 'run_cwd': 'code', 'run_argv': [str(exe)],
            'input_restart': str(restart.relative_to(root)), 'input_sha256': sha(restart),
            'source_sha256': {p.name: sha(p) for p in sources + sorted(code.glob('*.h'))},
            'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}
(data / (name + '_setparameters.txt')).write_text((code/'setparameters.c').read_text(), encoding='utf-8')
start = time.perf_counter()
with (data / (name + '_build.txt')).open('w') as log:
    subprocess.run(command, cwd=code, stdout=log, stderr=subprocess.STDOUT, check=True)
metadata['compile_seconds'] = time.perf_counter()-start
meta_path = data / (name+'_meta.json')
meta_path.write_text(json.dumps(metadata, indent=2), encoding='utf-8')
start = time.perf_counter()
print(f'Starting {name}', flush=True)
with output.open('w', encoding='utf-8') as stdout, (data / (name+'_stderr.txt')).open('w') as stderr:
    result = subprocess.run([str(exe)], cwd=code, stdout=stdout, stderr=stderr)
metadata.update(wall_seconds=time.perf_counter()-start, returncode=result.returncode,
                output_sha256=sha(output))
if protocol != 1:
    final = data / (name+'_restart.dat')
    if final.exists(): metadata['final_restart_sha256'] = sha(final)
meta_path.write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print(f'Finished {name}: exit {result.returncode}, {metadata["wall_seconds"]:.1f} s', flush=True)
raise SystemExit(result.returncode)
