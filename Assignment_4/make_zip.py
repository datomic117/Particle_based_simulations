import os
import zipfile

OUT = '8_A4.zip'
if os.path.exists(OUT):
    os.remove(OUT)
SKIP_DIRS = {'__pycache__', 'html', '.vscode'}      # html/ = generated documentation of the starting code
n = 0
with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED) as z:
    z.write('A4.ipynb', 'A4.ipynb')
    z.write('A4.html', 'A4.html')
    for top in ('code', 'data'):
        for root, dirs, files in os.walk(top):
            dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
            for f in sorted(files):
                if f.lower().endswith(('.exe', '.o', '.pyc', '.xyz')) or f in ('restart.dat', 'trajectories.xyz'):
                    print('skipped', root, f)
                    continue
                p = os.path.join(root, f).replace(os.sep, '/')
                z.write(p, p)
                n += 1
print('files from code/ and data/:', n)
with zipfile.ZipFile(OUT) as z:
    names = z.namelist()
    print(len(names), 'members; backslashes:', sum('\\' in x for x in names))
    print('size', os.path.getsize(OUT))
    for nm in names:
        if nm.startswith('code/'):
            print(' ', nm)
    print('data files:', sum(x.startswith('data/') for x in names))
