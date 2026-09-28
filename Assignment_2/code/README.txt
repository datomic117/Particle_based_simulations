MD_student_v7.0.0
=================

Starting code for 6EMA02 Particle-based Simulations.

Build and run (Linux, macOS, or MSYS2 on Windows):

    gcc -O3 *.c -o md -lm                      production
    gcc -g -Wall -Wextra *.c -o md_debug -lm   debugging
    ./md > run.csv

On Windows the supplied .vscode configuration builds and debugs this for you
through MSYS2; open this folder in VS Code and use Terminal > Run Build Task.

Your task list is the generated documentation. Open

    html/index.html

in a browser and follow the "Related Pages > Todo List" link: every \todo in
it is one piece of code you have to write. The same markers are in the source.

Run configuration lives in setparameters.c only. There are no input files:
change a value there and recompile.

Everything here is generated from the course source tree; do not expect the
Doxygen site to update until you rebuild it with "doxygen Doxyfile".
