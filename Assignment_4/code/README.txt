DEM_student_v7.1.0
==================

Starting code for 6EMA02 Particle-based Simulations.

Build and run (Linux, macOS, or MSYS2 on Windows):

    gcc -O3 *.c -o dem -lm                      production
    gcc -g -Wall -Wextra *.c -o dem_debug -lm   debugging
    ./dem > run.dat

Choose the run with the macro RUN at the top of setparameters.c, or on the
command line: -DRUN=1 a test collision (Part A), -DRUN=2 the tilted bed
(Part B), -DRUN=3 the silo (Parts C and D). For example

    gcc -O3 -DRUN=3 -DD_ORIFICE=8 *.c -o dem -lm

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
