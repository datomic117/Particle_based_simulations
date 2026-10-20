#!/bin/bash
# Part D3: collect the measured wall-clock times of the silo runs (the `time` output kept in runs/<name>/time.txt) and the
# simulated time into one small table, data/D_walltimes.csv. Run from Assignment_4/:  bash code/collect_times.sh
# The simulated time is the last time in the table of the run (runs/<name>/run.dat); steps = last Step of that table.
echo "run,D_over_dp,wall_s,t_end_s,steps" > data/D_walltimes.csv
for spec in D3:3 D6:6 D7:7 D8:8 D9:9; do
  name=${spec%%:*}; n=${spec##*:}
  real=$(grep real runs/$name/time.txt | awk '{print $2}')            # e.g. 3m54.010s
  wall=$(echo "$real" | sed 's/m/ /; s/s//' | awk '{printf "%.1f", $1*60+$2}')
  last=$(tail -1 runs/$name/run.dat | awk '{print $2","$1}')            # time,step
  echo "$name,$n,$wall,$last" >> data/D_walltimes.csv
done
cat data/D_walltimes.csv
