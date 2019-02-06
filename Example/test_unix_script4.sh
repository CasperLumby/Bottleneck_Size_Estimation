mkdir Transmission1_bottleneck
for t in `seq 0 7`; do
   cd Transmission1_bottleneck
   mkdir segment_$t
   cd ..
   cat /rds/user/mg878/hpc-work/package_lab/Transmission1_qStarStar/*/qStarStar_$t.txt > /rds/user/mg878/hpc-work/package_lab/Transmission1_bottleneck/segment_$t/test_$t.txt
done 