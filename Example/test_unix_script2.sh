C=660 #change this to your inferred noise parameter 
for s2 in `seq 1 100`; do
## Update path to .run file
/rds/user/mg878/hpc-work/package_lab/Codes/Codes_for_xStar/SimRealObsCount.run -if /rds/user/mg878/hpc-work/package_lab/Transmission1 -of /rds/user/mg878/hpc-work/package_lab/Transmission1_xStar -simOnly true -importHapsMahan /rds/user/mg878/hpc-work/package_lab/Transmission1_qStar -format Mahan -rmMonoSim false -C $C -ss $s2 -filterData False 
done 