# DOCKER_ID=simic_run_normal
# docker run -itd --name $DOCKER_ID guisesanz/simic:1.1.0 bash

# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Normal.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Normal_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Normal_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Normal_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Normal_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# 					                /root/SimiC/Data_MDS/POST_Run_Normal_1000.TF.pickle

# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/15.2_LaunchSimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3 -u /root/SimiC/pickle_helper.py


# nohup docker exec $DOCKER_ID python3 -u /root/SimiC/15.2_LaunchSimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &

# DOCKER_ID=simic_run_del1
# docker run -itd --name $DOCKER_ID guisesanz/simic:1.1.0 bash

# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del1.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del1_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del1_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del1_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del1_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/15.2_LaunchSimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3 -u /root/SimiC/pickle_helper.py


# nohup docker exec $DOCKER_ID python3 -u /root/SimiC/15.2_LaunchSimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &


# DOCKER_ID=simic_run_del2
# docker run -itd --name $DOCKER_ID guisesanz/simic:1.1.0 bash

# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del2.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del2_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del2_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del2_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del2_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/15.2_LaunchSimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3 -u /root/SimiC/pickle_helper.py


# nohup docker exec $DOCKER_ID python3 -u /root/SimiC/15.2_LaunchSimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &


# DOCKER_ID=simic_run_del3
# docker run -itd --name $DOCKER_ID guisesanz/simic:1.1.0 bash

# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del3.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del3_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del3_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del3_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del3_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/15.2_LaunchSimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3 -u /root/SimiC/pickle_helper.py


# nohup docker exec $DOCKER_ID python3 -u /root/SimiC/15.2_LaunchSimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &



DOCKER_ID=simic_run_prepost
docker run -itd --name $DOCKER_ID guisesanz/simic:1.1.0 bash

docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# Copy the data
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del_PREPOST.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del_PREPOST_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del_PREPOST_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del_PREPOST_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/POST_Run_Del_PREPOST_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# Copy the script
docker cp /home/tereshkova/data/gserranos/MDS/15.2_LaunchSimiC.py $DOCKER_ID:/root/SimiC/
docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
#### WARNING PICKLE NOT WORKING
docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

docker exec $DOCKER_ID python3 -u /root/SimiC/pickle_helper.py


nohup docker exec $DOCKER_ID python3 -u /root/SimiC/15.2_LaunchSimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &

