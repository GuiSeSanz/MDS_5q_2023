


# DOCKER_IMAGE='guisesanz/simic:1.1.0'
# eval "DOCKER_ID=$( docker run -d -t $DOCKER_IMAGE )";
# echo "The id of the running container is: $DOCKER_ID"


# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/5qVsElder.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/5qVsElder_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/5qVsElder_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS

# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/06.15_Launch_SimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3  /root/SimiC/pickle_helper.py
# docker exec $DOCKER_ID python3  /root/SimiC/06.15_Launch_SimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt

# # docker exec -it $DOCKER_ID bash
# # docker stop $DOCKER_ID && docker rm $DOCKER_ID




# DOCKER_IMAGE='guisesanz/simic:1.1.0'
# eval "DOCKER_ID=$( docker run -d -t $DOCKER_IMAGE )";
# echo "The id of the running container is: $DOCKER_ID"

# DOCKER_ID=boring_knuth

# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/5qVsnon5q.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/5qVsnon5q_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/5qVsnon5q_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS

# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/06.15_Launch_SimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3  /root/SimiC/pickle_helper.py
# docker exec $DOCKER_ID python3 -u /root/SimiC/06.15_Launch_SimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt

# docker exec -it $DOCKER_ID bash
# docker stop $DOCKER_ID && docker rm $DOCKER_ID


# DOCKER_ID=simicprepost

# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/simicprepost.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/06.15_Launch_SimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3  /root/SimiC/pickle_helper.py
# nohup docker exec $DOCKER_ID python3 -u /root/SimiC/06.15_Launch_SimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &




# docker run -itd --name Non5qVsElder guisesanz/simic:1.1.0 bash
# DOCKER_ID=Non5qVsElder

# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/Non5qVsElder.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/Non5qVsElder_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/Non5qVsElder_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/Non5qVsElder_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/Non5qVsElder_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/06.15_Launch_SimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3  /root/SimiC/pickle_helper.py
# nohup docker exec $DOCKER_ID python3 -u /root/SimiC/06.15_Launch_SimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &




# DOCKER_ID=prepost
# docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# # Copy the data
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/simicprepost.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
# docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/PreVsPost_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# # Copy the script
# docker cp /home/tereshkova/data/gserranos/MDS/06.15_Launch_SimiC.py $DOCKER_ID:/root/SimiC/
# docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
# #### WARNING PICKLE NOT WORKING
# docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

# docker exec $DOCKER_ID python3  /root/SimiC/pickle_helper.py
# nohup docker exec $DOCKER_ID python3 -u /root/SimiC/06.15_Launch_SimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &




docker run -itd --name simicthreeway guisesanz/simic:1.1.0 bash
DOCKER_ID=simicthreeway

docker exec $DOCKER_ID mkdir -p /root/SimiC/Data_MDS

# Copy the data
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/ThreeWayComparison.clustAssign.txt $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/ThreeWayComparison_1000.DF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/ThreeWayComparison_1000.TF.csv  $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/ThreeWayComparison_1000.DF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS
docker cp /home/tereshkova/data/gserranos/MDS/Data/SimiC/ThreeWayComparison_1000.TF.pickle  $DOCKER_ID:/root/SimiC/Data_MDS


# Copy the script
docker cp /home/tereshkova/data/gserranos/MDS/06.15_Launch_SimiC.py $DOCKER_ID:/root/SimiC/
docker cp /home/tereshkova/data/gserranos/MDS/Filter_weigths.R  $DOCKER_ID:/root/SimiC/
#### WARNING PICKLE NOT WORKING
docker cp /home/tereshkova/data/gserranos/MDS/pickle_helper.py  $DOCKER_ID:/root/SimiC/

docker exec $DOCKER_ID python3  /root/SimiC/pickle_helper.py


nohup docker exec $DOCKER_ID python3 -u /root/SimiC/06.15_Launch_SimiC.py > ./Data/SimiC/Results$DOCKER_ID.txt &

