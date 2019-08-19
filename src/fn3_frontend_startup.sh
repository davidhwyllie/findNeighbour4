dockerImageName=triendo/fnmonitor
dockerpid=`docker ps -a | grep $dockerImageName | grep "Up" | awk -F " " '{ print $1 }'`
if [[ $dockerpid != "" ]];then
   docker kill $dockerpid
   docker rm $dockerpid
fi
docker pull $dockerImageName
docker run -d -p 3020:80 --rm $dockerImageName

