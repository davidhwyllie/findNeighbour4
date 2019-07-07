# Frontend
There is an interactive front end to findNeighbour3, written in AngularJS.    
A short guide to how to use the front end is [here](https://davidhwyllie.github.io/FNMFINDNEIGHBOUR3/fnm3.pdf).  Source code is [here](https://gitlab.com/ModernisingMedicalMicrobiology/fnmonitor).   

Installation is possible using [Docker]((https://docs.docker.com/install/linux/docker-ce/ubuntu/):  

```
dockerImageName=triendo/fnmonitor
dockerpid=`docker ps -a | grep $dockerImageName | grep "Up" | awk -F " " '{ print $1 }'`
if [[ $dockerpid != "" ]];then
   docker kill $dockerpid
   docker rm $dockerpid
fi
docker pull $dockerImageName
docker run -d -p 3020:80 --rm $dockerImageName
```

A script launching the frontend is [here](src/frontend.sh).  To run:
```
cd src
chmod +x frontend.sh
sudo ./frontend.sh   # may be able to run without root access depending how docker is set up
```

This has only been tested on Linux.  Depending on how Docker is configured, these commands (which can be put in a shell script) may need to be run as root.
After execution, an Angular server should be running.

If accessing from the server, this should be accessible at http://localhost:3020.  
The page will ask you for the URL of the findNeighbour3 server.

When it is running, you can enter its url, such as http://127.0.0.1:5020.
If you are accessing the server from a different computer, this works fine but note:
* Access is not restricted or authenticated as configured at present
* The LISTENTO parameter in the configuration file has to be set to allow the server to respond to external connections.



  

