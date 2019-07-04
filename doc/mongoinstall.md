# Installing mongodb
We have tested this application with Mongo 3.6.1 (Windows, and MongoDb atlas), and Mongo 4 (linux, windows).
It does not work with 2.6.  

## Remote server
If you are using a remote mongodb server, you do not need to install mongodb locally.  We have tested findNeighbour3 with [MongoDB Atlas](https://www.mongodb.com/cloud/atlas).

## Windows
MSI installers are [available](https://www.mongodb.com/download-center).
The server can be started/stopped using the Services application, available via the Control Panel.

## Linux
Please read the [documentation](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/).
The below installs MongoDB 4 on Ubuntu 16.04 LTS:
```
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 9DA31620334BD75D9DCB49F368818C72E52529D4
echo "deb [ arch=amd64,arm64 ] https://repo.mongodb.org/apt/ubuntu xenial/mongodb-org/4.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb.list
sudo apt update
sudo apt install mongodb-org
```

It is automatically started but can be started with
```
sudo systemctl enable mongod
sudo systemctl start mongod
```

It can be stopped/restarted with
```
sudo systemctl stop mongod
sudo systemctl restart mongod
```  

If a sharded cluster is considered necessary, please see [here](mongosharding.md).

