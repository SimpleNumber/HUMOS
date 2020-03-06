**HUMOS** is a web application that requires a server to run. There are two ways how you can get **HUMOS**: either pull the docker container or install an instance of it manually. The first option is easier, while the second one allows better tuning **HUMOS** for your needs. Please see the corresponding sections below for more information.

# Manual installation of HUMOS

## 1. Requirements 

**HUMOS** has very moderate requirements regarding computational power, however, if you want to set it up for multiple concurent users, for example, to be used during the workshop, you should consider using server with several available processors. More detailed discussion can be found below (#4).

Use your favorite package manager to install `python3` and `pip`. You will need *sudo rights* to perform the installation.
Depending on the server configuration you might need *sudo rights* to use `pip` as well. Please, ask your system adminstrator to help you, if necessary.
    
## 2. Copying files

Download the latest vesion of repository from Github either manually using [Clone or download](https://github.com/SimpleNumber/HUMOS/archive/master.zip), or use `wget` (for example).
Unpack the zip into some folder, for example `HUMOS`

```shell
    wget https://github.com/SimpleNumber/HUMOS/archive/master.zip #downloads the last version
    unzip master.zip #this creates folder called HUMOS-master, 
    mkdir HUMOS #let's move it to folder HUMOS for simplicity
    mv HUMOS-master/* HUMOS/
    rmdir HUMOS-master
    rm master.zip
```

You can clone the complete repository, if you like

```
    git clone https://github.com/SimpleNumber/HUMOS.git
```
All files necessary for the installation are located in *installation* folder. 

## 3. Python packages

You will need to have the following python packages installed:

```
    numpy
    pandas
    IsoSpecPy
    pyteomics
    dash
    dash_table
    plotly
```

They are listed in `requirements.txt` that is downloaded as a part of repository and you can install all of them using the following command:

```shell
    pip install -r requirements.txt
```

To serve the app using several processors we will use `gunicorn`, if you want to run with single processor you can skip this step.
Install `gunicorn` using `pip`:

```shell
    pip install gunicorn
```

## 4a. Setting up Green Unicorn

`Green Unicorn` allows you to serve the same app using several processes on Linux. In the simplest case you can invoke it like this from the `HUMOS` folder.
```shell
    gunicorn -b 0.0.0.0:[port] -w [number_of_workers] --preload app:server
```
The description of the parameters:

`-b 0.0.0.0:[port]` binds the server to all network interfaces of the current machine, optinally using the specific port.

**NOTE**: The default port is 8000. If you want to deploy on standard HTTP port (80) you will need to setup `authbind` to run `gnicorn`
(Check this [guide](https://mutelight.org/authbind), for example). Alternatively, you can use `nginx`. Refer to [gunicorn documentation](http://docs.gunicorn.org/en/stable/deploy.html)
for details.

`-w` sets the number of worker processes to spawn. This number depends on the number of **concurrent** users that will use the app. In our hands, 5 workers were able to serve a group
of 10 users (for example, 10 students are working with the app at the same time during the workshop) without any visible lag. Thus, a rule of thumb is \[numer_of_concurent_users\]/2.
You might need to optimize it for your needs. For more detailed discussion on the selection of the number of workers, refer to [gunicorn documentation](http://docs.gunicorn.org/en/stable/design.html#how-many-workers).

`--preload` is used to set up the same peptide set for all workers. Since **HUMOS** initializes the set of peptides used for modelling at the start, it is important to fork the app process only after the initialization, otherwise
workers will receive different set of peptides for modelling.

Slightly more elaborated invokation is the following:
```shell
    gunicorn -b 0.0.0.0:[port] -w [number_of_workers] -D --logfile /.gunicorn-log --pid /.gunicorn.pid --preload app:server
```

This also demonizes `gunicorn` master process, set logging into `.gunicorn-log` and writes PID of the master process to `.gunicorn.pid`.

In the repository you can find shell script `humos.sh`, that simplifies running HUMOS. You might need to adjust the `HUMOSPATH` and `NUMWORKERS` variables in the script to
match your requirements. The usage is:
```shell
    humos.sh [start|restart|stop]
```

## 4b. Starting single instance

If you do not want to use `gunicorn` you can start the app with `flask` server provided with `dash`:
```shell
    python app.py
```

## 5. Enjoy your fresh HUMOS

# Docker container

If you are unfamiliar with Docker containers, check [Docker website](https://www.docker.com/get-started) to learn more. You will need to install Docker Engine suitable for your system before you can use any containers.

When you have Docker installed, run the following command to pull the latest build of the container from Docker Hub

```shell
    docker pull simplenumber/humos_dock
```

Container expose port 8000, that can be rebound to any other port when running container. Environmental variable NUMWORKERS allows changing number of concurent workers serving user requests (default is 2). Later you can run HUMOS container using the following command.

```shell
    docker run -d -p [port]:8000 --rm --env NUMWORKERS=[N] simplenumber/humos_dock
```

The arguments are:

   `-d` runs container in detached mode, i.e. the container will run in background. If you prefer to run in foreground substitute this key with `-it`, to attach STDIN and assign TTY for it.

   `-p [port]:8000` binds the exposed port 8000 to [port] on container host, for example, you can use `-p 80:8000` to run **HUMOS** on a standard port.
   
   `--rm` instruct Docker Engine to clean up the container and its filesystem when it quits (so next time you will start with fresh **HUMOS** again)
   
   `--env` allows sending environmental variables to container, for example, `--env NUMWORKERS=8` will set NUMWORKERS variable to 8 and, thus, 8 worker processes will be created. Section [4a](#4a-setting-up-green-unicorn) have a discussion on how many workers you need. If you are happy with the default value (2) just ommit this key.

To stop the container running in detached mode you should first find the Container ID by running the command. The ID is in the very first column.

```shell
    docker container ls | grep simplenumber/humos_dock
```

Later you can stop the container by running the following with the Container ID

```shell
    docker container stop [ContainerID]
```
