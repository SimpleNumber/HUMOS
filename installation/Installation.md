#Installation of HUMOS#

**HUMOS** is a web application that requires a server to run. In this guide we will use Linux to set it up . 

##1. Requirements##
**HUMOS** has very moderate requirements regarding computational power, however, if you want to set it up for multiple concurent users, for example,
to be used during the workshop, you should consider using server with several available processors. More detailed discussion can be found below (#4).

Use your favorite package manager to install `python3` and `pip`. You will need *sudo rights* to perform the installation.
Depending on the server configuration you might need *sudo rights* to use `pip` as well. Please, ask your system adminstrator to help you, if necessary.
    
##2. Copying files##
Download the latest vesion of repository from BitBucket either manually using [Downloads section](https://bitbucket.org/J_Bale/humos/downloads/), or use `curl` (for example).
Unpack the zip into some foler, for example `HUMOS`

```shell
    curl -L https://bitbucket.org/J_Bale/humos/get/tip.zip -o humos.zip #downloads the last version
    unzip humos.zip #there is a folder with not that pretty name, like J_Bale-humos-xxxx
    mkdir HUMOS #let's make it prettier, obviosly you have to adjust all pathes to your system
    cp -r J_Bale-humos-xxxx/* HUMOS/
    rm -r J_Bale-humos-xxxx
```

You can clone the complete repository, if you like

```
    hg clone https://bitbucket.org/J_Bale/humos
```

##3. Python packages##
You will need to have the following python packages installed

```
    numpy
    pandas
    IsoSpecPy
    pyteomics
    dash
    dash_table
    plotly
```

They are listed in `requirements.txt` that is downloaded as a part of repository and you can install all of them using the following command

```shell
    pip instlall -r requirements.txt
```

To serve the app using several processors we will use `gunicorn`, if you want to run with single processor you can skip this step.
Install `gunicorn` using `pip`

```shell
    pip install gunicorn
```

##4a. Setting up Green Unicorn##
`Green Unicorn` allows you to serve the same app using several processes on Linux. In the simplest case you can invoke it like this from the `HUMOS` folder.
```shell
    gunicorn -b 0.0.0.0:[port] -w [number_of_workers] --preload app:server
```
The description of the parameters:

`-b 0.0.0.0:[port]` binds the server to all network interfaces of the current machine, optinally using the specific port

**NOTE** The default port is 8000. If you want to deploy on standard HTTP port (80) you will need either setup `authbind` to run `gnicorn`
(Check this [guide](https://mutelight.org/authbind), for example). Alternatively, you can use `nginx`. Refer to [gunicorn documentation](http://docs.gunicorn.org/en/stable/deploy.html)
for details.

`-w` sets the number of worker processes to spawn. This number depends on the number of **concurrent** users that will use the app. In our hands, 5 workers was able to serve a group
of 10 users (for example, 10 students working with the app at the same time during the workshop) without any visible lag. Thus, a rule of thumb is \[numer_of_concurent_users\]/2.
You might need to optimize it for your needs. For more detailed discussion on the selection of the number of workers, refer to [gunicorn documentation](http://docs.gunicorn.org/en/stable/design.html#how-many-workers)

`--preload` since **HUMOS** initializes the set of peptides used for modelling at the start, it is importnat to fork the app process only after the initialization, otherwise
workers will receive different set of peptides for modelling

Slightly more elaborated invokation is the following
```shell
    gunicorn -b 0.0.0.0:[port] -w [number_of_workers] -D --logfile /.gunicorn-log --pid /.gunicorn.pid --preload app:server
```

This also demonizes `gunicorn` master process, set logging into `.gunicorn-log` and writes PID of the master process to `.gunicorn.pid`

In the repository you can find shell script `humos.sh`, that simplifies running HUMOS. You might need to adjust the `HUMOSPATH` and `NUMWORKERS` variables in the script to
match your requirements. The usage is:
```
    humos.sh [start|restart|stop]
```

##4b. Starting single instance##
If you do not want to use `gunicorn` you can start the app with `flask` server provided with `dash`.
```
    python app.py
```

##5. Enjoy your fresh HUMOS##