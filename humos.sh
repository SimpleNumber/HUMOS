#!/bin/sh
HUMOSPATH="./HUMOS"
NUMWORKERS=2

if [ $1 = "start" ]
then
   cd $HUMOSPATH
   truncate -s 0 ~/.gunicorn-log
   gunicorn -b 0.0.0.0 -w $NUMWORKERS -D --log-file ~/.gunicorn-log --pid ~/.gunicorn.pid --preload app:server
elif [ $1 = "restart" ]
then
   kill -s TERM $(cat ~/.gunicorn.pid)
   cd $HUMOSPATH
   truncate -s 0 ~/.gunicorn-log
   gunicorn -b 0.0.0.0 -w $NUMWORKERS -D --log-file ~/.gunicorn-log --pid ~/.gunicorn.pid --preload app:server
elif [ $1 = "stop" ]
then
   kill -s TERM $(cat ~/.gunicorn.pid)
fi
