#!/bin/bash

touch /tmp/foobar

ulimit -c unlimited
echo "/home/debian/cores/core.%p" > /proc/sys/kernel/core_pattern
DIR="$( cd "$( dirname $0 )" && pwd )"
$DIR/bringup-can.sh
sleep 1 # Shouldn't be necessary, but give can interface time to come up
$DIR/logger_main -logfilename /home/debian/logfile-$(date +%s) &
$DIR/can-dump &
$DIR/server_main &
$DIR/simple_control_main &
$DIR/state_estimator_main &
$DIR/scamp_main &
$DIR/sbus-test-run &
