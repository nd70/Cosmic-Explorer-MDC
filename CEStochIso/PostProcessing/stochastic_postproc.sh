#!/bin/sh -e

disp=$$
nohup /usr/X11R6/bin/Xvfb :${disp} 1>/dev/null &
pid=$!
env DISPLAY=:${disp}.0 ./postProcessScriptFull $*
kill $pid

