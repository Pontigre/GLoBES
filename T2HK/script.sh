#!/bin/bash

echo 'Turn off screensaver'
sleep 10 &&
xset s 0 0
xset s off

if [[ -e th13delta ]] ; then
  echo 'The file exists.'
  ./th13delta
else
  echo 'The file does not exist.';  sleep 3; exit 0
fi

#Shutdown and poweroff
sleep 3
sudo shutdown -P now
