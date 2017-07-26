#!/bin/bash

echo ''
echo 'Checking Platypus installation'

# Test running from the directory in which the script was invoked
platypus --help &> /dev/null

if [ $? != 0 ]
then
    echo 'Error: Platypus installation failed. Cannot run platypus from directory ' $(pwd)
    exit 1
fi

# Test running from directory that contains the script
ABSOLUTE_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd $ABSOLUTE_PATH

platypus --help &> /dev/null

if [ $? != 0 ]
then
    echo 'Error: Platypus installation failed. Cannot run platypus from directory ' $ABSOLUTE_PATH
    exit 1
fi


echo 'Platypus installation succeeded'
echo ''
