#!/bin/bash

echo "Eseguire suite di test? [Y/n]"
read a

if [ $a == "n" ]
then
    sed -i 's/main/main_iplib/g' Makefile
    mv main_iplib main
else
    sed -i 's/main_iplib/main/g' Makefile
fi

make