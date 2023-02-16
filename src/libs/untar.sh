#!/bin/bash

item=$1

if [ ! -d $item ];
then
    tar xzf $item.tgz
fi
