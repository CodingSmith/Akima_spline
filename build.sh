#!/bin/sh

dire="build"

if [ -d "$dire" ]; then
rmdir "$dire"
mkdir "$dire"
else
mkdir "$dire"
fi

cd "$dire"

cmake .. && make install
