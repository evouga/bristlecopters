#!/bin/bash

./mattia2poly < "$1" > tmp.poly
./triangle/triangle -pDBPNEga"$2" tmp.poly
./off2obj < tmp.1.off
rm tmp.poly
rm tmp.1.off
