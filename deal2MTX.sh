#!/bin/bash
sed -i -e 's/[\(\)\,]/ /g' matrix.out
./deal2MTX
