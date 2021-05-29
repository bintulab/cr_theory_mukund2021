#!/bin/bash

while IFS= read -r line; do
    tube="../data/${line}.fcs"
    echo "Copying: ${tube}"
    cp ${tube} .
done < tubes.txt
