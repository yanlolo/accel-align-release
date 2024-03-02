#!/bin/bash

while true; do
    read -ep $'\033[1;33m [rmi.sh] \033[0mwhich index would you like to train? [1-10] ' choice

    # Check if the input is a number
    if [[ "$choice" =~ ^[1-9]$|^10$ ]]; then
        break  # Break out of the loop if it's a valid number in the range [1-10]
    else
        echo "            Please enter a number between 1 and 10."
    fi
done

choice=$((choice + 1))

read type branching size avg_err max_err b_time <<< $(sed -n "${choice}p" optimizer.out)
echo $type $branching $size $avg_err $max_err $b_time > rmi_type.txt