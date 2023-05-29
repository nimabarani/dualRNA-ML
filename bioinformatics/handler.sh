#!/bin/bash

json_data=$(cat setup.json)

for data in $(echo $json_data | jq -c '.[]')
do
    sbatch featureCounts.sh $data
done
