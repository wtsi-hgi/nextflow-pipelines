#!/bin/bash

jq --arg study_id $study_id -n '{avus: [
       {attribute: "sample", value: "INT_RNA7427205", o: "="}, 
      {attribute: "target", value: "1", o: "="}]}' |\
singularity exec /software/hgi/containers/singularity-baton/baton.simg baton-metaquery \
	    --zone seq --obj --avu

       #{attribute: "study_id", value: $study_id, o: "="}, 
#//    |\
#jq '.[] as $a| 
#"\($a.avus | .[] | select(.attribute == "sample") | .value)____\($a.avus | .[] | select(.attribute == "sample_supplier_name") | .value)____\($a.avus | .[] | select(.attribute == "id_run") | .value)____\($a.avus | .[] | select(.attribute == "is_paired_read") | .value)____\($a.avus | .[] | select(.attribute == "study_id") | .value)____\($a.avus | .[] | select(.attribute == "study") | .value)"' |\
#    sed s"/$(printf '\t')//"g |\
#    sed s"/\"//"g |\
#    sed s"/____/$(printf '\t')/"g |\
#sort | uniq >> samples.tsv




# jq '.[] as $a | $a.avus | .[] as $item | "\($a) \($item)"' #  | .attribute as $key | "\($key)_____\(.value)"' | sed s'/\"//'g > study_id.txt

#cat study_id.txt |\
#    grep -P 'sample' > sample_i


#jq -n '{avus: [
#       {attribute: "study_id", value: "5643", o: "="}, 
#       {attribute: "target", value: "1", o: "="}]}' | \
#    singularity exec /software/hgi/containers/singularity-baton/baton.simg baton-metaquery \
#		--zone seq --obj --avu \
#    | jq '.[] | .avus | .[] | 
#        select(.attribute == "sample", 
#               .attribute == "sample_supplier_name")' 

#jq -n '{avus: [
#       {attribute: "study_id", value: "5643", o: "="}, 
#       {attribute: "target", value: "1", o: "="}]}' | \
#    singularity exec /software/hgi/containers/singularity-baton/baton.simg baton-metaquery \
#		--zone seq --obj --avu | jq -r '.[] | .avus | .[] | select(.attribute == "sample") | .value' \
#    | sort | uniq > samples.txt
