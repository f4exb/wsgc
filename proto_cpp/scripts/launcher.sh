#!/bin/bash
# arg #1 : number of CPUs
# arg #2 : number of runs
# arg #3 : prefix tag i.e. _RES
# arg #4 : suffix for result file i.e. wsgc_result${suffix}
nb_cpus=$1
shift
nb_runs=$1
shift
tag=$1
shift
suffix=$1
shift
command=""

whitespace="[[:space:]]"
for i in "$@"
do
    if [[ $i =~ $whitespace ]]
    then
        i=\"$i\"
    fi
    command="${command} ${i}"
done

function run()
{
    sleep 1
    
    cpu_id=$1
    shift
    nb_runs=$1
    shift
    tag=$1
    shift
    suffix=$1
    shift
    command=$*
    
    #source /shared/softs/cuda-5.0/setenv
    
    for (( run_i=1; run_i<=${nb_runs}; run_i++))
    do
        echo ${cpu_id} ${run_i} ${command}
        eval ${command} | grep ${tag} | cut -d ' ' -f 2 >> "wsgc_result${suffix}-${cpu_id}"
        #${command} | grep "_SOR" >> "wsgc_result${suffix}-${cpu_id}"
    done
}


function ChildReturned()
{
    echo "A child returned"
    jobs -l
}


rm -f wsgc_result${suffix}*
ts_start=$(date +%s.%N)

# Launch sub-processes
for (( cpu_i=0; cpu_i<nb_cpus; cpu_i++ ))
do
    (run ${cpu_i} ${nb_runs} ${tag} ${suffix} ${command} -y ${cpu_i}) &
    #taskset -pc ${cpu_i},$((${cpu_i}+4)) $!
    taskset -pc ${cpu_i} $!
done

# Wait for all sub-processes to finish
for i in $(jobs -p)
do
    wait $i
done

ts_stop=$(date +%s.%N)
echo "${ts_stop}-${ts_start}" | bc

# Concatenate the results
for (( cpu_i=0; cpu_i<nb_cpus; cpu_i++ ))
do
    cat wsgc_result${suffix}-${cpu_i} >> wsgc_result${suffix}
    rm wsgc_result${suffix}-${cpu_i}
done
