#!/bin/bash

proc_all=$1; shift
proc_use=$1; shift
ptile=$1; shift
cmd=$@
host_file=hosts.`date +%s`

cat <<HERE
#!/bin/bash

# BSUB -I
# BSUB -q hpc_linux
# BSUB -n ${proc_all}
# BSUB -R "span[ptile=12]"

rm -f ${host_file}
for host in \`echo \${LSB_HOSTS}\`
do
    echo \${host} >> ${host_file}
done
n_hosts=\`cat ${host_file} | sort | uniq | wc -l\`

mpdcleanup -f ${host_file} -r ssh -u xuewei
mpdboot -f ${host_file} -n \${n_hosts} -r ssh

run_cmd="mpiexec"
if [ ${ptile} -eq 0 ]
then
    run_cmd="\${run_cmd} -n ${proc_use} ${cmd}"
else
    for host in \`cat ${host_file} | sort | uniq | head -q -n ${proc_use}\`
    do
        run_cmd="\${run_cmd} -n ${ptile} -host \${host} ${cmd} :"
    done
    run_cmd=\`echo \${run_cmd} | sed "s/..$//"\`
fi

echo \${run_cmd}
\${run_cmd}

mpdcleanup -f ${host_file} -r ssh -u xuewei
rm -f ${host_file}
HERE
