#!/bin/bash

proc_all=$1; shift
proc_use=$1; shift
ptile=$1; shift
res_file=$1; shift
cmd=$@

cat <<HERE
#!/bin/bash

# BSUB -q hpc_linux
# BSUB -n ${proc_all}
# BSUB -R "span[ptile=12]"

rm -f mpd.hosts
for host in \`echo \${LSB_HOSTS}\`
do
    echo \${host} >> mpd.hosts
done
n_hosts=\`cat mpd.hosts | sort | uniq | wc -l\`

mpdcleanup -f mpd.hosts -r ssh -u xuewei
mpdboot -n \${n_hosts} -r ssh

run_cmd="mpiexec"
if [ ${ptile} -eq 0 ]
then
    run_cmd="\${run_cmd} -n ${proc_use} ${cmd}"
else
    for host in \`cat mpd.hosts | sort | uniq | head -q -n ${proc_use}\`
    do
        run_cmd="\${run_cmd} -n ${ptile} -host \${host} ${cmd} :"
    done
    run_cmd=\`echo \${run_cmd} | sed "s/..$//"\`
fi

echo \${run_cmd}
\${run_cmd} &> \${res_file}

mpdcleanup -f mpd.hosts -r ssh -u xuewei
rm -f mpd.hosts
HERE
