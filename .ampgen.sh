#lhcb
export env="AmpGen"
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_98python3 x86_64-centos7-gcc10-opt
#ampgen
export PATH=/besfs5/groups/psipp/psippgroup/public/liukun/qmi_ampgen-master/build/bin:$PATH
#export PATH=/publicfs2/ucas/user/zengshh/LHCb/gammacombo/build/bin:$PATH
export PYTHONPATH=""
echo "start AmpGen environment"

export PS1="\[\e[32m\][\[\e[36m\]${env} \[\e[33m\]\t \[\e[32m\]\u@\h \[\e[35m\]\w\[\e[32m\]]$\[\e[0m\] "
