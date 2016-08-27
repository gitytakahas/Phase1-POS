#!/bin/sh
DTSTAMP=`awk "BEGIN{print strftime(\"%y%b%d_%H%M\");}" `
LOGFN=Logs/log.$DTSTAMP
rm -f lastlog
ln -s $LOGFN lastlog
echo $LOGFN
echo
#echo JMT turned off INFO level messages because the TTCciSupervisor is way too chatty...
#echo
nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1953 -e ${BUILD_HOME}/pixel/XDAQConfiguration/Profile.xml -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelSuper &
nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1954 -e ${BUILD_HOME}/pixel/XDAQConfiguration/Profile.xml -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelTKFEC &
nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1955 -e ${BUILD_HOME}/pixel/XDAQConfiguration/Profile.xml -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelFEC &
nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1956 -e ${BUILD_HOME}/pixel/XDAQConfiguration/Profile.xml -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelFED &

# nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1953  -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelSuper &
# nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1954  -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelTKFEC &
# nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1955  -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelFEC &
# nohup ${XDAQ_ROOT}/bin/xdaq.sh -p 1956  -c $1 2>&1 -z pixel | tee ${LOGFN}_PixelFED &