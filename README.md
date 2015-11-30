# SherpaHggInterference


    cd hintmcnlo
    autoreconf -fi
    svn co http://sherpa.hepforge.org/svn/branches/hintmcnlo
    ./configure <your preferred flags> -->

    ./configure --enable-analysis --disable-silent-rules --enable-hepmc2=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/h
    epmc/2.06.07-cms4 --with-sqlite3=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/sqlite/3.7.17
    see here the option neeed in cmssw:https://github.com/cms-sw/cmsdist/blob/IB/CMSSW_7_2_X/stable/sherpa.spec
    make install 