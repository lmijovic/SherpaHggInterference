# SherpaHggInterference

Fork with updates to get HiggsLHCWG HiggsInterference code running locally

## Prerequisites 

* HEPMC 

```

git clone https://github.com/hep-mirrors/hepmc.git
./bootstrap
./configure --prefix=$PWD/install --with-momentum=MEV --with-length=MM
make
make install
/home/lm/disk/files/sw/mc/hepmc/install 

```

## Installation w/o blackhat: 

This should work fine to get the interference running: 

```

./configure --with-sqlite3=install --enable-hepmc2=/home/lm/disk/files/sw/mc/hepmc/install --enable-analysis --disable-silent-rules 

make

make install

```

## Installation with blackhat: 

Install blackhat; 


```
https://blackhat.hepforge.org/trac/wiki/BlackHatInstallation

download qd
/home/lm/disk/files/sw/mc/qd-2.3.22
./configure --enable-shared --prefix=$PWD/install

wget http://www.hepforge.org/archive/blackhat/blackhat-0.9.9.tar.gz
tar -xzf blackhat-0.9.9.tar.gz
cd blackhat-0.9.9
./configure --with-QDpath=/home/lm/disk/files/sw/mc/qd-2.3.22/install --enable-pythoninterface  --prefix=$PWD/install

nb: do not do https://sherpa.hepforge.org/trac/ticket/330
--enable-sherpaplugin=/home/lm/disk/files/sw/mc/sherpa
```

Following this a bunch of files need to be edited for C++11/C++17 and modern gcc compliance 

```
files modified: 
/home/lm/disk/files/sw/mc/blackhat-0.9.9/src/Interface/OptionsHandler.hpp --< make_pair -- pair 
 /home/lm/disk/files/sw/mc/blackhat-0.9.9/src/tree1.cc --> remove default arguements in template functions when needed 
 /home/lm/disk/files/sw/mc/blackhat-0.9.9/src/tree2.cc
 /home/lm/disk/files/sw/mc/blackhat-0.9.9/src/tree3.cc 
 /home/lm/disk/files/sw/mc/blackhat-0.9.9/src/util.cc
 /home/lm/disk/files/sw/mc/blackhat-0.9.9/src/treebase.cc
/home/lm/disk/files/sw/mc/blackhat-0.9.9/src/cut_Darren.h
<-- makes up it's own static_assert, which is a thing in C++11 ->
luckily this isn't used anywhere apart from in the .h
changed static_assert -> extatic_assert
<--- has a bunch of pure virtual functions in BH::cut::Darren::triangle_Darren< , but then wants to use it explicitly as a return type 
==> remove '=0' in these pure virtuals  witd empyy returns


```

Then, Sherpa: 


```

./configure --with-sqlite3=install --enable-hepmc2=/home/lm/disk/files/sw/mc/hepmc/install --enable-analysis --disable-silent-rules --enable-blackhat=/home/lm/disk/files/sw/mc/blackhat-0.9.9/install

make

make install

```
