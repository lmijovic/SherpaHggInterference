# SherpaHggInterference

Fork with updates to get HiggsLHCWG HiggsInterference code running locally

## Prerequisites 

* HEPMC 

```

git clone https://github.com/hep-mirrors/hepmc.git
cd hepmc/
./bootstrap
./configure --prefix=$PWD/install --with-momentum=MEV --with-length=MM
make
make install


```

## Installation w/o blackhat: 

Installing sherpa without blackhat should work fine in order to produce the interference samples: 

replace $PATHTOHEPMC/install with your hepmcdir/install as per installation above

```
git clone https://github.com/lmijovic/SherpaHggInterference.git
cd SherpaHggInterference/hintmcnlo/
./configure --with-sqlite3=install --enable-hepmc2=$PATHTOHEPMC/install --enable-analysis --disable-silent-rules 

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

# Interference runs

The interference cards alla LHCHiggsWG, but with relaxed selectors are in: hintmcnlo/Examples/H_in_GluonFusion/ERF/ , see Run.dat - s for info on what gets generated in each directory 

```
cd hintmcnlo/Examples/H_in_GluonFusion/ERF/

```

The Selector sets: 
* HiggsFinder y1.pt,y2.pt,y.eta,myy-low,myy-high 
* IsolationCut photon Frixione isolation

```
(selector){
  HiggsFinder 40 30 2.5 100 150;
  IsolationCut 22 0.1 2 0.025;
}(selector);
```

# Run 

Replace $YOURSHERPADIR with where sherpa is , run with random seed 10


```
$YOURSHERPADIR/SherpaHggInterference/hintmcnlo/bin/Sherpa -R 10 | tee sh.out
```

The output is set to HepMC format hard process (no hadronization). 

Note: weighted interference has huge spread of weights (7 ord of mag). Thus, prefer unweighted production (INT_UNW/SIG_UNW) which is slower, but yields better stats coverage for differentail distibutions and ML training.
 

