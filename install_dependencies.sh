cd
git clone https://github.com/susannvorberg/CCmpredPy
cd CCmpredPy
python3 setup.py install
cd
git clone https://github.com/soedinglab/hh-suite
mkdir -p hh-suite/build && cd hh-suite/build
yes | apt-get install cmake
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j 4 && make install
echo "export PATH="/root/hh-suite/build/bin:/root/hh-suite/build/scripts:$PATH"" > path_tmp
cat ~/.bashrc path_tmp > tmp
rm path_tmp
mv tmp ~/.bashrc
cd
git clone https://github.com/scapella/trimal
cd trimal/source
make
echo "export PATH="/root/trimal/source:$PATH"" > path_tmp
cat ~/.bashrc path_tmp > tmp
mv tmp ~/.bashrc
rm path_tmp
cd
. ~/.bashrc
