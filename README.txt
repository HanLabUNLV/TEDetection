Dependencies
------------
- Python 2.6+
- Pysam
- NumPy
- SciPy
- Bamcollate2
- bwa

Install Instructions
--------------------
Pysam:
  To install pysam, you will need pip installed. If you have python 2.7.9+, pip
  will be automatically installed. If not, install pip by:
    wget https://bootstrap.pypa.io/get-pip.py
    sudo python get-pip.py

  Then install pysam by:
    sudo pip install pysam

NumPy/SciPy:
Use your distributions package manager to install these two:
  Ubuntu/Debian:
    sudo apt-get install python-numpy python-scipy
  Fedora/RHEL:
    sudo yum install numpy scipy

Bamcollate2:
Bamcollate2 is part of a software suite called biobambam. To install
biobambam, you need to first install libmaus. Install libmaus by:
  git clone https://github.com/gt1/libmaus.git
  cd libmaus
  autoreconf -i -f
  ./configure
  make
  sudo make install
Then install biobambam by:
  git clone https://github.com/gt1/biobambam.git
  cd biobambam
  autoreconf -i -f
  ./configure --with-libmaus=/usr/local (if you installed libmaus with default configuration, as above)
  make
  sudo make install

bwa:
Install bwa by:
  git clone https://github.com/lh3/bwa.git
  cd bwa
  make
Then you can either add this folder to your PATH environment variable by
opening ~/.bash_profile or ~/.bashrc and adding the lines:
  PATH=$HOME/bwa:$PATH
  export PATH
Or you can copy the bwa executable to a folder already in the path:
  sudo cp bwa /usr/local/bin
