Dependencies
------------
- Python 2.6+
- Samtools
- Pysam
- NumPy
- SciPy
- Bamcollate2
- bwa

Install Instructions
--------------------
Samtools:
  Samtools requires ncurses to be installed for some of its features (although
  not used by this program). To install ncurses, use your distributions
  package manager:
    Ubuntu/Debian:
      sudo apt-get install ncurses-dev
    Fedora/RHEL:
      sudo yum install ncurses-devel
  Then install samtools by:
    wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
    tar -xjf samtools-1.2.tar.bz2
    cd samtools-1.2
    make
    sudo make install
  If you get an error similar to this:
    bam_tview_curses.o: undefined reference to symbol 'keypad'
  Open up samtools-1.2/Makefile, find the line that says:
    LDFLAGS = 
  and change it to:
    LDFLAGS = -ltinfo

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
        - If you get any errors, make sure these packages are installed according to your distribution:
            Ubuntu/Debian:
                sudo apt-get install autoconf automake libtool
            Fedora/RHEL:
                sudo yum install autoconf automake libtool
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
    sudo cp bwa /usr/local/bin/

Running the Program
-------------------
Only 3 of the arguments are necessary to run the program:
  -sd Path to the directory of source scripts (included with the program in SourceScripts/)
  -tr Path to the TE reference fasta (included with the program in TERef/)
  -bf Path to the .bam file to run the analysis on
All other arguments are optional or have predefined default values the program
will use. The bamfile must be sorted and indexed. If it is not already, you
can do so using samtools:
  samtools sort -O bam -T sorted file.bam
  samtools index file.bam
The TE reference must also be index by bwa. The included reference is indexed,
but if you would like to use your own, index it by:
  bwa index TEreference.fasta
A basic run of the program would look like the following:
./main.py -sd SourceScripts -tr TERef/repBase_Homo_sapiens_all_TE.fasta -bf cancer_sample.bam
