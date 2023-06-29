# \SCRIPT\-------------------------------------------------------------------------
#
#  CONTENTS      : download, build and install tools
#
#  DESCRIPTION   : none
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------
import os, sys, site

# defaults
if not 'threads_build' in config:
    config['threads_build'] = 1

rule default:
    shell : ""

rule all:
    input:
        # "bin/minimap2",
        # "bin/samtools",
        # "bin/NGmerge",
        "lib/python3.9/C3POa.py",
        "lib/python3.9/site-packages/conk/conk.cpython-39-x86_64-linux-gnu.so",
        "lib/python3.9/site-packages/medaka/maple_smolecule.py"

rule minimap2:
    output:
        bin = "bin/minimap2"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d minimap2 ]; then
            git clone https://github.com/lh3/minimap2 --branch v2.14 --depth=1 && cd minimap2
        else
            cd minimap2 && git fetch --all --tags --prune && git checkout tags/v2.14
        fi
        make clean && make -j{threads}
        cp minimap2 ../../{output.bin}
        """

rule htslib:
    output:
        src = directory("src/htslib")
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d htslib ]; then
            git clone https://github.com/samtools/htslib --branch 1.9 --depth=1 && cd htslib
        else
            cd htslib && git fetch --all --tags --prune && git checkout tags/1.9
        fi
        autoheader && autoconf && ./configure && make -j{threads}
        """

rule samtools:
    input:
        rules.htslib.output.src
    output:
        bin = "bin/samtools"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d samtools ]; then
            git clone https://github.com/samtools/samtools --branch 1.9 --depth=1 && cd samtools
        else
            cd samtools && git fetch --all --tags --prune && git checkout tags/1.9
        fi
        autoheader --warning=none && autoconf -Wno-syntax && ./configure && make -j{threads}
        cp samtools ../../{output.bin}
        """

rule NGmerge:
    output:
        bin = "bin/NGmerge"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d NGmerge ]; then
            git clone https://github.com/harvardinformatics/NGmerge --branch v0.3 --depth=1 && cd NGmerge
        else
            cd NGmerge && git fetch --all --tags --prune && git checkout tags/v0.3
        fi
        make clean && make -j{threads}
        cp NGmerge ../../{output.bin}
        """

rule C3POa:
    output:
        C3POa = "lib/python3.9/C3POa.py",
        conk = "lib/python3.9/site-packages/conk/conk.cpython-39-x86_64-linux-gnu.so"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ -d C3POa ]; then
            rm -r -f C3POa
        fi
        git clone https://github.com/christopher-vollmers/C3POa.git

        if [ -d conk ]; then
            rm -r -f conk
        fi
        git clone https://github.com/rvolden/conk && cd conk
        python setup.py sdist bdist_wheel
        python -m pip install dist/conk*whl --force-reinstall
        cd ../..
        
        if [ -d lib/python3.9/bin ]; then
            mv src/C3POa/bin/* lib/python3.9/bin/
        else
            mv src/C3POa/bin lib/python3.9
        fi
        mv src/C3POa/C3POa.py src/C3POa/C3POa_postprocessing.py lib/python3.9
        rm -r -f src/C3POa
        """

rule maple_medaka:
    output:
        ms = "lib/python3.9/site-packages/medaka/maple_smolecule.py"
    shell:
        """
        mkdir -p src && cd src
        if [ -d maple ]; then
            rm -r -f maple
        fi
        git clone https://github.com/gordonrix/maple.git && cd ..
        mv src/maple/rules/utils/maple_smolecule.py src/maple/rules/utils/medaka.py lib/python3.9/site-packages/medaka
        rm -r -f src/maple
        """
