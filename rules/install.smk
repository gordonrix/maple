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
import os, sys, site, platform

def get_arch():
    arch = platform.machine()
    if arch in ["x86_64", "i686"]:
        return "x86"
    elif arch == "arm64":
        return "arm64"
    else:
        raise ValueError(f"Unsupported architecture: {arch}")

# defaults
if not 'threads_build' in config:
    config['threads_build'] = 1

rule default:
    shell : ""

rule all:
    input:
        "bin/minimap2",
        "bin/samtools",
        "bin/NGmerge",
        "bin/umicollapse.jar",
        "bin/medaka"

rule minimap2:
    output:
        bin = "bin/minimap2"
    threads: config['threads_build']
    params:
        arch = get_arch()
    shell:
        """
        mkdir -p src && cd src
        if [ -d minimap2 ]; then
            rm -rf minimap2
        fi
        git clone https://github.com/lh3/minimap2 --depth=1 && cd minimap2

        if [ "{params.arch}" = "x86" ]; then
            make clean && make -j{threads}
        elif [ "{params.arch}" = "arm64" ]; then
            make clean && make -j{threads} arm_neon=1 aarch64=1
        fi

        cp minimap2 ../../{output.bin}
        cd ../
        rm -rf minimap2
        """

rule htslib:
    output:
        src = directory("src/htslib")
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src

        if [ -d htslib ]; then
            rm -rf htslib
        fi
        git clone https://github.com/samtools/htslib --depth=1 && cd htslib

        git submodule update --init --recursive
        autoreconf -i && ./configure && make -j{threads}

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

        if [ -d samtools ]; then
            rm -rf samtools
        fi
        git clone https://github.com/samtools/samtools --depth=1 && cd samtools
        
        autoreconf -i && ./configure && make -j{threads}

        cp samtools ../../{output.bin}
        cd ../
        rm -rf samtools
        """

rule NGmerge:
    output:
        bin = "bin/NGmerge"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src

        if [ -d NGmerge ]; then
            rm -rf NGmerge
        fi

        git clone https://github.com/harvardinformatics/NGmerge --depth=1 && cd NGmerge

        # Detect operating system and set CFLAGS for macOS
        if [[ "$OSTYPE" == "darwin"* ]]; then
            echo "Modifying CFLAGS for macOS"
            sed -i '' 's/^CFLAGS?.*$/CFLAGS?=-Xpreprocessor -fopenmp -lomp/' Makefile
        fi
        
        make clean && make -j{threads}
        cp NGmerge ../../{output.bin}
        cd ../
        rm -rf NGmerge
        """

rule C3POa:
    output:
        C3POa = "lib/python3.10/C3POa.py",
        conk = "lib/python3.10/site-packages/conk/conk.cpython-39-x86_64-linux-gnu.so"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ -d C3POa ]; then
            rm -r -f C3POa
        fi
        git clone -b peakFinderCustomSettings https://github.com/gordonrix/C3POa.git

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

# umicollapse can be installed via conda but uses an outdated version of snappy that results in an error. See https://github.com/Daniel-Liu-c0deb0t/UMICollapse/issues/27
# note
rule umicollapse:
    output:
        jar = "bin/umicollapse.jar"
    threads: config['threads_build']
    shell:
        r"""
        mkdir -p src && cd src
        if [ -d UMICollapse ]; then
            rm -rf UMICollapse
        fi
        git clone https://github.com/Daniel-Liu-c0deb0t/UMICollapse.git --depth=1 && cd UMICollapse

        # Create the lib directory and download dependencies
        mkdir -p lib && cd lib
        curl -O -L https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar
        # Download the updated snappy-java jar (version 1.1.10.7 instead of 1.1.7.3)
        curl -O -L https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.10.7/snappy-java-1.1.10.7.jar
        cd ..

        # Patch build.sh to reference the new snappy version
        sed -i.bak 's/snappy-java-1\.1\.7\.3.jar/snappy-java-1.1.10.7.jar/g' build.sh

        chmod +x build.sh && ./build.sh

        # Copy the resulting jar to the bin directory
        cp umicollapse.jar ../../{output.jar}
        cd ../
        rm -rf src
        """

# rule maple_medaka:
#     output:
#         ms = "lib/python3.9/site-packages/medaka/maple_smolecule.py"
#     shell:
#         """
#         mkdir -p src && cd src
#         if [ -d maple ]; then
#             rm -r -f maple
#         fi
#         git clone https://github.com/gordonrix/maple.git && cd ..
#         mv src/maple/rules/utils/maple_smolecule.py src/maple/rules/utils/medaka.py lib/python3.9/site-packages/medaka
#         rm -r -f src/maple
#         """


rule medaka_wrapper:
    output:
        wrapper="bin/medaka"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ -d medaka ]; then
            rm -rf medaka
        fi
        git clone https://github.com/nanoporetech/medaka.git --depth=1 && cd medaka
        make install

        # modify the medaka to use code from maple
        cd ..
        if [ -d maple ]; then
            rm -rf maple
        fi
        git clone https://github.com/gordonrix/maple.git
        mv maple/rules/utils/maple_smolecule.py maple/rules/utils/medaka.py medaka/venv/lib/python3.10/site-packages/medaka
        rm -rf src/maple

        # Create a directory to hold the Medaka environment.
        mkdir -p ../bin/medaka_env
        cp -r medaka/venv ../bin/medaka_env/
        # Create a wrapper that opens the env with the medaka command
        cat << 'EOF' > ../bin/medaka
#!/bin/bash
# Activate the Medaka virtual environment and execute medaka,
# passing along all provided command-line flags.
source "$(dirname "$0")/medaka_env/venv/bin/activate"
exec medaka "$@"
EOF
        chmod +x ../bin/medaka
        cd ..
        rm -rf medaka
        """
