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
import os, sys, site, platform, sysconfig

def get_arch():
    arch = platform.machine()
    if arch in ["x86_64", "i686"]:
        return "x86"
    elif arch == "arm64":
        return "arm64"
    else:
        raise ValueError(f"Unsupported architecture: {arch}")

PYVER = "{0}.{1}".format(*__import__("sys").version_info[:2])
EXT_SUFFIX = sysconfig.get_config_var("EXT_SUFFIX")

print(f"Performing install for:\n    Python version: {PYVER}\n    architecture: {get_arch()}\n    extension suffix: {EXT_SUFFIX}")

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
        "bin/medaka",
        f"lib/python{PYVER}/C3POa.py"

rule minimal:
    input:
        "bin/minimap2",
        "bin/samtools",
        "bin/NGmerge",
        "bin/umicollapse.jar",
        f"lib/python{PYVER}/C3POa.py"

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

rule abPOA:
    output:
        bin = "bin/abPOA"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        rm -rf abPOA
        git clone --recursive https://github.com/yangao07/abPOA.git
        cd abPOA; make; ./bin/abpoa ./test_data/seq.fa > cons.fa
        cp ./bin/abpoa ../../{output.bin}
        """


    # shell:
    #     """
    #     mkdir -p src && cd src
    #     rm -rf abPOA
    #     wget https://github.com/yangao07/abPOA/releases/download/v1.5.3/abPOA-v1.5.3.tar.gz
    #     tar -zxvf abPOA-v1.5.3.tar.gz && cd abPOA-v1.5.3
    #     make; ./bin/abpoa ./test_data/seq.fa > cons.fa
    #     cp ./bin/abPOA ../../{output.bin}
    #     """

rule C3POa:
    input:
        "bin/abPOA"
    output:
        C3POa = f"lib/python{PYVER}/C3POa.py",
        conk  = f"lib/python{PYVER}/site-packages/conk/conk{EXT_SUFFIX}"
    params:
        pyver      = PYVER,
        ext_suffix = EXT_SUFFIX
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        rm -rf C3POa conk
        
        # clone & build
        git clone https://github.com/gordonrix/C3POa.git
        git clone https://github.com/rvolden/conk && cd conk
        python setup.py sdist bdist_wheel
        python -m pip install dist/conk*whl --force-reinstall
        cd ../..
        
        if [ -d lib/python{params.pyver}/bin ]; then
            mv src/C3POa/bin/* lib/python{params.pyver}/bin/
        else
            mv src/C3POa/bin lib/python{params.pyver}/
        fi
        mv src/C3POa/C3POa.py src/C3POa/C3POa_postprocessing.py lib/python{params.pyver}
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

# perform a build from source as that is required for the maple overrides, then construct a wrapper and move it to bin/
# also avoids git LFS because that is dependent on nanoporetech account not going over the limit
rule medaka_wrapper:
    output:
        wrapper="bin/medaka"
    log: "logs/medaka_wrapper.log"
    threads: config['threads_build']
    shell:
        """
        mkdir -p logs    # make sure the logs/ dir exists
        exec > {log} 2>&1

        echo ">>> Working dir: $(pwd)"
        mkdir -p src && cd src

        echo ">>> Skipping LFS smudge"
        # 1) Clone the repo without pulling down LFS objects
        export GIT_LFS_SKIP_SMUDGE=1
        rm -rf medaka
        git clone https://github.com/nanoporetech/medaka.git --depth=1 --branch=v2.1.1 medaka && cd medaka

        echo ">>> Detecting Medaka version"
        VERSION_TAG=$(git describe --tags --abbrev=0)
        echo ">>> VERSION_TAG=$VERSION_TAG"
        VERSION=${{VERSION_TAG#v}}
        echo ">>> VERSION (no leading v) = $VERSION"

        echo ">>> Downloading the Medaka wheel (no install)…"
        TMP=$(mktemp -d)
        pip download medaka==$VERSION --no-deps -d "$TMP" \
        2>&1 | sed 's/^/   pip: /'

        echo ">>> Unpacking only the data/ folder from the wheel"
        WHEEL=$(ls "$TMP"/medaka-"$VERSION"-*.whl)
        unzip -q "$WHEEL" -d "$TMP/pkg"

        echo ">>> Grafting real model blobs into the Python package data/ folder"
        # make sure the target folder exists
        mkdir -p medaka/data
        cp -v "$TMP/pkg/medaka/data/"* medaka/data/ | sed 's/^/   cp: /'

        echo ">>> Removing any remaining LFS‑pointer stubs (tiny text files)"
        # LFS pointers start with "version https://git-lfs.github.com/spec/v1"
        grep -rl '^version https://git-lfs.github.com/spec/v1' medaka/data/ \
        | xargs rm -v | sed 's/^/   rm: /' || true

        echo ">>> medaka/data/ now contains only real .tar.gz blobs:"
        ls -1 medaka/data/ | sed 's/^/   /'

        echo ">>> Cleaning up temp wheel dir"
        rm -rf "$TMP"

        echo ">>> Finally, running make install (check_lfs will now pass)"
        make install

        # apply maple overrides
        cd ..
        rm -rf maple
        git clone --branch main --single-branch https://github.com/gordonrix/maple.git
        mv \
          maple/rules/utils/maple_smolecule.py \
          maple/rules/utils/medaka.py \
          medaka/venv/lib/python3.11/site-packages/medaka/
        rm -rf maple

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
