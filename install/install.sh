#!/bin/bash
# Maple installation script
#
# This script automates the full installation of the Maple pipeline:
#   1. Creates a conda environment with all Python dependencies
#   (from pinned requirements in a lock file or requirements.yaml)
#   2. Activates the environment
#   3. Builds and installs additional tools (minimap2, samtools, medaka, etc.) via snakemake
#
# Usage: source install/install.sh [OPTIONS]
#
# Options:
#   -n, --name NAME    Environment name (default: from requirements.yaml)
#   -u, --unpinned     Use unpinned requirements.yaml instead of lock file
#   -f, --force        Remove existing environment before installing
#   -j, --jobs N       Number of parallel jobs for tool installation (default: 4)
#   -h, --help         Show this help message

# Get script directory (works in both bash and zsh, when sourced)
if [ -n "${BASH_SOURCE[0]}" ] 2>/dev/null; then
    # Bash
    _MAPLE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
elif [ -n "$0" ] && [ "$0" != "-zsh" ] && [ "$0" != "zsh" ] && [ "$0" != "-bash" ] && [ "$0" != "bash" ]; then
    # Zsh (when sourced, $0 is the script path)
    _MAPLE_SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
else
    # Fallback: try to find install directory relative to pwd
    if [ -f "install/requirements.yaml" ]; then
        _MAPLE_SCRIPT_DIR="$(pwd)/install"
    elif [ -f "requirements.yaml" ]; then
        _MAPLE_SCRIPT_DIR="$(pwd)"
    else
        echo "Error: Could not determine script directory."
        echo "Please run from the maple directory: source install/install.sh"
        return 1 2>/dev/null || exit 1
    fi
fi

# Validate that we found the right directory
if [ ! -f "${_MAPLE_SCRIPT_DIR}/requirements.yaml" ]; then
    echo "Error: requirements.yaml not found in ${_MAPLE_SCRIPT_DIR}"
    echo "Please run this script from the maple directory using: source install/install.sh"
    return 1 2>/dev/null || exit 1
fi

# Parse environment name from requirements.yaml
_MAPLE_DEFAULT_ENV_NAME=$(grep "^name:" "${_MAPLE_SCRIPT_DIR}/requirements.yaml" 2>/dev/null | head -1 | sed 's/name: *//')
if [ -z "${_MAPLE_DEFAULT_ENV_NAME}" ]; then
    _MAPLE_DEFAULT_ENV_NAME="maple"
fi

# Defaults
_MAPLE_ENV_NAME="${_MAPLE_DEFAULT_ENV_NAME}"
_MAPLE_USE_UNPINNED=false
_MAPLE_FORCE=false
_MAPLE_JOBS=4

# Parse arguments
while [ $# -gt 0 ]; do
    case $1 in
        -n|--name)
            _MAPLE_ENV_NAME="$2"
            shift 2
            ;;
        -u|--unpinned)
            _MAPLE_USE_UNPINNED=true
            shift
            ;;
        -f|--force)
            _MAPLE_FORCE=true
            shift
            ;;
        -j|--jobs)
            _MAPLE_JOBS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Maple installation script"
            echo ""
            echo "Usage: source install/install.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -n, --name NAME    Environment name (default: ${_MAPLE_DEFAULT_ENV_NAME})"
            echo "  -u, --unpinned     Use unpinned requirements.yaml instead of lock files"
            echo "  -f, --force        Remove existing environment before installing"
            echo "  -j, --jobs N       Number of parallel jobs for tool installation (default: 4)"
            echo "  -h, --help         Show this help message"
            unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS
            return 0 2>/dev/null || exit 0
            ;;
        *)
            echo "Unknown option: $1"
            unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS
            return 1 2>/dev/null || exit 1
            ;;
    esac
done

# Check conda is available
if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda not found. Please install conda first."
    unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS
    return 1 2>/dev/null || exit 1
fi

_MAPLE_CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "${_MAPLE_CONDA_BASE}" ]; then
    echo "Error: Could not determine conda base directory."
    unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE
    return 1 2>/dev/null || exit 1
fi

_MAPLE_ENV_PATH="${_MAPLE_CONDA_BASE}/envs/${_MAPLE_ENV_NAME}"

# Check if environment already exists
if [ -d "${_MAPLE_ENV_PATH}" ]; then
    if [ "$_MAPLE_FORCE" = true ]; then
        echo "Removing existing environment '${_MAPLE_ENV_NAME}'..."
        if ! conda env remove -n "${_MAPLE_ENV_NAME}" -y; then
            echo "Error: Failed to remove existing environment."
            unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH
            return 1 2>/dev/null || exit 1
        fi
    else
        echo "Error: Environment '${_MAPLE_ENV_NAME}' already exists at ${_MAPLE_ENV_PATH}"
        echo "Use --force to remove and reinstall, or --name to use a different name."
        unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH
        return 1 2>/dev/null || exit 1
    fi
fi

# Check for mamba, bootstrap if needed
_MAPLE_MAMBA_BOOTSTRAP="${_MAPLE_CONDA_BASE}/envs/mamba-bootstrap"
if command -v mamba >/dev/null 2>&1; then
    _MAPLE_CONDA_CMD="mamba"
elif [ -x "${_MAPLE_MAMBA_BOOTSTRAP}/bin/mamba" ]; then
    echo "Using mamba from bootstrap environment..."
    _MAPLE_CONDA_CMD="${_MAPLE_MAMBA_BOOTSTRAP}/bin/mamba"
else
    echo "mamba not found. Installing mamba-bootstrap environment..."
    if ! conda create -n mamba-bootstrap -c conda-forge mamba -y; then
        echo "Error: Failed to create mamba-bootstrap environment."
        unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP
        return 1 2>/dev/null || exit 1
    fi
    _MAPLE_CONDA_CMD="${_MAPLE_MAMBA_BOOTSTRAP}/bin/mamba"
fi

echo "Installing maple environment '${_MAPLE_ENV_NAME}'..."

if [ "$_MAPLE_USE_UNPINNED" = true ]; then
    echo "Using unpinned requirements.yaml..."
    if ! ${_MAPLE_CONDA_CMD} env create -f "${_MAPLE_SCRIPT_DIR}/requirements.yaml" -n "${_MAPLE_ENV_NAME}" -y; then
        echo "Error: Failed to create environment from requirements.yaml"
        unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP _MAPLE_CONDA_CMD
        return 1 2>/dev/null || exit 1
    fi
else
    # Detect platform and select lock file
    _MAPLE_PLATFORM="$(uname -s)-$(uname -m)"
    case "${_MAPLE_PLATFORM}" in
        Linux-x86_64)  _MAPLE_LOCK="conda-linux-64.lock" ;;
        Darwin-x86_64) _MAPLE_LOCK="conda-osx-64.lock" ;;
        Darwin-arm64)  _MAPLE_LOCK="conda-osx-arm64.lock" ;;
        *)
            echo "Unsupported platform: ${_MAPLE_PLATFORM}"
            unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP _MAPLE_CONDA_CMD _MAPLE_PLATFORM
            return 1 2>/dev/null || exit 1
            ;;
    esac

    if [ ! -f "${_MAPLE_SCRIPT_DIR}/${_MAPLE_LOCK}" ]; then
        echo "Error: Lock file ${_MAPLE_LOCK} not found in ${_MAPLE_SCRIPT_DIR}"
        echo "Use --unpinned to install from requirements.yaml instead."
        unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP _MAPLE_CONDA_CMD _MAPLE_PLATFORM _MAPLE_LOCK
        return 1 2>/dev/null || exit 1
    fi

    echo "Using pinned ${_MAPLE_LOCK}..."
    if ! ${_MAPLE_CONDA_CMD} create -n "${_MAPLE_ENV_NAME}" --file "${_MAPLE_SCRIPT_DIR}/${_MAPLE_LOCK}" -y; then
        echo "Error: Failed to create environment from lock file"
        unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP _MAPLE_CONDA_CMD _MAPLE_PLATFORM _MAPLE_LOCK
        return 1 2>/dev/null || exit 1
    fi
fi

echo "Activating environment..."
eval "$(conda shell.bash hook)" 2>/dev/null || eval "$(conda shell.zsh hook)" 2>/dev/null
if ! conda activate "${_MAPLE_ENV_NAME}"; then
    echo "Error: Failed to activate environment"
    unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP _MAPLE_CONDA_CMD _MAPLE_PLATFORM _MAPLE_LOCK
    return 1 2>/dev/null || exit 1
fi

echo "Installing additional tools (minimap2, samtools, medaka, etc.)..."
if ! snakemake --snakefile "${_MAPLE_SCRIPT_DIR}/../rules/install.smk" --directory "${_MAPLE_ENV_PATH}" -j "${_MAPLE_JOBS}" all; then
    echo "Error: Failed to install additional tools"
    unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP _MAPLE_CONDA_CMD _MAPLE_PLATFORM _MAPLE_LOCK
    return 1 2>/dev/null || exit 1
fi

echo ""
echo "Installation complete! Environment '${_MAPLE_ENV_NAME}' is now active."
echo "For future sessions, activate with: conda activate ${_MAPLE_ENV_NAME}"

# Clean up variables
unset _MAPLE_SCRIPT_DIR _MAPLE_DEFAULT_ENV_NAME _MAPLE_ENV_NAME _MAPLE_USE_UNPINNED _MAPLE_FORCE _MAPLE_JOBS _MAPLE_CONDA_BASE _MAPLE_ENV_PATH _MAPLE_MAMBA_BOOTSTRAP _MAPLE_CONDA_CMD _MAPLE_PLATFORM _MAPLE_LOCK
