#!/bin/bash
echo "Installing third party libraries"
# List of usefull colors
COLOR_RESET="\033[0m"
COLOR_INFO="\033[0;32m"
COLOR_ITEM="\033[1;34m"
COLOR_QUES="\033[1;32m"
COLOR_WARN="\033[0;33m"
COLOR_BOLD="\033[1m"
COLOR_UNDE="\033[4m"

FWP_INSTALL_PREFIX=/usr/local
COMMON_INSTALL_PREFIX=/usr

# Getting the current directory of this script
CURRENT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


## This function detects the current os and distro
CURRENT_OS="UBUNTU" #CENTOS, UBUNTU are other valid options
echo CURRENT_OS
function install_glpk
{
	# Remove old folder (sanity procedure)
	cd $CURRENT_DIR/thirdparty
	sudo rm -rf glpk
	sudo rm -rf glpk-4.42.tar.gz

	if  [ "$CURRENT_OS" == "UBUNTU" ]; then
		# Getting Eigen 3.2.10
		echo "install GLPK"
		wget http://ftp.gnu.org/gnu/glpk/glpk-4.42.tar.gz
		tar -zxvf glpk-4.42.tar.gz
		rm -rf glpk-4.42.tar.gz
		cd glpk-4.42
		./configure --enable-dl --enable-shared
		sudo make -j
		sudo make install
    fi
}

## install_glpk

function install_politopix
{
	# Remove old folder (sanity procedure)
	cd $CURRENT_DIR/thirdparty/politopix
	rm -rf build

	if  [ "$CURRENT_OS" == "UBUNTU" ]; then
		# Getting Eigen 3.2.10
		echo "install Politopix"
		mkdir build && cd build
		cmake ..
		sudo make
    fi
}

install_politopix