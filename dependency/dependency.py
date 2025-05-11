#!/usr/bin/env python3

import subprocess
import platform
import sys
import os

# Define system package requirements
COMMON_PACKAGES = [
    "build-essential", "libcurl4-openssl-dev", "libssl-dev",
    "libxml2-dev", "libgit2-dev", "libharfbuzz-dev", "libfribidi-dev",
    "libfontconfig1-dev", "libfreetype6-dev", "libpng-dev", "libtiff5-dev",
    "libjpeg-dev", "gfortran", "curl", "git", "plink2", "bcftools", "tabix"
]

def run_cmd(cmd):
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"Failed to run: {' '.join(cmd)}")
        sys.exit(1)

def install_on_debian():
    print("Detected Debian/Ubuntu system")
    run_cmd(["sudo", "apt-get", "update"])
    run_cmd(["sudo", "apt-get", "install", "-y"] + COMMON_PACKAGES)

def install_on_redhat():
    print("Detected RedHat/CentOS/Fedora system")
    run_cmd(["sudo", "dnf", "install", "-y"] + COMMON_PACKAGES)

def install_on_arch():
    print("Detected Arch Linux system")
    run_cmd(["sudo", "pacman", "-Sy", "--noconfirm"] + COMMON_PACKAGES)

def detect_and_install():
    distro = platform.freedesktop_os_release().get("ID", "").lower()

    if distro in ["ubuntu", "debian"]:
        install_on_debian()
    elif distro in ["fedora", "centos", "rhel"]:
        install_on_redhat()
    elif distro in ["arch", "manjaro"]:
        install_on_arch()
    else:
        print(f"Unsupported distro: {distro}")
        sys.exit(1)

if __name__ == "__main__":
    detect_and_install()
    os.system("Rscript dependency.R")
    print("Dependency installation complete!")