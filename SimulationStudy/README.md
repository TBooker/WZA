# Simulations to validate the WZA

This repository contains the code and SLiM configuration files to run and analyse the simulations we devised to validate the WZA.

The SLiM config files are stored in [slim_configs/](slim_configs/).

These configuration files are designed to be run from the command line. I use servers that require the use of scheduling software, so it's a little different to a situation where you have full control of a server or machine.

The [bin/](bin/) directory contains all the code (mostly Python, with some R) to perform the analyses we conducted on the simulation data.

The shell scripts in [shell_scripts/](shell_scripts/) run the analyse on all the simulation replicates we had data for.  
