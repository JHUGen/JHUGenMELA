# JHUGenMELA

## Instruction for using JHUGenMELA/MELA

Compile MELA using

``` bash
./setup.sh -j N # (N being the number of cores to compile, or blank for max. allowed)
```

* The script decides whether to use integration with the experiments' software.
* Either the one in this folder, or that in the main JHUGenMELA folder works; they do the same thing.
* Feedback on how this script works in different environments is appreciated.

Once you compile, make sure to follow the on-screen instructions to set up all extra environment variables, i.e., run

``` bash
eval $(./setup.sh env) # or equivalents in non-bash shells
```

## Previous Tutorials

There are a few tutorials lying around for MELA. They are listed below:

* [LPC EFT Workshop Notre Dame](https://indico.cern.ch/event/1378665/timetable/?view=standard#30-mela-tools)
* [LPC Offshell Workshop at Fermilab](https://indico.cern.ch/event/1375252/timetable/#8-mc-generators-2-mela-tools)
