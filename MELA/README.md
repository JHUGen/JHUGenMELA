# Instruction for using JHUGenMELA/MELA

Compile using

```
./setup.sh -j N # (N being the number of cores to compile, or blank for max. allowed)
```

* The script decides whether to use integration with the experiments' software.
* Either the one in this folder, or that in the main JHUGenMELA folder works; they do the same thing.
* Feedback on how this script works in different environments is appreciated.

Once you compile, make sure to follow the on-screen instructions to set up all extra environment variables, i.e., run

```
eval $(./setup.sh env) # or equivalents in non-bash shells
```
