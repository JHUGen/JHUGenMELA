# JHUGenMELA          {#mainpage}

@authors I. Anderson, S. Bolognesi, F. Caola, J. Davis, Y. Gao, A. V. Gritsan, 
@authors L. S. Mandacaru Guerra, Z. Guo, L. Kang, S. Kyriacou, C. B. Martin, T. Martini, 
@authors K. Melnikov, R. Pan, M. Panagiotou, R. Rontsch, J. Roskes, U. Sarica, 
@authors M. Schulze, M. V. Srivastav, N. V. Tran, A. Whitbeck, M. Xiao, Y. Zhou

This is the documentation for JHUGenMELA, a package that produces Matrix Element Calculations as used in JHUGen and MCFM-JHUGen.
The main user-interface in C++ can be found at the documentation for the Mela class, while the Python can be found in the @ref PyMela section.

Documentation written by Mohit Srivastav in 2024.

## Installation ##  {#install_sec}

Compile using

```{.sh}
./setup.sh -j N # (N being the number of cores to compile, or blank for max. allowed)
```

* The script decides whether to use integration with the experiments' software.
* Either the one in this folder, or that in the main JHUGenMELA folder works; they do the same thing.
* Feedback on how this script works in different environments is appreciated.

Once you compile, make sure to follow the on-screen instructions to set up all extra environment variables, i.e., run

```{.sh}
eval $(./setup.sh env) # or equivalents in non-bash shells
```

## Previous Tutorials  {#prev_tutorials}

There are a few tutorials lying around for MELA. They are listed below:

* [LPC EFT Workshop Notre Dame](https://indico.cern.ch/event/1378665/timetable/?view=standard#30-mela-tools)
* [LPC Offshell Workshop at Fermilab](https://indico.cern.ch/event/1375252/timetable/#8-mc-generators-2-mela-tools)