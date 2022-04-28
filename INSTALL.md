# Installing Cantera

## Basic Instructions

To compile using the default options, run `scons build` followed by
`scons install`.

Configuration options are specified with `name=value` on the command line, e.g.:
`scons build optimize=n prefix=/home/$USER/cantera`

The full list of configuration options and their default values can be shown by running
`scons help --options`. The list of available `scons` commands (e.g. `build`) can be
shown by running `scons help`.


## Detailed Instructions

See the instructions available [online](https://cantera.org/install/compiling-install.html).