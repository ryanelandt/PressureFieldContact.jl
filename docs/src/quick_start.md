# Quick-Start

This guide assumes that you are using Ubuntu Linux.

## Download Julia

Download the current stable release for your operating system (e.g., Linux x86 64-bit) from [this website](https://julialang.org/downloads/).
When the download is complete, copy the file into your home folder; and then extract the file.

## Open Julia

Extracting the file should have created a folder with a name similar to "julia-1.1.1" in your home folder.
Open this folder, and then open the "bin" folder inside of it.
Right click inside the "bin" folder and select "Open in Terminal".
Type "./julia" into the terminal and then press enter; this is the command line interface to Julia.

## Install PressureFieldContact

Install `PressureFieldContact` by typing the commands below into Julia's command line interface.

```julia
using Pkg
Pkg.develop("GenericLinearAlgebra")  # Do this even though the package is registered
Pkg.add("PressureFieldContact")
```

## Verify Installation

You can verify your installation of `PressureFieldContact` by typing the commands below into Julia's command line interface.
The test takes about 5 minutes.
If you have waited 5 minutes and the test hangs, make sure that you "developed" `GenericLinearAlgebra` as shown in the "Install PressureFieldContact" section of this guide.

```julia
using Pkg
Pkg.test("PressureFieldContact")
```

## Run example

Examples are currently in the test folder.
Run the boxes example first to get an idea of how geometry is created.
Then run the pencil example for a dexterous manipulation example.

## IDE (ATOM)

Atom is an IDE for languages including Julia.
Install Atom by following the instructions [here](http://docs.junolab.org/latest/man/installation/).
Be sure to install Atom after you completed the "Install PressureFieldContact" section of this quick-start guide.
If you do not do this, there may/will be a problem with the version of a package called `WebIO`, and `PressureFieldContact` will not install.
