# Installation

## Install Julia

The recommended way to install and manage Julia is [`juliaup`](https://github.com/JuliaLang/juliaup), the official Julia version manager. It allows different Julia versions to be installed side by side and makes it easy to update Julia to the latest stable release.

### Linux

On Linux, install `juliaup` from a terminal using:

```bash
curl -fsSL https://install.julialang.org | sh
```

After the installation, restart the terminal and verify the Julia installation using:

```bash
julia --version
```

The installed Julia versions and available channels can be inspected with:

```bash
juliaup status
```

The current stable release can be updated using:

```bash
juliaup update release
```

### Windows

On Windows, Julia and `juliaup` can be installed from the Microsoft Store or from a terminal using:

```powershell
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
```

After the installation, open a new terminal and verify the installation using:

```powershell
julia --version
```

The installed Julia versions can be inspected and updated using:

```powershell
juliaup status
juliaup update release
```

For additional installation options, see the [official Julia installation instructions](https://docs.julialang.org/en/v1/manual/installation/).

## Install GeoModBox.jl

For using `GeoModBox.jl` as a Julia package, installation through the Julia package registry is recommended.
For working with the examples, exercises, and course material, clone the complete GitHub repository and activate the corresponding Julia environment.

### Install from the Julia package registry

Once `GeoModBox.jl` is registered in the Julia General Registry, the package can be installed directly from the Julia package manager.

Start Julia and enter the package manager by pressing `]`:

```julia
pkg> add GeoModBox
```

Alternatively, from the Julia REPL:

```julia
using Pkg
Pkg.add("GeoModBox")
```

The package can then be loaded using:

```julia
using GeoModBox
```

Individual submodules can be loaded as required, for example:

```julia
using GeoModBox.HeatEquation.TwoD
```

This installation method is recommended when only the `GeoModBox.jl` package and its numerical routines are required.

### Install from GitHub

To obtain the complete repository, including the examples, exercises, and documentation, clone the GitHub repository:

```bash
git clone https://github.com/GeoSci-FFM/GeoModBox.jl.git
cd GeoModBox.jl
```

Activate the main `GeoModBox.jl` environment and install the required package dependencies using:

```bash
julia --project=.
```

and then, from the Julia REPL:

```julia
using Pkg
Pkg.instantiate()
```

Alternatively, the unregistered development version of `GeoModBox.jl` can be installed directly from GitHub using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/GeoSci-FFM/GeoModBox.jl")
```

For working with the examples and exercises, cloning the complete repository is recommended.

> **User and developer workflows:**  
> For regular users and students, the environments described below can be instantiated directly and will use the registered version of `GeoModBox.jl`. Developers who modify the `GeoModBox.jl` source code should instead link the `examples` and `exercises` environments to their local checkout using `Pkg.develop`. This ensures that examples and exercises are executed with the locally modified version of `GeoModBox.jl` rather than the registered release.

## Julia environments

`GeoModBox.jl` uses separate Julia environments for the core package, examples, exercises, and documentation. This keeps the dependencies of the main package small while providing the additional packages required for plotting, analysis, exercises, and documentation.

### Users and students

When running an example, activate the `examples` environment from the repository root:

```julia
using Pkg
Pkg.activate("examples")
Pkg.instantiate()
```

or start Julia directly with:

```bash
julia --project=examples
```

For the exercises, use the corresponding `exercises` environment:

```julia
using Pkg
Pkg.activate("exercises")
Pkg.instantiate()
```

or:

```bash
julia --project=exercises
```

### Developers

When developing `GeoModBox.jl`, the example and exercise environments should use the local version of the package. From the repository root, activate the respective environment and link the local package using:

```julia
using Pkg
Pkg.activate("examples")
Pkg.develop(path=pwd())
Pkg.instantiate()
```

For the exercises, use:

```julia
using Pkg
Pkg.activate("exercises")
Pkg.develop(path=pwd())
Pkg.instantiate()
```

This links the respective environment to the `GeoModBox.jl` source code in the repository root. Changes made to the local package can therefore be tested directly with the examples and exercises.

### VS Code

In VS Code, the active Julia environment can be selected from the Julia environment indicator in the status bar. Select the `examples` or `exercises` environment, respectively, **before starting the Julia REPL**.

The documentation uses a separate environment located in `docs/`.