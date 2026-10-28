# Using Conda Environments

Conda is a clean way to manage dependencies for software projects.

## General Files

The simplest way to set up the dependencies for DGSWEM is to load the necessary compilers and then add any desired features with the update command.

```
conda env create -f environment-gcc.yml
conda env update -f feat-netcdf.yml
```

## Pinned Files

Pinned environment files reproduce a particular environment down to the build hash. These are specified by compiler, feature, and computer architecture.

```
conda env create -n dgswem-env -f environment-gcc-linux-64.pinned.yml
```

## Locked Files

Conda lock files are the most exacting, fully bypassing the conda solver. These require `conda-lock` to load.

```
conda-lock install --name dgswem-env envs/environment-gcc.lock.yml
```