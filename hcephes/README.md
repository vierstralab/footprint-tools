# hcephes

[![Travis](https://travis-ci.com/limix/hcephes.svg?branch=master)](https://travis-ci.com/limix/hcephes)

The [Netlib Cephes library](https://www.netlib.org/cephes/) reformatted for the 21st
century.

## Install

The easiest way of installing it is via [conda](https://conda.io/)

```bash
conda install -c conda-forge hcephes
```

Alternatively, one can compile and install it.
From Linux, MacOS, or Windows (bash terminal) systems, enter

```bash
# DO_CMD=sudo
curl -fsSL https://git.io/JerYI | GITHUB_USER=limix GITHUB_PROJECT=hcephes bash
```

## Usage

It requires you to provide the library name `hcephes` to your linker and the path
inclusion of its `hcephes.h` header to your compiler.
For example, suppose you on MacOS and you are using [gcc](https://www.gnu.org/software/gcc/).
A C file like

```c
/* example.c */
#include "hcephes.h"

#include <stdio.h>
#include <stdlib.h>

int main() {
    printf("%f\n", hcephes_bdtr(4, 5, 0.25));
    return 0;
}
```

might require the following command to create an executable file:

```bash
gcc example.c -lhcephes -I/usr/local/include -o example
```

For the complete list of the available functions, we refer the reader to
[include/hcephes.h](include/hcephes.h) file itself and to the [cephes library](https://www.netlib.org/cephes/)
documentation.

## CMake

Add the following to your `CMakeLists.txt`:

```
find_package(hcephes REQUIRED)

target_link_libraries(mylib PRIVATE HCEPHES::hcephes)
```

## Authors

* [Danilo Horta](https://github.com/horta)

## License

This project is licensed under the [MIT License](https://raw.githubusercontent.com/limix/hcephes/master/LICENSE.md).
