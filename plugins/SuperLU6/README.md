# SuperLU plugin

**UG4-Plugin** implementing an interface for the [SuperLU solver](https://github.com/xiaoyeli/superlu)

Copyright 2014-2023 Goethe University Frankfurt

## Install
Please install/clone this repository through UG4's package manager[ughub](https://github.com/UG4/ughub).

The default is the *internal* SuperLU v6.0.0, which is included as a submodule. Switching to an *external* version can be achieved via cmake:

 ```
 cmake -DSuperLU6=ON [-DUSE_INTERNAL_SUPERLU=ON|OFF] ..
 ```


## Dependency
* [SuperLU solver (v6.0.0)](https://github.com/xiaoyeli/superlu) is included as a submodule. It can be initiated  and updated via git:
```
  git submodule init 
  git submodule update
 ```
* If you are using SuperLU, remember to cite it! Bibtex information is [here](http://crd-legacy.lbl.gov/~xiaoye/SuperLU/referencing.html).

## License
The plugin can be used in accordance with SuperLU license. 
