# BMI203 HW3 Skeleton

[![Build
Status](https://travis-ci.org/cjmathy/bmi203_hw3.svg?branch=master)](https://travis-ci.org/cjmathy/bmi203_hw3)



## structure

CHANGE THIS The main file that you will need to modify is `cluster.py` and the corresponding `test_cluster.py`. `utils.py` contains helpful classes that you can use to represent Active Sites. `io.py` contains some reading and writing files for interacting with PDB files and writing out cluster info.

```
.
├── README.md
│   ...
├── bmi203_hw2
│   ├── __init__.py
│   ├── __main__.py
│   ├── cluster.py
│   ├── io.py
│   └── utils.py
└── test
    ├── test_cluster.py
    └── test_io.py
```

## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. Then the package's
main function (located in `bmi203_hw3/__main__.py`) can be run as
follows

```
python -m bmi203_hw3
```

## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.


## contributors

Original design by Scott Pegg. Refactored and updated by Tamas Nagy. Repurposed for HW3 by Chris Mathy.
