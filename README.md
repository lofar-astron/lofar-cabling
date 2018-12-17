# Lofar cabling

A Python 3 project that generates cable layouts for future international LOFAR fields.
Documentation for this program is currently only available in the form of docstrings.

A very basic usage example is as follows:
```
import lofarcabling
layout = lofarcabling.go(90, (40, -40), "points.csv", "layout.csv", True)
````
90 is the amount of degrees the field is rotated clockwise \
(40, -40) are the x and y coordinates of the startpoint \
"points.csv": The filelocation where you want to save a .csv with the coordinates of all antennas \
"layout.csv": The filelocation where you want to save a .csv with the coordinates of the points and lines \
True / False: True if you want the "improved" version of the field, False if you just want a quick version (True takes ~1:30 hour and False takes ~ minute)

### Prerequisites

To run this project, you will need python 3 and some packages listed in [requirements.txt](https://github.com/lofar-astron/lofar-cabling/blob/master/requirements.txt) \
\
Python version 3.6 or greater is required

### Installing

To install, use `python setup.py install` (or see `python setup.py --help` for options).

## Running tests

Tests are in development still.

## Built with

[Anaconda](https://anaconda.org/anaconda/anaconda-navigator) - GUI used to launch applications \
[Jupyter Notebook](https://jupyter.org/) - The web application used

## Authors

* **Thomas Tijsma** - Realized the project 
* **Tammo Jan Dijkema** - Initial setup and support - [ASTRON](http://astron.nl/)

## License

See [LICENSE](https://github.com/lofar-astron/lofar-cabling/blob/master/LICENSE) for details