Environmentally-driven population dynamics model of Aedes albopictus
====================================================================

This is a python (v2.7) package implementing the environmentally-driven population dynamics model of Aedes albopictus. 

Contents
--------

1) Prerequisites
2) Linux installation
3) Usage

1) Prerequisites
----------------

The model requires the following packages, which are not included in this package:

	numpy
    
	pkg_resources

2) Linux installation
---------------------

1) Easy way:

If you have pip installed, you can use the following command to download and install the package.
	pip install albopictus

Alternatively, you can download the source code from PyPI and run pip on the latest version xxx.
	pip install albopictus-xxx.tar.gz

2) Hard way:

If pip is not available, you can unpack the package contents and perform a manual install.
	tar -xvzf albopictus-xxx.tar.gz
    
	cd albopictus-xxx
    
	python setup.py install

This will install the package in the site-packages directory of your python distribution. If you do not have root privileges or you wish to install to a different directory, you can use the --prefix argument.

	python setup.py install --prefix=<dir>

In this case, please make sure that <dir> is in your PYTHONPATH, or you can add it with the following command.

In bash shell:
	export PYTHONPATH=<dir>:$PYTHONPATH
In c shell:
	setenv PYTHONPATH <dir>:$PYTHONPATH

Usage
-----

Information for usage and contents are documented in the package, and can be accessed with the help utility.

	import albopictus
    
	help(albopictus)

Credits
-------

'modern-package-template' - http://pypi.python.org/pypi/modern-package-template
