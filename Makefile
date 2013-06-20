msg:
	@echo run python setup.py [install-directory]
	@echo or make install to install the ase - espresso interface

install:
	python setup.py install

EXTRAFILES = espsite.py.example

py = python
esp = $(instdir)/espresso
doinstall:
	@echo -e "\n\nInstalling ase - espresso interface to directory\n$(esp)\n\n"
	mkdir -p $(esp)
	cp `ls *.py|sed 's/setup.py\|__init__.py//'` $(EXTRAFILES) $(esp)
	sed s/SVNVERSION/$(shell svnversion 2>/dev/null|| echo no_version_found)/g <__init__.py >$(esp)/__init__.py
	$(py) -m compileall $(esp)
	$(CC) -s -O2 -o $(esp)/espfilter c-src/espfilter.c
