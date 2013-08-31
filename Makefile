msg:
	@echo run python setup.py [--prefix=install-directory]
	@echo or make install to install the ase - espresso interface

install:
	python setup.py install

EXTRAFILES = espsite.py.example.*

py = python
esp = $(pylib)/espresso
bindir = $(instdir)/bin

define installscript
echo "#!$(py)" >$(bindir)/$(1)
tail -n +2 $(1) >>$(bindir)/$(1)
chmod 755 $(bindir)/$(1)
endef

doinstall:
	@echo -e "\n\nInstalling ase - espresso interface to directory\n$(instdir)\n\n"
	mkdir -p $(esp)
	mkdir -p $(bindir)
	cp `ls *.py|sed 's/setup.py\|__init__.py//'` $(EXTRAFILES) $(esp)
	sed s/GITVERSION/$(shell git describe --always 2>/dev/null|| echo no_version_found)-git/g <__init__.py >$(esp)/__init__.py
	$(py) -m compileall $(esp)
	$(CC) -s -O2 -o $(esp)/espfilter c-src/espfilter.c
	$(CC) -s -O2 -o $(bindir)/cubecutperiodic c-src/cubecutperiodic.c
	$(call installscript,pwlog2trajectory)
