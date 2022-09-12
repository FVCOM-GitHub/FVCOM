# $Id$
# $Name$
# $Revision$
#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh

SUBDIRS = nctools parallel interp


all:
	for dir in $(SUBDIRS); do (cd $$dir && make) || exit 1; done

clean:
	for dir in $(SUBDIRS); do (cd $$dir && make clean) || exit 1; done

clobber:
	for dir in $(SUBDIRS); do (cd $$dir && make clobber) || exit 1; done

