#!/bin/sh
rm -f configure
aclocal
autoconf
touch AUTHORS NEWS README ChangeLog
automake --add-missing

