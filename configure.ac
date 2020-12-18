# Initialize AC: ACINIT( package_name, version, bug_report_address)
AC_INIT([bcfcpp], [0.1])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# AM stands for Automake commands, AC stands for Autoconf commands
# We use libraries
# Let's check for a C++ compiler
AC_LANG(C++)
# Let's set a C++ compiler
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_17

# This suppresses the error we get for also having the cc version of helper library.

# Let's specify where the Makefiles should be produced.
# These are the same locations as your Makefile.in's, but named as Makefile only
# We need to do this because both Autoconf and Automake read this file and
# produce Makefiles from this list.
PKG_CHECK_MODULES([HTS], [htslib >= 1.11])
AC_CONFIG_FILES([Makefile src/Makefile])
# Finally produce "configure" script
AC_OUTPUT
