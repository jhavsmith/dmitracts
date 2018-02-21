# makefile for dfuse

CC = /usr/bin/gcc		# use to compile with gcc


OPTIONS =   -O2 -DGSL_RANGE_CHECK_OFF #-g -pg -O2	
#g gives debugging code, pg adds code for profiling with gprof (Saturn is front-end on OS X
# -DGSL_RANGE_CHECK_OFF
CFLAGS = -c $(OPTIONS) -Wall  #Wall turns on all are warning messages


# This is a rule for making .o files from .c files of the same name
.c.o:
	$(CC) $(CFLAGS) $*.c

# This rule says that .o files will be compiled whenever the .c or
# .h of the same name has been modified more recently than the .o has.
# This is a shortcut to typing the names of the .c and .h files in the
# dependency lists below.
.SUFFIXES: .o .c .h

# This is a list of all object files necessary for linking. Any new
# source files that are created should have their corresponding .o
# added to this list.
OBJECTS = main_dfuse.o func_dfuse.o 

# This is the default compilation "target" and its dependencies. Here
# the only dependencies are the OBJECTS listed above. The second line
# is the command for compilation. (Note that there MUST be a tab
# before any command.) This command is issued whenever one of the
# dependencies has been modified more recently than the target.
# The macro $@ refers to the name of the target.
dfuse: $(OBJECTS)
	$(CC) $(OPTIONS) $(OBJECTS) -o $@ -lgsl  -lm -lgslcblas

# swapped "-framework Accelerate" for "-lgslcblas", 21 Aug 06

# This is a list of all the .o files and their dependencies. Note that
# the .c and .h files with the same name as the .o are not included
# because of the .SUFFIXES rule above. Any file #included by the
# corresponding .c or .h files should be placed in the dependency list.
# This does not include standard include files.
# The command in the .c.o rule will be issued whenever one of the
# dependencies has been modified more recently than the target. So if
# simhybrid.h is modified, then all four object files will be
# remade. But if only tree.c is modified, then only tree.o will be
# remade.
# Any new source files that are created should have their corresponding
# .o and dependency lists added here.
main_dfuse.o: head_dfuse.h
func_dfuse.o: head_dfuse.h

# This is just a target to make it easier to clean up the directory.
# To invoke this command, type "make clean" at the prompt. make with
# no arguments defaults to the first target in the makefile, but any
# target can be forced to be made by giving that target as a command
# line argument.
clean:
	rm -f *.o *.bak *~ dfuse
