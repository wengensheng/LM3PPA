# Makefile created by mkmfweng $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 


.DEFAULT:
	-echo $@ does not exist.
all: a.out
SRC =
OBJ =
clean: neat
	-rm -f .a.out.cppdefs $(OBJ) a.out
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
a.out: $(OBJ) 
	$(LD) $(OBJ) -o a.out  $(LDFLAGS)
