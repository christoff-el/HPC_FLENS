PC = $(shell uname -n)

ifeq ($(PC),pacioli)
  $(shell module load gcc/4.7.2)
  $(shell module load sge/6.2u5)
  $(shell module load openmpi/gcc/64/1.4.2)
endif

subdirs = Flens_supl Fem examples
.PHONY: all clean
all:
	for i in $(subdirs); do (cd $$i && $(MAKE)) || exit 1; done

clean:
	for i in $(subdirs); do (cd $$i && $(MAKE) clean) || exit 1; done
	cd examples/output; rm *.dat; cd ../..
