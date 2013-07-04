subdirs = LinearAlgebra Flens_supl Fem examples
.PHONY: all clean
all:
	for i in $(subdirs); do (cd $$i && $(MAKE)) || exit 1; done

clean:
	for i in $(subdirs); do (cd $$i && $(MAKE) clean) || exit 1; done
	cd examples/output; rm *.dat; cd ../..
