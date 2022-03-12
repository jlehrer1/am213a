OBJECTS = LinAl.o Driver_LinAl.o
CHECKS = -fimplicit-none  -Wall  -Wextra -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all -pedantic -fbacktrace
MODULES = LinAl.mod
CONTAINER = jmlehrer/fortran-test

.PHONY: clean test parta bone btwo 

output.txt: LinAl.exe
	./LinAl.exe > output.txt

test:
	gfortran $(CHECKS) -llapack LinAl.f90 minimal.f90 -o test

parta:
	gfortran $(CHECKS) -llapack LinAl.f90 parta.f90 -o parta

bone:
	gfortran $(CHECKS) LinAl.f90 partb_1.f90 -o bone
	 
btwo:
	gfortran $(CHECKS) LinAl.f90 partb_2.f90 -o btwo

%.o: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8 $(CHECKS) $<

%.mod: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8 $(CHECKS) $<

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl.exe

exec:
	docker exec -it $(CONTAINER) /bin/bash

build:
	docker build -t $(CONTAINER) .

run:
	docker run -it $(CONTAINER) /bin/bash

push:
	docker push $(CONTAINER)
