OBJECTS = LinAl.o Driver_LinAl.o
CHECKS = -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all -pedantic  -fbacktrace
MODULES = LinAl.mod

.PHONY: clean

output.txt: LinAl.exe
	./LinAl.exe > output.txt

LinAl.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) $(CHECKS) -o LinAl.exe

Q2: 
	gfortran $(CHECKS) LinAl.f90 question_2.f90 -o question2

Q1:
	gfortran $(CHECKS) LinAl.f90 question_1.f90 -o question1

%.o: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8 $(CHECKS) $<

%.mod: %.f90
	gfortran -c -fdefault-real-8 -fdefault-double-8 $(CHECKS) $<

clean:
	rm -f $(OBJECTS) $(MODULES) LinAl.exe
