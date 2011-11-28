
include ../makefile.def


OBJS = qlc3d.o mesh.o simu.o electrodes.o lc.o box.o alignment.o readsettings.o ReadGiDMesh3D.o reordernodes.o\
		solutionvector.o setvolumeq.o setsurfacesq.o sparsematrix.o createsparsematrix.o calcpot3d.o writeLCD.o\
		tensortovector.o calcQ3d.o  assembleQ.o solve_pcg.o eventlist.o settings.o geometry.o energy.o reader.o\
		inits.o refinement.o refinement2.o line.o meshrefinement.o autorefinement.o filesysfun.o



qlc3d.exe:$(OBJS)
	g++  -fopenmp -o qlc3d.exe  $(OBJS) $(LDFLAGS) -lgomp -lpthread

%.o:%.cpp
	$(CXX) -c $(CFLAGS) $(CCCFLAGS) $(D_OF) $< $(INCLUDES)

clean:
	rm -f qlc3d.exe *.o
