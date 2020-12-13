INCLDIR	:= include
OBJDIR	:= obj
SRCDIR	:= src
BINDIR	:= bin


BLASINC= -DWITHBLAS #-framework Accelerate
UMFPACKINC= -I/usr/local/include
UMFPACKLIBS= -L/usr/local/lib -lumfpack  -lcholmod -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig $(BLASLIBS)
#  macos install a  avec MKL

DIR_UMFPACK= /usr/include/suitesparse
UMFPACKINC = -I$(DIR_UMFPACK)
UMFPACKLIBS = -L$(DIR_UMFPACK)/lib -lumfpack  -lcholmod  -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig $(BLASLIBS)


CC      := g++ -fopenmp
VPATH	:=
LDFLAGS :=
LIBRARY :=
CFLAGS  := -std=c++11 $(UMFPACKINC) -I /usr/include/ -I $(INCLDIR) -MMD -MP

#Source and object files (automatic)
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(subst $(SRCDIR)/,$(OBJDIR)/, $(subst .cpp,.o, $(SRCS)))

# Define here your main source files separated by spaces (without suffix!)
EXEC = mainlapP1 mainlapP2 SteadyStokes check421 TimeSteadyStokes NavierStokes

#Phony = do not represent a file
#.PHONY: all
all : makedir $(EXEC)

-include $(SRC:%.cpp=%.d)

# For multiple binaries
$(EXEC) : %: %.cpp $(OBJS)
	$(CC) -o $(BINDIR)/$@ $^ $(CFLAGS) $(UMFPACKINC) $(UMFPACKLIBS)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC)  -c -o $@ $< $(CFLAGS) $(UMFPACKINC) $(UMFPACKLIBS)

#Clean: delete every binaries and object files
.PHONY: clean
clean :
	rm -rf $(OBJDIR)/*
	rm -rf $(BINDIR)/*
#Building folders (-p : no error if folder do not exist)
.PHONY: makedir
makedir :
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)

#For some debug
.PHONY: print
print :
	echo $(SRCS)
	echo $(OBJS)

#Remarks:
# $@ : filename representing the target
# $< : filename of the first prerequisite
# $^ : filenames of all prerequisites, separated by spaces. Dupplicated are removed.
# $? : names of all prerequisites that are newer than the target, separated by spaces
# $+ : similar to $^ but include dupplicates
# $* : stem of target file (filename without suffix)
