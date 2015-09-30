SRC      = $(wildcard *.cpp)
RM       = rm -rf
#CC       = g++
#CFLAGS   = -g -O3 -fomit-frame-pointer -march=opteron -c -I./jsim -I./jsim/include
CC       = mpicxx
CFLAGS   = -O3 -DMPICH_IGNORE_CXX_SEEK -c -I./ -I../include
LDFLAGS  = -L../lib -lstdc++

OBJ      = $(SRC:%.cpp=%.o)
EXE      = $(SRC:%.cpp=%.exe)

all: $(OBJ) $(EXE)
	@echo 'Done making...'

%.exe: %.o
	@echo Making $< from $@
	$(CC) -o $@ $< ${LDFLAGS}

%.o: %.cpp
	@echo Making $< from $@
	$(CC) $(CFLAGS) $< -o $@

clean:
	@echo 'Cleaning...'
	-$(RM) ${OBJ}
	-$(RM) ${EXE}
	-$(RM) *~
	@echo 'Done cleaning...'
