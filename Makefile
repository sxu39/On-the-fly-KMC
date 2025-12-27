D_SRC=./src
D_OBJ=./obj

SRC_CPP = $(wildcard $(D_SRC)/*.cpp)
OBJ_CPP = $(addprefix $(D_OBJ)/, $(patsubst %.cpp, %.o, $(notdir $(SRC_CPP))))

TARGET=./bin/kmc

#DEEPMD     = /root/lbg_2_bk/Software/deepmd
#TENSORFLOW = /root/lbg_2_bk/Software/tf

DEEPMD = /opt/deepmd-kit-2.2.1
TENSORFLOW = /opt/deepmd-kit-2.2.1

#CC=icc -g
#CC=mpiicc 
#CC=mpicxx
#CC=CC
#HONG=-D__MPI
#CC=mpiicc
CC=g++
OPTION=-O3 -Wno-unused-result -D HIGH_PREC -L${DEEPMD}/lib -L${TENSORFLOW}/lib -I${DEEPMD}/include -Wl,--no-as-needed -ldeepmd_cc -lstdc++ -ltensorflow_cc -Wl,-rpath=${DEEPMD}/lib -Wl,-rpath=${TENSORFLOW}/lib
#OPTION=-O3 -Wno-unused-result

${TARGET}:$(OBJ_CPP)
	$(CC) $(OPTION) -o  $@ $^

$(D_OBJ)/%.o: $(D_SRC)/%.cpp
	@if [ ! -d $(D_OBJ) ]; then mkdir $(D_OBJ); fi
	@if [ ! -d ./bin ]; then mkdir ./bin; fi
	$(CC) -c  $(OPTION) $< -o $@

.PHONY: clean

clean:
	rm -f $(D_OBJ)/* $(TARGET) 

