CXX=mpic++ -Wall -O3 -fopenmp

CPPFLAGS=-Iinc 

LDFLAGS=-lboost_program_options -lm 

FILES=export forward init memory parse_args shalw
OBJECTS=$(addsuffix .o, $(FILES))
BIN=bin/shalw

all : $(BIN)

$(BIN) : $(addprefix obj/, $(OBJECTS))
	$(CXX) -o $@ $^ $(LDFLAGS)

obj/%.o : src/%.cpp
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

obj/%.o : src/%.c
	$(CXX) -c -o $@ $^ $(CPPFLAGS)

clean :
	rm -f bin/* obj/*

exec:
	mpirun -n 2 -bynode ./bin/shalw --export --export-path /tmp/3302011
	(export N = 40 ; for e in $(BIN) ; do echo "Running $$e ..." ; echo "" >> $(OUT) ; echo "Running $$e ..." >> $(OUT) ; ./$$e $$N >> $(OUT); done)
