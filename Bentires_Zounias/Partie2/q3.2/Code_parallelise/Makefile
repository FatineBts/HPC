CXX=mpic++ -Wall -O3 -mavx

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