TARGET := libSRobot.so

CXX := g++
CFLAGS :=  -Wreturn-type -std=c++11 -DCALCULATE_MODE=1 -DPLATFORM=1 -ISType -I./ -fPIC -shared -O3

LDFLAGS := -pthread 

SRCS := $(wildcard *.cpp SType/*.cpp)

OBJS := $(patsubst %cpp,%o,$(SRCS))  

all: $(OBJS) 
	$(CXX) $^ $(CFLAGS) $(LDFLAGS) -o $(TARGET)
clean:
	rm $(TARGET) SType/*.o ./*.o
clean_all:  
	rm $(TARGET) $(OBJS) 
install:
	sudo mv $(TARGET) /lib
	
%.o:%.cpp  
	$(CXX) $< $(CFLAGS) $(LDFLAGS) -c -o $@ 
