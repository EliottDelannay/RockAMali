#run
## ushort = 2uchar: 4096*2 = 8192BoF
## uint   = 4uchar: 2048*2 = 8192BoF
FRAME_SIZE=2048
NP=4
USE_GPU=--use-GPU --GPU-factory program
#USE_GPU=
DO_CHECK=--do-check
#DO_CHECK=

DATA=./
DATA=/media/temp/
DIN=samples/
DOUT=results/
FIN=sample.cimg
FOUT=$(FIN)

#compiler options
LIB_XWINDOWS=-I/usr/X11R6/include -L/usr/X11R6/lib -lX11
LIB_CIMG=-I../CImg -Wall -W -ansi -pedantic -Dcimg_use_vt100 -lpthread -lm -fopenmp
LIB_BOOST_ASIO=-lboost_system
LIB_BOOST_COMPUTE=-lMali -L/usr/lib/aarch64-linux-gnu/ -DBOOST_COMPUTE_MAX_CL_VERSION=102

DO_GPU=-DDO_GPU $(LIB_BOOST_COMPUTE)
#DO_GPU=

#source package
SRC_DATA_BUFFER=thread_lock.hpp CDataAccess.hpp CDataBuffer.hpp
SRC_DATA_PROCESS=CDataProcessor.hpp CDataProcessorGPU.hpp CDataProcessorGPUfactory.hpp
HELP_OUTPUT=process_sequential_help.output process_help.output send_help.output receive_help.output store_help.output

all: process_sequential process send receive doc
#all: process process_sequential doc

#all: time_copy
time_copy: time_copy.cpp
	g++ -O0 -o time_copy time_copy.cpp $(DO_GPU) && ./time_copy

gui: main.cpp
	g++ -O0 -o generate.X main.cpp -I../CImg -Wall -W -ansi -pedantic -Dcimg_use_vt100 -lpthread -lm -fopenmp -lboost_system $(LIB_XWINDOWS) && ./generate.X -h -I && ./generate.X -v > VERSION
	./generate.X -h 2> generateX_help.output

process: process.cpp $(SRC_DATA_BUFFER) CDataGenerator.hpp $(SRC_DATA_PROCESS) CDataStore.hpp
	g++ -O0 -o process   process.cpp $(LIB_CIMG) -Dcimg_display=0 $(DO_GPU) && ./process -h -I && ./process -v > VERSION
	./process -h 2> process_help.output

#SEQ_GPU=
#SEQ_GPU=-DDO_GPU_SEQ_QUEUE
SEQ_GPU=-DDO_GPU_NO_QUEUE
process_sequential: process_sequential.cpp $(SRC_DATA_BUFFER) CDataGenerator.hpp $(SRC_DATA_PROCESS) CDataStore.hpp
	g++ $(SEQ_GPU) -O0 -o process_sequential   process_sequential.cpp $(LIB_CIMG) -Dcimg_display=0 $(DO_GPU) && ./process_sequential -h -I && ./process_sequential -v > VERSION
	./process_sequential -h 2> process_sequential_help.output

send: send.cpp $(SRC_DATA_BUFFER) CDataGenerator.hpp CDataSend.hpp
	g++ -O0 -o send   send.cpp  $(LIB_CIMG) $(LIB_BOOST_ASIO) -Dcimg_display=0 && ./send -h -I && ./send -v > VERSION
	./send -h 2> send_help.output

receive: receive.cpp $(SRC_DATA_BUFFER) CDataReceive.hpp $(SRC_DATA_PROCESS) CDataStore.hpp
	g++ -O0 -o receive receive.cpp  $(LIB_CIMG) $(LIB_BOOST_ASIO) -Dcimg_display=0 $(DO_GPU) && ./receive -h -I && ./receive -v > VERSION
	./receive -h 2> receive_help.output

doc: doxygen.cpp VERSION VERSIONS $(HELP_OUTPUT) process.cpp process_sequential.cpp send.cpp receive.cpp  $(SRC_DATA_BUFFER) CDataReceive.hpp $(SRC_DATA_PROCESS) CDataStore.hpp
	./doxygen.sh

#NP=4
NT=`echo $(NP)+2   | bc`
NB=`echo $(NP)*4096| bc`
NS=`echo $(NP)*8192| bc`
process_run:
#	./process -c 4 -s $(FRAME_SIZE) -o $(DATA)$(DIN)$(FIN) -r $(DATA)$(DOUT)$(FOUT) -b 8 -n 16 --use-GPU --do-check #2>/dev/null | grep -e info -e test
	./process -c $(NT) -s $(FRAME_SIZE) -o $(DATA)$(DIN)$(FIN) -r $(DATA)$(DOUT)$(FOUT) -b $(NB) -n $(NS) $(USE_GPU) $(DO_CHECK) 2>&1 | grep -e info -e test -e failed -e double -e fault --color

process_sequential_run:
	./process_sequential -s $(FRAME_SIZE) -o $(DATA)$(DIN)$(FIN) -r $(DATA)$(DOUT)$(FOUT) -n 123 $(USE_GPU) $(DO_CHECK)

NS=123456
send_run:
	./send    -c 2 -s $(FRAME_SIZE) -b  8 -n $(NS) -w 2345678

#NP=4
NT=`echo $(NP)+3   | bc`
NB=`echo $(NP)*4096| bc`
receive_run: clear
#	./receive -c 2 -s $(FRAME_SIZE) -b 128 -n 12345 -o $(DATA)$(DIN)$(FIN) -r $(DATA)$(DOUT)$(FOUT) -C -E -W
#	./receive -c 3 -s $(FRAME_SIZE) -b 16 -n 123 -o $(DATA)$(DIN)$(FIN) -r $(DATA)$(DOUT)$(FOUT) $(USE_GPU) $(DO_CHECK) -E -W
#	./receive -c 4 -s $(FRAME_SIZE) -b 16 -n 123 -o $(DATA)$(DIN)$(FIN) -r $(DATA)$(DOUT)$(FOUT) $(USE_GPU) $(DO_CHECK) -E -W
	./receive -c $(NT) -s $(FRAME_SIZE) -b $(NB) -n $(NS) -o $(DATA)$(DIN)$(FIN) -r $(DATA)$(DOUT)$(FOUT) $(USE_GPU) $(DO_CHECK) -W

clear:
	rm -fr $(DATA)/samples/ $(DATA)/results/
	mkdir  $(DATA)/samples/ $(DATA)/results/
	rm -f sample_??????.cimg
	sync

clean: clear
	rm -f send.X    send
	rm -f receive.X receive
	rm -f process.X process
	rm -f process_sequential.X process_sequential

display:
	convert -append $(DATA)/samples/sample*.png $(DATA)/samples.png && display samples.png &
	convert -append $(DATA)/results/sample*.png $(DATA)/results.png && display results.png

