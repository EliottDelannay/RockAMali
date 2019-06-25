all: gui nogui
	
gui: main.cpp
	g++ -O0 -o store.X main.cpp -I../CImg -Wall -W -ansi -pedantic -Dcimg_use_vt100 -I/usr/X11R6/include  -lm -L/usr/X11R6/lib -lpthread -lX11 && ./store.X -h -I && ./store.X -v > VERSION
nogui: main.cpp
	g++ -O0 -o store   main.cpp -I../CImg -Wall -W -ansi -pedantic -Dcimg_use_vt100 -lpthread -Dcimg_display=0 && ./store -h -I && ./store -v > VERSION

