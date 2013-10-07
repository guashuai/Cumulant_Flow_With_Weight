cal_cumulant_new : cal_cumulant_new.o
	g++ -O3 -m64 cal_cumulant_new.o -L/usr/local/root_v5.27/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lfreetype -pthread -Wl,-rpath,/usr/local/root_v5.27/lib -lm -ldl  -o cal_cumulant_new
cal_cumulant_new.o : cal_cumulant_new.C
	g++  -O3 -pipe -Wall -W -Woverloaded-virtual -D_REENTRANT -pthread -m64 -I/usr/local/root_v5.27/include -c cal_cumulant_new.C
