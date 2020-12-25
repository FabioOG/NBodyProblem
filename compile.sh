g++ -Wall -fexceptions -O2  -c main.cpp -o obj/main.o
g++ -Wall -fexceptions -O2  -c src/body.cpp -o obj/body.o
g++ -Wall -fexceptions -O2  -c src/NBody.cpp -o obj/NBody.o
g++ -Wall -fexceptions -O2  -c src/vector3.cpp -o obj/vector3.o
g++  -o bin/nBodyProblem obj/main.o obj/body.o obj/NBody.o obj/vector3.o -s 
