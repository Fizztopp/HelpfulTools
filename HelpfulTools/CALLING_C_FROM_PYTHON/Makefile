all: pybind libmatmul

pybind:
	g++ -O3 -Wall -shared -std=c++11 -fPIC -Ieigen3/ `python3-config --includes` -Ipybind11/include `python3 -m pybind11 --includes` pybind11.cpp -o matmul_pybind`python3-config --extension-suffix`

libmatmul:
	gcc -c -fPIC -O3 matmul.c -o matmul.o
	gcc -shared -o libmatmul.so matmul.o
