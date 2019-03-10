@echo off

rem Build object files
g++ -c src\process.cc -o src\process.o
g++ -c src\cpu_core.cc -o src\cpu_core.o
g++ -c src\main.cc -o src\main.o

rem Link object files
g++ src\process.o src\cpu_core.o src\main.o -o bin\lab1_part1.exe

rem Remove object files
if exist "src\*.o" del "src\*.o" /q
