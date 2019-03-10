@echo off

rem Build object files
g++ -c src\main.cc -o src\main.o

rem Link object files
g++ src\main.o -o bin\lab1_part2.exe

rem Remove object files
if exist "src\*.o" del "src\*.o" /q
