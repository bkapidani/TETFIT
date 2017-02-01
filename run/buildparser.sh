cd ../src/FIT/
g++ -O3 -std=gnu++14 -fdiagnostics-color=always -DENABLE_CYGWIN_MATH_FIX -fpermissive -o parser.exe parser.cpp
cd ../../run/
