cd ../src/FIT/
./parser ../../run/input.txt
cd ../../build_folder
cmake ../src/ && make -j7
cd ../output/ && ../build_folder/fit_uniud
