cd ../src/mesher/
./parser ../../run/input.txt
cd ../../build_folder
cmake ../src/ && make -j7
cd ../run/ && ../build_folder/fdtd_uniud
