bash scripts

compile.sh compiles the code
run.sh runs the compiled code
both.sh compiles then runs the code
yeet.sh compiles (and runs) in 03 and march=native (basically very fast) 


parameter file

number of x cells  (the simulation time scales with ~ sqrt( n_x^3 * n_y^3 ) )
number of y cells
spatial accuracy of derivatives, 1,2,3 (3 most accurate but slower, should probably use 3)
physical length scale (ie total size)
inverse sigma squared for initial gaussian profile
total simulation time (try 0.005)



compiling and running by hand

g++ ScalarField.hpp ScalarField.cpp ErwinSchrodinger.cpp -o prog
./prog
python plot.py

note, in bashscripts i put g++-8 as g++ defaults to clang on mac and the compiling with multiple files breaks
depending on your machine may want to replace g++-8 -> g++, gcc, or other 



plot.py

edit variable max_plot_value to change the colorbar limits 


