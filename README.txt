bash scripts

yeet.sh compiles (and runs) in 02 and march=native (basically very fast) 
clear.sh clears up all dat files and all saved png's from the python script

parameter file

number of x cells  (the simulation time scales with ~ sqrt( n_x^4 * n_y^4 ) )
number of y cells
spatial accuracy of derivatives, 1,2,3 (3 most accurate but slower, should probably use 3)
physical length scale (ie total size)
inverse sigma squared for initial gaussian profile
total simulation time (try 0.03)


note, in bashscripts i put g++-8 as g++ defaults to clang on mac and the compiling with multiple files breaks
depending on your machine may want to replace g++-8 -> g++, gcc, or other 


to change frequency of verbosity and .dat saving modify main(), specifically the loop containing rk4




to turn multiple png's into an mp4, use

ffmpeg -r 15 -pattern_type glob -i '*.png' -c:v libx264 output.mp4

where 15 is framerate (can change this obv) and output.mp4 is the file name. other flags should not be changed
NB: this requires the download of ffmpeg


