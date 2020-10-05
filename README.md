# ornamental tool

## dependencies:
- g++-7 (>= 7.5.0)
- boost-1.74 (>= boost-1.72)
- gsl (>= 2.1)
- gnuplot (>= 4.6)
- ncurses (>= libncurses5.so)

## compilation:
- g++-7 -c -O2 -std=c++17 -I./include -o main.o main.cpp
- g++-7 -o ornamental main.o -L /usr/local/lib -lgsl -lgslcblas -lcurses

## running:
./ornamental [ornamental.init]

## list of commands (v.2):
```
help|h                                      help
environment|env [<var> <value>]             list of environment variables
window|w                                    open window to view the plots
files|f                                     view status of the data files
files|f <fileID> R|W|RW                     set status of the data file
fload|f< <filename> [<fileID>]              load data file
fsave|f> <fileID> <filename>                save the file as
fadd|f+ <fileID>                            create new data file
fdel|f- <fileID>                            delete the file 
festimate|f^ MLM|random <fileID>            estimate of missing data
models|m                                    view status of all the models
models|m <modelID>                          set the model as active
madd|m+ <fileID:sproutID> [<modelID>]       build a model of the sprout
mdel|m- <modelID>                           delete the model 
plots|p                                     show IDs of the all plots
plot|p <plotID>                             show the plot
padd|p+ <plotID> [width value]              add the plot with smoooth width
pdel|p-                                     delete the current plot
pdel|p- <plotID>                            delete the plot
psave|p> <filename>                         save the current plot (as pdf)
scheme <modelID> [<date|day=0>]             show the scheme of inflorescence
scheme <fileID>:<sproutID> [<date|day=0>]   show the scheme of inflorescence
test|t <plotID> [nsamples]                  statistics and normality test for plot
generate|g <modelID> <N>[<n1>:<n2>]         generate N model sprouts
quit|exit                                   exit the program
```
