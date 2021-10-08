# scale for fig 3.4
#set xrange [-4:4]
#set yrange [-2:4]
#set zrange [2:-6]
# scale for fig 3.5
#set xrange [-2:2]
#set yrange [-1:4]
#set zrange [3:0]
# scale for fig 3.6
#set xrange [-2:2]
#set yrange [5:10]
#set zrange [-3:-8]

# scale for three figs at once
set xrange [-4:4]
set yrange [-2:22]
set zrange [3:-42]

#set autoscale
unset log
unset label
set xtic auto
set ytic auto

splot "./data_02.txt" u 1:2:3 with lines, \
"./data_02.txt" u 4:5:6 with lines, \
"./data_02.txt" u 7:8:9 with lines

set terminal postscript eps color
set output "trajectory.eps"
replot
save "plotTrajectory.plt"
set term pop

