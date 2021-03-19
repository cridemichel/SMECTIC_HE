set size ratio -1
unset key
set term pngcairo size 1500, 1500 font ",30"
set xlabel 'y'
set ylabel 'z' rotate by 0
set cbrange [:5]
set palette defined (-1 "black", -0.5 "dark-violet", 0 "violet", 1 "white")
set output 'gyz2.png'
plot 'gyz.txt' with image
