set term x11
set xrange [0.9:1.2]
set yrange [-0.01:0.2]
set xlabel "FX Rate"
set ylabel "Option Value ($)"
set title "Value Functions for Options at Different Time to Expiry\n**From Regressor Pricer**"
plot "tmp/BucketFit0.dat","tmp/BucketFit1.dat","tmp/BucketFit2.dat","tmp/BucketFit3.dat","tmp/BucketFit4.dat","tmp/BucketFit5.dat","tmp/BucketFit6.dat","tmp/BucketFit7.dat","tmp/BucketFit8.dat","tmp/BucketFit9.dat","tmp/BucketFit10.dat","tmp/BucketFit11.dat","tmp/BucketFit12.dat",

# 