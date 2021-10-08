for i in `seq 1 21`;
do 

    n=`printf %0.2d $i`;
    j=`echo $i - 1.0 | bc -l`;
    k=`echo $i - 0.01 | bc -l`;

    echo "~/tweezers/standalone_double/double_driver -l ${j}e-6 -r ${k}e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >> cpu_double_prob_z_${n}.txt";

    ./../standalone_double/double_driver -l ${j}e-6 -r ${k}e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100 >> cpu_double_prob_z_${n}.txt

done;


# a sample output:
# -l 0e-6 -r .99e-6 -t -20.0e-6 -b 8.1e-6 -d 1.0 -n 100