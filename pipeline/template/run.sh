module load python/3.5.1

#########Submission########
n=1
mem=200
queue="--partition=general"
days=7
###########################

large=DECREASE
de=$(($large / 100))

####Simulation Settings#####
loci=LOCI
eloci=2
reps=30
decrease=$de
dispersal=0.001
mu=0.0001
selection=SELECTION
minfreq=0.04
path=/pine/scr/b/s/bsr/quantinemo/
###########################

#to make the naming convention for decrease to fitness
name=$loci"l_"$large"df_"$selection"s"


##All parameters have default settings, to specify settings you want to change, change them above, and then add the name of the variable to the string below
##In the example below I have changed loci and fecundity, so they are included and everything else is default
id=$(sbatch $dep --ntasks=1 --cpus-per-task=1 -N 1 --job-name=${name}.gen --output=${name}_gen.out --mem=1000 $queue --time=${days}-0 --wrap="python3 genEpi.py -n $name -l $eloci -d $decrease")
echo $id
id=${id##* }
dep="--dependency=afterok:$id"


id=$(sbatch $dep --ntasks=1 --cpus-per-task=1 -N 1 --job-name=${name}.gen --output=${name}_gen.out --mem=1000 $queue --time=${days}-0 --wrap="python3 GenSettings.py -s $selection -n $name -l $loci -le $eloci -r $reps -u $mu -t template")
echo $id
id=${id##* }
dep="--dependency=afterok:$id"

id=$(sbatch $dep --ntasks=1 --cpus-per-task=$n -N 1 --job-name=${name}.qntnm --output=${name}_qtnm.out --mem=${mem}000 $queue --time=${days}-0 --wrap="./quantinemo ${name}.ini")
echo $id
id=${id##* }
dep="--dependency=afterok:$id"


id=$(sbatch $dep --ntasks=1 --cpus-per-task=1 -N 1 --job-name=${name}.anlys --output=${name}_anlys.out --mem=${mem}000 $queue --time=${days}-0 --wrap="python3 genLag.py -r $reps -n $name")
echo $id
id=${id##* }
dep="--dependency=afterok:$id"

id=$(sbatch $dep --ntasks=1 --cpus-per-task=1 -N 1 --job-name=${name}.gen --output=${name}_gen.out --mem=${mem}00 $queue --time=${days}-0 --wrap="python3 genoTime.py -n $name -le $eloci -l $loci -r $reps")
echo $id
id=${id##* }
dep="--dependency=afterok:$id"

id=$(sbatch $dep --ntasks=1 --cpus-per-task=1 -N 1 --job-name=${name}.gen --output=${name}_poptime.out --mem=${mem}00 $queue --time=${days}-0 --wrap="python3 poptime.py -l $i -s $p -d $k")
echo $id
id=${id##* }
dep="--dependency=afterok:$id"

