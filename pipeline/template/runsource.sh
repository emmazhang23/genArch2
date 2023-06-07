module load python/3.5.1
#WIP
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
path=./
###########################

#to make the naming convention for decrease to fitness
name=$loci"l_"$large"df_"$selection"s"

id=$(sbatch $dep --ntasks=1 --cpus-per-task=1 -N 1 --job-name=${name}.gen --output=${name}._genarch_analyses_source_out --mem=${mem}00 $queue --time=${days}-0 --wrap="python3 genarch_analyses_source.py -n $name -p $path -l $loci -r $reps -f 2")
echo $id