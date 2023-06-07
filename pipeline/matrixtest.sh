#to make deleting the tests easier lol
mkdir test
cd test

#array for decrease in fitness k (1-effect size)*100
FITNESS=(0 95 100)

#array for selection intensity p
SELECTION=(1 3 1000)


#makes template file to copy into each smaller directory
mkdir template
cd ..
cp -r template/. test/template/

cd test

#creates directories for each 1, 10, 100, 1k loci in one directory titled "test"
#creates directories for each combination of decrease in fitness and selection we want using the arrays 

for i in  1 10 100 1000
do
   mkdir $i"loci"
   cd $i"loci"
   for k in ${FITNESS[@]}; do
   for p in ${SELECTION[@]}; do
   mkdir $k"df_"$p"s"
   cd ..
   cp -r template/. $i"loci"/$k"df_"$p"s"
   cd $i"loci"/$k"df_"$p"s"
   chmod +x quantinemo
   sed -i "s/DECREASE/"$k"/" run.sh
   sed -i "s/SELECTION/"$p"/" run.sh
   sed -i "s/LOCI/"$i"/" run.sh
   bash run.sh
   cd ..
   done
   done
   cd ..
done
cd 
