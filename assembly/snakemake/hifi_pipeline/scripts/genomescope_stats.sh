## get useful values
cd $1
VAR="$(grep -n "Genome Haploid Length" summary.txt | cut -f1 -d:)"
sed -n $VAR\p summary.txt | sed -e 's/  \+/\t/g' | cut -f3 | sed -e 's/,//g' | sed -e 's/ bp//g' > estimated_genome_size.txt
VAR="$(grep -n "kmercov " model.txt | cut -f1 -d:)"
KCOV="$(printf "%.2f\n" $(sed -n $VAR\p model.txt | sed -e 's/ \+/\t/g' | cut -f2))"
printf "%.0f\n" $(echo "$KCOV * 1.5" | bc) > transition_parameter.txt
printf "%.0f\n" $(echo ""$(cat transition_parameter.txt)" * 3" | bc) > maximum_depth.txt

