
##########################  MERGE FILES

# You need the list of Fulani populations in one file called "List_of_Pops.txt". One row per population

## First prepare the new files
cat List_of_pops.txt | while read line; do
echo ${line}
for i in {1..22}; do
sed -n '/BEGIN GENOTYPES/,$p' Populations_inputs/${line}/${line}_${i}_hapguess_switch.out | tail -n +2 | sed '$ d' > Populations_inputs/${line}/NEW_${line}_${i}_hapguess_switch.out
done; done

## Second merge the files for each group. Please replace the word "Group" with the name of the group

Group="Eastern" # Please indicate here the name of the group
mkdir Populations_inputs/${Group}
for i in {1..22}; do
	for Pop in Mali_Fulani_InnerDelta Mali_Fulani_Diafarabe BurkinaFaso_Fulani_Tindangou BurkinaFaso_Fulani_Banfora Cameroon_Fulani_Tcheboua BurkinaFaso_Fulani_Ziniare; do 
	cat Populations_inputs/${Pop}/NEW_${Pop}_${i}_hapguess_switch.out >> Populations_inputs/${Group}/${Group}_${i}_hapguess_switch.out; done
done

# Append fastPHASE intro and END to the file.
Group="Eastern"
for i in {1..22}; do
ed -s Populations_inputs/${Group}/${Group}_${i}_hapguess_switch.out << 'EOF'
0a
********************************************
*                                          *
*      Output from fastPHASE 1.4.0         *
*      Code by P Scheet                    *
*      pscheet@alum.wustl.edu              *
*                                          *
********************************************
BEGIN COMMAND_LINE
fastPHASE -S-77 -I_chr1_allstarts -XI -T25 -C0 -K25 -ooutput_final
END COMMAND_LINE

BEGIN COMMAND_EXPLAIN
 K no. clusters (chosen or supplied): 25
 S seed for random numbers (chosen or supplied): -77
END COMMAND_EXPLAIN

BEGIN DESCRIBE_TASKS
minimize switch error
END DESCRIBE_TASKS

BEGIN GENOTYPES
.
$a
END GENOTYPES
.
w
EOF
echo "${Group}_${i} hapguess done"
done

