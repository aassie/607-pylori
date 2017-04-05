#Getting the files from NCBI
rsync --copy-links --recursive --include "*.gbff.gz" --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Helicobacter_pylori/latest_assembly_versions/ ./

################################
####RUNNING DNA AND PLASMIDS####
################################

#Using the fnas
ls Pylori/*/*genomic.fna.gz | grep -v "cds" | grep -v "rna" > Original.names.w.folder.tab
Oname=$(<Original.names.w.folder.tab)
for i in $Oname; do cp $i fasta/; done
cd fasta
gunzip -d *

#Renaming the files and their contigs with a shorter name
mkdir Renamed

a=1; 
for i in *.fna; 
	do old=$(printf $i); 
	new=$(printf "HP-%04d.fna" "$a"); 
	echo -e $old ' \t ' $new >> Old.and.new.names.tab; 
	awk -v A=$new '/^>/{print ">"A"_" ++i; next}{print}' < $i > Renamed/$new; let a=a+1; 
done

#Running CheckM on lineage mode

mkdir checkm
checkm lineage_wf -x fna --pplacer_threads 10 -t 10 -f checkm/Checkm.results.log ./ checkm/

#Making a tree of the Checkm alignment On geneious

#Getting some meta data out of the Genbank file

for i in *gbff; do printf "%s\t %s\t %s\t %s\t %s\t %s\n" "$i" "$(grep "TITLE" $i | head -n 1)" "$(grep "AUTHORS" $i | head -n 1)" "$(grep "BioProject" $i | head -n 1)" "$(grep "country" $i | head -n 1)" "$(grep "LOCUS" $i | head -n 1)"; done > Meta.data.tab


#Prokka
#Run on quiet mode with single echo of analysis milestone

ls *.fna | sed 's/.fna//g' > name
cd ../..
mv fasta/Renamed/name ./New.names.tab
Nname=$(<New.names.tab)
for i in $Nname; 
	do prokka --cpus 10 --outdir prokka/$i/ --force --quiet --prefix $i fasta/Renamed/$i.fna && echo "$i annotation done";
done

#single line for multiline fasta file (found here https://www.biostars.org/p/9262/)
for i in {001..607}; 	
	do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < prokka/HP-0$i/HP-0$i.faa > prokka/HP-0$i/HP-0$i.sl.faa ; 
done

#rename genes by replacing the Prokka locus tag by the HP one, but keeping the numbering

for i in {001..607};
	do prefi=$(awk -F "_" '/^>/ {print $1}' prokka/HP-0$i/HP-0$i.sl.faa | head -n 1) ;
	sed "s/$prefi/>HP-0$i/g" prokka/HP-0$i/HP-0$i.sl.faa > prokka/HP-0$i/HP-0$i.sl.rn.faa;
done

#Same BS but for nucleotide

for i in {001..607}; 	
	do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < prokka/HP-0$i/HP-0$i.ffn > prokka/HP-0$i/HP-0$i.sl.ffn;
	prefi=$(awk -F "_" '/^>/ {print $1}' prokka/HP-0$i/HP-0$i.sl.ffn | head -n 1) ;
	sed "s/$prefi/>HP-0$i/g" prokka/HP-0$i/HP-0$i.sl.ffn > prokka/HP-0$i/HP-0$i.sl.rn.ffn;
done

#Extracting vacA, -h remove file path and sed remove the "--" between grep
grep -h -A 1 "Vacuolating cytotoxin autotransporter" prokka/HP-*/*.rn.faa | sed 's/--//g' > vacA.gene.collection.fasta

#Extracting gene name from tree file to edit color information
grep description cag.gene.collection.tree | cut -d "[" -f 1 > tmp && cut -d "_" -f 1 tmp > tmp1 && cut -d "_" -f 2 tmp > tmp2 && cut -d "_" -f 3 tmp > tmp3 && paste -d "_" tmp1 tmp2 > tmp4 && paste tmp4 tmp3 > cag.tree.name 

#Renaming Cag gene in nex file

grep "description" cag.gene.collection.tree.nex | cut -d "[" -f1 | sed 's/\t//g' > tmp && cut -d "_" -f 2 tmp > tmp2 && readarray OldHP < ./tmp && readarray NewHP < ./tmp2 
cp cag.gene.collection.tree.nex cag.gene.collection.tree.2.nex

#-2 in the condition because the array will grep a "tree tree=1" value and the last one is always null

for ((i=0;i<${#OldHP[@]}-2;i++)); 
	do a=$(printf '%s' "${OldHP[i]}"); 
	b=$(printf '%s' "${NewHP[i]}");
	sed "s/$a/HP_$b/g" cag.gene.collection.tree.2.nex > renamed.cag.tree.nex;
	mv renamed.cag.tree.nex cag.gene.collection.tree.2.nex;
done

#Checking tree with Dendroscope

##########################
#### WITHOUT plasmisd ####
##########################

mkdir Renamed_wo_plasmid
a=1; 
for i in Original/*.fna; 
	do old=$(printf $i); 
	new=$(printf "HP-%04d.wop.fna" "$a");
	nname=$(printf "HPwop-%04d" "$a"); 
	echo -e $old ' \t ' $new >> Old.and.new.wo.plasmid.tab;
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $i > Renamed_wo_plasmid/tmp.$new; let a=a+1;
	sed '/lasmid/,+1 d' Renamed_wo_plasmid/tmp.$new | awk -v A=$nname '/^>/{print ">"A"_" ++i; next}{print}' > Renamed_wo_plasmid/$new;
	rm Renamed_wo_plasmid/tmp.$new;
done

#Running CheckM on lineage mode
cd Renamed_wo_plasmid
mkdir checkm
checkm lineage_wf -x fna --pplacer_threads 10 -t 10 -f checkm/Checkm.results.log ./ checkm/


###########################
#### JUST the plasmisd ####
###########################

mkdir Renamed_plasmid
a=1; 
for i in Original/*.fna; 
	do old=$(printf $i); 
	new=$(printf "HP-%04d.p.fna" "$a");
	nname=$(printf "HPp-%04d" "$a"); 
	echo -e $old ' \t ' $new >> Old.and.new.plasmid.tab;
	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $i > Renamed_plasmid/tmp.$new; let a=a+1;
	grep -A 1 "lasmid" Renamed_plasmid/tmp.$new | awk -v A=$nname '/^>/{print ">"A"_" ++i; next}{print}' > Renamed_plasmid/$new;
	gcount=$(grep -c ">" Renamed_plasmid/$new);
	if [ $gcount == 0 ];
	then rm Renamed_plasmid/$new;
	fi; 	
	rm Renamed_plasmid/tmp.$new;
done

##############################################
####Reannotation of plasmid vs non plasmid####
##############################################

cd fasta/Renamed_wo_plasmid
ls *.fna | sed 's/.fna//g' > namewop
cd ../..
cd fasta/Renamed_plasmid
ls *fna | sed 's/.fna//g' > namep
cd ../..
mv fasta/Renamed_wo_plasmid/namewop ./New.names.wop.tab
mv fasta/Renamed_plasmid/namep ./New.names.p.tab
mkdir prokka/Wo_Plasmid prokka/Plasmid

#Apparently at one point a space introduced himself in the fasta file to correct this:
sed  -i '/^ *$/d' fasta/Renamed_wo_plasmid/HP-*.fna
sed  -i '/^ *$/d' fasta/Renamed_plasmid/HP-*.fna

Nnamewop=$(<New.names.wop.tab)
for i in $Nnamewop; 
	do prokka --cpus 10 --outdir prokka/Wo_Plasmid/$i/ --force --quiet --prefix $i fasta/Renamed_wo_plasmid/$i.fna && echo "$i annotation done";
done
Nnamep=$(<New.names.p.tab)
for i in $Nnamep; 
	do prokka --cpus 10 --outdir prokka/Plasmid/$i/ --force --quiet --prefix $i fasta/Renamed/$i.fna && echo "$i annotation done";
done

####################################
####Ortholog analysis with roary####
####################################

mkdir roary
cd roary
mkdir gff
cp ../prokka/Wo_Plasmid/*/*gff ./gff/

roary -f woplasmid/ -p 10 -e -n -v ./gff/*gff

###########################################
####Second Ortholog analysis with roary####
###########################################

#used -s option for not looking at paralogs

cd roary
mkdir wo_plasmid_gff
cp ../prokka/Wo_Plasmid/*/*gff ./wo_plasmid_gff/

roary -f woplasmid/ -s -p 10 -e -n -v ./wo_plasmid_gff/*gff


#Gene blast DB

#Blast on server

#FastOrtho

#ANI test
mkdir ani
cd ani/
Nname=$(<../New.names.tab)
for i in $Nname; do for j in $Nname; do /tool/enveomics/Scripts/ani.rb -1 $i.fasta -2 $j.fasta -r $i.vs.$j.tab; done; done

#To get two way ANI value from save results
for i in $Sname; do for j in $Sname; do printf "%s \t%s \n" "$i" " $j" "$(grep "Two-way" $i.vs.$j.tab | cut -d " " -f 3)"; done ; done
