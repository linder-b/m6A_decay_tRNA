#!/bin/bash
#!/bin/awk -f


awk 'BEGIN{OFS="\t";FS="\t";lib=0}ARGIND==1{
	len[$1]=$2;
	hash[$1"\t"$2]=1;
	codon[$1"\t"$2]=$3
}ARGIND==2{
	a=$2/3+1;
	lib=lib+$3;
	if(hash[$1"\t"a]){
		hash1[$1]=hash1[$1]+$3;
		hash2[$1"\t"a]=$3
	}
}ARGIND==3{
	if(hash1[$1]!=""){
		print  $0"\t"hash2[$1"\t"$2]*1000000/lib;
		if(count[$1]!=1){
			sum++;
			count[$1]=1
		}
	}
}END{
	print "'$1'\t"sum >> "./sum_0.5.log"
}'  $2 $1  $2 | \

awk '{
	hash[$1]=hash[$1]""$4","; 
	len[$1]=$2
}END{
	for(i in hash){
		print i"\t"len[i]"\t"hash[i]
	}
}' > $3