#!/bin/bash
#!/bin/awk -f

fasta=$1
output=$2
extend=$3
exclude=$4


 awk 'BEGIN{OFS="\t"}{

        if($1~/>/){

          split($1,a,">")
          hash[NR]=a[2]

        }else if(hash[NR-1]){
          
          len=length($0)
          
          for(i="'$3'"+1+3*"'$exclude'";i<=len-"'$extend'"-3*"'$exclude'";i=i+3){
            
            codon=substr($0,i,3)
            gsub(/t/,"T",codon)
            gsub(/g/,"G",codon)
            gsub(/c/,"C",codon)
            gsub(/a/,"A",codon)
            print hash[NR-1],(i-"'$extend'"-1)/3+1,codon
          
          }
        
        }
      
      }' $fasta > $output # a file that contains the codons and their locations with both end excluded.