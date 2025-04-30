#!/bin/bash
#!/bin/awk -f



help_info(){
  echo "
 **********************************************************************
*   _____    _____      _____     _____      _____   __    ___   ____  *
*  |  __  \ |  __  \   |  __  \  /  __  \   /  ___| |  |  /  /  /  __| *
*  | |__\ | | |__\ |   | |__\ | |  /  \  | |  /     |  | /  /  |  /    *
*  |   ___| |  ____|   |   ___| | |    | | |  |     |  |/  /    \  \   *
*  |   \    | |        |   \    | |    | | |  |     |      \     \  \  *
*  | |\ \   | |        | |\ \   |  \__/  | |  |___  |  |\   \   __/  / *
*  |_| \_\  |_|        |_| \_\   \ ____ /   \ ____| |__| \___\ |____/  *
*                                                                      *
*                                                                      *
*          Copyright by Jie Wu   2017.01.23                            *
*                                                                      *
 **********************************************************************
  "
  echo "./rp_analysis_for_short_reads.sh : tool for analyzing Ribosome profiling data and make Codon plot and/or Wave plot"
  echo "Usage: ./rp_analysis_for_short_reads.sh -M [Mode] [Options]"

  echo "
             [ -M|--mode ]                = Analyzation mode, choose one from \"Mapping\", \"CodonPlot\" or \"CodonPlot,WavePlot\". 
             [ -H|--help|-h  ]               = Help information.

      1) When [Mode] is \"Mapping\", Options: 
      
             [ -Q|--fastq ]               = Input fastq file.
             [ -T|--transcriptom-index ]  = Prefix of bowtie index of cds transcriptome. 
             [ -O|--output ]              = Name of output folder. 
             [ -X|--max_length ]          = Maximum read length to check the frame information.
             [ -I|--min_length ]          = Minimum read length to check the frame information.
             [ -E|--extend ]              = Length extended flanking annotated CDS, which was built as the transcriptome index.
      2) When [Mode] is \"CodonPlot\", Options:

             [ -A|--fasta ]               = Fasta file that were used to build the transcriptome index.
             [ -D|--index ]               = Index file (7 columns seperated by tab including sam/bed6 file name, sample name, replicates, read length, frame, A site, class). eg. \"CodonPlot.index\".
             [ -O|--output ]              = Name of output folder.
             [ -E|--extend ]              = Length extended flanking annotated CDS, which was built as the transcriptome index.       
             [ -F|--fexcluded ]           = Codon length excluded at the 5' end of each transcript when plotting the codon plot. 
             [ -B|--bexcluded ]           = Codon length excluded at the 3' end of each transcript when plotting the codon plot. 

      3) When [Mode] is \"WavePlot\" or \"CodonPlot,WavePlot\", Options:
             [ -A|--fasta ]               = Fasta file that were used to build the transcriptome index.
             [ -D|--index ]               = If \"WavePlot\", use index file like \"WavePlot.index\" (6 columns seperated by tab).
                                            If \"CodonPlot,WavePlot\", use index file like \"CodonPlot_WavePlot.index\" (7 columns seperated by tab).
             [ -O|--output ]              = Name of output folder.
             [ -E|--extend ]              = Length extended flanking annotated CDS, which was built as the transcriptome index.       
             [ -F|--fexcluded ]           = Codon length excluded at the 5' end of each transcript when plotting the codon plot and wave plot. 
             [ -B|--bexcluded ]           = Codon length excluded at the 3' end of each transcript when plotting the codon plot and wave plot. 
             [ -C|--codons ]              = codons seperated by "," to make the wave plot, eg. \"CGA,CAA,ACC\".
             [ -S|--offset ]              = offset file. Only used when offset file is available and only -S, -O, -C is needed. eg. \"all_codon_Asite_offset.txt\".
  " 
  echo "
  eg. 1) ./rp_analysis_for_short_reads.sh -M Mapping -Q fastq_file -T index_prefix -O output -X 31 -I 28 -E 18
      2) ./rp_analysis_for_short_reads.sh -M CodonPlot -A transcriptome_fasta_file -D index_file -O output -E 18 -F 15 -B 15 
      3) ./rp_analysis_for_short_reads.sh -M WavePlot -A  transcriptome_fasta_file -D index_file -O output -E 18 -F 15 -B 15 -C CGA,CGG,CAA
      4) ./rp_analysis_for_short_reads.sh -M CodonPlot,WavePlot -A  transcriptome_fasta_file -D index_file -O output -E 18 -F 15 -B 15 -C CGA,CGG,CAA
      5) ./rp_analysis_for_short_reads.sh -S offset_file -C CGA,CGG,CAA -O output
" 
  if ! [ -x "$(command -v R)" ];then
     echo -e "Warnings: \n\tR is not installed in this server so pdf files can't be made. Please install it or use the files from the results such as \"all_codon_Asite_nor_wt.txt\", \"all_codon_Asite_offset.txt\", to make the plot in Rstudio."

  fi

}


if [ $# -eq 0 ]
then
  help_info
fi

TEMP=`getopt -o M:Q:T:O:L:A:X:I:D:E:S:C:F:B:f:a:Hh --long mode:,fastq:,transcriptome_index:,output:,length:,fasta:,max_length:,min_length:,index:,extend:,fexcluded:,bexcluded:,offset:,codons:,frame:,Asite:,help  -- "$@"`
eval set -- "$TEMP"


while true
do
        case "$1" in
    -M|--mode)
          ARG_mode=$2
          shift 2
        ;;
    -Q|--fastq)
          ARG_fastq=$2
          shift 2
        ;;
    -T|--transcriptom-index)
          ARG_transcriptome_index=$2
          shift 2
        ;;
    -O|--output)
          ARG_output=$2
          shift 2
        ;;
    -A|--fasta)
          ARG_fasta=$2
          shift 2
        ;;
    -X|--max_length)
          ARG_max_length=$2
          shift 2
        ;;
    -I|--min_length)
          ARG_min_length=$2
          shift 2
        ;;
    -D|--index)
         ARG_index=$2
         shift 2
        ;;
    -E|--extend)
          ARG_extend=$2
          shift 2
        ;;
    -F|--fexcluded)
          ARG_exclude_5end=$2
          shift 2
        ;;
    -B|--bexcluded)
          ARG_exclude_3end=$2
          shift 2
        ;;
    -S|--offset)
          ARG_offsetfile=$2
          shift 2
        ;;
    -C|--codons)
          ARG_codons=$2
          shift 2
        ;;
    -f|--frame)
          ARG_frame=$2
          shift 2
        ;;
    -a|--asite)
          ARG_frame=$2
          shift 2
        ;;
    -H|-h|--help)
          help_info
          shift
        ;;

    --)
          shift
          break
        ;;
    *)
          echo "Internal error!"
          #help_info
          exit 1
        ;;
        esac
done
#echo $ARG_fastq
if [[ $ARG_output != "" ]]; then

   if [[ ! -d $ARG_output ]]; then
  
      mkdir $ARG_output

   fi
fi




if [[ $ARG_mode == "Mapping" ]]; then
    
    Start_mapping=`date +%s`
    
    echo -e "\nMapping......"
    name=${ARG_fastq/\.fastq}
    name=${name/*\//}

   bowtie -p 16 -v 1 -m 1 --no-unal --norc  --best --strata -t $ARG_transcriptome_index -q $ARG_fastq -S $ARG_output/$name"_mapped.sam"  2>$ARG_output/$name"_mapped.log"
   

   awk  'BEGIN{FS="\t";OFS="\t"} $2==0{
   
     sub(/M/,"",$6)
     if($0~/MD:Z:0[A-Z]/){              
         p=substr($10,2,length($10))
         print $3,$4-1+1,$4-1+$6,"@"$1"|"p,"0","+"    # reads with the first unmapped nucleotide.
       }else{
         print $3,$4-1,$4-1+$6,"@"$1"|"$10,"0","+"
     }
   
   }' $ARG_output/$name"_mapped.sam" > $ARG_output/$name".bed6" #produce bed6 file for Codon plot and Wave plot. 
    
    echo -e "\nProducing $name frame file..."
    awk 'BEGIN{OFS="\t";FS="\t"}$2==0{                                                                                                                                                                                                                                                                                                                                                                                                                                                    

      sub(/M/,"",$6)
      if($6>="'$ARG_min_length'" && $6<="'$ARG_max_length'"){
          print $4-"'$ARG_extend'"-1,$6
      }
    
    }' $ARG_output/$name"_mapped.sam" | awk 'BEGIN{OFS="\t";FS="\t"}{hash[$0]++}END{for(i in hash){print i,hash[i]}}'  > $ARG_frame # read counts in frame
    
    if ! [ -x "$(command -v R)" ];then
       echo -e "Warnings: \n\tR is not installed in this server so pdf files can't be made. Please install it or use the files from the results such as files with suffix \"_frame.txt\" to make the plot in Rstudio."

    else
       dir=${0/rp_analysis.sh/}
       echo $dir
       chmod +x $dir/frame_plot.r
       R --vanilla --slave --args $ARG_frame $ARG_max_length $ARG_min_length < $dir/frame_plot.r  # making the frame plot. see the pdf file. 
    fi
    End_mapping=`date +%s`
    Time_mapping=$((End_mapping-Start_mapping))
    echo "Took $((Time_mapping/60)) mins to map and produce frame plot"
fi




if [[ $ARG_mode =~ "WavePlot" ]]; then

   Start_plotting=`date +%s`

    
      
      
      echo -e "\n Assigning ribosome occupancy of A site...... "
      IFS=$'\n'
      fasta_filename=${ARG_fasta/\.fa*/}
      fasta_filename=${fasta_filename/*\//}
     
     
      for i in `cat $ARG_index`;do
     
        OLD_IFS="$IFS"
        IFS=$'\t'
        array=($i)
        IFS="$OLD_IFS"
        if [[ ${array[0]} =~ ".sam" ]]; then
              #statements
           file_name=${array[0]/\.sam/}
           file_name1=${file_name/*\//}
           awk  'BEGIN{OFS="\t";FS="\t"} $2==0{
 
             sub(/M/,"",$6)
             if($0~/MD:Z:0[A-Z]/){            
               p=substr($10,2,length($10))
               print $3,$4-1+1,$4-1+$6,"@"$1"|"p,"0","+"    # reads with the first unmapped nucleotide. 
             }else{
               print $3,$4-1,$4-1+$6,"@"$1"|"$10,"0","+"
             } 
           
           }' $file_name".sam" > $file_name".bed6" # producing bed file.

        elif [[ ${array[0]} =~ ".bed6" ]]; then
          #statements
            file_name=${array[0]/\.bed6/}  
            file_name1=${file_name/*\//}
     
        fi
        echo -e "\n\tProcessing ${array[0]} ......"
         awk 'BEGIN{
           OFS="\t"
          
           split("'${array[1]}'",length_array,",")
          
           for(len=1;len<=length(length_array);len++){
             hash_len[length_array[len]]=1
           
           }
           split("'${array[2]}'",frame_array,",")
           split("'${array[3]}'",A_site_array,",")
           for(i=1;i<=length(frame_array);i++){
             split(frame_array[i],frame_length,":")
             split(A_site_array[i],A_site_length,":")
             for(j=1;j<=length(frame_length);j++ ){
               frame_length_hash[length_array[i]"\t"frame_length[j]]=A_site_length[j]
             }
           }
         }NR==FNR{
               if($1~/>/){
                 sub(/>/,"",$1)
                 line[NR]=$1
               }else{
                 if(line[NR-1]){
                   length_gene[line[NR-1]]=length($0)
                 }
               }
                
         }NR!=FNR{
         
           if(hash_len[$3-$2]){
     
                 x=$2-"'$ARG_extend'"
                 if(x<0){
         
                   if(x%3+3==3){frame=0}else if(x%3+3==1){frame=1}else if(x%3+3==2){frame=2}
                 }else if(x>=0){
                   if(x%3==0){frame=0}else if(x%3==1){frame=1}else if(x%3==2){frame=2}
                 }
                 split($4,a,"|")
                 if($3-$2<60){
                    for(j=$3+1;j<=$2+60;j++){
                      a[2]=a[2]""hash_codon[$1"\t"j]
                    }
                 }
                 if(frame==0 && frame_length_hash[$3-$2"\t0"]!=""){
                       if(x+frame_length_hash[$3-$2"\t0"]-1>= 3*"'$ARG_exclude_5end'" && x+frame_length_hash[$3-$2"\t0"]-1<=length_gene[$1]-2*"'$ARG_extend'"-1-3*"'$ARG_exclude_3end'"){
                           print $0,substr(a[2],frame_length_hash[$3-$2"\t0"],3),substr(a[2],frame_length_hash[$3-$2"\t0"]+24,3),substr(a[2],frame_length_hash[$3-$2"\t0"]+27,3),substr(a[2],frame_length_hash[$3-$2"\t0"]+30,3),x+frame_length_hash[$3-$2"\t0"]-1
                        }
                 }else if(frame==2 && frame_length_hash[$3-$2"\t2"]!=""){
                       if(x+frame_length_hash[$3-$2"\t2"]-1>=3*"'$ARG_exclude_5end'" && x+frame_length_hash[$3-$2"\t2"]-1<=length_gene[$1]-2*"'$ARG_extend'"-1-3*"'$ARG_exclude_3end'"){
                           print $0,substr(a[2],frame_length_hash[$3-$2"\t2"],3),substr(a[2],frame_length_hash[$3-$2"\t2"]+24,3),substr(a[2],frame_length_hash[$3-$2"\t2"]+27,3),substr(a[2],frame_length_hash[$3-$2"\t2"]+30,3),x+frame_length_hash[$3-$2"\t2"]-1
                       }
                }else if(frame==1 && frame_length_hash[$3-$2"\t1"]!=""){
                      if(x+frame_length_hash[$3-$2"\t1"]-1>=3*"'$ARG_exclude_5end'" && x+frame_length_hash[$3-$2"\t1"]-1<=length_gene[$1]-2*"'$ARG_extend'"-1-3*"'$ARG_exclude_3end'"){  
                            print $0,substr(a[2],frame_length_hash[$3-$2"\t1"],3),substr(a[2],frame_length_hash[$3-$2"\t1"]+24,3),substr(a[2],frame_length_hash[$3-$2"\t1"]+27,3),substr(a[2],frame_length_hash[$3-$2"\t1"]+30,3),x+frame_length_hash[$3-$2"\t1"]-1  # only useful when 5 end locates at -11 in the frame plot.   
                      }
                }
            }
         
         }' $ARG_fasta $file_name".bed6" > $ARG_output/$file_name1"_codon_Asite.txt"
       
  
        awk 'BEGIN{OFS="\t";FS="\t"}{
           hash[$1"\t"$11]++
         }END{for(i in hash){
                 print i,hash[i]
               }
        }' $ARG_output/$file_name1"_codon_Asite.txt" > $ARG_output/$file_name1"_codon_Asite_count.txt"
     done
  End_plotting=`date +%s` 
  Time_plotting=$((End_plotting-Start_plotting))
  echo "Took $((Time_plotting/60)) mins to make codon plot and wave plot"
fi
