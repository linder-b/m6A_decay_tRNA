#!/bin/bash
#!/bin/awk -f


frame_file=$1
offset_file=$2
extend=$3
bed_file=${frame_file/_frame.txt/.bed6}


fp_length=`awk 'BEGIN{FS=OFS="\t"}{

        sum=sum+$3

        hash[$2]=hash[$2]+$3

    }END{

        for(i in hash){

          print i"\t"hash[i]"\t"hash[i]/sum
        }

    }' $frame_file | sort -k2,2rn | \

    awk '{

      if(sum>0.75){

        exit

      }else{

        sum=sum+$3
        print
      }
    }' | cut -f 1 | sort -k1,1n |  paste -s -d ','` 


    if [[ ! $fp_length == "" ]]; then
      
      OLD_IFS="$IFS"

      IFS=$','

      fp_length_array=($fp_length)

      IFS="$OLD_IFS"

      fp_frames=""

      fp_offset_all=""

      for i in "${fp_length_array[@]}"

      do

       fp_frame=`awk '$2=="'$i'"' $frame_file | \

       sort -k1,1n | \

       awk '{

        x=$1

        if(x>=0){

          frame=x%3

        }else if(x%3+3==3){

          frame=0

        }else if(x%3+3==1){

          frame=1

        }else if(x%3+3==2){

          frame=2

        }

        print $0"\t"frame
      }' | \

      awk '{

        hash[$4]=hash[$4]+$3

        sum=sum+$3

      }END{

        for(i in hash){

          if(hash[i]/sum>0.3){

            print i"\t"hash[i]"\t"hash[i]/sum
          }

        }

      }' | cut -f 1 | sort -k1,1n |  paste -s -d ':'`


      if [[ $fp_frames == "" ]]; then



        fp_frames=$fp_frame

      else

        fp_frames=$fp_frames","$fp_frame

      fi
      

      OLD_IFS="$IFS"

      IFS=$':'

      fp_offset_array=($fp_frame)
      
      IFS="$OLD_IFS"

      fp_offsets=""

      for j in "${fp_offset_array[@]}"

      do 

        fp_offset=`awk '$2=="'$i'"'  $frame_file | \

        sort -k1,1n | \

        awk 'BEGIN{FS=OFS="\t"}$1<0{
        
          $3=$3+1
          hash[$1]=$0
        
        }END{
        
          for(i=-"'$extend'";i<=-1;i++){
          
            if(hash[i]){
            
              print hash[i]
              
            }else{
            
              print i"\t'$i'\t1"
            }
          }
        }' | \

        awk '{

          x=$1

          if(x>=0){

            frame=x%3

          }else if(x%3+3==3){

            frame=0

          }else if(x%3+3==1){

            frame=1

          }else if(x%3+3==2){

            frame=2

          }

          print $0"\t"frame

        }' | \

        awk '$4=="'$j'" && $1<0' | \

        awk '{

          hash[NR]=$3

          if(NR==1){

            print $0"\t1"

          }else if(NR>1){

            print $0"\t"$3/hash[NR-1]
          }
        }' | sort -k5,5rn | head -n 1 | awk '{print 4-$1}'`

        if [[ $fp_offsets == ""  ]]; then

          fp_offsets=$fp_offset

        else

          fp_offsets=$fp_offsets":"$fp_offset

          
        fi

      done


       if [[ $fp_offset_all == ""  ]]; then

          fp_offset_all=$fp_offsets

        else

          fp_offset_all=$fp_offset_all","$fp_offsets
          
        fi      

      done


    fi




  echo -e $bed_file"\t"$fp_length"\t"$fp_frames"\t"$fp_offset_all > $offset_file
     
