c=1

g++ dreyfus-wagner.cpp -O3 -o dreyfus-wagner

for x in `ls testcases`
do  
    # if c < 14 continue
    if [ $c -lt 14 ]
    then
        c=$((c + 1))
        continue
    fi
    echo $c
    start=$(date +%s%N)
    p=`./dreyfus-wagner < "testcases/$x"`
    end=$(date +%s%N)
    input_file="testcases/$x"
    output_file="output_dreyfus.csv"
    echo -n "$c," >> "$output_file"
    awk 'NR==2 {printf "%s,", $2} NR==3 {printf "%s,", $2} /Terminals / {printf "%s,", $2}' "$input_file" >> "$output_file"
    echo "$(($(($end-$start))/1000000)) ms" >> "$output_file"
    c=$((c + 1))
done 