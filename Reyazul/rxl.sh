c=1

g++ minpath+local_best.cpp -o minpath+local_best

for x in `ls ../testcases`
do  
    # if c < 14 continue
    if [ $c -lt 1 ]
    then
        c=$((c + 1))
        continue
    fi

    echo $c
    start=$(date +%s%N)
    p=$(./minpath+local_best < "../testcases/$x")  # Capture the program output
    end=$(date +%s%N)

    # Extract the value after 'VALUE' from the program output
    value=$(echo "$p" | awk '/^VALUE/ {print $2}')

    # Extract the corresponding opt value from track1.csv
    testcase_name=$(basename "$x")  # Get the test case file name
    opt=$(awk -F, -v t="$testcase_name" '$1 ~ t {gsub(/[[:space:]]/, "", $2); print $2}' ../ground_truth/track1.csv)

    # Calculate the approximation ratio
    if [ -n "$opt" ] && [ "$opt" -ne 0 ]; then
        ratio=$(echo "scale=6; $value / $opt" | bc)
    else
        ratio="N/A"  # Handle cases where opt is missing or zero
    fi

    input_file="../testcases/$x"
    output_file="output_minpath+local_best.csv"

    echo -n "$c," >> "$output_file"
    awk 'NR==2 {printf "%s,", $2} NR==3 {printf "%s,", $2} /Terminals / {printf "%s,", $2}' "$input_file" >> "$output_file"
    echo -n "$value," >> "$output_file"  # Append the extracted VALUE
    echo -n "$ratio," >> "$output_file"  # Append the approximation ratio
    echo "$(($(($end-$start))/1000000)) ms" >> "$output_file"

    c=$((c + 1))
done
