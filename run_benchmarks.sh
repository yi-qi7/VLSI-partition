#!/bin/bash

BENCHMARK_DIR="ISPD_benchmark"
RESULT_FILE="benchmark_results.txt"

echo "========================================" | tee $RESULT_FILE
echo "Benchmark Test Results" | tee -a $RESULT_FILE
echo "Date: $(date)" | tee -a $RESULT_FILE
echo "========================================" | tee -a $RESULT_FILE
echo "" | tee -a $RESULT_FILE

printf "%-15s %-15s %-15s\n" "Benchmark" "Cut Size" "Partition Ratio" | tee -a $RESULT_FILE
echo "----------------------------------------" | tee -a $RESULT_FILE

for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18; do
    benchmark_file="${BENCHMARK_DIR}/ibm${i}.hgr"
    
    if [ -f "$benchmark_file" ]; then
        echo ""
        echo "Processing ibm${i}.hgr..."
        
        output=$(./main "$benchmark_file" 2>&1)
        
        cut_size=$(echo "$output" | grep "Cut size:" | awk '{print $3}')
        ratio_x=$(echo "$output" | grep "Ratio X:" | awk '{print $3}')
        ratio_y=$(echo "$output" | grep "Ratio Y:" | awk '{print $3}')
        
        if [ -n "$cut_size" ]; then
            printf "%-15s %-15s X: %-6s Y: %-6s\n" "ibm${i}.hgr" "$cut_size" "$ratio_x" "$ratio_y" | tee -a $RESULT_FILE
        else
            printf "%-15s %-15s\n" "ibm${i}.hgr" "ERROR" | tee -a $RESULT_FILE
        fi
    else
        printf "%-15s %-15s\n" "ibm${i}.hgr" "NOT FOUND" | tee -a $RESULT_FILE
    fi
done

echo "" | tee -a $RESULT_FILE
echo "========================================" | tee -a $RESULT_FILE
echo "Test completed. Results saved to $RESULT_FILE" | tee -a $RESULT_FILE
echo "========================================" | tee -a $RESULT_FILE
