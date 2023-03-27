#!/bin/bash
filename='samples_SE.txt'
echo Start
while read p; do
        echo $p
        alignment_rate="$(grep -oPm1 'Mapping rate = \K\w+' "./counts/SE/${p}/logs/salmon_quant.log")";
        libr_type="$(grep -oPm1 'most likely library type as \K\w+' "./counts/SE/${p}/logs/salmon_quant.log")";
        #lib_type="$(grep -oPm1 '"expected_format": \K\w+' "./counts_SE/${p}/lib_format_counts.json")";
        num_seq="$(grep -oPm1 'Total Sequences</td><td>\K\w+' "./fastqc_files/SE/${p}_trimmed_fastqc.html")";
        printf "%s\t%s\t%s\t%s\n" "${p}" "${num_seq}" "${alignment_rate}" "${libr_type}" >> ./alignment_info_SE.txt
done < $filename  
