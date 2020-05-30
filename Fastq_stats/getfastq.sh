##
zcat /data_NCCL_02/IlluminaHiSeq/Illumina_pt2_B1700_R1.fq.gz  | head -90000 | gzip > sequence_1.fastq.gz
zcat /data_NCCL_02/IlluminaHiSeq/Illumina_pt2_B1700_R2.fq.gz  | head -90000 | gzip > sequence_2.fastq.gz


##awk
awk '{if(length($0)<100)print}' sequence_1_trimmed.fq | less
awk '{if(length($0)<100)print(length($0))}' sequence_1_trimmed.fq | sort

grep -B 3 "TTGATAGGCATGGACTCGCTGCTCATCCTTGCCTCTAGCTGCCTATTGCG" sequence_1_trimmed.fq
