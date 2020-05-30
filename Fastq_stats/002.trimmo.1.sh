TRIMMOMATIC_JAR=/thinker/net/raomingming/Software/Trimmomatic-0.36/trimmomatic-0.36.jar

java -jar $TRIMMOMATIC_JAR PE -phred64 \
sequence_1.fastq.gz sequence_2.fastq.gz \
sequence_1_trimmed.fq sequence_1_unpair_trimmed.fq \
sequence_2_trimmed.fq sequence_2_unpair_trimmed.fq \
HEADCROP:0 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#sequence_2_trimmed.fq sequence_2_unpair_trimmed.fq \
#HEADCROP:15 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
