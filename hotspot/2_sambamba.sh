bam=$1
ref=$2
tmpbam=$(mktemp --suffix '.bam')
tmpgbbam=$(mktemp --suffix '.bam')

fastq=$(echo $bam | cut -d . -f 1 | awk '{print $1".fastq.gz"}')
echo $(zcat $fastq | wc -l)/4 | bc 

sambamba view -q -c $bam
sambamba markdup -q -r $bam $tmpbam &> /dev/null 
sambamba view -q -c $tmpbam 
sambamba view -q -c -F 'mapping_quality>=20' $tmpbam
sambamba view -q -c -F 'mapping_quality>=20 and ref_name=~/^chr(\d+|X|Y)$/' $tmpbam
sambamba view -q -f bam -F 'mapping_quality>=20 and ref_name=~/^chr(\d+|X|Y)$/ and sequence_length>=21 and sequence_length <=31' -o $tmpgbbam $tmpbam 
sambamba view -q -c $tmpgbbam
sambamba slice -q -L $2 $tmpgbbam | wc -l


rm -f $tmpbam* $tmpgbbam*

