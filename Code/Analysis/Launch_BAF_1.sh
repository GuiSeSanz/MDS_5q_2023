# nohup samtools view /home/tereshkova/data/gserranos/MDS/Data/BAM_files/SMD132114579_possorted_genome_bam.bam | /opt/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin /home/tereshkova/data/gserranos/MDS/Data/Annotation/hg38.list /home/tereshkova/data/gserranos/MDS/Data/BAM_files/SMD132114579 50 0 > nohup_SMD132114579.out &
# nohup samtools view /home/tereshkova/data/gserranos/MDS/Data/BAM_files/SMD211420_possorted_genome_bam.bam | /opt/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin /home/tereshkova/data/gserranos/MDS/Data/Annotation/hg38.list /home/tereshkova/data/gserranos/MDS/Data/BAM_files/SMD211420 50 0 > nohup_SMD211420.out &
nohup samtools view /home/tereshkova/data/gserranos/MDS/Data/FS-0634-post/outs/possorted_genome_bam.bam | /opt/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin /home/tereshkova/data/gserranos/MDS/Data/Annotation/hg38.list /home/tereshkova/data/gserranos/MDS/Data/BAM_files/FS-0634-post 50 0 > nohup_FS-0634-post.out &
nohup samtools view /home/tereshkova/data/gserranos/MDS/Data/FS-0406-post/outs/possorted_genome_bam.bam | /opt/BAFExtract/bin/BAFExtract -generate_compressed_pileup_per_SAM stdin /home/tereshkova/data/gserranos/MDS/Data/Annotation/hg38.list /home/tereshkova/data/gserranos/MDS/Data/BAM_files/FS-0406-post 50 0 > nohup_FS-0406-post.out &