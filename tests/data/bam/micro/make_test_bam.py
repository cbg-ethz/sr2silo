# This python script generates the test BAM file from
# the SAM / BAM specifications: https://samtools.github.io/hts-specs/SAMv1.pdf

"""
Coor  12345678901234 *5678901234567890123456789012345
ref   AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT

+r001/1     TTAGATAAAGGATA*CTG
+r002      aaaAGATAA*GGATA
+r003    gcctaAGCTAA
+r004                  ATAGCT..............TCAGC
-r003                         ttagctTAGGC
-r001/2                                     CAGCGGCAT


leads to a SAM file with the following content:


@HD VN:1.6 SO:coordinate
@SQ SN:ref LN:45
r001    99 ref  7 30 8M2I4M1D3M * 37  39 TTAGATAAAGGATACTG *
r002     0 ref  9 30 3S6M1P1I4M *  0   0 AAAAGATAAGGATA     * SA:Z:ref,29,-,6H5M,17,0;
r003     0 ref  9 30 5S6M       *  0   0 GCCTAAGCTAA        *
r004  2064 ref 16 30 6M14N5M    *  0   0 ATAGCTTCAGC        *
r003  2064 ref 29 17 6H5M       *  0   0 TAGGC              * SA:Z:ref,9,+,5S6M,30,1;
r001   147 ref 37 30 9M         *  7 -39 CAGCGGCAT          * NM:i:1
"""

from __future__ import annotations

import pysam

# lets create a sam file out of the above
with open("test.sam", "w") as f:
    f.write(
        """@HD VN:1.6 SO:coordinate
@SQ SN:ref LN:45
r001    99  ref 7 30    8M2I4M1D3M = 37  39 TTAGATAAAGGATACTG   *
r002    0   ref 9 30    3S6M1P1I4M *  0   0 AAAAGATAAGGATA  * SA:Z:ref,29,-,6H5M,17,0;
r003    0   ref 9 30    5S6M       *  0   0 GCCTAAGCTAA *
r004    2064    ref 16  30 6M14N5M    *  0   0 ATAGCTTCAGC  *
r003    2064    ref 29  17 6H5M       *  0   0 TAGGC    * SA:Z:ref,9,+,5S6M,30,1;
r001    147 ref 37  30 9M         =  7 -39 CAGCGGCAT    * NM:i:1
"""
    )

tmpfilename = "test.bam"

header = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "ref", "LN": 45}],
}

with pysam.AlignmentFile(tmpfilename, "wb", header=header) as outf:
    a = pysam.AlignedSegment()
    a.query_name = "r001"
    a.query_sequence = "TTAGATAAAGGATACTG"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 7
    a.mapping_quality = 30
    a.cigar = [(0, 8), (1, 2), (0, 4), (2, 1), (0, 3)]  # Fixed CIGAR: 8M2I4M1D3M
    a.next_reference_id = 0
    a.next_reference_start = 37
    a.template_length = 39
    a.query_qualities = [0] * 17
    a.tags = []
    outf.write(a)

    b = pysam.AlignedSegment()
    b.query_name = "r002"
    b.query_sequence = "AAAAGATAAGGATA"
    b.flag = 0
    b.reference_id = 0
    b.reference_start = 9
    b.mapping_quality = 30
    b.cigar = [(4, 3), (0, 6), (6, 1), (1, 1), (0, 4)]
    b.next_reference_id = 0
    b.next_reference_start = 0
    b.template_length = 0
    b.query_qualities = [0] * 14
    outf.write(b)

    c = pysam.AlignedSegment()
    c.query_name = "r003"
    c.query_sequence = "GCCTAAGCTAA"
    c.flag = 0
    c.reference_id = 0
    c.reference_start = 9
    c.mapping_quality = 30
    c.cigar = [(4, 5), (0, 6)]
    c.next_reference_id = 0
    c.next_reference_start = 0
    c.template_length = 0
    c.query_qualities = [0] * 11
    outf.write(c)

    d = pysam.AlignedSegment()
    d.query_name = "r004"
    d.query_sequence = "ATAGCTTCAGC"
    d.flag = 2064
    d.reference_id = 0
    d.reference_start = 16
    d.mapping_quality = 30
    d.cigar = [(0, 6), (3, 14), (0, 5)]
    d.next_reference_id = 0
    d.next_reference_start = 0
    d.template_length = 0
    d.query_qualities = [0] * 11
    outf.write(d)

    e = pysam.AlignedSegment()
    e.query_name = "r003"
    e.query_sequence = "TAGGC"
    e.flag = 2064
    e.reference_id = 0
    e.reference_start = 29
    e.mapping_quality = 17
    e.cigar = [(5, 6), (0, 5)]
    e.next_reference_id = 0
    e.next_reference_start = 0
    e.template_length = 0
    e.query_qualities = [0] * 5
    e.tags = [("SA", "Z:ref,9,+,5S6M,30,1")]
    outf.write(e)

    f = pysam.AlignedSegment()
    f.query_name = "r001"
    f.query_sequence = "CAGCGGCAT"
    f.flag = 147
    f.reference_id = 0
    f.reference_start = 37
    f.mapping_quality = 30
    f.cigar = [(0, 9)]
    f.next_reference_id = 0
    f.next_reference_start = 7
    f.template_length = -39
    f.query_qualities = [0] * 9
    f.tags = [("NM", "i:1")]
    outf.write(f)

# close the file
outf.close()


# now read the bam file and print the reads

index = pysam.index(tmpfilename)


# read in bam and print out the reads
with pysam.AlignmentFile(tmpfilename, "rb") as f:
    for read in f.fetch():
        print(
            read.query_name,
            read.query_sequence,
            read.query_qualities,
            read.flag,
            read.reference_id,
            read.reference_start,
            read.mapping_quality,
            read.cigarstring,
            read.next_reference_id,
            read.next_reference_start,
            read.template_length,
        )
        print("tags:", read.tags)
        print("query_qualities:", read.query_qualities)
        print("query_length:", read.query_length)
        print("query_alignment_length:", read.query_alignment_length)
        print("query_alignment_start:", read.query_alignment_start)
        print("query_alignment_end:", read.query_alignment_end)
        print("reference_length:", read.reference_length)
        print("reference_start:", read.reference_start)
        print("reference_end:", read.reference_end)
