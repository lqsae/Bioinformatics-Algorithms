import pysam
from pathlib import Path

def sumMatchesAndMismatches(segment):
    """
    Get total matches/mismatches from CIGAR string (M field)
    Code dictionary:
    M    BAM_CMATCH    0
    I    BAM_CINS      1
    D    BAM_CDEL      2
    N    BAM_CREF_SKIP 3
    S    BAM_CSOFT_CLIP 4
    H    BAM_CHARD_CLIP 5
    P    BAM_CPAD      6
    =    BAM_CEQUAL    7
    X    BAM_CDIFF     8
    B    BAM_CBACK     9
    """
    return sum([value for (code, value) in segment.cigartuples if code == 0])

def getNumberOfMatches(segment):
    """
    Get the number of matches from alignment
    Do not consider insertion/deletion as mismatches
    """
    parsed_MD = segment.get_aligned_pairs(with_seq=True)
    return len([
        base for (read_pos, ref_pos, base) in parsed_MD 
        if ((base is not None and base.isupper()) and read_pos is not None)
    ])

def percent_identity(segment):
    """
    Compute percent identity from MD tag of aligned segment.
    segment: pysam AlignedSegment object.
    """
    return 100 * (getNumberOfMatches(segment) / sumMatchesAndMismatches(segment))

def calculate_identity(bam_file):
    # 打开BAM文件
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    identities = []

    # 遍历每个比对
    for segment in bam:
        if segment.is_unmapped:
            continue
        
        identity = percent_identity(segment)
        identities.append(identity)
    
    bam.close()
    return identities

# 示例使用
bam_file_path = "your_alignment.bam"
identities = calculate_identity(bam_file_path)

# 输出每个比对的identity
for idx, identity in enumerate(identities):
    print(f"Alignment {idx+1}: Identity = {identity:.2f}%")
