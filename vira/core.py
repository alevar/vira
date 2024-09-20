import os
import re
import argparse
import subprocess

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
# from Bio import substitution_matrices

from utils.common import *

from classes.txgroup import Transcriptome, Gene, Bundle
from classes.transcript import Transcript

class Vira:
    def __init__(self, args):
        if gtf_or_gff(args.annotation) is None:
            raise ValueError(f"{args.annotation} is not a valid GTF/GFF file.")
        
        if not os.path.exists(args.genome):
            raise FileNotFoundError(f"Genome file {args.genome} not found.")
        
        if not os.path.exists(args.target):
            raise FileNotFoundError(f"Input file {args.target} not found.")
        
        if args.guide and not os.path.exists(args.guide):
            raise FileNotFoundError(f"Guide annotation file {args.guide} not found.")
        
        # INPUT FILES
        self.annotation = args.annotation # reference annotation
        self.genome = args.genome # reference genome
        self.target = args.target # target genome
        self.guide = args.guide # guide annotation of the target genome - used to verify the inferred target annotation
        self.output = args.output
        self.gtf = gtf_or_gff(args.annotation)

        # OPTIONS
        self.keep_tmp = args.keep_tmp
        self.tmp_dir = standard_path(args.tmp_dir)
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        # TOOLS
        self.minimap2 = args.minimap2
        self.gffread = args.gffread
        self.sam2gtf = args.sam2gtf
        self.miniprot = args.miniprot
        
        # TMP FILES
        self.cds_nt_fasta_fname = self.tmp_dir+"cds_nt.fasta"
        self.exon_nt_fasta_fname = self.tmp_dir+"exon_nt.fasta"
        self.cds_sam_fname = self.tmp_dir+"cds_nt.sam"
        self.exon_sam_fname = self.tmp_dir+"exon_nt.sam"
        self.exon_sam2gtf_fname = self.tmp_dir+"exon_nt.sam2gtf.gtf"
        self.cds_sam2gtf_fname = self.tmp_dir+"cds_nt.sam2gtf.gtf"

        self.check_tools()

    def check_tools(self):
        if subprocess.call(f"command -v {self.minimap2}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"minimap2 is not installed or not available in PATH. Please install minimap2 before running vira. Installation instructions can be found at: https://github.com/lh3/minimap2")
        
        if subprocess.call(f"command -v {self.gffread}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"gffread is not installed or not available in PATH. Please install gffread before running vira. Installation instructions can be found at: https://github.com/gpertea/gffread")
        
        if subprocess.call(f"command -v {self.sam2gtf}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"sam2gtf is not installed or not available in PATH. Please install sam2gtf before running vira. Installation instructions can be found at: https://github.com/alevar/sam2gtf")
        
        # if subprocess.call(f"command -v {self.miniprot}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        #     raise EnvironmentError(f"miniprot is not installed or not available in PATH. Please install miniprot before running vira. Installation instructions can be found at: https://github.com/lh3/miniprot")

    def write_output(self):
        print("Writing output file...")

    def run(self):
        # extract the reference transcript and protein sequences
        cmd = [self.gffread,"-g",self.genome,"-x",self.cds_nt_fasta_fname,"-w",self.exon_nt_fasta_fname,self.annotation]
        print(f"Extracting reference transcript and protein sequences: {' '.join(cmd)}")
        subprocess.call(cmd)

        # run minimap of transcripts and proteins to the CRF genome
        cmd = [self.minimap2,"--for-only","-ax","splice",self.target,self.cds_nt_fasta_fname]
        print(" ".join(cmd)+" > "+self.cds_sam_fname)
        with open(self.cds_sam_fname,"w+") as outFP:
            subprocess.call(cmd,stdout=outFP)

        cmd = [self.minimap2,"--for-only","-ax","splice",self.target,self.exon_nt_fasta_fname]
        print(" ".join(cmd)+" > "+self.exon_sam_fname)
        with open(self.exon_sam_fname,"w+") as outFP:
            subprocess.call(cmd,stdout=outFP)

        # run sam2gtf
        cmd = [self.sam2gtf,
            "-i",self.exon_sam_fname,
            "-o",self.exon_sam2gtf_fname]
        print(" ".join(cmd))
        subprocess.call(cmd)

        cmd = [self.sam2gtf,
            "-i",self.cds_sam_fname,
            "-o",self.cds_sam2gtf_fname]
        print(" ".join(cmd))
        subprocess.call(cmd)

        # combine annotated transcripts, CDSs and guide annotation together
        # for each transcript/cds annotate any differences
        self.build()
        
    def compare_intron_chains(self, ref_tome, target_tome):
        donor_map = {} # holds the mapping between reference and target donor sites
        acceptor_map = {} # holds the mapping between reference and target acceptor sites
        for ref_tx in ref_tome.transcript_it():
            for ref_i,ref_intron in enumerate(ref_tx.introns_intron()):
                # find position of the intron in the target genome
                
                read_tid = ref_tx.get_attr("read_name")
                target_tx = target_tome.get_by_tid(read_tid)
                
                
        
        return

    def process_cigar(self, cigar_string, ref_start, tx_start):
        """
        Process CIGAR string to create a mapping from reference genome positions
        to target (transcript) positions.

        :param cigar_string: CIGAR string from alignment
        :param ref_start: Start position on the reference genome (1-based)
        :param tx_start: Start position on the transcript (1-based)
        :return: A dictionary mapping reference positions to transcript positions
        """
        # CIGAR operation regex
        cigar_operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
        
        # Initialize positions
        ref_pos = ref_start
        tx_pos = tx_start

        # Maps to track reference to query and query to reference
        ref_to_tx_map = {}
        tx_to_ref_map = {}

        # Process each CIGAR operation
        for length, op in cigar_operations:
            length = int(length)

            if op == 'M' or op == '=' or op == 'X':  # Match/Mismatch
                for _ in range(length):
                    # Map reference to transcript and transcript to reference
                    ref_to_tx_map[ref_pos] = tx_pos
                    tx_to_ref_map[tx_pos] = ref_pos
                    ref_pos += 1
                    tx_pos += 1
            elif op == 'I':  # Insertion in query (relative to the reference)
                # Insertions affect only the transcript position
                for _ in range(length):
                    tx_to_ref_map[tx_pos] = None  # No reference mapping for inserted bases
                    tx_pos += 1
            elif op == 'D' or op == 'N':  # Deletion or skipped region in the reference
                # Deletions affect only the reference position
                for _ in range(length):
                    ref_to_tx_map[ref_pos] = None  # No transcript mapping for deleted bases
                    ref_pos += 1
            elif op == 'S':  # Soft clipping (not aligned, still present in the transcript)
                # Soft clipping affects only the transcript position
                tx_pos += length
            elif op == 'H':  # Hard clipping (not aligned and not present in the transcript)
                # Hard clipping affects neither reference nor transcript positions
                continue

        return ref_to_tx_map, tx_to_ref_map

    def build(self):
        # start by building transcriptomes for reference and target
        ref_tome = Transcriptome()
        ref_tome.load_genome(self.genome)
        ref_tome.build_from_file(self.annotation)
        ref_tome.extract_introns()

        target_tome = Transcriptome()
        target_tome.load_genome(self.target)
        target_tome.build_from_file(self.exon_sam2gtf_fname)
        target_tome.extract_introns()
        
        # extract ref and target donor/acceptor pairs into dictionaries.
        # Validate that all are conserved between the two genomes

        guide_tome = Transcriptome()
        if self.guide:
            guide_tome.build_from_file(self.guide)
            
            # verify compatibility of guide transcripts and CDSs with the target annotation
            
        # iterate over target transcripts
        for target_tx in target_tome:
            target_read_tid = target_tx.get_attr("read_name")
            target_tx.data = {"seq":"", "cds":"", "ref2tx_map":None, "tx2ref_map":None}
            print(f"Target transcript {target_read_tid}")

            # pull the corresponding transcript from reference
            ref_tx = ref_tome.get_by_tid(target_read_tid)
            ref_tx.data = {"seq":"", "cds":""}
            
            # pull the corresponding transcript from guide
            guide_tx = None
            if self.guide:
                guide_tx = guide_tome.get_by_tid(target_read_tid)
                guide_tx.data = {"seq":"", "cds":""}

            target_tx.data["seq"] = target_tx.get_sequence(target_tome.genome)
            ref_tx.data["seq"] = ref_tx.get_sequence(ref_tome.genome)
            
            target_tx.data["ref2tx_map"], target_tx.data["tx2ref_map"] = self.process_cigar(target_tx.get_attr("cigar"), ref_tx.start, target_tx.start)
            
            # check all donor and acceptor sites noting whether they are conserved or not
            
            continue
            
        # iterate over reference transcripts and check if any were not annotated in the target
        for ref_tx in ref_tome:
            if ref_tx.tid not in target_tome.tid_map:
                print(f"Reference transcript {ref_tx.tid} not annotated in the target genome")
        
def main():
    parser = argparse.ArgumentParser(description="Tool for HIV-1 genome annotation")

    parser.add_argument('-a', '--annotation', required=True, type=str, help='Path to the reference GTF/GFF annotation file')
    parser.add_argument('-g', '--genome', required=True, type=str, help='Path to the reference genome FASTA file')
    parser.add_argument('-t', '--target', required=True, type=str, help='Path to the target genome FASTA file')
    parser.add_argument('-q', '--guide', type=str, help='Optional path to the guide annotation file for the target genome. Transcripts and CDS from the guide will be used to validate the annotation')
    parser.add_argument('-o', '--output', type=str, help='Path to the output GTF file')

    parser.add_argument('--gffread', type=str, default='gffread', help='Path to the gffread executable')
    parser.add_argument('--minimap2', type=str, default='minimap2', help='Path to the minimap2 executable')
    parser.add_argument('--sam2gtf', type=str, default='sam2gtf', help='Path to the sam2gtf executable')
    parser.add_argument('--miniprot', type=str, default='miniprot', help='Path to the miniprot executable')

    parser.add_argument('--keep-tmp', action='store_true', help='Keep temporary files')
    parser.add_argument('--tmp-dir', type=str, default='./tmp', help='Directory to store temporary files')

    args = parser.parse_args()

    vira = Vira(args)
    vira.run()
    vira.write_output()

if __name__ == "__main__":
    main()
