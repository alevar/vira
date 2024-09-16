import argparse
import subprocess
import os

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

    def build(self):
        # start by building transcriptomes for reference and target
        ref_tome = Transcriptome()
        ref_tome.load_genome(self.genome)
        ref_tome.build_from_file(self.annotation)

        target_tome = Transcriptome()
        target_tome.load_genome(self.target)
        target_tome.build_from_file(self.exon_sam2gtf_fname)

        guide_tome = Transcriptome()
        if self.guide:
            guide_tome.build_from_file(self.guide)
            
        # iterate over target transcripts
        for target_tx in target_tome:
            target_read_tid = target_tx.get_attr("read_name")
            print(f"Target transcript {target_read_tid}")

            # pull the corresponding transcript from reference
            ref_tx = ref_tome.get_by_tid(target_read_tid)
            
            # pull the corresponding transcript from guide
            guide_tx = None
            if self.guide:
                guide_tx = guide_tome.get_by_tid(target_read_tid)

            target_tx_seq = target_tx.get_sequence(target_tome.genome)
            ref_tx_seq = ref_tx.get_sequence(ref_tome.genome)
            
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
