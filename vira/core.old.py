import os
import re
import sys
import copy
import pysam
import shutil
import argparse
import subprocess

from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

from .utils.common import *

from .classes.txgroup import Transcriptome, Gene, Bundle
from .classes.transcript import Transcript, Object
from .classes.splicegraph import SpliceGraph

class Vira:
    def __init__(self, args):
        
        # OPTIONS
        self.force_cds = args.force_cds
        self.keep_tmp = args.keep_tmp
        self.tmp_dir = standard_path(args.tmp_dir)
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        # TOOLS
        self.minimap2 = args.minimap2
        self.sam2gtf = args.sam2gtf
        self.miniprot = args.miniprot
        self.gffread = args.gffread
        self.snapper = args.snapper
        self.check_tools()
        
        
        if gtf_or_gff(args.annotation) is None:
            raise ValueError(f"{args.annotation} is not a valid GTF/GFF file.")
        
        if not os.path.exists(args.genome):
            raise FileNotFoundError(f"Genome file {args.genome} not found.")
        
        if not os.path.exists(args.target):
            raise FileNotFoundError(f"Input file {args.target} not found.")
        
        if args.guide and not os.path.exists(args.guide):
            raise FileNotFoundError(f"Guide annotation file {args.guide} not found.")           
        
        
        # INPUT FILES
        # create copies of files in tmp directory for use in the pipeline
        self.annotation = self.tmp_dir+"reference.gtf"
        self.genome = self.tmp_dir+"reference.fasta"
        self.target = self.tmp_dir+"target.fasta"
        shutil.copyfile(args.annotation, self.annotation)
        shutil.copyfile(args.genome, self.genome)
        shutil.copyfile(args.target, self.target)
        self.output = args.output
        self.gtf = gtf_or_gff(args.annotation)
        self.guide = None # guide annotation of the target genome - used to verify the inferred target annotation
        if args.guide:
            self.guide = self.tmp_dir+"target.guide."+args.guide.rsplit(".",1)[-1]
            shutil.copyfile(args.guide, self.guide)
        
        # TMP FILES
        self.dedup_reference_gtf_fname = self.tmp_dir+"dedup_reference.gtf"
        self.dedup_reference_cds_id_map = {}
        self.cds_nt_fasta_fname = self.tmp_dir+"cds_nt.fasta"
        self.cds_aa_fasta_fname = self.tmp_dir+"cds_aa.fasta"
        self.exon_nt_fasta_fname = self.tmp_dir+"exon_nt.fasta"
        self.cds_sam_fname = self.tmp_dir+"cds_nt.sam"
        self.exon_sam_pass1_fname = self.tmp_dir+"exon_nt.pass1.sam"
        self.exon_sam_pass2_fname = self.tmp_dir+"exon_nt.pass2.sam"
        self.exon_sam_snapper_fname = self.tmp_dir+"exon_nt.snapper.sam"
        self.exon_sam_fname = self.tmp_dir+"exon_nt.sam"
        self.exon_sam2gtf_pass1_fname = self.tmp_dir+"exon_nt.pass1.sam2gtf.gtf"
        self.exon_sam2gtf_fname = self.tmp_dir+"exon_nt.sam2gtf.gtf"
        self.pass1_junc_bed_fname = self.tmp_dir+"exon_nt.pass1.sam2gtf.junc.bed"
        self.cds_gtf_fname = self.tmp_dir+"cds.miniprot.gtf"
        self.guide_junc_bed_fname = self.tmp_dir+"guide.junc.bed"
        
        # Alignment
        self.aligner = Align.PairwiseAligner()
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        self.aligner.substitution_matrix = extend_matrix_alphabet(
            self.aligner.substitution_matrix,
            codes='BXZJUO-.',
        )

    def check_tools(self):
        if subprocess.call(f"command -v {self.minimap2}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"minimap2 is not installed or not available in PATH. Please install minimap2 before running vira. Installation instructions can be found at: https://github.com/lh3/minimap2")
        
        if subprocess.call(f"command -v {self.gffread}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"gffread is not installed or not available in PATH. Please install gffread before running vira. Installation instructions can be found at: https://github.com/gpertea/gffread")
        
        if subprocess.call(f"command -v {self.sam2gtf}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"sam2gtf is not installed or not available in PATH. Please install sam2gtf before running vira. Installation instructions can be found at: https://github.com/alevar/sam2gtf")
        
        if subprocess.call(f"command -v {self.miniprot}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"miniprot is not installed or not available in PATH. Please install miniprot before running vira. Installation instructions can be found at: https://github.com/lh3/miniprot")
        
        if subprocess.call(f"command -v {self.snapper}", shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
            raise EnvironmentError(f"snapper is not installed or not available in PATH. Please install snapper before running vira. Installation instructions can be found at: https://github.com/alevar/snapper")
        
        
        # verify version of miniprot as well to be above 0.13-r248
        version = subprocess.check_output(f"{self.miniprot} --version", shell=True).decode("utf-8").strip()
        print(version)
        v0 = int(version.split(".")[0])
        v1 = int(version.split(".")[1].split("-")[0])
        v2 = int(version.split("-")[1].split("r")[1])
        if v0 < 0 or (v0 == 0 and v1 < 13) or (v0 == 0 and v1 == 13 and v2 < 248):
            raise EnvironmentError(f"miniprot version {version} is not supported. Please install miniprot version 0.13-r248 or later. The miniprot version used is required to have the --spsc=<fname> option for specifying a list of trusted junctions")

    def extract_gene_protein_gtf(self, in_fname, genome, out_fname):
        # takes a gtf file, and extract one transcript per gene_id with its protein annotated
        # asserts there is only one unique protein per gene_id
        
        tome = Transcriptome()
        tome.load_genome(genome)
        tome.build_from_file(in_fname)
        
        cds_map = {}
        cds_id_map = {} # maps each tid to the id of the cds that was chosen
        
        for tx in tome:
            tx.data = {"cds": ""}
            nt = tx.get_sequence(tome.genome,use_cds=True)
            tx.data["cds"] = translate(nt)
            
            tid = tx.attrs["transcript_id"]
            cds_map.setdefault(tx.data["cds"],tid)
            cds_id_map[tid] = cds_map[tx.data["cds"]]

        # write out the output
        with open(out_fname,"w+") as outFP:
            for cds, tid in cds_map.items():
                tx = tome.get_by_tid(tid)
                outFP.write(tx.to_gtf()+"\n")
                
        return cds_id_map
    
    def run_commands(self):
        # extract junctions from the guide if available
        if self.guide is not None:
            cmd = ["paftools.js","gff2bed","-j",self.guide]
            print(" ".join(cmd)+" > "+self.guide_junc_bed_fname)
            with open(self.guide_junc_bed_fname,"w+") as outFP:
                subprocess.call(cmd,stdout=outFP)

        # extract the reference transcript
        cmd = [self.gffread,
                "-g",self.genome,
                "-w",self.exon_nt_fasta_fname,
                self.annotation]
        print(f"Extracting reference transcript sequences: {' '.join(cmd)}")
        subprocess.call(cmd)

        # run minimap of transcript sequences to the target genome
        cmd = [self.minimap2,"--for-only","-a",
               "-k9","-w3","--splice","--splice-flank=no","-g2k","-G9k","-A1","-B2","-O2,32","-E1,0","-b0","-C4","-z200","-un","--cap-sw-mem=0"]
        if self.guide is not None:
            cmd.extend(["--junc-bed",self.guide_junc_bed_fname,
                        "--junc-bonus","100"]) # set high bonus for guide junctions
        cmd.extend([self.target,self.exon_nt_fasta_fname])
        print(" ".join(cmd)+" > "+self.exon_sam_pass1_fname)
        with open(self.exon_sam_pass1_fname,"w+") as outFP:
            subprocess.call(cmd,stdout=outFP)

        # run sam2gtf
        cmd = [self.sam2gtf,
            "-i",self.exon_sam_pass1_fname,
            "-o",self.exon_sam2gtf_pass1_fname,
            "-p","50"]
        print(" ".join(cmd))
        subprocess.call(cmd)

        # extract consensus junctions from the pass1 alignment
        self.extract_consensus_sjs(self.annotation, self.genome, self.exon_sam2gtf_pass1_fname, self.target, self.pass1_junc_bed_fname)
        
        # do 2nd pass using consensus junctions this time
        cmd = [self.minimap2,"--for-only","-a",
               "-k9","-w3","--splice","--splice-flank=no","-g2k","-G9k","-A1","-B2","-O2,32","-E1,0","-b0","-C4","-z200","-un","--cap-sw-mem=0"]
        cmd.extend(["--junc-bed",self.pass1_junc_bed_fname,
                    "--junc-bonus","100"])
        cmd.extend([self.target,self.exon_nt_fasta_fname])
        print(" ".join(cmd)+" > "+self.exon_sam_pass2_fname)
        with open(self.exon_sam_pass2_fname,"w+") as outFP:
            subprocess.call(cmd,stdout=outFP)
        
        # run snapper to align introns
        cmd = [self.snapper,"--reference", self.annotation, "--sam", self.exon_sam_pass2_fname, "--output", self.exon_sam_snapper_fname, "--qry_intron_match_score", "100"]
        print(" ".join(cmd))
        subprocess.call(cmd)


        # snapper only performs correction such that the donor/acceptor sites match with the reference
        # however, if the intron length was a bit off (due to an extra D or I operation in CIGAR) - snapper can't do anything about it
        # we can instead perform majority vote for a consensus donor and acceptor site (for each reference donor/acceptor site since we know the mapping now from snapper)
        # and then force the consensus positions for each

        # now that the alignment matches reference donor/acceptor positions - time to build ocnsensus intron length
        donor_map, acceptor_map = self.extract_intron_map(self.annotation, self.exon_sam_snapper_fname)
        # build consensus acceptor maps
        donor_consensus_map = {}
        acceptor_consensus_map = {}
        for dp, data in donor_map.items():
            # set the donor position to the position with the highest count or first in case of a tie
            max_freq_pos = max(data.items(), key=lambda x: x[1])[0]
            for pos, freq in data.items():
                donor_consensus_map[pos] = max_freq_pos
        for ap, data in acceptor_map.items():
            # set the acceptor position to the position with the highest count or first in case of a tie
            max_freq_pos = max(data.items(), key=lambda x: x[1])[0]
            for pos, freq in data.items():
                acceptor_consensus_map[pos] = max_freq_pos
            
        # now we need to apply these consensus positions to the alignment
        self.apply_consensus_introns(donor_consensus_map, acceptor_consensus_map, self.exon_sam_snapper_fname, self.exon_sam_fname)
        
        # run sam2gtf
        cmd = [self.sam2gtf,
            "-i",self.exon_sam_fname,
            "-o",self.exon_sam2gtf_fname,
            "-p","50"]
        print(" ".join(cmd))
        subprocess.call(cmd)

        # do the proteins
        # begin by extracting deduplicated CDSs from the target genome
        # and building a map of the cds duplicate tids to the gene id
        # make sure there is a single CDS per gene_id
        self.dedup_reference_cds_id_map = self.extract_gene_protein_gtf(self.annotation, self.genome, self.dedup_reference_gtf_fname)
        cmd = [self.gffread,
                "-g",self.genome,
                "-y",self.cds_aa_fasta_fname,
                "-x",self.cds_nt_fasta_fname,
                self.dedup_reference_gtf_fname]
        print(f"Extracting reference protein sequences: {' '.join(cmd)}")
        subprocess.call(cmd)

        # get transcript_id to gene_id mapping
        tid2gid = {}
        with open(self.dedup_reference_gtf_fname,"r") as inFP:
            for line in inFP:
                if line.startswith("#"):
                    continue
                lcs = line.strip().split("\t")
                if lcs[2] == "transcript":
                    attrs = extract_attributes(lcs[8])
                    tid = attrs["transcript_id"]
                    gid = attrs["gene_id"]
                    tid2gid[tid] = gid
        
        miniprot_gff_fname = self.tmp_dir+"miniprot.gff"
        cmd = [self.miniprot,
               "--gff",
                "--spsc="+self.pass1_junc_bed_fname,
               self.target, self.cds_aa_fasta_fname]
        print(" ".join(cmd)+" > "+miniprot_gff_fname)
        with open(miniprot_gff_fname,"w+") as outFP:
            subprocess.call(cmd,stdout=outFP)

        # need to standardize the miniprot output
        tome = Transcriptome()
        tome.build_from_file(miniprot_gff_fname)
        # use comments to extract the PAF alignment notes and append to the records
        with open(miniprot_gff_fname,"r") as inFP:
            cur_cigar = None
            cur_tid = None
            for line in inFP:
                if line.startswith("##PAF"):
                    cur_tid = line.split("\t")[1]
                    cur_cigar = line.split("cg:Z:",1)[1].split("\t",1)[0]
                else:
                    if cur_cigar is not None:
                        new_tid = line.split("\t")[8].split("ID=",1)[1].split(";",1)[0]
                        tx = tome.get_by_tid(new_tid)
                        tx.set_tid(cur_tid)
                        tx.add_attribute("cigar",cur_cigar)
                        tx.set_gid(tid2gid[cur_tid])

                        for e in tx.exons:
                            e[2].set_tid(cur_tid)
                        for c in tx.cds:
                            c[2].set_tid(cur_tid)
                        cur_cigar = None
                        cur_tid = None

        # write out the standardized file
        with open(self.cds_gtf_fname,"w+") as outFP:
            outFP.write(tome.to_gtf())

    def run(self):
        try:
            self.run_commands()
        except Exception as e:
            sys.stderr.write(f"Error running the pipeline: {str(e)}\n")
            sys.exit(1)
        # combine annotated transcripts, CDSs and guide annotation together
        # for each transcript/cds annotate any differences
        try:
            self.build()
        except Exception as e:
            sys.stderr.write(f"Error building the final output: {str(e)}\n")
            sys.exit(1)

    def extract_consensus_sjs(self, ref_gtf_fname, ref_fasta_fname, trg_gtf_fname, trg_fasta_fname, out_gtf_fname) -> None:
        # given a gtf file  produced by the sam2gtf tool
        # extracts a mapping of query junctions to target junctions
        # for each query junction, computes what the consensus position is
        # output the result in the format compatible with miniprot
        # expected format: ctg  offset  +|-  D|A  score

        # start by building transcriptomes for reference and target
        ref_tome = Transcriptome()
        ref_tome.load_genome(ref_fasta_fname)
        ref_tome.build_from_file(ref_gtf_fname)
        ref_tome.extract_introns()

        target_tome = Transcriptome()
        target_tome.load_genome(trg_fasta_fname)
        target_tome.build_from_file(trg_gtf_fname)
        target_tome.extract_introns()
        # deduplicate target transcripts and convert transcript_ids
        self.reassign_tids(target_tome)

        # iterate over target transcripts
        for target_tx in target_tome:
            target_tx.data = dict()
            target_tx.data["ref2trg_map"] = None
            target_tx.data["trg2ref_map"] = None

            # pull the corresponding transcript from reference
            ref_tx = ref_tome.get_by_tid(target_tx.get_tid())

            # assign gene_id based on the reference along with other attributes
            target_tx.set_gid(ref_tx.get_attr("gene_id"))
            for e in target_tx.get_exons():
                e[2].set_gid(ref_tx.get_attr("gene_id"))
            for c in target_tx.get_cds():
                c[2].set_gid(ref_tx.get_attr("gene_id"))

            target_tx.data["ref2trg_map"], target_tx.data["trg2ref_map"] = self.process_cigar(target_tx.get_attr("cigar"), ref_tx, target_tx.get_start())

        # extract junction mapping
        donor_map = {} # holds the mapping between reference and target donor sites
        acceptor_map = {} # holds the mapping between reference and target acceptor sites
        for ref_tx in ref_tome.transcript_it():
            for ref_i,ref_intron in enumerate(ref_tx.introns_it()):
                # find position of the intron in the target genome
                target_tx = target_tome.get_by_tid(ref_tx.get_tid())
                if target_tx is None:
                    continue

                if ref_intron[0]-1 not in target_tx.data["ref2trg_map"] or ref_intron[1] not in target_tx.data["ref2trg_map"]:
                    continue
                trg_donor_pos = target_tx.data["ref2trg_map"][ref_intron[0]-1]
                trg_acceptor_pos = target_tx.data["ref2trg_map"][ref_intron[1]]

                if trg_donor_pos is None or trg_acceptor_pos is None:
                    continue
                if trg_donor_pos[1] != "M" or trg_acceptor_pos[1] != "M":
                    continue

                donor_map.setdefault(ref_intron[0],[]).append(trg_donor_pos)
                acceptor_map.setdefault(ref_intron[1],[]).append(trg_acceptor_pos)

        # verify consistency
        for donor_pos in donor_map:
            # assign the most common mapping as the target site
            donor_map[donor_pos] = max(set(donor_map[donor_pos]), key=donor_map[donor_pos].count)[0]
        for acceptor_pos in acceptor_map:
            acceptor_map[acceptor_pos] = max(set(acceptor_map[acceptor_pos]), key=acceptor_map[acceptor_pos].count)[0]

        # write out the results
        with open(out_gtf_fname,"w+") as outFP:
            for donor_pos in donor_map:
                outFP.write(f"{target_tx.get_seqid()}\t{donor_map[donor_pos]}\t+\tD\t100\n")
            for acceptor_pos in acceptor_map:
                outFP.write(f"{target_tx.get_seqid()}\t{acceptor_map[acceptor_pos]-1}\t+\tA\t100\n")
        return None

    def extract_intron_map(self, qry_tome_fname: str, trg_aln_fname: str):
        # extract the mapping of introns from reference to target
        # for each donor and acceptor site in the reference - find the corresponding donor and acceptor site in the target
        donor_map = {} # holds the mapping between reference and target donor sites
        acceptor_map = {} # holds the mapping between reference and target acceptor sites

        # load the reference transcriptome
        qry_tome = Transcriptome()
        qry_tome.build_from_file(qry_tome_fname)

        # iterate over alignments
        for read in pysam.AlignmentFile(trg_aln_fname, "r"):
            if read.is_unmapped or not read.cigarstring or read.is_secondary:
                continue

            try:
                qry_tx = qry_tome.get_by_tid(read.query_name)

                # process the cigar to get the mapping
                trg_pos = read.reference_start
                ref2trg_map, _ = self.process_cigar(read.cigarstring, qry_tx, trg_pos)

                for it in qry_tx.introns_it():
                    donor_pos = it[0]-1
                    acceptor_pos = it[1]

                    if donor_pos not in ref2trg_map:
                        raise ValueError(f"Donor position {donor_pos} not found in ref2trg_map")
                    if acceptor_pos not in ref2trg_map:
                        raise ValueError(f"Acceptor position {acceptor_pos} not found in ref2trg_map")
                    
                    trg_donor_pos, _ = ref2trg_map[donor_pos]
                    trg_acceptor_pos, _ = ref2trg_map[acceptor_pos]

                    donor_map.setdefault(donor_pos,{}).setdefault(trg_donor_pos,0)
                    donor_map[donor_pos][trg_donor_pos] += 1
                    acceptor_map.setdefault(acceptor_pos,{}).setdefault(trg_acceptor_pos,0)
                    acceptor_map[acceptor_pos][trg_acceptor_pos] += 1

            except Exception as e:
                sys.stderr.write(f"Error processing {read.query_name}: {str(e)}\n")

        return donor_map, acceptor_map
    
    def process_cigar(self, cigar_string: str, qry_tx: Transcript, trg_start:int):
        """
        Process CIGAR string to create a mapping from query genome positions
        to target positions.

        :param cigar_string: CIGAR string from alignment
        :param qry_tx: query transcript
        :param trg_tx: target transcript
        :return: A dictionary mapping qry positions to target positions
        """
        # CIGAR operation regex
        cigar_operations = parse_cigar_into_tuples(cigar_string)
        
        # Initialize positions
        qry_pos = 0
        trg_pos = trg_start

        # Maps to track reference to query and query to reference
        qry_to_trg_map = {}
        trg_to_qry_map = {}

        # Process each CIGAR operation
        for length, op in cigar_operations:

            if op == 'M' or op == '=' or op == 'X':  # Match/Mismatch
                for _ in range(length):
                    qry_genome_pos = qry_tx.genome_coordinate(qry_pos)
                    # Map query to target and target to query
                    qry_to_trg_map[qry_genome_pos] = (trg_pos,op)
                    trg_to_qry_map[trg_pos] = (qry_genome_pos,op)
                    qry_pos += 1
                    trg_pos += 1
            elif op == 'I':  # Insertion in query (relative to the target)
                # Insertions affect only the query position
                for _ in range(length):
                    qry_genome_pos = qry_tx.genome_coordinate(qry_pos)
                    qry_to_trg_map[qry_genome_pos] = (None,op)  # No target mapping for inserted bases
                    qry_pos += 1
            elif op == 'D':
                prev_qry_pos = qry_pos - 1
                prev_qry_genome_pos = qry_tx.genome_coordinate(prev_qry_pos)
                if prev_qry_genome_pos not in qry_to_trg_map:
                    raise ValueError(f"Query position {prev_qry_genome_pos} not found in qry_to_trg_map")
                if qry_to_trg_map[prev_qry_genome_pos][0] is None:
                    raise ValueError(f"Previous target position at {prev_qry_genome_pos} is None")
                
                for _ in range(length):
                    qry_to_trg_map[prev_qry_genome_pos] = [trg_pos,op]  # Extend the previous match to cover the deletion
                    trg_pos += 1
            elif op == 'N':
                # we need to check here if it matches any of the query sites
                # this way we could convert novel introns to matches
                # alternatively, we could also do this 
                for _ in range(length):
                    trg_pos += 1
            elif op == 'S':  # Soft clipping (not aligned, still present in the query)
                # Soft clipping affects only the query position
                for _ in range(length):
                    qry_genome_pos = qry_tx.genome_coordinate(qry_pos)
                    qry_to_trg_map[qry_genome_pos] = (trg_pos,op)
                    qry_pos += 1
            elif op == 'H':  # Hard clipping (not aligned and not present in the target)
                # Hard clipping affects neither query nor target positions
                continue

        return qry_to_trg_map, trg_to_qry_map
    
    def adjust_cigar_introns(self, cigarstring, aln_start, donor_consensus_map, acceptor_consensus_map):
        # # adjust the positions of the introns in the cigar string based on the provided maps
        # map contains a dict of positions mapping them to consensus positions
        # for every donor/acceptor site in the alignment - the method checks the provided maps
        # and applies the consensus position if positions are matching
        # the alignment is modified by updating cigar with I or D operations at the earliest available site before donor or after acceptor
    
        ops = parse_cigar_into_tuples(cigarstring)

        pos = aln_start

        cigar_idx = 0
        while cigar_idx < len(ops):
            oplen, op = ops[cigar_idx]
            if op == 'N':
                # get donor and acceptor positions
                cur_donor_pos = pos - 1
                pos += oplen
                cur_acceptor_pos = pos

                # if we adjust donor - does that mean the donor_consensus map is invalidated?
                new_donor_pos = donor_consensus_map.get(cur_donor_pos, cur_donor_pos)
                new_acceptor_pos = acceptor_consensus_map.get(cur_acceptor_pos, cur_acceptor_pos)

                if new_donor_pos != cur_donor_pos:
                    sub_ops = ops[:cigar_idx]
                    if new_donor_pos < cur_donor_pos:
                        # insert an insertion at the earliest available site after acceptor
                        shorten_cigar_inplace(sub_ops, cur_donor_pos-new_donor_pos, from_end=True, offset=2)
                    else:
                        # insert a deletion at the earliest available site before donor
                        elongate_cigar_inplace(sub_ops, new_donor_pos-cur_donor_pos, from_end=True, offset=2)
                    ops = sub_ops + ops[cigar_idx:]
                    # account for possibly modified cigar length in the cigar_index
                    cigar_idx = len(sub_ops)

                # replace current RefSkip with updated length
                new_intron_length = (new_acceptor_pos - new_donor_pos)-1
                ops[cigar_idx] = (new_intron_length,'N')

                if new_acceptor_pos != cur_acceptor_pos:
                    sub_ops = ops[cigar_idx+1:]
                    if new_acceptor_pos > cur_acceptor_pos:
                        # insert a deletion at the earliest available site after acceptor
                        shorten_cigar_inplace(sub_ops, new_acceptor_pos-cur_acceptor_pos, from_end=False, offset=2)
                    else:
                        # insert an insertion at the earliest available site before donor
                        elongate_cigar_inplace(sub_ops, cur_acceptor_pos-new_acceptor_pos, from_end=False, offset=2)
                    ops = ops[:cigar_idx+1] + sub_ops

            else:
                if op in ['M', '=', 'X']:
                    pos += oplen
                elif op == 'I':
                    pass
                elif op == 'D':
                    pos += oplen
                elif op == 'S':
                    pass
                elif op == 'H':
                    pass
                else:
                    raise ValueError(f"Invalid CIGAR operation {op}")
                
            cigar_idx += 1

        return build_cigar_from_tuples(ops)

    def apply_consensus_introns(self, donor_consensus_map, acceptor_consensus_map, in_sam_fname, out_sam_fname):
        # adjust the positions of the introns in the alignment based on the provided maps
        # map contains a dict of positions mapping them to consensus positions
        # for every donor/acceptor site in the alignment - the method checks the provided maps
        # and applies the consensus position if positions are matching
        # the alignment is modified by updating cigar with I or D operations at the earliest available site before donor or after acceptor
        input_mode = 'rb' if in_sam_fname.endswith('.bam') else 'r'
        output_mode = 'wb' if out_sam_fname.endswith('.bam') else 'wh'

        with pysam.AlignmentFile(in_sam_fname, input_mode) as infile, \
            pysam.AlignmentFile(out_sam_fname, output_mode, template=infile) as outfile:

            for read in infile:
                if read.is_unmapped or not read.cigarstring:
                    outfile.write(read)
                    continue

                try:
                    new_cigar = self.adjust_cigar_introns(read.cigarstring, read.reference_start, donor_consensus_map, acceptor_consensus_map)
                    
                    # Create modified read
                    modified_read = pysam.AlignedSegment(outfile.header)
                    modified_read.set_tags(read.get_tags())
                    modified_read.query_name = read.query_name
                    modified_read.query_sequence = read.query_sequence
                    modified_read.flag = read.flag
                    modified_read.reference_id = read.reference_id
                    modified_read.reference_start = read.reference_start
                    modified_read.mapping_quality = read.mapping_quality
                    modified_read.cigarstring = new_cigar
                    modified_read.query_qualities = read.query_qualities
                    modified_read.next_reference_id = read.next_reference_id
                    modified_read.next_reference_start = read.next_reference_start
                    modified_read.template_length = read.template_length

                    outfile.write(modified_read)

                except Exception as e:
                    sys.stderr.write(f"Error processing {read.query_name}: {str(e)}\n")
                    outfile.write(read)

    def reassign_tids(self, tome: Transcriptome, attr: str = "read_name"):
        # assigns the specified attribute as the transcript id
        # checks there are no duplicates
        assigned = set()
        for tx in tome:
            cur_tid = tx.get_tid()
            tid = tx.get_attr(attr)
            if tid in assigned:
                raise ValueError(f"Duplicate transcript id {tid} found in the annotation")
            tx.set_tid(tid)
            for e in tx.get_exons():
                e[2].set_tid(tid)
            for c in tx.get_cds():
                c[2].set_tid(tid)

            # change mapping in tome
            tome.tid_map[tid] = tome.tid_map.pop(cur_tid)

    def extract_junction_seq(self, tx: Transcript, genome):
        # for each transcript extract donor and acceptor sites for each intron
        sjs = []
        if len(tx.exons) == 1:
            return sjs
        for i,e in enumerate(tx.get_exons()):
            if i != 0:
                # skip acceptor extraction for the first exon
                acceptor_seq = genome[e[2].seqid][e[2].start-1-2:e[2].start-1].seq
                sjs[-1][1] = acceptor_seq
                e[2].add_attribute("acceptor_seq",acceptor_seq,replace=True)
            if i != len(tx.exons)-1:
                # skip donor extraction for the last exon
                donor_seq = genome[e[2].seqid][e[2].end-1:e[2].end-1+2].seq
                sjs.append([donor_seq,None])
                e[2].add_attribute("donor_seq",donor_seq,replace=True)
        return sjs

    def compare_sj_seq(self, ref_sj_seq: str, target_sj_seq: str):
        # compare donor and acceptor sites
        sj_comp = []
        for i in range(len(ref_sj_seq)):
            ref_donor, ref_acceptor = ref_sj_seq[i]
            target_donor, target_acceptor = target_sj_seq[i]
            if ref_donor == target_donor:
                sj_comp.append("D")
            else:
                sj_comp.append("d")
            if ref_acceptor == target_acceptor:
                sj_comp.append("A")
            else:
                sj_comp.append("a")
        return sj_comp
    
    def get_first_cds(self, tx: Transcript, tome: Transcriptome):
        # for a given transcript - identify the first available ORF and produce a list of CDS for that transcript
        cds = []
        if tx.data == None or tx.data["seq"] == "":
            tx.data["seq"] = tx.get_sequence(tome.genome)
        orf = find_first_orf(tx.data["seq"])
        if len(orf) == 0:
            return cds
        ostart, oend = orf
        # translate to genomic coordinates
        ostart = tx.genome_coordinate(ostart)
        oend = tx.genome_coordinate(oend)
        tx_chain = tx.get_chain()
        cds_chain = cut_chain(tx_chain, ostart, oend)
        for c in cds_chain:
            obj = Object()
            obj.set_attributes({"transcript_id":tx.get_tid()})
            obj.set_start(c[0])
            obj.set_end(c[1])
            obj.set_seqid(tx.get_seqid())
            obj.set_strand(tx.get_strand())
            obj.obj_type = Types.CDS
            cds.append(obj)
        
        return cds
    
    def build_target2guide_map(self, guide_tome: Transcriptome, target_tome: Transcriptome, ref_tome: Transcriptome):
        target2guide_map = {}
        
        guide_cds_map = {}
        for tx in guide_tome:
            if tx.has_cds():
                aa = tx.data["cds"]
                guide_cds_map.setdefault(aa,tx.get_tid())
        
        # load a map of all transcripts for each cds chain in the reference
        ref_cds_map = {}
        for tx in ref_tome:
            if not tx.get_tid() in target_tome: # make sure the reference transcripts we are including are only those that were mapped over
                continue
            if tx.has_cds():
                aa = tx.data["cds"]
                ref_cds_map.setdefault(aa,[]).append(tx.get_tid())
                
        # for each reference protein - find the corresponding guide protein
        for aa, tids in ref_cds_map.items():
            # find matching guide protein by aligning against all guide proteins
            alignment, identity, guide_tid = find_best_alignment(self.aligner, aa, guide_cds_map)
            for tid in tids:
                target2guide_map[tid] = guide_tid

        return target2guide_map

    def build(self):
        # start by building transcriptomes for reference and target
        ref_tome = Transcriptome()
        ref_tome.load_genome(self.genome)
        ref_tome.build_from_file(self.annotation)
        ref_tome.extract_introns()
        for tx in ref_tome:
            tx.data = {"seq":"", "cds":""}
            tx.data["seq"] = tx.get_sequence(ref_tome.genome)
            nt = tx.get_sequence(ref_tome.genome,use_cds=True)
            tx.data["cds"] = translate(nt)

        target_tome = Transcriptome()
        target_tome.load_genome(self.target)
        target_tome.build_from_file(self.exon_sam2gtf_fname)
        target_tome.extract_introns()
        # deduplicate target transcripts and convert transcript_ids
        try:
            self.reassign_tids(target_tome)
        except ValueError as e:
            sys.stderr.write(f"Error: {str(e)}\n")
            return
        for tx in target_tome:
            tx.data = {"seq":"", "cds":""}
            tx.data["seq"] = tx.get_sequence(target_tome.genome)

        # load the cds results
        target_cds_tome = Transcriptome()
        target_cds_tome.load_genome(self.target)
        target_cds_tome.build_from_file(self.cds_gtf_fname)
        target_cds_tome.extract_introns()
        for tx in target_cds_tome:
            tx.data = {"seq":"", "cds":""}
            tx.data["seq"] = tx.get_sequence(target_cds_tome.genome)
            nt = tx.get_sequence(target_cds_tome.genome,use_cds=True)
            tx.data["cds"] = translate(nt)

        guide_tome = Transcriptome()
        target2guide_map = {}
        if self.guide is not None:
            guide_tome.load_genome(self.target)
            guide_tome.build_from_file(self.guide)
            guide_tome.extract_introns()
            # extract cds sequnces from the guide
            for tx in guide_tome:
                tx.data = {"cds": ""}
                nt = tx.get_sequence(guide_tome.genome,use_cds=True)
                tx.data["cds"] = translate(nt)
                tx.merge_cds("longest")

            target2guide_map = self.build_target2guide_map(guide_tome, target_tome, ref_tome)

        # iterate over reference transcripts and report any that were not annotated in the target
        for ref_tx in ref_tome:
            if ref_tx.tid not in target_tome.tid_map:
                print(f"Reference transcript {ref_tx.tid} not annotated in the target genome")
            
        # iterate over target transcripts
        for target_tx in target_tome:
            target_tx.data["ref2trg_map"] = None
            target_tx.data["trg2ref_map"] = None

            # pull the corresponding transcript from reference
            ref_tx = ref_tome.get_by_tid(target_tx.get_tid())

            # assign gene_id based on the reference along with other attributes
            target_tx.set_gid(ref_tx.get_attr("gene_id"))
            for e in target_tx.get_exons():
                e[2].set_gid(ref_tx.get_attr("gene_id"))
            for c in target_tx.get_cds():
                c[2].set_gid(ref_tx.get_attr("gene_id"))
            
            target_tx.data["ref2trg_map"], target_tx.data["trg2ref_map"] = self.process_cigar(target_tx.get_attr("cigar"), ref_tx, target_tx.get_start())
            
            self.fix_with_guide(target_tx, ref_tx, guide_tome)

            # # check all donor and acceptor sites noting whether they are conserved or not
            # ref_sj_seq = self.extract_junction_seq(ref_tx, ref_tome.genome)
            # target_sj_seq = self.extract_junction_seq(target_tx, target_tome.genome)
            # # compare donor acceptor pairs
            # # sj_comp = self.compare_sj_seq(ref_sj_seq, target_sj_seq)
            
        # check all donor and acceptor positions noting whether they are conserved or not
        donor_map, acceptor_map = self.extract_intron_map(self.annotation, self.exon_sam_fname)
        for donor_pos in donor_map:
            if len(set(donor_map[donor_pos])) > 1:
                raise ValueError(f"Multiple target donor sites found for reference donor site {donor_pos}: {donor_map[donor_pos]}")
        for acceptor_pos in acceptor_map:
            if len(set(acceptor_map[acceptor_pos])) > 1:
                raise ValueError(f"Multiple target acceptor sites found for reference acceptor site {acceptor_pos}: {acceptor_map[acceptor_pos]}")
            
        cds_choices = {"miniprot":{}, "guide":{}}
        
        #========================================================================
        #===========================   MINIPROT   ===============================
        #========================================================================
        # load the CDS for each transcript
        for target_tx in target_tome:
            tid = target_tx.get_tid()
            # get the tid of the transcript whose cds was used in the deduplicated reference
            cds_tid = self.dedup_reference_cds_id_map[tid]
            if not cds_tid in target_cds_tome: # skipped if not mapped over
                continue
            target_cds_tx = target_cds_tome.get_by_tid(cds_tid)

            # check compatibility of the CDS with the transcript
            target_chain = target_tx.get_chain()
            target_cds_chain = target_cds_tx.get_chain(use_cds=True)
            if not target_cds_chain == cut_chain(target_chain, target_cds_chain[0][0], target_cds_chain[-1][1]):
                continue
            # add the CDS to the transcript
            tmp_tx = copy.deepcopy(target_tx)
            for c in target_cds_tx.get_cds():
                tmp = copy.deepcopy(c[2])
                tmp.add_attribute("transcript_id",tid,replace=True)
                tmp_tx.add_cds(tmp)
            # get translated sequence
            nt = tmp_tx.get_sequence(target_tome.genome,use_cds=True)
            tmp_tx.data["cds"] = translate(nt)
            cds_choices["miniprot"][tid] = tmp_tx
            
        #========================================================================
        #============================   GUIDE   =================================
        #========================================================================
        # load the guide annotation where available
        if self.guide is not None:
            for tid, guide_tid in target2guide_map.items():
                if guide_tid is None:
                    continue
                target_tx = target_tome.get_by_tid(tid)
                if target_tx is None:
                    raise ValueError(f"Transcript {tid} not found in the target genome")
                guide_tx = guide_tome.get_by_tid(guide_tid)
                
                # check compatibility of the CDS with the transcript
                target_chain = target_tx.get_chain()
                guide_cds_chain = guide_tx.get_chain(use_cds=True)
                if not guide_cds_chain == cut_chain(target_chain, guide_cds_chain[0][0], guide_cds_chain[-1][1]):
                    continue
                # add the CDS to the transcript
                tmp_tx = copy.deepcopy(target_tx)
                for c in guide_tx.get_cds():
                    tmp = copy.deepcopy(c[2])
                    tmp.add_attribute("transcript_id",tid,replace=True)
                    tmp_tx.add_cds(tmp)
                # get translated sequence
                nt = tmp_tx.get_sequence(target_tome.genome,use_cds=True)
                tmp_tx.data["cds"] = translate(nt)
                cds_choices["guide"][tid] = tmp_tx
                
        # compare the CDS choices ensuring consistency
        # for each transcript compare choices
        # also ensure all agree between transcripts of the same gene
        for tx in target_tome:
            if tx.get_tid() in cds_choices["miniprot"] and not tx.get_tid() in cds_choices["guide"]:
                tx.cds = cds_choices["miniprot"][tx.get_tid()].cds
                tx.add_attribute("cds_source","miniprot")
            elif tx.get_tid() in cds_choices["guide"] and not tx.get_tid() in cds_choices["miniprot"]:
                tx.cds = cds_choices["guide"][tx.get_tid()].cds
                tx.add_attribute("cds_source","guide")
            elif tx.get_tid() in cds_choices["guide"] and tx.get_tid() in cds_choices["miniprot"]:
                tx.cds = cds_choices["guide"][tx.get_tid()].cds
                tx.add_attribute("cds_source","guide")
            else:
                continue

        # write out the final GTF file
        with open(self.output,"w+") as outFP:
            outFP.write(target_tome.to_gtf())

    def fix_with_guide(self, tx: Transcript, ref_tx: Transcript, guide_tome: Transcriptome):
        # uses the ref2trg_map and trg2ref_map to find which positions on the transcript 
        # should be contiguous based on the guide exons
        # for every guide interval, maps over to the reference transcript,
        # identifies intervals on reference where both reference transcript and guide are contiguous
        # and returns their positions relative to the guide transcript

        target_sub_chain = []
        for gi in guide_tome.intervals():
            # get the minimum gi start which is in the trg2ref_map
            keys_in_interval = [key for key in tx.data["trg2ref_map"] if gi[0] <= key <= gi[1]]

            min_gi_start = None
            max_gi_end = None
            if keys_in_interval:
                min_gi_start = min(keys_in_interval)
                max_gi_end = max(keys_in_interval)
            else: # gi not found in reference
                continue
            
            ref_start = tx.data["trg2ref_map"][min_gi_start][0]
            ref_end = tx.data["trg2ref_map"][max_gi_end][0]
            
            # now take the reference chain and cut it to the ref_start and ref_end
            ref_sub_chain = cut_chain(ref_tx.get_chain(), ref_start, ref_end)
            # convert each interval back into the target space
            for c in ref_sub_chain:
                # get the minimum gi start which is in the trg2ref_map
                keys_in_interval = [key for key in tx.data["ref2trg_map"] if c[0] <= key <= c[1]]

                min_c_start = None
                max_c_end = None
                if keys_in_interval:
                    min_c_start = min(keys_in_interval)
                    max_c_end = max(keys_in_interval)
                else: # c not found in reference
                    continue
                target_sub_chain.append([tx.data["ref2trg_map"][min_c_start][0],tx.data["ref2trg_map"][max_c_end][0]])

        # each interval should now be contiguous
        # apply these intervals to the target transcript to create contiguous blocks
        dummy_exon = tx.get_exons()[0][2]
        for c in target_sub_chain:
            if c[1]-c[0] < 1:
                continue
            obj = Object()
            obj.set_seqid(dummy_exon.get_seqid())
            obj.set_strand(dummy_exon.get_strand())
            obj.set_attributes({"transcript_id":tx.get_tid()})
            obj.set_start(c[0])
            obj.set_end(c[1])
            tx.add_exon(obj)
        tx.merge_exons()
        
    def fix_with_local(self, tx: Transcript, ref_tx: Transcript):
        # use the ref2trg_map and trg2ref_map to identify a windown within which to perform unspliced alignment
        # for each individual exon. This way, hopefully we can refine the exon boundaries better
        for re in ref_tx.exons():
            # get positions in the target transcript
            keys_in_interval = [key for key in tx.data["ref2trg_map"] if re[0] <= key <= re[1]]
            
            min_re_start = None
            max_re_end = None
            
        return
        
        
def main():
    parser = argparse.ArgumentParser(description="Tool for HIV-1 genome annotation")

    parser.add_argument('-a', '--annotation', required=True, type=str, help='Path to the reference GTF/GFF annotation file')
    parser.add_argument('-g', '--genome', required=True, type=str, help='Path to the reference genome FASTA file')
    parser.add_argument('-t', '--target', required=True, type=str, help='Path to the target genome FASTA file')
    parser.add_argument('-q', '--guide', type=str, help='Optional path to the guide annotation file for the target genome. Transcripts and CDS from the guide will be used to validate the annotation')
    parser.add_argument('-o', '--output', type=str, help='Path to the output GTF file')
    
    parser.add_argument('--force-cds', action='store_true', help='Force the CDS from the guide onto the transcript chain, even if that means merging adjacent exons together (can fix alignment artifacts such as spurious introns). If the CDS does not fit the transcript chain, the transcript will be skipped')

    parser.add_argument('--gffread', type=str, default='gffread', help='Path to the gffread executable')
    parser.add_argument('--minimap2', type=str, default='minimap2', help='Path to the minimap2 executable')
    parser.add_argument('--sam2gtf', type=str, default='sam2gtf', help='Path to the sam2gtf executable')
    parser.add_argument('--miniprot', type=str, default='miniprot', help='Path to the miniprot executable. If not set - minimap2 will be used to align nucleotide sequence of the CDS instead')
    parser.add_argument('--snapper', type=str, default='snapper', help='Path to the snapper executable')

    parser.add_argument('--keep-tmp', action='store_true', help='Keep temporary files')
    parser.add_argument('--tmp-dir', type=str, default='./tmp', help='Directory to store temporary files')

    args = parser.parse_args()

    try:
        vira = Vira(args)
        vira.run()
    except Exception as e:
        sys.stderr.write(f"Error running the pipeline: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()

# TODO:
# 1. Attributes
#   - whether guide or miniprot used for CDS
#   - whether guide or miniprot chains did not fit transcript chain



# What we can do to rescue some of the genomes that are currently failing
# Once we get consensus donor and acceptor positions
# we can reconstruct transcripts from those positions only
# Then we can check if the transcript matches the CIGAR and if not - we flag it as such
# for every genome, we can compute a score of how good annotation is
# this could factor in several things such as 
# how many inconsisten donor/acceptor sites there are
# how many introns are missing or are extra
# how many donor/acceptor sites match AG/GT expected pairing