import os
import sys
import copy
import pysam
import shutil
import argparse
import subprocess

from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

from .utils.common import *

from .classes.txgroup import Transcriptome
from .classes.transcript import Transcript

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
        self.annotation = self.tmp_dir+"query.gtf"
        self.genome = self.tmp_dir+"query.fasta"
        self.target = self.tmp_dir+"target.fasta"
        shutil.copyfile(args.annotation, self.annotation)
        shutil.copyfile(args.genome, self.genome)
        shutil.copyfile(args.target, self.target)
        
        # load seqid from the query genome and target genome
        self.qry_seqid = list(SeqIO.index(self.genome,"fasta").keys())[0]
        self.trg_seqid = list(SeqIO.index(self.target,"fasta").keys())[0]
        self.qry_strand = "+"
        self.trg_strand = "+"
        
        self.qry_tome = Transcriptome()
        self.qry_tome.load_genome(self.genome)
        self.qry_tome.build_from_file(self.annotation)
        self.qry_tome.extract_introns()
        for tx in self.qry_tome:
            tx.data = {"seq":"", "cds":""}
            tx.data["seq"] = tx.get_sequence(self.qry_tome.genome)
            nt = tx.get_sequence(self.qry_tome.genome,use_cds=True)
            tx.data["cds"] = translate(nt)


        self.output = args.output
        self.gtf = gtf_or_gff(args.annotation)
        self.guide = None # guide annotation of the target genome - used to verify the inferred target annotation
        self.guide_tome = None
        if args.guide:
            self.guide = self.tmp_dir+"target.guide."+args.guide.rsplit(".",1)[-1]
            shutil.copyfile(args.guide, self.guide)
            self.guide_tome = Transcriptome()
            self.guide_tome.build_from_file(self.guide)
            self.guide_tome.load_genome(self.target)
            # extract cds sequnces from the guide
            for tx in self.guide_tome:
                tx.data = {"cds": ""}
                nt = tx.get_sequence(self.guide_tome.genome,use_cds=True)
                tx.data["cds"] = translate(nt)
                tx.merge_cds("longest")
        
        # TMP FILES
        self.dedup_qry_gtf_fname = self.tmp_dir+"dedup_query.gtf"
        self.dedup_qry_cds_id_map = {}
        self.cds_nt_fasta_fname = self.tmp_dir+"cds_nt.fasta"
        self.cds_aa_fasta_fname = self.tmp_dir+"cds_aa.fasta"
        self.exon_nt_fasta_fname = self.tmp_dir+"exon_nt.fasta"
        self.cds_sam_fname = self.tmp_dir+"cds_nt.sam"
        self.exon_sam_mm2_fname = self.tmp_dir+"exon_nt.mm2.sam" # initial alignment with minimap2
        self.exon_sam_snapper_fname = self.tmp_dir+"exon_nt.snapper.sam" # minimap2 alignment adjusted by snapper
        self.consensus_gtf = self.tmp_dir+"consensus.gtf"
        self.cds_gtf_fname = self.tmp_dir+"cds.miniprot.gtf"
        self.guide_junc_bed_fname = self.tmp_dir+"guide.junc.bed"
        self.junc_bed_for_miniprot_fname = self.tmp_dir+"junc.for_miniprot.bed"
        
        # stats files
        self.junction_stats_fname = args.output.rsplit(".",1)[0]+".junction_stats.tsv"
        with open(self.junction_stats_fname,"w+") as stats_outFP:
            # write header
            stats_outFP.write("position\tquery_position\ttype\tmap_consistency\tsequence\tquery_sequence\n")

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

    def process_cigar(self, cigar_string: str, qry_tx: Transcript, trg_start:int):
        """
        Process CIGAR string to create a mapping from query genome positions
        to target positions.

        :param cigar_string: CIGAR string from alignment
        :param qry_tx: query transcript
        :param trg_tx: target transcript
        :return: A dictionary mapping qry positions to target positions
        """
        cigar_operations = parse_cigar_into_tuples(cigar_string)
        
        qry_pos = 0
        trg_pos = trg_start

        # Maps to track query to target and target to query
        qry_to_trg_map = {}
        trg_to_qry_map = {}

        for length, op in cigar_operations:

            if op == 'M' or op == '=' or op == 'X':
                for _ in range(length):
                    qry_genome_pos = qry_tx.genome_coordinate(qry_pos)
                    qry_to_trg_map[qry_genome_pos] = (trg_pos,op)
                    trg_to_qry_map[trg_pos] = (qry_genome_pos,op)
                    qry_pos += 1
                    trg_pos += 1
            elif op == 'I':
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
            elif op == 'S':
                # Soft clipping affects only the query position
                for _ in range(length):
                    qry_genome_pos = qry_tx.genome_coordinate(qry_pos)
                    qry_to_trg_map[qry_genome_pos] = (trg_pos,op)
                    qry_pos += 1
            elif op == 'H':
                # Hard clipping affects neither query nor target positions
                continue

        return qry_to_trg_map, trg_to_qry_map

    def extract_intron_map(self, qry_tome: str, trg_aln_fname: str):
        # extract the mapping of introns from query to target
        # for each donor and acceptor site in the query - find the corresponding donor and acceptor site in the target
        donor_map = {} # holds the mapping between query and target donor sites
        acceptor_map = {} # holds the mapping between query and target acceptor sites

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

    def build_consensus_junction_map(self, donor_consensus_map, acceptor_consensus_map, query_tome):        
        qry2trg_junc_map = dict()
        trg2qry_junc_map = dict()
        for qry_tx in query_tome.transcript_it():
            for qry_it in qry_tx.introns_it():
                qry_donor = qry_it[0]-1
                qry_acceptor = qry_it[1]
                
                trg_donor = donor_consensus_map.get(qry_donor)+1
                trg_acceptor = acceptor_consensus_map.get(qry_acceptor)+1
                
                qry2trg_junc_map[(qry_donor, qry_acceptor)] = {
                    "donor": trg_donor,
                    "acceptor": trg_acceptor,
                    "seqid": qry_tx.seqid,
                    "strand": qry_tx.strand,
                    "tid": qry_tx.tid,
                    "gid": qry_tx.get_attr("gene_id")
                }
                trg2qry_junc_map[(trg_donor, trg_acceptor)] = {
                    "donor": qry_donor,
                    "acceptor": qry_acceptor,
                    "seqid": qry_tx.seqid,
                    "strand": qry_tx.strand,
                    "tid": qry_tx.tid,
                    "gid": qry_tx.get_attr("gene_id")
                }

        return qry2trg_junc_map, trg2qry_junc_map
    
    def write_consensus_gtf(self, qry_tome, aln_sam_fname, qry2trg_junc_map, out_gtf_fname):
        # using information from query transcriptome and alignment - construct a GTF file for the target genome
        
        # preload the alignments into a dictionary organized by query name
        aln_dict = {}
        for read in pysam.AlignmentFile(aln_sam_fname, "r"):
            if read.is_unmapped or not read.cigarstring or read.is_secondary:
                continue
            aln_dict[read.query_name] = read
        
        with open(out_gtf_fname,"w+") as outFP:
            # iterate over qry_tome transcripts
            for qry_tx in qry_tome.transcript_it():
                if qry_tx.tid not in aln_dict:
                    raise ValueError(f"ALIGNMENT: Transcript {qry_tx.tid} not found in alignment")
                
                read = aln_dict[qry_tx.tid]
                trg_tx_start = read.reference_start+1
                trg_tx_end = read.reference_end
                
                # iterate over the introns of the qry_transcript and construct the target transcript
                trg_exon_start = trg_tx_start
                trg_exon_end = None
                trg_exons = []
                for qry_it in qry_tx.introns_it():
                    qry_donor = qry_it[0]-1
                    qry_acceptor = qry_it[1]
                    trg_donor = qry2trg_junc_map[(qry_donor, qry_acceptor)]["donor"]
                    trg_acceptor = qry2trg_junc_map[(qry_donor, qry_acceptor)]["acceptor"]
                    
                    trg_exon_end = trg_donor
                    trg_exons.append((trg_exon_start, trg_exon_end))
 
                    trg_exon_start = trg_acceptor
                    
                trg_exons.append((trg_exon_start, trg_tx_end))

                # compute transcript alignment score
                cigar_exons = cigar2exons(read.reference_start,read.cigarstring)
                tx_aln_score = self.compare_chains(cigar_exons, trg_exons, mismatch_penalty=1, da_penalty=5)

                outFP.write(f"{read.reference_name}\t{'VIRA'}\ttranscript\t{trg_tx_start}\t{trg_tx_end}\t.\t{qry_tx.strand}\t.\ttranscript_id \"{qry_tx.tid}\"; gene_id \"{qry_tx.get_attr('gene_id')}\"; cigar \"{read.cigarstring}\"; aln_start \"{read.reference_start}\"; vira_tx_aln_score \"{tx_aln_score['score']}\";\n")
                
                # write out exons
                for exon in trg_exons:
                    outFP.write(f"{read.reference_name}\t{'VIRA'}\texon\t{exon[0]}\t{exon[1]}\t.\t{qry_tx.strand}\t.\ttranscript_id \"{qry_tx.tid}\"; gene_id \"{qry_tx.get_attr('gene_id')}\";\n")

    def compare_chains(self, chain1, chain2, mismatch_penalty=1, da_penalty=5):
        """
        Compare two chains of intervals and compute matches, mismatches, and penalties.

        Args:
            chain1: List of tuples [(start, end), ...] for the first chain.
            chain2: List of tuples [(start, end), ...] for the second chain.
            mismatch_penalty: Penalty for mismatched positions (general).
            da_penalty: Penalty for donor site mismatches.

        Returns:
            A dictionary with keys:
            - matches: Total number of matched positions.
            - mismatches: Total number of mismatched positions.
            - donor_mismatches: Total number of donor site mismatches.
            - acceptor_mismatches: Total number of acceptor site mismatches.
            - total_penalty: Sum of all mismatch penalties.
            - score: noralized score relative to the maximum possible penalty [0, 1].
        """
        def flatten_chain(chain):
            """Flatten intervals into a set of positions."""
            return set(pos for start, end in chain for pos in range(start, end + 1))
        
        def get_donor_acceptor_sites(chain):
            """Get donor (end) and acceptor (start) sites separately."""
            donors = {interval[1] for interval in chain[:-1]}  # End of exons except last
            acceptors = {interval[0] for interval in chain[1:]}  # Start of exons except first
            return donors, acceptors
        
        # Flatten the chains into sets of positions
        positions1 = flatten_chain(chain1)
        positions2 = flatten_chain(chain2)
        
        matches = positions1 & positions2
        mismatches = (positions1 - positions2) | (positions2 - positions1)
        
        # Identify donor and acceptor mismatches
        donors1, acceptors1 = get_donor_acceptor_sites(chain1)
        donors2, acceptors2 = get_donor_acceptor_sites(chain2)
        
        donor_mismatches = (donors1 - donors2) | (donors2 - donors1)
        acceptor_mismatches = (acceptors1 - acceptors2) | (acceptors2 - acceptors1)
        
        # Calculate penalties
        match_count = len(matches)
        mismatch_count = len(mismatches)
        donor_mismatch_count = len(donor_mismatches)
        acceptor_mismatch_count = len(acceptor_mismatches)
        
        total_penalty = (
            mismatch_count * mismatch_penalty +
            donor_mismatch_count * da_penalty +
            acceptor_mismatch_count * da_penalty
        )
        
        # Calculate maximum possible penalty
        total_positions = len(positions1 | positions2)
        max_penalty = (
            total_positions * mismatch_penalty +
            len(donors1 | donors2) * da_penalty +
            len(acceptors1 | acceptors2) * da_penalty
        )
        
        # Normalize penalty to [0, 1]
        normalized_penalty = total_penalty / max_penalty if max_penalty > 0 else 0

        return {
            "matches": match_count,
            "mismatches": mismatch_count,
            "donor_mismatches": donor_mismatch_count,
            "acceptor_mismatches": acceptor_mismatch_count,
            "total_penalty": total_penalty,
            "score": 1.0 - normalized_penalty,
        }
        
    def write_junction_stats(self, qry_tome, trg_tome, site_map, site_type, seq_offset, stats_outFP):
        for qry_site_pos, trg_site_data in site_map.items():
            trg_site_pos = max(trg_site_data.items(), key=lambda x: x[1])[0]
            
            total_occurences = sum(trg_site_data.values())
            consensus_occurences = trg_site_data[trg_site_pos]
            
            # extract query sequence around the site
            qry_site_start_pos = qry_site_pos if site_type == "donor" else qry_site_pos-1
            qry_site_end_pos = qry_site_pos+1 if site_type == "donor" else qry_site_pos
            
            trg_site_start_pos = trg_site_pos+1 if site_type == "donor" else trg_site_pos+1-1
            trg_site_end_pos = trg_site_pos+1+1 if site_type == "donor" else trg_site_pos+1
            
            # add +2 to the site offset to account for the 2bp length of the site itself and then add 2 positions extra on each side
            qry_site_seq = qry_tome.genome[self.qry_seqid][qry_site_start_pos-seq_offset:qry_site_end_pos+seq_offset+1].seq
            trg_site_seq = trg_tome.genome[self.trg_seqid][trg_site_start_pos-seq_offset:trg_site_end_pos+seq_offset+1].seq
            
            stats_lcs = [trg_site_pos+1, qry_site_pos, site_type, consensus_occurences/total_occurences, trg_site_seq, qry_site_seq]
            stats_outFP.write("\t".join(map(str,stats_lcs))+"\n")
        return
    
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


    def write_junction_bed_for_miniprot(self, donor_map, acceptor_map, out_bed_fname):
        """
        example with acceptor where CTA is exon and ACGAAG is intron
        ACGAAG|CTA
        123456 789
        chr  6  +  A  10
        """
        with open(out_bed_fname,"w+") as outFP:
            # process donors
            for _, trg_pos in donor_map.items():
                outFP.write(f"{self.trg_seqid}\t{trg_pos+1}\t{self.trg_strand}\tD\t100\n")
                
            # process acceptors
            for _, trg_pos in acceptor_map.items():
                outFP.write(f"{self.trg_seqid}\t{trg_pos}\t{self.trg_strand}\tA\t100\n")
    
    def build(self):
        # 1. get junc-bed from the guide annotation for the initial minimap2 alignment
        if self.guide is not None:
            cmd = ["paftools.js","gff2bed","-j",self.guide]
            print("Extracting guide junctions into BED for minimap: "+" ".join(cmd)+" > "+self.guide_junc_bed_fname)
            with open(self.guide_junc_bed_fname,"w+") as outFP:
                subprocess.call(cmd,stdout=outFP)

        # extract the reference transcript
        cmd = [self.gffread,
                "-g",self.genome,
                "-w",self.exon_nt_fasta_fname,
                self.annotation]
        print(f"Extracting reference transcript sequences: {' '.join(cmd)}")
        subprocess.call(cmd)

        # 2. align with minimap2
        cmd = [self.minimap2,"--for-only","-a",
               "-k9","-w3","--splice","--splice-flank=no","-g2k","-G9k","-A1","-B2","-O2,32","-E1,0","-b0","-C4","-z200","-un","--cap-sw-mem=0"]
        if self.guide is not None:
            cmd.extend(["--junc-bed",self.guide_junc_bed_fname,
                        "--junc-bonus","15"]) # set high bonus for guide junctions
        cmd.extend([self.target,self.exon_nt_fasta_fname])
        print("Mapping transcripts: "+" ".join(cmd)+" > "+self.exon_sam_mm2_fname)
        with open(self.exon_sam_mm2_fname,"w+") as outFP:
            subprocess.call(cmd,stdout=outFP)
        
        # 3. correct with snapper
        cmd = [self.snapper,"--reference", self.annotation, "--sam", self.exon_sam_mm2_fname, "--output", self.exon_sam_snapper_fname, "--qry_intron_match_score", "100"]
        print("Executing intron correction with snapper: " + " ".join(cmd))
        subprocess.call(cmd)

        # NOTE: snapper only performs correction such that the donor/acceptor sites match with the reference
        # however, if the intron length was a bit off (due to an extra D or I operation in CIGAR) - snapper can't do anything about it
        # we can instead perform majority vote for a consensus donor and acceptor site (for each reference donor/acceptor site since we know the mapping now from snapper)
        # and then force the consensus positions for each

        # 4. extract consensus coordinates for all donors and acceptors
        donor_map, acceptor_map = self.extract_intron_map(self.qry_tome, self.exon_sam_snapper_fname)
        # build consensus maps (qry donor/acceptor -> consensus trg donor/acceptor)
        donor_consensus_map = {site: max(counts.items(), key=lambda x: x[1])[0] for site, counts in donor_map.items()}
        acceptor_consensus_map = {site: max(counts.items(), key=lambda x: x[1])[0] for site, counts in acceptor_map.items()}
        qry2trg_junc_map, trg2qry_junc_map = self.build_consensus_junction_map(donor_consensus_map, acceptor_consensus_map, self.qry_tome)
        
        # 5. verify consensus coordinates with donors/acceptors available from the guide annotation
        if self.guide_tome is not None:
            for guide_tx in self.guide_tome.transcript_it():
                for guide_it in guide_tx.introns_it():
                    guide_donor = guide_it[0]-1
                    guide_acceptor = guide_it[1]
                    if (guide_donor,guide_acceptor) not in trg2qry_junc_map:
                        raise ValueError(f"GUIDE: Guide intron {guide_it} not found in trg2qry_junc_map")
                    
        # 6. verify that consensus junctions are indeed junctions (there is an intron > 100bp for example)
        for _, trg_data in qry2trg_junc_map.items():
            if trg_data["acceptor"] - trg_data["donor"] < 2:
                raise ValueError(f"ALIGNMENT: Consensus junction {trg_data} is not a valid junction")

        # 7. write out the transcripts into the 1st stage GTF file (before the CDSs are introduced based on miniprot and guide)
        self.write_consensus_gtf(self.qry_tome, self.exon_sam_snapper_fname, qry2trg_junc_map, self.consensus_gtf)

        # compute junction metrics
        # 1. does the alignment match consensus site?
        # 2. extract nucleotides around each donor/acceptor site on query and target genomes
        
        # collect junction sequences
        # position query_position type map_consistency sequence query_sequence
        with open(self.junction_stats_fname,"a") as stats_outFP:
            # load target tome
            trg_tome = Transcriptome()
            trg_tome.build_from_file(self.consensus_gtf)
            trg_tome.load_genome(self.target)
            
            self.write_junction_stats(self.qry_tome, trg_tome, donor_map, "donor", 5, stats_outFP)
            self.write_junction_stats(self.qry_tome, trg_tome, acceptor_map, "acceptor", 5, stats_outFP)
        
        
        """
        ===================================================================================================
        PROTEINS
        ====================================================================================================
        """
        # extract guide junctions for miniprot
        self.write_junction_bed_for_miniprot(donor_consensus_map, acceptor_consensus_map, self.junc_bed_for_miniprot_fname)
        
        # begin by extracting deduplicated CDSs from the target genome
        # and building a map of the cds duplicate tids to the gene id
        # make sure there is a single CDS per gene_id
        self.dedup_qry_cds_id_map = self.extract_gene_protein_gtf(self.annotation, self.genome, self.dedup_qry_gtf_fname)
        cmd = [self.gffread,
                "-g",self.genome,
                "-y",self.cds_aa_fasta_fname,
                "-x",self.cds_nt_fasta_fname,
                self.dedup_qry_gtf_fname]
        print(f"Extracting query protein sequences: {' '.join(cmd)}")
        subprocess.call(cmd)
        
        # get transcript_id to gene_id mapping
        tid2gid = {}
        with open(self.dedup_qry_gtf_fname,"r") as inFP:
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
                "--spsc="+self.junc_bed_for_miniprot_fname,
                "-j","0", # disable splicing model
               self.target, self.cds_aa_fasta_fname]
        print("Mapping reference proteins: " + " ".join(cmd)+" > "+miniprot_gff_fname)
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

        self.merge_exons_cds()

        return
    
    def merge_exons_cds(self):
        """
        ===============================================================
        Combine Protein annotations with Transcript annotations
        use the outputs of the previous steps to create a finalized annotation
        ===============================================================
        """

        # consensus transcript annotation without CDS
        trg_tome = Transcriptome()
        trg_tome.load_genome(self.target)
        trg_tome.build_from_file(self.consensus_gtf)
        for tx in trg_tome:
            tx.data = {"seq":"", "cds":""}
            tx.data["seq"] = tx.get_sequence(trg_tome.genome)

        # miniprot mappings
        trg_cds_tome = Transcriptome()
        trg_cds_tome.load_genome(self.target)
        trg_cds_tome.build_from_file(self.cds_gtf_fname)
        for tx in trg_cds_tome:
            tx.data = {"seq":"", "cds":""}
            tx.data["seq"] = tx.get_sequence(trg_cds_tome.genome)
            nt = tx.get_sequence(trg_cds_tome.genome,use_cds=True)
            tx.data["cds"] = translate(nt)

        # identify which guide CDS corresponds to which reference protein and which transcripts
        qry2guide_map = {}
        if self.guide is not None:
            qry2guide_map = self.build_query2guide_map(self.guide_tome, trg_tome, self.qry_tome)
            
        # iterate over target transcripts
        for target_tx in trg_tome:
            # pull the corresponding transcript from reference
            ref_tx = self.qry_tome.get_by_tid(target_tx.get_tid())

            # assign gene_id based on the reference along with other attributes
            target_tx.set_gid(ref_tx.get_attr("gene_id"))
            for e in target_tx.get_exons():
                e[2].set_gid(ref_tx.get_attr("gene_id"))
            for c in target_tx.get_cds():
                c[2].set_gid(ref_tx.get_attr("gene_id"))

        cds_choices = {"miniprot":{}, "guide":{}}
        
        #========================================================================
        #===========================   MINIPROT   ===============================
        #========================================================================
        # load the CDS for each transcript
        for target_tx in trg_tome:
            tid = target_tx.get_tid()
            # get the tid of the transcript whose cds was used in the deduplicated query
            cds_tid = self.dedup_qry_cds_id_map[tid]
            if not cds_tid in trg_cds_tome: # skipped if not mapped over
                tx.add_attribute("miniprot_mapped","0")
                continue
            else:
                tx.add_attribute("miniprot_mapped","1")
            target_cds_tx = trg_cds_tome.get_by_tid(cds_tid)

            # check compatibility of the CDS with the transcript
            target_chain = target_tx.get_chain()
            target_cds_chain = target_cds_tx.get_chain(use_cds=True)
            if not target_cds_chain == cut_chain(target_chain, target_cds_chain[0][0], target_cds_chain[-1][1]):
                tx.add_attribute("miniprot_compatibility","0")
                continue
            else:
                tx.add_attribute("miniprot_compatibility","1")
            # add the CDS to the transcript
            tmp_tx = copy.deepcopy(target_cds_tx)
            for c in target_cds_tx.get_cds():
                c[2].add_attribute("transcript_id",tid,replace=True)
                c[2].add_attribute("gene_id",target_tx.get_attr("gene_id"),replace=True)
            # get translated sequence
            nt = tmp_tx.get_sequence(trg_tome.genome,use_cds=True)
            tmp_tx.data["cds"] = translate(nt)
            cds_choices["miniprot"][tid] = tmp_tx
            
        #========================================================================
        #============================   GUIDE   =================================
        #========================================================================
        # load the guide annotation where available
        if self.guide is not None:
            for tid, (guide_tid, alignment, identity) in qry2guide_map.items():
                if guide_tid is None:
                    continue
                target_tx = copy.deepcopy(trg_tome.get_by_tid(tid))
                if target_tx is None:
                    raise ValueError(f"Transcript {tid} not found in the target genome")
                guide_tx = self.guide_tome.get_by_tid(guide_tid)
                
                # check compatibility of the CDS with the transcript
                target_chain = target_tx.get_chain()
                guide_cds_chain = guide_tx.get_chain(use_cds=True)
                if not guide_cds_chain == cut_chain(target_chain, guide_cds_chain[0][0], guide_cds_chain[-1][1]):
                    target_tx.add_attribute("guide_compatibility","0")
                    continue
                else:
                    target_tx.add_attribute("guide_compatibility","1")
                # add the CDS to the transcript
                tmp_tx = copy.deepcopy(target_tx)
                for c in guide_tx.get_cds():
                    tmp = copy.deepcopy(c[2])
                    tmp.add_attribute("transcript_id",tid,replace=True)
                    tmp_tx.add_cds(tmp)
                # get translated sequence
                nt = tmp_tx.get_sequence(trg_tome.genome,use_cds=True)
                tmp_tx.data["cds"] = translate(nt)
                cds_choices["guide"][tid] = tmp_tx
                
        # compare the CDS choices ensuring consistency
        # for each transcript compare choices
        # also ensure all agree between transcripts of the same gene
        for tx in trg_tome:
            if tx.get_tid() in cds_choices["miniprot"] and not tx.get_tid() in cds_choices["guide"]:
                miniprot_tx = cds_choices["miniprot"][tx.get_tid()]
                tx.cds = miniprot_tx.cds
                tx.add_attribute("cds_source","miniprot")
                # extract the identity attribute and add it to the transcript
                miniprot_identity = miniprot_tx.get_attr("Identity")
                if miniprot_identity is not None:
                    tx.add_attribute("identity_miniprot",miniprot_identity)
                else:
                    raise ValueError(f"Identity not found for transcript {tx.get_tid()} in miniprot")

            elif tx.get_tid() in cds_choices["guide"] and not tx.get_tid() in cds_choices["miniprot"]:
                tx.cds = cds_choices["guide"][tx.get_tid()].cds
                tx.add_attribute("cds_source","guide")
                guide_identity = None
                if tx.get_tid() in cds_choices["guide"]:
                    guide_identity = qry2guide_map[tx.get_tid()][2]
                tx.add_attribute("identity_guide",guide_identity)
            elif tx.get_tid() in cds_choices["guide"] and tx.get_tid() in cds_choices["miniprot"]:
                miniprot_tx = cds_choices["miniprot"][tx.get_tid()]
                tx.cds = cds_choices["guide"][tx.get_tid()].cds
                tx.add_attribute("cds_source","guide")
                # extract the identity attribute and add it to the transcript
                miniprot_identity = miniprot_tx.get_attr("Identity")
                if miniprot_identity is not None:
                    tx.add_attribute("identity_miniprot",miniprot_identity)
                else:
                    raise ValueError(f"Identity not found for transcript {tx.get_tid()} in miniprot")
                guide_identity = None
                if tx.get_tid() in cds_choices["guide"]:
                    guide_identity = qry2guide_map[tx.get_tid()][2]
                tx.add_attribute("identity_guide",guide_identity)
            else:
                tx.add_attribute("cds_source","none")

        # write out the final GTF file
        with open(self.output,"w+") as outFP:
            outFP.write(trg_tome.to_gtf())
        
    def build_query2guide_map(self, guide_tome: Transcriptome, trg_tome: Transcriptome, qry_tome: Transcriptome):
        qry2guide_map = {}
        
        guide_cds_map = {}
        for tx in guide_tome:
            if tx.has_cds():
                aa = tx.data["cds"]
                guide_cds_map.setdefault(aa,tx.get_tid())
        
        # load a map of all transcripts for each cds chain in the query
        qry_cds_map = {}
        for tx in qry_tome:
            if not tx.get_tid() in trg_tome: # make sure the query transcripts we are including are only those that were mapped over
                continue
            if tx.has_cds():
                aa = tx.data["cds"]
                qry_cds_map.setdefault(aa,[]).append(tx.get_tid())
                
        # for each query protein - find the corresponding guide protein
        for aa, tids in qry_cds_map.items():
            # find matching guide protein by aligning against all guide proteins
            alignment, identity, guide_tid = find_best_alignment(self.aligner, aa, guide_cds_map)
            for tid in tids:
                qry2guide_map[tid] = (guide_tid,alignment,identity)

        return qry2guide_map

    def run(self):
        try:
            self.build()
        except Exception as e:
            sys.stderr.write(f"Error running the pipeline: {str(e)}\n")
            
    def cleanup(self):
        # run cleanup routine to remove any tmp data if necessary
        try:
            if not self.keep_tmp:
                shutil.rmtree(self.tmp_dir)
        except Exception as e:
            sys.stderr.write(f"Error cleaning up: {str(e)}\n")
        return
        
        
def main():
    parser = argparse.ArgumentParser(description="By-Reference Exon and CDS Viral Genome Annotation.")

    parser.add_argument('-a', '--annotation', required=True, type=str, help='Path to the query GTF/GFF annotation file')
    parser.add_argument('-g', '--genome', required=True, type=str, help='Path to the query genome FASTA file')
    parser.add_argument('-t', '--target', required=True, type=str, help='Path to the target genome FASTA file')
    parser.add_argument('-q', '--guide', type=str, help='Optional path to the guide annotation file for the target genome. Transcripts and CDS from the guide will be used to validate the annotation')
    parser.add_argument('-o', '--output', type=str, help='Path to the output GTF file')
    
    parser.add_argument('--force-cds', action='store_true', help='Force the CDS from the guide onto the transcript chain, even if that means merging adjacent exons together (can fix alignment artifacts such as spurious introns). If the CDS does not fit the transcript chain, the transcript will be skipped')

    parser.add_argument('--gffread', type=str, default='gffread', help='Path to the gffread executable')
    parser.add_argument('--minimap2', type=str, default='minimap2', help='Path to the minimap2 executable')
    parser.add_argument('--miniprot', type=str, default='miniprot', help='Path to the miniprot executable. If not set - minimap2 will be used to align nucleotide sequence of the CDS instead')
    parser.add_argument('--snapper', type=str, default='snapper', help='Path to the snapper executable')

    parser.add_argument('--keep-tmp', action='store_true', help='Keep temporary files')
    parser.add_argument('--tmp-dir', type=str, default='./tmp', help='Directory to store temporary files')

    args = parser.parse_args()
    
    vira = Vira(args)

    try:
        vira.run()
        vira.cleanup()
    except Exception as e:
        vira.cleanup()
        sys.stderr.write(f"Error running the pipeline: {str(e)}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()