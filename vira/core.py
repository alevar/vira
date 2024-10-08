import os
import re
import copy
import shutil
import argparse
import subprocess

from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

from .utils.common import *

from .classes.txgroup import Transcriptome, Gene, Bundle
from .classes.transcript import Transcript, Object

class Vira:
    def __init__(self, args):
        
        # OPTIONS
        self.keep_tmp = args.keep_tmp
        self.tmp_dir = standard_path(args.tmp_dir)
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        # TOOLS
        self.minimap2 = args.minimap2
        self.sam2gtf = args.sam2gtf
        self.miniprot = args.miniprot
        self.gffread = args.gffread
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
        # create cop;ies of files in tmp directory for use in the pipeline
        self.annotation = self.tmp_dir+"reference.gtf"
        shutil.copyfile(args.annotation, self.annotation)
        self.genome = self.tmp_dir+"reference.fasta"
        shutil.copyfile(args.genome, self.genome)
        self.target = self.tmp_dir+"target.fasta"
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

            target_tx.data["ref2trg_map"], target_tx.data["trg2ref_map"] = self.process_cigar(target_tx.get_attr("cigar"), ref_tx, target_tx)

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

    def run(self):
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
        print(" ".join(cmd)+" > "+self.exon_sam_fname)
        with open(self.exon_sam_fname,"w+") as outFP:
            subprocess.call(cmd,stdout=outFP)
        
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

        # combine annotated transcripts, CDSs and guide annotation together
        # for each transcript/cds annotate any differences
        self.build()
        
    def compare_intron_sets(self, ref_tome: Transcriptome, target_tome: Transcriptome):
        # verifies consistency of intron mapping between reference and target genomes
        # for every reference donor/acceptor - make sure there is only one corresponding target donor/acceptor
        # raise issues otherwise

        donor_map = {} # holds the mapping between reference and target donor sites
        acceptor_map = {} # holds the mapping between reference and target acceptor sites
        for ref_tx in ref_tome.transcript_it():
            for ref_i,ref_intron in enumerate(ref_tx.introns_it()):
                # find position of the intron in the target genome
                target_tx = target_tome.get_by_tid(ref_tx.get_tid())
                
                trg_donor_pos = target_tx.data["ref2trg_map"][ref_intron[0]-1]
                trg_acceptor_pos = target_tx.data["ref2trg_map"][ref_intron[1]]

                assert trg_donor_pos is not None and trg_donor_pos[1]=="M", f"Target donor site not found for reference donor site {ref_intron[0]}"
                assert trg_acceptor_pos is not None and trg_acceptor_pos[1]=="M", f"Target acceptor site not found for reference acceptor site {ref_intron[1]}"

                donor_map.setdefault(ref_intron[0],[]).append(trg_donor_pos)
                acceptor_map.setdefault(ref_intron[1],[]).append(trg_acceptor_pos)

        # verify consistency
        for donor_pos in donor_map:
            assert len(set(donor_map[donor_pos])) == 1, f"Multiple target donor sites found for reference donor site {donor_pos}: {donor_map[donor_pos]}"
            # set the target donor site to the first element in the list
            donor_map[donor_pos] = donor_map[donor_pos][0][0]
        for acceptor_pos in acceptor_map:
            assert len(set(acceptor_map[acceptor_pos])) == 1, f"Multiple target acceptor sites found for reference acceptor site {acceptor_pos}: {acceptor_map[acceptor_pos]}"
            # set the target acceptor site to the first element in the list
            acceptor_map[acceptor_pos] = acceptor_map[acceptor_pos][0][0]

        return donor_map, acceptor_map

    def process_cigar(self, cigar_string: str, qry_tx: Transcript, trg_tx: Transcript):
        """
        Process CIGAR string to create a mapping from query genome positions
        to target positions.

        :param cigar_string: CIGAR string from alignment
        :param qry_tx: query transcript
        :param trg_tx: target transcript
        :return: A dictionary mapping qry positions to target positions
        """
        # CIGAR operation regex
        cigar_operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
        
        # Initialize positions
        qry_pos = 0
        trg_pos = trg_tx.start

        # Maps to track reference to query and query to reference
        qry_to_trg_map = {}
        trg_to_qry_map = {}

        # Process each CIGAR operation
        for length, op in cigar_operations:
            length = int(length)

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
                assert prev_qry_genome_pos in qry_to_trg_map, f"Query position {prev_qry_genome_pos} not found in qry_to_trg_map"
                assert qry_to_trg_map[prev_qry_genome_pos][1] in ['M', '=', 'X'], f"Previous operation at {prev_qry_genome_pos} is not a match"

                for _ in range(length):
                    qry_to_trg_map[prev_qry_genome_pos] = [trg_pos,op]  # Extend the previous match to cover the deletion
                    trg_pos += 1
            elif op == 'N':
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

    def reassign_tids(self, tome: Transcriptome, attr: str = "read_name"):
        # assigns the specified attribute as the transcript id
        # checks there are no duplicates
        assigned = set()
        for tx in tome:
            cur_tid = tx.get_tid()
            tid = tx.get_attr(attr)
            assert tid not in assigned, f"Duplicate transcript id {tid} found in the annotation"
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
        self.reassign_tids(target_tome)
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
            
            # target_tx.data["ref2trg_map"], target_tx.data["trg2ref_map"] = self.process_cigar(target_tx.get_attr("cigar"), ref_tx, target_tx)
            
            # check all donor and acceptor sites noting whether they are conserved or not
            ref_sj_seq = self.extract_junction_seq(ref_tx, ref_tome.genome)
            target_sj_seq = self.extract_junction_seq(target_tx, target_tome.genome)
            # compare donor acceptor pairs
            # sj_comp = self.compare_sj_seq(ref_sj_seq, target_sj_seq)
            
        # check all donor and acceptor positions noting whether they are conserved or not
        # donor_map, acceptor_map = self.compare_intron_sets(ref_tome, target_tome)
          
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
            # load a map of guide proteins against which we will be searching 
            guide_cds_map = {}
            for tx in guide_tome:
                if tx.has_cds():
                    aa = tx.data["cds"]
                    guide_cds_map.setdefault(aa,tx.get_tid())
            
            # load a map of all transcripts for each cds chain in the reference
            ref_cds_map = {}
            for tx in ref_tome:
                if not tid in target_tome: # make sure the reference transcripts we are including are only those that were mapped over
                    continue
                if tx.has_cds():
                    aa = tx.data["cds"]
                    ref_cds_map.setdefault(aa,[]).append(tx.get_tid())
                    
            # for each reference protein - find the corresponding guide protein
            for aa, tids in ref_cds_map.items():
                # find matching guide protein by aligning against all guide proteins
                alignment, identity, guide_tid = find_best_alignment(self.aligner, aa, guide_cds_map)
                if guide_tid is None:
                    continue
                # assign the guide protein to the reference proteins
                for tid in tids:
                    target_tx = target_tome.get_by_tid(tid)
                    assert target_tx is not None, f"Transcript {tid} not found in the target genome"
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
            if tx.get_tid in cds_choices["miniprot"] and not tx.get_tid() in cds_choices["guide"]:
                tx.cds = cds_choices["miniprot"][tx.get_tid()].cds
            elif tx.get_tid in cds_choices["guide"] and not tx.get_tid() in cds_choices["miniprot"]:
                tx.cds = cds_choices["guide"][tx.get_tid()].cds
            elif tx.get_tid() in cds_choices["guide"] and tx.get_tid() in cds_choices["miniprot"]:
                tx.cds = cds_choices["guide"][tx.get_tid()].cds
            else:
                continue

        # write out the final GTF file
        with open(self.output,"w+") as outFP:
            outFP.write(target_tome.to_gtf())
        
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
    parser.add_argument('--miniprot', type=str, default='miniprot', help='Path to the miniprot executable. If not set - minimap2 will be used to align nucleotide sequence of the CDS instead')

    parser.add_argument('--keep-tmp', action='store_true', help='Keep temporary files')
    parser.add_argument('--tmp-dir', type=str, default='./tmp', help='Directory to store temporary files')

    args = parser.parse_args()

    vira = Vira(args)
    vira.run()

if __name__ == "__main__":
    main()

# TODO:
# 1. Attributes
#   - whether guide or miniprot used for CDS
#   - whether guide or miniprot chains did not fit transcript chain