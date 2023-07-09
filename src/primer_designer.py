from cmath import inf
from functools import cmp_to_key
import subprocess
from primerdesigner.database import Database
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from primerdesigner.utils import parseNCBIFastaDescription
from primer3 import bindings
import re
from primerdesigner.interface import MSA
import time
from math import ceil
from tabulate import tabulate
from pyfaidx import Fasta
import platform
from functools import cmp_to_key

"""
    Primer Design tool
    APIs:
        find(database, include, exclude, output_dir)
            Pre-built pipeline for primer design

        PrimerDesigner(database, include, exclude, output_dir)
            Initialize primer design tool

        PrimerDesigner.selectRepresentative()
            select an annotated assembly with best quality
            return assembly accession

        PrimerDesigner.findHomologousGroups(reference_id)
            find homologous genes for each gene in the reference assembly
            return {<id>:[SeqRecord_1, SeqRecord_2,...]}

        PrimerDesigner.alignHomologousGroup(template_gene_id, unaligned_group)
            Use mafft to align list of unaligned SeqRecord
            return list of aligned SeqRecord

        PrimerDesigner.calculateMask(template_gene_id, aligned_group)
            calculate mask for region with low consistency
            return calculated mask(str)

        PrimerDesigner.designPrimers(aligned_group,masks)
            design primers for the aligned_group, using primer3
            return designed primer and other information in a dictionary

        PrimerDesigner.checkPrimers(primer_id, left_primer, right_primer)
            Use blast to check the amplicons of the primer set against all genomes in the database
            And determine whether the primer is valid for included and excluded organism
            return whether the primer set passed specificity check and amplicons

"""


def find(db: Database, include, exclude, workers, pick_probe=False, reference_id=None, output_dir=None):
    """
    Pre-built pipeline of primer Design tool
    """

    print("\n\nStage 0: Initalizing")
    primer_designer = PrimerDesigner(db, include, exclude, workers, pick_probe, output_dir)

    representative_id = None
    print("\n\nStage 1: Selecting one representative genome for includes")
    if reference_id is not None:
        if reference_id not in primer_designer.include_assembly_id_list:
            print("Warning: reference id not found in included organisms, assigning reference automatically")
        else:
            print("Representative genome accession already given")
            print(f"Choosing {reference_id} as reference genome")
            representative_id = reference_id

    if representative_id is None:
        representative_id = primer_designer.selectRepresentative()

    print("\n\nStage 2: Find homologous sequences for each gene of the representative assembly")
    unaligned_groups = primer_designer.findHomologousGroups(reference_id=representative_id)

    print("\n\nStage 3: Iterate through all qualified homologous group to design primers")
    primer_designer.makeBriefReport()
    for template_id, unaligned_group in unaligned_groups.items():
        aligned_group = primer_designer.alignHomologousGroup(template_id, unaligned_group)
        mask = primer_designer.calculateMask(template_id, aligned_group)
        primers = primer_designer.designPrimers(aligned_group, mask)
        primers_amplicons = {}
        primers_status = {}
        for primer_id, primer in primers.items():
            if pick_probe is False:
                status, amplicons = primer_designer.checkPrimers(primer_id, primer["PRIMER_LEFT_SEQUENCE"], primer["PRIMER_RIGHT_SEQUENCE"])
            else:
                status, amplicons = primer_designer.checkPrimers(primer_id, primer["PRIMER_LEFT_SEQUENCE"], primer["PRIMER_RIGHT_SEQUENCE"], primer["PRIMER_INTERNAL_SEQUENCE"])
            primers_amplicons[primer_id] = amplicons
            primers_status[primer_id] = status
            primer_designer.makePrimerReport(primer_id, status, amplicons, primer)
            primer_designer.addBriefReport(status, primer)
            print("")
        primer_designer.makeHomologousGroupReport(template_id, aligned_group, mask, primers_status)
        print("")


class PrimerDesigner:
    # MACROS
    TEMP_BLAST_INPUT_DIR = os.path.join("temp", "BLAST", "input")
    TEMP_BLAST_INPUT_FILE = os.path.join("temp", "BLAST", "input", "query.fna")
    TEMP_BLAST_OUTPUT_DIR = os.path.join("temp", "BLAST", "output")
    TEMP_BLAST_OUTPUT_FILE = os.path.join("temp", "BLAST", "output", "blast_results.tsv")

    TEMP_MSA_INPUT_DIR = os.path.join("temp", "MSA", "input")
    TEMP_MSA_INPUT_FILE_TEMPLATE = os.path.join("temp", "MSA", "input", "{0}_unaligned.fna")
    TEMP_MSA_OUTPUT_DIR = os.path.join("temp", "MSA", "output")
    TEMP_MSA_OUTPUT_FILE_TEMPLATE = os.path.join("temp", "MSA", "output", "{0}_aligned.fna")

    TEMP_SPECIFICITY_INPUT_DIR = os.path.join("temp", "specificity", "input")
    TEMP_SPECIFICITY_INPUT_FILE_TEMPLATE = os.path.join("temp", "specificity", "input", "primer_set_{0}.fasta")
    TEMP_SPECIFICITY_OUTPUT_DIR = os.path.join("temp", "specificity", "output")
    TEMP_SPECIFICITY_OUTPUT_FILE_TEMPLATE = os.path.join("temp", "specificity", "output", "primer_set_{0}.tsv")

    PRIMER_REPORT_DIR = os.path.join("reports", "PRIMER_REPORT")
    PRIMER_REPORT_FILE_TEMPLATE = os.path.join("reports", "PRIMER_REPORT", "{0}_primer_report.md")
    HOMOLOGOUS_GROUP_REPORT_DIR = os.path.join("reports", "HOMOLOGOUS_GROUP_REPORT")
    HOMOLOGOUS_GROUP_REPORT_FILE_TEMPLATE = os.path.join("reports", "HOMOLOGOUS_GROUP_REPORT", "{0}_homologous_group_report.md")
    BRIEF_REPORT_FILE = os.path.join("reports", "brief_report.tsv")

    BLAST_FILTER_MIN_RECORD_COVERAGE = 30
    BLAST_FILTER_MIN_RECORD_LENGTH = 300
    BLAST_FILTER_MIN_TARGET_LENGTH = 300
    BLAST_FILTER_MIN_OVERLAP = 300
    MASK_MAX_MISMATCH = 4
    MASK_CHECK_RADIUS = 20
    PRIMER_CHECK_FILTER_MAX_PRODUCT_LENGTH = 4000
    PRIMER_CHECK_FILTER_MIN_PRODUCT_LENGTH = 40
    PRIMER_CHECK_FILTER_MAX_MISMATCH = 11
    PRIMER_CHECK_FILTER_MIN_COVERAGE = 0.8
    REPORT_ALIGNMENT_ROW_LENGTH = 240
    MAX_PENALTY_FOR_UNINTENDED_PRIMER_BINDING = 6
    MAX_PENALTY_FOR_UNINTENDED_PROBE_BINDING = 4
    MAX_PENALTY_FOR_INTENDED_PRIMER_BINDING = 2
    MAX_PENALTY_FOR_INTENDED_PROBE_BINDING = 2
    QPCR_PRIMER_DESIGN_PARAMETERS = {
        "PRIMER_FIRST_BASE_INDEX": 0,
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_DNTP_CONC": 0.6,
        "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_GC_CLAMP": 1,
        "PRIMER_INSIDE_PENALTY": -1.0,
        "PRIMER_LIBERAL_BASE": 1,
        "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS": 0,
        "PRIMER_LOWERCASE_MASKING": 0,
        "PRIMER_MAX_END_GC": 3,
        "PRIMER_MAX_END_STABILITY": 9.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_MAX_HAIRPIN_TH": 47.00,
        "PRIMER_MAX_LIBRARY_MISPRIMING": 12.00,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_POLY_X": 5,
        "PRIMER_MAX_SELF_ANY": 8.00,
        "PRIMER_MAX_SELF_ANY_TH": 47.00,
        "PRIMER_MAX_SELF_END": 3.00,
        "PRIMER_MAX_SELF_END_TH": 47.00,
        "PRIMER_MAX_SIZE": 23,
        "PRIMER_MAX_TEMPLATE_MISPRIMING": 24.00,
        "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": 47.00,
        "PRIMER_MAX_TM": 62.0,
        "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4,
        "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7,
        "PRIMER_MIN_END_QUALITY": 0,
        "PRIMER_MIN_GC": 30.0,
        "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": 3,
        "PRIMER_MIN_QUALITY": 0,
        "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": 3,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MIN_TM": 58.0,
        "PRIMER_NUM_RETURN": 5,
        "PRIMER_OPT_GC_PERCENT": 50.0,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_OUTSIDE_PENALTY": 0.0,
        "PRIMER_PAIR_MAX_COMPL_ANY": 8.00,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 47.00,
        "PRIMER_PAIR_MAX_COMPL_END": 3.00,
        "PRIMER_PAIR_MAX_COMPL_END_TH": 47.00,
        "PRIMER_PAIR_MAX_DIFF_TM": 3,
        "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING": 12.00,
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING": 12.00,
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH": 47.00,
        "PRIMER_PAIR_WT_COMPL_ANY": 0.1,
        "PRIMER_PAIR_WT_COMPL_ANY_TH": 0.1,
        "PRIMER_PAIR_WT_COMPL_END": 0.1,
        "PRIMER_PAIR_WT_COMPL_END_TH": 0.1,
        "PRIMER_PAIR_WT_DIFF_TM": 0.1,
        "PRIMER_PAIR_WT_IO_PENALTY": 0.0,
        "PRIMER_PAIR_WT_LIBRARY_MISPRIMING": 0.0,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.0,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.0,
        "PRIMER_PAIR_WT_PRODUCT_TM_GT": 0.0,
        "PRIMER_PAIR_WT_PRODUCT_TM_LT": 0.0,
        "PRIMER_PAIR_WT_PR_PENALTY": 1.0,
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING": 0.0,
        "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH": 0.0,
        "PRIMER_PICK_ANYWAY": 0,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PRODUCT_SIZE_RANGE": [70, 150],
        "PRIMER_QUALITY_RANGE_MAX": 100,
        "PRIMER_QUALITY_RANGE_MIN": 0,
        "PRIMER_SALT_CORRECTIONS": 1,
        "PRIMER_SALT_DIVALENT": 1.5,
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_SEQUENCING_ACCURACY": 20,
        "PRIMER_SEQUENCING_INTERVAL": 250,
        "PRIMER_SEQUENCING_LEAD": 50,
        "PRIMER_SEQUENCING_SPACING": 500,
        "PRIMER_TASK": "generic",
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": 0,
        "PRIMER_TM_FORMULA": 1,
        "PRIMER_WT_END_QUAL": 0.0,
        "PRIMER_WT_END_STABILITY": 0.0,
        "PRIMER_WT_GC_PERCENT_GT": 0.1,
        "PRIMER_WT_GC_PERCENT_LT": 0.1,
        "PRIMER_WT_HAIRPIN_TH": 0.0,
        "PRIMER_WT_LIBRARY_MISPRIMING": 0.0,
        "PRIMER_WT_NUM_NS": 0.0,
        "PRIMER_WT_POS_PENALTY": 0.0,
        "PRIMER_WT_SELF_ANY": 0.0,
        "PRIMER_WT_SELF_ANY_TH": 0.0,
        "PRIMER_WT_SELF_END": 0.0,
        "PRIMER_WT_SELF_END_TH": 0.0,
        "PRIMER_WT_SEQ_QUAL": 0.0,
        "PRIMER_WT_SIZE_GT": 0.5,
        "PRIMER_WT_SIZE_LT": 0.5,
        "PRIMER_WT_TEMPLATE_MISPRIMING": 0.0,
        "PRIMER_WT_TEMPLATE_MISPRIMING_TH": 0.0,
        "PRIMER_WT_TM_GT": 0.5,
        "PRIMER_WT_TM_LT": 0.5,
    }

    def makeBriefReport(self):
        with open(self.brief_report_file, "w") as f:
            f.truncate()
            f.write("# Genome-wide target-specific primer searching tool, version 1.1\n")
            f.write(f"# Primers are supposed to specific to {', '.join(self.raw_include)}\n")
            f.write(f"# Using {', '.join(self.raw_exclude)} for misprime examination\n")
            f.write(f"# Database config path: {self.db.DATABASE_CONFIG_FILE}\n")
            if self.pick_probe is False:
                f.write(
                    "\t".join(
                        [
                            "template sequence accession",
                            "template protein function",
                            "template assembly id",
                            "template gene id",
                            "template organism",
                            "pair penalty",
                            "left penalty",
                            "right penalty",
                            "left sequence",
                            "right sequence",
                            "left GC%",
                            "right GC%",
                            "left TM",
                            "right TM",
                            "product size",
                        ]
                    )
                )
            else:
                f.write(
                    "\t".join(
                        [
                            "template sequence accession",
                            "template protein function",
                            "template assembly id",
                            "template gene id",
                            "template organism",
                            "pair penalty",
                            "left penalty",
                            "right penalty",
                            "internal penalty",
                            "left sequence",
                            "right sequence",
                            "internal sequence",
                            "left GC%",
                            "right GC%",
                            "internal GC%",
                            "left TM",
                            "right TM",
                            "internal TM",
                            "product size",
                        ]
                    )
                )
            f.write("\n")

    def addBriefReport(self, status, primer):
        if status:
            with open(self.brief_report_file, "a") as f:
                if self.pick_probe is False:
                    f.write(
                        "\t".join(
                            [
                                primer["TEMPLATE_SEQUENCE_ACCESSION"],
                                primer["TEMPLATE_PRODUCT"],
                                primer["TEMPLATE_ASSEMBLY_ID"],
                                primer["TEMPLATE_GENE_ID"],
                                primer["TEMPLATE_ORGANISM"],
                                primer["PRIMER_PAIR_PENALTY"],
                                primer["PRIMER_LEFT_PENALTY"],
                                primer["PRIMER_RIGHT_PENALTY"],
                                primer["PRIMER_LEFT_SEQUENCE"],
                                primer["PRIMER_RIGHT_SEQUENCE"],
                                primer["PRIMER_LEFT_GC_PERCENT"],
                                primer["PRIMER_RIGHT_GC_PERCENT"],
                                primer["PRIMER_LEFT_TM"],
                                primer["PRIMER_RIGHT_TM"],
                                primer["PRIMER_PAIR_PRODUCT_SIZE"],
                            ]
                        )
                    )
                else:
                    f.write(
                        "\t".join(
                            [
                                primer["TEMPLATE_SEQUENCE_ACCESSION"],
                                primer["TEMPLATE_PRODUCT"],
                                primer["TEMPLATE_ASSEMBLY_ID"],
                                primer["TEMPLATE_GENE_ID"],
                                primer["TEMPLATE_ORGANISM"],
                                primer["PRIMER_PAIR_PENALTY"],
                                primer["PRIMER_LEFT_PENALTY"],
                                primer["PRIMER_RIGHT_PENALTY"],
                                primer["PRIMER_INTERNAL_PENALTY"],
                                primer["PRIMER_LEFT_SEQUENCE"],
                                primer["PRIMER_RIGHT_SEQUENCE"],
                                primer["PRIMER_INTERNAL_SEQUENCE"],
                                primer["PRIMER_LEFT_GC_PERCENT"],
                                primer["PRIMER_RIGHT_GC_PERCENT"],
                                primer["PRIMER_INTERNAL_GC_PERCENT"],
                                primer["PRIMER_LEFT_TM"],
                                primer["PRIMER_RIGHT_TM"],
                                primer["PRIMER_INTERNAL_TM"],
                                primer["PRIMER_PAIR_PRODUCT_SIZE"],
                            ]
                        )
                    )
                f.write("\n")

    def checkDependency(self):
        sys = platform.system()
        if sys != "Linux":
            raise RuntimeError(f"This program is not supported on this platform: {sys}")
        dependency_valid = 1
        try:
            mafft = subprocess.run(["which", "mafft"], capture_output=True, encoding="utf-8")
            assert mafft.stdout != "" and mafft.stderr == ""
            print(f"Found mafft at {mafft.stdout[:-1]}")
        except BaseException:
            print("Mafft is not installed, this program is required for multiple sequence alignment, please check their official website: https://mafft.cbrc.jp/alignment/software/")
            dependency_valid = 0

        try:
            blastn = subprocess.run(["which", "blastn"], capture_output=True, encoding="utf-8")
            assert blastn.stdout != "" and blastn.stderr == ""
            print(f"Found blastn at {blastn.stdout[:-1]}")
        except BaseException:
            print("BLAST is not installed, this program is required for homologous group identification, please check their official website: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")

        if not dependency_valid:
            raise RuntimeError("Dependency invalid")

    # init method is called if seperate steps need to be debugged
    # If you want to find candidates, the find method is recommended
    def __init__(self, db: Database, include, exclude, workers, pick_probe=False, output_dir=None):
        time_stamp = time.strftime("%Y%m%d%H%M%S", time.localtime())
        if output_dir is None:
            output_dir = os.path.join(os.path.abspath("."), f"results_{time_stamp}")

        # check if output_dir is a valid path
        try:
            output_dir = os.path.abspath(output_dir)
        except BaseException:
            raise ValueError(f"Output path invalid at: {output_dir} ")

        # check if output_dir already exists, if not then create the directory
        try:
            assert not os.path.isfile(output_dir)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            elif len(os.listdir(output_dir)) != 0:
                print(f"Found existing file at output directory: {output_dir}")
                c = input("Delete existing files? (Y/n): ")
                while c != "Y" and c != "n":
                    print(f"Can't recognize input character: {c}")
                    c = input("Delete existing files? (Y/n): ")
                if c == "Y":
                    print("Deleting existing files ...", end=" ")
                    subprocess.run(["rm", "-rf", output_dir])
                    os.makedirs(output_dir)
                    print("Done")
                elif c == "n":
                    print("Continue anyway")

        except BaseException:
            raise RuntimeError(f"Can not create output directory at: {output_dir}")

        self.checkDependency()

        self.output_dir = os.path.abspath(output_dir)
        self.db = db
        self.assembly_data = {}
        self.raw_include = include
        self.raw_exclude = exclude
        self.workers = workers
        self.pick_probe = pick_probe

        # other directory
        self.temp_blast_input_dir = os.path.join(self.output_dir, PrimerDesigner.TEMP_BLAST_INPUT_DIR)
        self.temp_blast_output_dir = os.path.join(self.output_dir, PrimerDesigner.TEMP_BLAST_OUTPUT_DIR)
        self.temp_msa_input_dir = os.path.join(self.output_dir, PrimerDesigner.TEMP_MSA_INPUT_DIR)
        self.temp_msa_output_dir = os.path.join(self.output_dir, PrimerDesigner.TEMP_MSA_OUTPUT_DIR)
        self.temp_specificity_input_dir = os.path.join(self.output_dir, PrimerDesigner.TEMP_SPECIFICITY_INPUT_DIR)
        self.temp_specificity_output_dir = os.path.join(self.output_dir, PrimerDesigner.TEMP_SPECIFICITY_OUTPUT_DIR)

        self.temp_blast_input_file = os.path.join(self.output_dir, PrimerDesigner.TEMP_BLAST_INPUT_FILE)
        self.temp_blast_output_file = os.path.join(self.output_dir, PrimerDesigner.TEMP_BLAST_OUTPUT_FILE)
        self.temp_msa_input_file_template = os.path.join(self.output_dir, PrimerDesigner.TEMP_MSA_INPUT_FILE_TEMPLATE)
        self.temp_msa_output_file_template = os.path.join(self.output_dir, PrimerDesigner.TEMP_MSA_OUTPUT_FILE_TEMPLATE)
        self.temp_specificity_input_file_template = os.path.join(self.output_dir, PrimerDesigner.TEMP_SPECIFICITY_INPUT_FILE_TEMPLATE)
        self.temp_specificity_output_file_template = os.path.join(self.output_dir, PrimerDesigner.TEMP_SPECIFICITY_OUTPUT_FILE_TEMPLATE)

        self.primer_report_dir = os.path.join(self.output_dir, PrimerDesigner.PRIMER_REPORT_DIR)
        self.primer_report_file_template = os.path.join(self.output_dir, PrimerDesigner.PRIMER_REPORT_FILE_TEMPLATE)
        self.homologous_group_report_dir = os.path.join(self.output_dir, PrimerDesigner.HOMOLOGOUS_GROUP_REPORT_DIR)
        self.homologous_group_report_file_template = os.path.join(self.output_dir, PrimerDesigner.HOMOLOGOUS_GROUP_REPORT_FILE_TEMPLATE)
        self.brief_report_file = os.path.join(self.output_dir, PrimerDesigner.BRIEF_REPORT_FILE)

        if not os.path.exists(self.temp_blast_input_dir):
            os.makedirs(self.temp_blast_input_dir)
        if not os.path.exists(self.temp_blast_output_dir):
            os.makedirs(self.temp_blast_output_dir)
        if not os.path.exists(self.temp_msa_input_dir):
            os.makedirs(self.temp_msa_input_dir)
        if not os.path.exists(self.temp_msa_output_dir):
            os.makedirs(self.temp_msa_output_dir)
        if not os.path.exists(self.temp_specificity_input_dir):
            os.makedirs(self.temp_specificity_input_dir)
        if not os.path.exists(self.temp_specificity_output_dir):
            os.makedirs(self.temp_specificity_output_dir)
        if not os.path.exists(self.primer_report_dir):
            os.makedirs(self.primer_report_dir)
        if not os.path.exists(self.homologous_group_report_dir):
            os.makedirs(self.homologous_group_report_dir)

        self.PRIMER_GLOBAL_ARGS = PrimerDesigner.QPCR_PRIMER_DESIGN_PARAMETERS
        if pick_probe is True:
            self.PRIMER_GLOBAL_ARGS.update(
                {
                    "PRIMER_PICK_INTERNAL_OLIGO": 1,
                    "PRIMER_PAIR_WT_IO_PENALTY": 1.0,
                    "PRIMER_INTERNAL_DNA_CONC": 50.0,
                    "PRIMER_INTERNAL_DNTP_CONC": 0.0,
                    "PRIMER_INTERNAL_MAX_GC": 80.0,
                    "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 47.00,
                    "PRIMER_INTERNAL_MAX_LIBRARY_MISHYB": 12.00,
                    "PRIMER_INTERNAL_MAX_NS_ACCEPTED": 0,
                    "PRIMER_INTERNAL_MAX_POLY_X": 5,
                    "PRIMER_INTERNAL_MAX_SELF_ANY": 12.00,
                    "PRIMER_INTERNAL_MAX_SELF_ANY_TH": 47.00,
                    "PRIMER_INTERNAL_MAX_SELF_END": 12.00,
                    "PRIMER_INTERNAL_MAX_SELF_END_TH": 47.00,
                    "PRIMER_INTERNAL_MAX_SIZE": 30,
                    "PRIMER_INTERNAL_MAX_TM": 72.0,
                    "PRIMER_INTERNAL_MIN_GC": 30.0,
                    "PRIMER_INTERNAL_MIN_QUALITY": 0,
                    "PRIMER_INTERNAL_MIN_SIZE": 13,
                    "PRIMER_INTERNAL_MIN_TM": 68.0,
                    "PRIMER_INTERNAL_OPT_GC_PERCENT": 50.0,
                    "PRIMER_INTERNAL_OPT_SIZE": 25,
                    "PRIMER_INTERNAL_OPT_TM": 70.0,
                    "PRIMER_INTERNAL_SALT_DIVALENT": 0.0,
                    "PRIMER_INTERNAL_SALT_MONOVALENT": 50.0,
                    "PRIMER_INTERNAL_WT_END_QUAL": 0.0,
                    "PRIMER_INTERNAL_WT_GC_PERCENT_GT": 0.1,
                    "PRIMER_INTERNAL_WT_GC_PERCENT_LT": 0.1,
                    "PRIMER_INTERNAL_WT_HAIRPIN_TH": 0.0,
                    "PRIMER_INTERNAL_WT_LIBRARY_MISHYB": 0.0,
                    "PRIMER_INTERNAL_WT_NUM_NS": 0.0,
                    "PRIMER_INTERNAL_WT_SELF_ANY": 0.1,
                    "PRIMER_INTERNAL_WT_SELF_ANY_TH": 0.0,
                    "PRIMER_INTERNAL_WT_SELF_END": 0.1,
                    "PRIMER_INTERNAL_WT_SELF_END_TH": 0.0,
                    "PRIMER_INTERNAL_WT_SEQ_QUAL": 0.0,
                    "PRIMER_INTERNAL_WT_SIZE_GT": 1.0,
                    "PRIMER_INTERNAL_WT_SIZE_LT": 1.0,
                    "PRIMER_INTERNAL_WT_TM_GT": 1.0,
                    "PRIMER_INTERNAL_WT_TM_LT": 1.0,
                }
            )

        # check if input is sane
        # map identifier to id

        self.exclude_assembly_id_list = []
        exclude_not_mapped = []
        for identifier in exclude:
            cursor = self.db.cursor.execute(
                f'SELECT ID FROM ASSEMBLY WHERE \
                ORGANISM_NAME="{identifier}" OR \
                ORGANISM_ID="{identifier}" OR \
                ORGANISM_GENUS="{identifier}" OR \
                ORGANISM_SPECIES="{identifier}" OR \
                SCIENTIFIC_NAME="{identifier}"'
            )
            unmapped = 1
            for row in cursor:
                self.exclude_assembly_id_list.append(row[0])
                unmapped = 0
            if unmapped:
                exclude_not_mapped.append(identifier)

        self.include_assembly_id_list = []
        include_not_mapped = []
        for identifier in include:
            cursor = self.db.cursor.execute(
                f'SELECT ID FROM ASSEMBLY WHERE \
                ORGANISM_NAME="{identifier}" OR \
                ORGANISM_ID="{identifier}" OR \
                ORGANISM_GENUS="{identifier}" OR \
                ORGANISM_SPECIES="{identifier}" OR \
                SCIENTIFIC_NAME="{identifier}"'
            )
            unmapped = 1
            for row in cursor:
                if row[0] not in self.exclude_assembly_id_list:
                    self.include_assembly_id_list.append(row[0])
                unmapped = 0
            if unmapped:
                include_not_mapped.append(identifier)

        self.include_assembly_id_list = list(set(self.include_assembly_id_list))
        self.exclude_assembly_id_list = list(set(self.exclude_assembly_id_list))

        # output unmapped and mapped organism
        for unmapped in include_not_mapped:
            print(f"Warning: No genomic information about include identifier: {unmapped}")

        for unmapped in exclude_not_mapped:
            print(f"Warning: No genomic information about exclude identifier: {unmapped}")

        print("Include:")
        for mapped in self.include_assembly_id_list:
            organism_name = self.db.cursor.execute(f'SELECT ORGANISM_NAME FROM ASSEMBLY WHERE ID="{mapped}"').fetchone()[0]
            print(f"{mapped}: {organism_name}")

        print("Exclude:")
        for mapped in self.exclude_assembly_id_list:
            organism_name = self.db.cursor.execute(f'SELECT ORGANISM_NAME FROM ASSEMBLY WHERE ID="{mapped}"').fetchone()[0]
            print(f"{mapped}: {organism_name}")

        temp_query_string = '("' + '","'.join(self.include_assembly_id_list) + '")'
        include_annotation_count = self.db.cursor.execute(f'SELECT COUNT(*) FROM FILE WHERE ASSEMBLY IN {temp_query_string} AND TYPE="GFF3"').fetchone()[0]
        if include_annotation_count <= 0:
            raise ValueError("At least one record should have annotation in include group")

    def selectRepresentative(self):
        assemblies_data = {}  # select assemblyies with parsable annotation
        assemblies_id = []
        # select
        cursor = self.db.cursor.execute('SELECT ID, LEVEL, ORGANISM_NAME ,SUBMISSION_DATE, CONTIG_N50, TOTAL_UNGAPPED_LENGTH FROM ASSEMBLY WHERE ID IN (SELECT ASSEMBLY FROM FILE WHERE TYPE = "GFF3")')
        for row in cursor:
            if row[0] not in self.include_assembly_id_list:
                continue
            assemblies_data[row[0]] = {
                "level": row[1],
                "organism_name": row[2],
                "submission_date": row[3],
                "contig_n50": int(row[4]),
                "total_ungapped_length": int(row[5]),
            }
            assemblies_id.append(row[0])

        def __cmp(accession1: str, accession2: str) -> int:
            assembly_data_1 = assemblies_data[accession1]
            assembly_data_2 = assemblies_data[accession2]
            score = {"Contig": 0, "Scaffold": 1, "Chromosome": 2, "Complete Genome": 3}
            if score[assembly_data_1["level"]] > score[assembly_data_2["level"]]:
                return 1
            elif score[assembly_data_1["level"]] < score[assembly_data_2["level"]]:
                return -1
            else:
                if assembly_data_1["contig_n50"] > assembly_data_2["contig_n50"]:
                    return 1
                elif assembly_data_1["contig_n50"] < assembly_data_2["contig_n50"]:
                    return -1
            return 0

        if len(assemblies_id) == 0:
            raise ValueError("At least one of the includes should have annotation")
        sorted_assemblies_id = sorted(assemblies_id, key=cmp_to_key(__cmp), reverse=True)
        print(f"Select {sorted_assemblies_id[0]} as representative genome")
        print(f'Name: {assemblies_data[sorted_assemblies_id[0]]["organism_name"]}')
        print(f'Level: {assemblies_data[sorted_assemblies_id[0]]["level"]}')
        print(f'N50: {assemblies_data[sorted_assemblies_id[0]]["contig_n50"]}')
        print(f'Total ungapped length: {assemblies_data[sorted_assemblies_id[0]]["total_ungapped_length"]}')
        print(f'Submission date: {assemblies_data[sorted_assemblies_id[0]]["submission_date"]}')
        return sorted_assemblies_id[0]

    def getSequenceFromID(self, assembly_id, gene_id):
        genes = []
        cursor = self.db.cursor.execute(f'SELECT ID,CHROMOSOME,START,END,STRAND,GENE_BIOTYPE,PRODUCT FROM GENE WHERE ASSEMBLY = "{assembly_id}" AND ID = "{gene_id}"')
        for row in cursor:
            genes.append({"id": row[0], "chromosome": row[1], "start": row[2], "end": row[3], "strand": row[4], "gene_biotype": row[5], "product": row[6]})
        assert len(genes) == 1
        gene = genes[0]
        chromosome = gene["chromosome"]

        nucleotide_fasta_files = []
        cursor = self.db.cursor.execute(f'SELECT PATH FROM FILE WHERE TYPE = "GENOMIC_NUCLEOTIDE_FASTA" AND ASSEMBLY = "{assembly_id}"')
        for row in cursor:
            nucleotide_fasta_files.append(Fasta(row[0]))

        sequence = None
        for file in nucleotide_fasta_files:
            if chromosome in file.keys():
                sequence = file[chromosome][int(gene["start"]) - 1 : int(gene["end"])]
                if gene["strand"] == "-":
                    sequence = sequence.reverse.complement
                break
        try:
            assert sequence is not None
        except AssertionError:
            raise KeyError(f"Sequence not found, assembly: {assembly_id}, gene: {gene_id}")

        return SeqRecord(
            seq=Seq(str(sequence)),
            id=f'{assembly_id}-{gene["id"]}',
            description=f'[assembly_id={assembly_id}] [id={gene["id"]}] [gene_biotype={gene["gene_biotype"]}] [product={gene["product"]}] [chromosome={gene["chromosome"]}] [start={gene["start"]}] [end={gene["end"]}] [strand={gene["strand"]}]',
        )

    def findHomologousGroups(self, reference_id):
        # write gene fasta
        genes = []
        # get gene data from database
        cursor = self.db.cursor.execute(f'SELECT ID,CHROMOSOME,START,END,STRAND,GENE_BIOTYPE,PRODUCT FROM GENE WHERE ASSEMBLY = "{reference_id}"')
        for row in cursor:
            genes.append({"id": row[0], "chromosome": row[1], "start": row[2], "end": row[3], "strand": row[4], "gene_biotype": row[5], "product": row[6]})

        # merge to file
        reference_nucleotide_fasta_files = []
        cursor = self.db.cursor.execute(f'SELECT PATH FROM FILE WHERE TYPE = "GENOMIC_NUCLEOTIDE_FASTA" AND ASSEMBLY = "{reference_id}"')
        for row in cursor:
            reference_nucleotide_fasta_files.append(Fasta(row[0]))

        with open(self.temp_blast_input_file, "w") as f:
            f.truncate()
        for gene in genes:
            chromosome = gene["chromosome"]
            for file in reference_nucleotide_fasta_files:
                if chromosome in file.keys():
                    sequence = file[chromosome][int(gene["start"]) - 1 : int(gene["end"])]
                    if gene["strand"] == "-":
                        sequence = sequence.reverse.complement
                    with open(self.temp_blast_input_file, "a") as f:
                        f.write(
                            f'>{reference_id}-{gene["id"]} sequence of {gene["id"]} from assembly {reference_id} [assembly_id={reference_id}] [id={gene["id"]}] [gene_biotype={gene["gene_biotype"]}] [product={gene["product"]}] [chromosome={gene["chromosome"]}] [start={gene["start"]}] [end={gene["end"]}] [strand={gene["strand"]}]\n'
                        )
                        f.write(str(sequence) + "\n")
                    break

        # run blast
        BLAST_COMMAND = [
            "blastn",
            "-outfmt",
            "6 qacc saccver qcovhsp evalue length qstart qend sstart send sseq sstrand",
            "-num_threads",
            str(self.workers),
            "-query",
            self.temp_blast_input_file,
            "-db",
            self.db.blast_db_prefix,
            "-out",
            self.temp_blast_output_file,
        ]

        print("Running BLAST")
        subprocess.run(BLAST_COMMAND)
        # parse blast results

        print("Parsing BLAST results")

        hit_count = {}
        organism_name = {}
        cursor = self.db.cursor.execute("SELECT ID, SCIENTIFIC_NAME FROM ASSEMBLY")
        for row in cursor:
            hit_count[row[0]] = 0
            organism_name[row[0]] = row[1]

        total_genes_num = len(genes)

        homologous_groups = {}
        total_sequence = 0
        with open(self.temp_blast_output_file) as f:
            for line in f.readlines():
                if line[0] == "#":
                    continue
                # Parse hit
                hit = dict(
                    zip(
                        ["template_accession", "sequence_accession", "coverage", "evalue", "alignment_length", "qstart", "qend", "sstart", "send", "sequence", "strand"],
                        line.replace("\n", "").split("\t"),
                    )
                )
                hit["assembly_id"] = hit["sequence_accession"].split("-")[0]

                hit_count[hit["assembly_id"]] += 1

                # Filter hit
                if float(hit["coverage"]) < PrimerDesigner.BLAST_FILTER_MIN_RECORD_COVERAGE or abs(int(hit["send"]) - int(hit["sstart"])) < PrimerDesigner.BLAST_FILTER_MIN_RECORD_LENGTH:  # low length
                    continue
                # Add hit
                if hit["template_accession"] not in homologous_groups.keys():
                    # if the template of this group hasn't been processed
                    homologous_groups[hit["template_accession"]] = []

                    """
                        if hit["assembly_id"] in templates[hit["template_accession"]].keys():
                        # if not the first hit in the homologous group
                        previous: SeqRecord = templates[hit["template_accession"]][hit["assembly_id"]]
                        description_dict: Dict[str, str] = parseNCBIFastaDescription(previous.description)
                        if float(hit["evalue"]) > float(description_dict["evalue"]):
                            continue
                        else:
                            # passed filtering
                            templates[hit["template_accession"]][hit["assembly_id"]] = SeqRecord(
                                seq=Seq(hit["sequence"].replace("-", "").upper()),
                                id=templates[hit["template_accession"]][hit["assembly_id"]].id,
                                description=f"[organism={hit_assembly.data['sciname']}] [sequence_accession={hit['sequence_accession']}] [assembly_id={hit_assembly.data['accession']}] [evalue={hit['evalue']}] [qstart={hit['qstart']}] [qend={hit['qend']}] genomic coding sequence, partial, homologous sequence of template {hit['template_accession']}",)
                    """

                total_sequence += 1
                homologous_groups[hit["template_accession"]].append(
                    SeqRecord(
                        seq=Seq(hit["sequence"].replace("-", "").upper()),
                        id="Segment_" + str(total_sequence),
                        description=f"[sequence_accession={hit['sequence_accession']}] [assembly_id={hit['assembly_id']}] [evalue={hit['evalue']}] [qstart={hit['qstart']}] [qend={hit['qend']}] [template_accession={hit['template_accession']}]genomic nucleotide sequence, partial, homologous sequence of template {hit['template_accession']}",
                    )
                )

        print(f"Found {len(homologous_groups.keys())} homologous groups")

        include_hit_count = {}
        exclude_hit_count = {}
        others_hit_count = {}
        for key, value in hit_count.items():
            if key in self.exclude_assembly_id_list:
                exclude_hit_count[key] = value
            elif key in self.include_assembly_id_list:
                include_hit_count[key] = value
            else:
                others_hit_count[key] = value

        sorted_include_hit_count = sorted(include_hit_count.items(), key=lambda x: x[1], reverse=True)
        sorted_exclude_hit_count = sorted(exclude_hit_count.items(), key=lambda x: x[1], reverse=True)
        sorted_others_hit_count = sorted(others_hit_count.items(), key=lambda x: x[1], reverse=True)

        print("Type\tAccession\tOrganism Name\tCoverage Rate%")
        for key, value in sorted_include_hit_count:
            print(f"Include\t{key}\t{organism_name[key]}\t{round(hit_count[key]/total_genes_num*100,2)}%")

        for key, value in sorted_exclude_hit_count:
            print(f"Exclude\t{key}\t{organism_name[key]}\t{round(hit_count[key]/total_genes_num*100,2)}%")

        for key, value in sorted_others_hit_count:
            print(f"OTHERS\t{key}\t{organism_name[key]}\t{round(hit_count[key]/total_genes_num*100,2)}%")

        print("Filtering BLAST results")
        # filter based on whether the sequence covers all assemblies needed to be included
        for template_id, homologous_group in list(homologous_groups.items()):
            covered_assembly_id_list = set()
            not_covered_assembly_id_list = set()
            for record in homologous_group:
                parsed_description = parseNCBIFastaDescription(record.description)
                covered_assembly_id_list.add(parsed_description["assembly_id"])
            span_include = True
            for accession in self.include_assembly_id_list:
                if accession not in covered_assembly_id_list:
                    not_covered_assembly_id_list.add(accession)
                    span_include = False
            if not span_include:
                # print(f'{template_id} failed because homologous sequence can not cover {", ".join(not_covered_assembly_id_list)}')
                del homologous_groups[template_id]

        # filter based on minimal overlap range
        for template_id, homologous_group in list(homologous_groups.items()):
            max_qstart = -1
            min_qend = inf
            for record in homologous_group:
                parsed_description = parseNCBIFastaDescription(record.description)
                assembly_id = parsed_description["assembly_id"]
                if assembly_id not in self.include_assembly_id_list:
                    continue
                qstart = min(int(parsed_description["qstart"]), int(parsed_description["qend"]))
                qend = max(int(parsed_description["qstart"]), int(parsed_description["qend"]))
                max_qstart = max(qstart, max_qstart)
                min_qend = min(qend, min_qend)
            if min_qend - max_qstart + 1 < self.BLAST_FILTER_MIN_OVERLAP:
                del homologous_groups[template_id]

        print(f"Remain {len(homologous_groups)} homologous groups")
        return homologous_groups

    def alignHomologousGroup(self, template_id, unaligned_group):
        print(f"Aligning homologous genes of template: {template_id}")
        input_path = self.temp_msa_input_file_template.format(template_id)
        output_path = self.temp_msa_output_file_template.format(template_id)
        with open(input_path, "w") as f:
            f.truncate()
        with open(output_path, "w") as f:
            f.truncate()

        template_assembly_id = template_id.split("-")[0]
        template_gene_id = template_id.replace(template_assembly_id + "-", "")
        template_record = self.getSequenceFromID(template_assembly_id, template_gene_id)
        template_record.id = "TEMPLATE_" + template_record.id
        unaligned_group.append(template_record)
        SeqIO.write(unaligned_group, input_path, "fasta")
        try:
            MSA(input_fasta_path=input_path, output_fasta_path=output_path, workers=self.workers)
        except BaseException:
            raise RuntimeError(f"Multiple sequence alignment failed for file {input_path}, please recheck the input file")

        aligned_group = [record for record in SeqIO.parse(output_path, "fasta")]
        return aligned_group

    def calculateMask(self, template_id, aligned_group):
        print(f"Calculating mask for {template_id}")
        mask = ["t" for _ in range(len(str(aligned_group[0].seq)))]
        template = None
        for record in aligned_group:
            if record.id[0:8] == "TEMPLATE":
                template = record
        assert template is not None

        for index in range(len(str(aligned_group[0].seq))):
            start = max(index - self.MASK_CHECK_RADIUS, 0)
            end = min(index + self.MASK_CHECK_RADIUS, len(str(aligned_group[0].seq)))
            template_segment = str(template.seq)[start : end + 1]
            mismatch_score = {id: inf for id in self.include_assembly_id_list}
            for record in aligned_group:
                parsed_description = parseNCBIFastaDescription(record.description)
                if parsed_description["assembly_id"] not in self.include_assembly_id_list:
                    continue
                segment = str(record.seq)[start : end + 1]
                assert len(template_segment) == len(segment)
                temp_mismatch_score = 0
                for i in range(len(template_segment)):
                    if template_segment[i] != segment[i]:
                        if template_segment[i] == "-" or segment[i] == "-":
                            temp_mismatch_score = inf
                        else:
                            temp_mismatch_score += 1
                mismatch_score[parsed_description["assembly_id"]] = min(mismatch_score[parsed_description["assembly_id"]], temp_mismatch_score)

            if max(mismatch_score.values()) > self.MASK_MAX_MISMATCH:
                mask[index] = "f"
            elif template[index] == "-":
                mask[index] = "-"
        return mask

    def designPrimers(self, aligned_group, mask):
        assert len(mask) == len(str(aligned_group[0].seq))
        length = len(mask)
        misprime_lib = {}
        template = None
        for index, record in enumerate(aligned_group):
            parsed_description = parseNCBIFastaDescription(record.description)
            if parsed_description["assembly_id"] in self.exclude_assembly_id_list:
                misprime_lib[record.id] = str(record.seq).replace("-", "").upper()
            if record.id[0:8] == "TEMPLATE":
                template = record
        assert template is not None
        template_parsed_description = parseNCBIFastaDescription(template.description)

        sequence = template.seq
        sequence_for_design = ""
        primer_site = ""

        for i in range(length):
            base = sequence[i].upper()
            value = mask[i]
            # 如果这个位置在所有目标序列都是空位那可以跳过
            # 如果这个位置只在这个目标序列不是空位那不可以跳过，要增加占位符
            if base == "-" and value == "-":  # empty value span all include
                continue
            elif base == "-" and value == "f":  # empty value not span all include
                sequence_for_design += "N"
                primer_site += "f"
            elif base != "-" and value == "t":
                sequence_for_design += base
                primer_site += "t"
            elif base != "-" and value == "f":
                sequence_for_design += base
                primer_site += "f"

        template_organism = None
        cursor = self.db.cursor.execute(f'SELECT ORGANISM_NAME FROM ASSEMBLY WHERE ID = "{template_parsed_description["assembly_id"]}"')
        for row in cursor:
            template_organism = str(row[0])
        assert template_organism is not None
        excluded_region = [[match.span()[0], match.span()[1] - match.span()[0]] for match in re.finditer("f+", primer_site)]
        print(f"Designing primer for {template.id}")
        try:
            if self.pick_probe is False:
                designed_primers = bindings.designPrimers(
                    seq_args={"SEQUENCE_ID": template.id, "SEQUENCE_TEMPLATE": sequence_for_design, "SEQUENCE_EXCLUDED_REGION": excluded_region}, global_args=self.PRIMER_GLOBAL_ARGS, misprime_lib=misprime_lib
                )
            else:
                designed_primers = bindings.designPrimers(
                    seq_args={"SEQUENCE_ID": template.id, "SEQUENCE_TEMPLATE": sequence_for_design, "SEQUENCE_EXCLUDED_REGION": excluded_region}, global_args=self.PRIMER_GLOBAL_ARGS, mishyb_lib=misprime_lib, misprime_lib=misprime_lib
                )
            if self.pick_probe is False:
                print(f"left: {designed_primers['PRIMER_LEFT_EXPLAIN']}, right: {designed_primers['PRIMER_RIGHT_EXPLAIN']}, pair: {designed_primers['PRIMER_PAIR_EXPLAIN']}")
            else:
                print(f"left: {designed_primers['PRIMER_LEFT_EXPLAIN']}, internal: {designed_primers['PRIMER_INTERNAL_EXPLAIN']}, right: {designed_primers['PRIMER_RIGHT_EXPLAIN']}, pair: {designed_primers['PRIMER_PAIR_EXPLAIN']}")
            answer = {}
            index = 0
            while f"PRIMER_PAIR_{index}_PENALTY" in designed_primers.keys():
                if self.pick_probe is False:
                    new_primer_set = {
                        "TEMPLATE_SEQUENCE": str(sequence),
                        "CONSISTENCY_MASK": str("".join(mask)),
                        "TEMPLATE_SEQUENCE_ACCESSION": str(template_parsed_description["chromosome"]),
                        "TEMPLATE_PRODUCT": str(template_parsed_description["product"]),
                        "TEMPLATE_ASSEMBLY_ID": str(template_parsed_description["assembly_id"]),
                        "TEMPLATE_GENE_ID": str(template_parsed_description["id"]),
                        "TEMPLATE_ORGANISM": template_organism,
                        "PRIMER_PAIR_PENALTY": str(designed_primers[f"PRIMER_PAIR_{index}_PENALTY"]),
                        "PRIMER_LEFT_PENALTY": str(designed_primers[f"PRIMER_LEFT_{index}_PENALTY"]),
                        "PRIMER_RIGHT_PENALTY": str(designed_primers[f"PRIMER_RIGHT_{index}_PENALTY"]),
                        "PRIMER_LEFT_SEQUENCE": str(designed_primers[f"PRIMER_LEFT_{index}_SEQUENCE"]),
                        "PRIMER_RIGHT_SEQUENCE": str(designed_primers[f"PRIMER_RIGHT_{index}_SEQUENCE"]),
                        "PRIMER_LEFT": designed_primers[f"PRIMER_LEFT_{index}"],
                        "PRIMER_RIGHT": designed_primers[f"PRIMER_RIGHT_{index}"],
                        "PRIMER_LEFT_TM": str(designed_primers[f"PRIMER_LEFT_{index}_TM"]),
                        "PRIMER_RIGHT_TM": str(designed_primers[f"PRIMER_RIGHT_{index}_TM"]),
                        "PRIMER_LEFT_GC_PERCENT": str(designed_primers[f"PRIMER_LEFT_{index}_GC_PERCENT"]),
                        "PRIMER_RIGHT_GC_PERCENT": str(designed_primers[f"PRIMER_RIGHT_{index}_GC_PERCENT"]),
                        "PRIMER_LEFT_SELF_ANY_TH": str(designed_primers[f"PRIMER_LEFT_{index}_SELF_ANY_TH"]),
                        "PRIMER_RIGHT_SELF_ANY_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_SELF_ANY_TH"]),
                        "PRIMER_LEFT_SELF_END_TH": str(designed_primers[f"PRIMER_LEFT_{index}_SELF_ANY_TH"]),
                        "PRIMER_RIGHT_SELF_END_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_SELF_ANY_TH"]),
                        "PRIMER_LEFT_HAIRPIN_TH": str(designed_primers[f"PRIMER_LEFT_{index}_SELF_ANY_TH"]),
                        "PRIMER_RIGHT_HAIRPIN_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_SELF_ANY_TH"]),
                        # "PRIMER_LEFT_LIBRARY_MISPRIMING": str(designed_primers[f"PRIMER_LEFT_{index}_LIBRARY_MISPRIMING"]),
                        # "PRIMER_RIGHT_LIBRARY_MISPRIMING": str(designed_primers[f"PRIMER_RIGHT_{index}_LIBRARY_MISPRIMING"]),
                        # "PRIMER_PAIR_LIBRARY_MISPRIMING": str(designed_primers[f"PRIMER_PAIR_{index}_LIBRARY_MISPRIMING"]),
                        "PRIMER_LEFT_END_STABILITY": str(designed_primers[f"PRIMER_LEFT_{index}_END_STABILITY"]),
                        "PRIMER_RIGHT_END_STABILITY": str(designed_primers[f"PRIMER_RIGHT_{index}_END_STABILITY"]),
                        "PRIMER_LEFT_TEMPLATE_MISPRIMING": str(designed_primers[f"PRIMER_LEFT_{index}_TEMPLATE_MISPRIMING"]),
                        "PRIMER_RIGHT_TEMPLATE_MISPRIMING": str(designed_primers[f"PRIMER_RIGHT_{index}_TEMPLATE_MISPRIMING"]),
                        "PRIMER_PAIR_TEMPLATE_MISPRIMING": str(designed_primers[f"PRIMER_PAIR_{index}_TEMPLATE_MISPRIMING"]),
                        "PRIMER_LEFT_TEMPLATE_MISPRIMING_TH": str(designed_primers[f"PRIMER_LEFT_{index}_TEMPLATE_MISPRIMING_TH"]),
                        "PRIMER_RIGHT_TEMPLATE_MISPRIMING_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_TEMPLATE_MISPRIMING_TH"]),
                        "PRIMER_PAIR_COMPL_ANY_TH": str(designed_primers[f"PRIMER_PAIR_{index}_COMPL_ANY_TH"]),
                        "PRIMER_PAIR_COMPL_END_TH": str(designed_primers[f"PRIMER_PAIR_{index}_COMPL_END_TH"]),
                        "PRIMER_PAIR_PRODUCT_SIZE": str(designed_primers[f"PRIMER_PAIR_{index}_PRODUCT_SIZE"]),
                    }
                else:
                    new_primer_set = {
                        "TEMPLATE_SEQUENCE": str(sequence),
                        "CONSISTENCY_MASK": str("".join(mask)),
                        "TEMPLATE_SEQUENCE_ACCESSION": str(template_parsed_description["chromosome"]),
                        "TEMPLATE_PRODUCT": str(template_parsed_description["product"]),
                        "TEMPLATE_ASSEMBLY_ID": str(template_parsed_description["assembly_id"]),
                        "TEMPLATE_GENE_ID": str(template_parsed_description["id"]),
                        "TEMPLATE_ORGANISM": template_organism,
                        "PRIMER_PAIR_PENALTY": str(designed_primers[f"PRIMER_PAIR_{index}_PENALTY"]),
                        "PRIMER_LEFT_PENALTY": str(designed_primers[f"PRIMER_LEFT_{index}_PENALTY"]),
                        "PRIMER_RIGHT_PENALTY": str(designed_primers[f"PRIMER_RIGHT_{index}_PENALTY"]),
                        "PRIMER_INTERNAL_PENALTY": str(designed_primers[f"PRIMER_INTERNAL_{index}_PENALTY"]),
                        "PRIMER_LEFT_SEQUENCE": str(designed_primers[f"PRIMER_LEFT_{index}_SEQUENCE"]),
                        "PRIMER_RIGHT_SEQUENCE": str(designed_primers[f"PRIMER_RIGHT_{index}_SEQUENCE"]),
                        "PRIMER_INTERNAL_SEQUENCE": str(designed_primers[f"PRIMER_INTERNAL_{index}_SEQUENCE"]),
                        "PRIMER_LEFT": designed_primers[f"PRIMER_LEFT_{index}"],
                        "PRIMER_RIGHT": designed_primers[f"PRIMER_RIGHT_{index}"],
                        "PRIMER_INTERNAL": designed_primers[f"PRIMER_INTERNAL_{index}"],
                        "PRIMER_LEFT_TM": str(designed_primers[f"PRIMER_LEFT_{index}_TM"]),
                        "PRIMER_RIGHT_TM": str(designed_primers[f"PRIMER_RIGHT_{index}_TM"]),
                        "PRIMER_INTERNAL_TM": str(designed_primers[f"PRIMER_INTERNAL_{index}_TM"]),
                        "PRIMER_LEFT_GC_PERCENT": str(designed_primers[f"PRIMER_LEFT_{index}_GC_PERCENT"]),
                        "PRIMER_RIGHT_GC_PERCENT": str(designed_primers[f"PRIMER_RIGHT_{index}_GC_PERCENT"]),
                        "PRIMER_INTERNAL_GC_PERCENT": str(designed_primers[f"PRIMER_INTERNAL_{index}_GC_PERCENT"]),
                        "PRIMER_LEFT_SELF_ANY_TH": str(designed_primers[f"PRIMER_LEFT_{index}_SELF_ANY_TH"]),
                        "PRIMER_RIGHT_SELF_ANY_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_SELF_ANY_TH"]),
                        "PRIMER_INTERNAL_SELF_ANY_TH": str(designed_primers[f"PRIMER_INTERNAL_{index}_SELF_ANY_TH"]),
                        "PRIMER_LEFT_SELF_END_TH": str(designed_primers[f"PRIMER_LEFT_{index}_SELF_ANY_TH"]),
                        "PRIMER_RIGHT_SELF_END_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_SELF_ANY_TH"]),
                        "PRIMER_INTERNAL_SELF_END_TH": str(designed_primers[f"PRIMER_INTERNAL_{index}_SELF_ANY_TH"]),
                        "PRIMER_LEFT_HAIRPIN_TH": str(designed_primers[f"PRIMER_LEFT_{index}_SELF_ANY_TH"]),
                        "PRIMER_RIGHT_HAIRPIN_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_SELF_ANY_TH"]),
                        "PRIMER_INTERNAL_HAIRPIN_TH": str(designed_primers[f"PRIMER_INTERNAL_{index}_SELF_ANY_TH"]),
                        # "PRIMER_LEFT_LIBRARY_MISPRIMING": str(designed_primers[f"PRIMER_LEFT_{index}_LIBRARY_MISPRIMING"]),
                        # "PRIMER_RIGHT_LIBRARY_MISPRIMING": str(designed_primers[f"PRIMER_RIGHT_{index}_LIBRARY_MISPRIMING"]),
                        # "PRIMER_PAIR_LIBRARY_MISPRIMING": str(designed_primers[f"PRIMER_PAIR_{index}_LIBRARY_MISPRIMING"]),
                        "PRIMER_LEFT_END_STABILITY": str(designed_primers[f"PRIMER_LEFT_{index}_END_STABILITY"]),
                        "PRIMER_RIGHT_END_STABILITY": str(designed_primers[f"PRIMER_RIGHT_{index}_END_STABILITY"]),
                        "PRIMER_LEFT_TEMPLATE_MISPRIMING": str(designed_primers[f"PRIMER_LEFT_{index}_TEMPLATE_MISPRIMING"]),
                        "PRIMER_RIGHT_TEMPLATE_MISPRIMING": str(designed_primers[f"PRIMER_RIGHT_{index}_TEMPLATE_MISPRIMING"]),
                        "PRIMER_PAIR_TEMPLATE_MISPRIMING": str(designed_primers[f"PRIMER_PAIR_{index}_TEMPLATE_MISPRIMING"]),
                        "PRIMER_LEFT_TEMPLATE_MISPRIMING_TH": str(designed_primers[f"PRIMER_LEFT_{index}_TEMPLATE_MISPRIMING_TH"]),
                        "PRIMER_RIGHT_TEMPLATE_MISPRIMING_TH": str(designed_primers[f"PRIMER_RIGHT_{index}_TEMPLATE_MISPRIMING_TH"]),
                        "PRIMER_PAIR_COMPL_ANY_TH": str(designed_primers[f"PRIMER_PAIR_{index}_COMPL_ANY_TH"]),
                        "PRIMER_PAIR_COMPL_END_TH": str(designed_primers[f"PRIMER_PAIR_{index}_COMPL_END_TH"]),
                        "PRIMER_PAIR_PRODUCT_SIZE": str(designed_primers[f"PRIMER_PAIR_{index}_PRODUCT_SIZE"]),
                    }

                answer[f"PRIMER_{index}_" + template.id] = new_primer_set
                index += 1
            return answer
        except OSError as e:
            print(e)
            return {}

    def checkPrimers(self, identifier, left_primer: str, right_primer: str, internal_primer=None):

        print(f"Checking the specificity of designed primer {identifier} in other organism")
        input_path = self.temp_specificity_input_file_template.format(identifier)
        output_path = self.temp_specificity_output_file_template.format(identifier)
        with open(input_path, "w") as f:
            f.write(f">{identifier}_LEFT\n")
            f.write(left_primer + "\n")
            f.write(f">{identifier}_RGHT\n")
            f.write(right_primer + "\n")
            if self.pick_probe is True:
                if internal_primer is None:
                    raise ValueError("Please provide an internal probe for checking when passing pick_probe = True")
                f.write(f">{identifier}_INTL\n")
                f.write(internal_primer + "\n")

        BLAST_COMMAND_FOR_PRIMER_CHECK = [
            "blastn",
            "-task",
            "blastn-short",
            "-word_size",
            "6",
            "-num_threads",
            str(self.workers),
            "-outfmt",
            "6 qacc saccver qstart qend sstart send qseq sseq sstrand evalue mismatch qlen",
            "-query",
            input_path,
            "-db",
            self.db.blast_db_prefix,
            "-out",
            output_path,
        ]

        print(f"Running BLAST for primer {identifier}")
        subprocess.run(BLAST_COMMAND_FOR_PRIMER_CHECK)
        print(f"Parsing BLAST results for primer {identifier}")
        primer_hits = {}
        probe_hits = {}
        # key: sequence accession, value: list of location ,strand, etc.
        with open(output_path, "r") as f:
            # process and filter blast result here
            for row in f.readlines():
                if row[0] == "#":
                    continue
                else:
                    keys = ["query_accession", "subject_accession", "qstart", "qend", "sstart", "send", "query_sequence", "subject_sequence", "strand", "evalue", "mismatch", "query_length"]
                    values = row.replace("\n", "").split("\t")
                    hit = dict(zip(keys, values))

                    real_qend = None
                    real_qstart = None
                    real_send = None
                    real_sstart = None
                    real_query_sequence = None
                    real_subject_sequence = None

                    if hit["query_accession"][-4:] == "LEFT":
                        real_query_sequence = left_primer
                    elif hit["query_accession"][-4:] == "RGHT":
                        real_query_sequence = right_primer
                    else:
                        real_query_sequence = internal_primer

                    if hit["query_sequence"].find("-") != -1 or hit["subject_sequence"].find("-") != -1:
                        print("Found Insertion in Primer")
                        # continue

                    if hit["strand"] == "plus":
                        real_sstart = int(hit["sstart"]) - int(hit["qstart"]) + 1
                        real_send = int(hit["send"]) + int(hit["query_length"]) - int(hit["qend"])
                        real_qstart = 1
                        real_qend = int(hit["query_length"])
                        real_subject_id = hit["subject_accession"].split("-")[0]
                        real_sequence_id = hit["subject_accession"].replace(real_subject_id + "-", "")
                        cursor = self.db.cursor.execute(f'SELECT PATH FROM FILE WHERE TYPE = "GENOMIC_NUCLEOTIDE_FASTA" AND ASSEMBLY = "{real_subject_id}"')
                        for row in cursor:
                            file = Fasta(row[0])
                            if real_sequence_id in file:
                                real_subject_sequence = str(file[real_sequence_id][real_sstart - 1 : real_send]).upper()

                    if hit["strand"] == "minus":
                        real_sstart = int(hit["sstart"]) + int(hit["qstart"]) - 1
                        real_send = int(hit["send"]) - int(hit["query_length"]) + int(hit["qend"])
                        real_qstart = 1
                        real_qend = int(hit["query_length"])
                        real_subject_id = hit["subject_accession"].split("-")[0]
                        real_sequence_id = hit["subject_accession"].replace(real_subject_id + "-", "")
                        cursor = self.db.cursor.execute(f'SELECT PATH FROM FILE WHERE TYPE = "GENOMIC_NUCLEOTIDE_FASTA" AND ASSEMBLY = "{real_subject_id}"')
                        for row in cursor:
                            file = Fasta(row[0])
                            if real_sequence_id in file:
                                real_subject_sequence = str(file[real_sequence_id][real_send - 1 : real_sstart].reverse.complement).upper()

                    if len(real_subject_sequence) != len(real_query_sequence):
                        continue

                    hit["sstart"] = real_sstart
                    hit["send"] = real_send
                    hit["qstart"] = real_qstart
                    hit["qend"] = real_qend
                    hit["subject_sequence"] = real_subject_sequence
                    hit["query_sequence"] = real_query_sequence

                    subject_accession = hit["subject_accession"]
                    penalty = 0
                    if hit["query_accession"][-4:] == "LEFT" or hit["query_accession"][-4:] == "RGHT":
                        for i in range(len(hit["query_sequence"])):
                            if hit["query_sequence"][i] == "-" or hit["subject_sequence"][i] == "-":
                                penalty += 5
                            elif hit["query_sequence"][i] != hit["subject_sequence"][i]:
                                if i >= len(hit["query_sequence"]) - 3:
                                    penalty += 3
                                else:
                                    penalty += 1
                        if penalty > self.MAX_PENALTY_FOR_UNINTENDED_PRIMER_BINDING:
                            continue
                        hit["penalty"] = penalty
                        if subject_accession not in primer_hits.keys():
                            primer_hits[subject_accession] = []
                        primer_hits[subject_accession].append(hit)
                    else:
                        for i in range(len(hit["query_sequence"])):
                            if hit["query_sequence"][i] == "-" or hit["subject_sequence"][i] == "-":
                                penalty += 7
                            elif hit["query_sequence"][i] != hit["subject_sequence"][i]:
                                if i >= len(hit["query_sequence"]) - 3:
                                    penalty += 3
                                else:
                                    penalty += 1
                        if penalty > self.MAX_PENALTY_FOR_UNINTENDED_PROBE_BINDING:
                            continue
                        hit["penalty"] = penalty
                        if subject_accession not in probe_hits.keys():
                            probe_hits[subject_accession] = []
                        probe_hits[subject_accession].append(hit)

        effective_pairs = []
        for accession, hits in primer_hits.items():
            if len(hits) < 2:
                continue

            for i in range(len(hits)):
                for j in range(i + 1, len(hits)):
                    hit1 = hits[i]
                    hit2 = hits[j]
                    if hit1["strand"] == hit2["strand"] or hit1["subject_accession"] != hit2["subject_accession"]:
                        continue
                    if hit1["strand"] == "plus":
                        forward = hit1
                        reverse = hit2
                    else:
                        forward = hit2
                        reverse = hit1

                    assert int(forward["sstart"]) < int(forward["send"])
                    assert int(reverse["send"]) < int(reverse["sstart"])
                    start = int(forward["sstart"])
                    end = int(reverse["sstart"])
                    if start >= end:
                        continue
                    amplicon_assembly_id = forward["subject_accession"].split("-")[0]
                    amplicon_organism = None
                    cursor = self.db.cursor.execute(f'SELECT ORGANISM_NAME FROM ASSEMBLY WHERE ID = "{amplicon_assembly_id}"')
                    for row in cursor:
                        amplicon_organism = row[0]
                    product_length = end - start + 1
                    if product_length < self.PRIMER_CHECK_FILTER_MIN_PRODUCT_LENGTH or product_length > self.PRIMER_CHECK_FILTER_MAX_PRODUCT_LENGTH:
                        continue

                    if self.pick_probe is False:
                        effective_pair = {
                            "LEFT_PRIMER_IDENTIFIER": forward["query_accession"],
                            "LEFT_PRIMER_MATCHED": forward["query_sequence"],
                            "LEFT_PRIMER_TEMPLATE_MATCHED": forward["subject_sequence"],
                            "LEFT_PRIMER_TEMPLATE_START": int(forward["sstart"]),
                            "LEFT_PRIMER_TEMPLATE_END": int(forward["send"]),
                            "LEFT_PRIMER_PENALTY": int(forward["penalty"]),
                            "LEFT_PRIMER_LENGTH": int(forward["qend"]) - int(forward["qstart"]) + 1,
                            "LEFT_PRIMER_TEMPLATE_LENGTH": int(forward["send"]) - int(forward["sstart"]) + 1,
                            "RIGHT_PRIMER_IDENTIFIER": reverse["query_accession"],
                            "RIGHT_PRIMER_MATCHED": reverse["query_sequence"],
                            "RIGHT_PRIMER_TEMPLATE_MATCHED": reverse["subject_sequence"],
                            "RIGHT_PRIMER_TEMPLATE_START": int(reverse["sstart"]),
                            "RIGHT_PRIMER_TEMPLATE_END": int(reverse["send"]),
                            "RIGHT_PRIMER_PENALTY": int(reverse["penalty"]),
                            "RIGHT_PRIMER_LENGTH": int(reverse["qend"]) - int(reverse["qstart"]) + 1,
                            "RIGHT_PRIMER_TEMPLATE_LENGTH": int(reverse["sstart"]) - int(reverse["send"]) + 1,
                            "AMPLICON_SEQUENCE_ACCESSION": forward["subject_accession"],
                            "AMPLICON_ASSEMBLY_ID": amplicon_assembly_id,
                            "AMPLICON_ORGANISM": amplicon_organism,
                            "AMPLICON_SIZE": int(reverse["sstart"]) - int(forward["sstart"]) + 1,
                        }
                        effective_pairs.append(effective_pair)

                    elif forward["subject_accession"] in probe_hits:
                        # enumerate through all probe
                        for probe in probe_hits[forward["subject_accession"]]:
                            if int(forward["sstart"]) > int(probe["sstart"]) or int(forward["sstart"]) > int(probe["send"]) or int(probe["send"]) > int(reverse["sstart"]) or int(probe["send"]) > int(reverse["send"]):
                                continue
                            effective_pair = {
                                "LEFT_PRIMER_IDENTIFIER": forward["query_accession"],
                                "LEFT_PRIMER_MATCHED": forward["query_sequence"],
                                "LEFT_PRIMER_TEMPLATE_MATCHED": forward["subject_sequence"],
                                "LEFT_PRIMER_TEMPLATE_START": int(forward["sstart"]),
                                "LEFT_PRIMER_TEMPLATE_END": int(forward["send"]),
                                "LEFT_PRIMER_PENALTY": int(forward["penalty"]),
                                "LEFT_PRIMER_LENGTH": int(forward["qend"]) - int(forward["qstart"]) + 1,
                                "LEFT_PRIMER_TEMPLATE_LENGTH": int(forward["send"]) - int(forward["sstart"]) + 1,
                                "RIGHT_PRIMER_IDENTIFIER": reverse["query_accession"],
                                "RIGHT_PRIMER_MATCHED": reverse["query_sequence"],
                                "RIGHT_PRIMER_TEMPLATE_MATCHED": reverse["subject_sequence"],
                                "RIGHT_PRIMER_TEMPLATE_START": int(reverse["sstart"]),
                                "RIGHT_PRIMER_TEMPLATE_END": int(reverse["send"]),
                                "RIGHT_PRIMER_PENALTY": int(reverse["penalty"]),
                                "RIGHT_PRIMER_LENGTH": int(reverse["qend"]) - int(reverse["qstart"]) + 1,
                                "RIGHT_PRIMER_TEMPLATE_LENGTH": int(reverse["sstart"]) - int(reverse["send"]) + 1,
                                "INTERNAL_PRIMER_IDENTIFIER": probe["query_accession"],
                                "INTERNAL_PRIMER_MATCHED": probe["query_sequence"],
                                "INTERNAL_PRIMER_TEMPLATE_MATCHED": probe["subject_sequence"],
                                "INTERNAL_PRIMER_TEMPLATE_START": int(probe["sstart"]),
                                "INTERNAL_PRIMER_TEMPLATE_END": int(probe["send"]),
                                "INTERNAL_PRIMER_PENALTY": int(probe["penalty"]),
                                "INTERNAL_PRIMER_LENGTH": int(probe["qend"]) - int(probe["qstart"]) + 1,
                                "INTERNAL_PRIMER_TEMPLATE_LENGTH": int(probe["sstart"]) - int(probe["send"]) + 1,
                                "AMPLICON_SEQUENCE_ACCESSION": forward["subject_accession"],
                                "AMPLICON_ASSEMBLY_ID": amplicon_assembly_id,
                                "AMPLICON_ORGANISM": amplicon_organism,
                                "AMPLICON_SIZE": int(reverse["sstart"]) - int(forward["sstart"]) + 1,
                            }
                            effective_pairs.append(effective_pair)

        include_covered = []
        exclude_not_excluded = []
        for effective_pair in effective_pairs:
            if effective_pair["AMPLICON_ASSEMBLY_ID"] not in self.include_assembly_id_list:
                continue
            if effective_pair["LEFT_PRIMER_PENALTY"] > self.MAX_PENALTY_FOR_INTENDED_PRIMER_BINDING:
                continue
            if effective_pair["RIGHT_PRIMER_PENALTY"] > self.MAX_PENALTY_FOR_INTENDED_PRIMER_BINDING:
                continue
            if self.pick_probe is True:
                if effective_pair["INTERNAL_PRIMER_PENALTY"] > self.MAX_PENALTY_FOR_INTENDED_PROBE_BINDING:
                    continue
            if effective_pair["AMPLICON_ASSEMBLY_ID"] not in include_covered:
                include_covered.append(effective_pair["AMPLICON_ASSEMBLY_ID"])

        for effective_pair in effective_pairs:
            if effective_pair["AMPLICON_ASSEMBLY_ID"] not in self.exclude_assembly_id_list:
                continue
            if effective_pair["LEFT_PRIMER_PENALTY"] > self.MAX_PENALTY_FOR_UNINTENDED_PRIMER_BINDING:
                continue
            if effective_pair["RIGHT_PRIMER_PENALTY"] > self.MAX_PENALTY_FOR_UNINTENDED_PRIMER_BINDING:
                continue
            if self.pick_probe is True:
                if effective_pair["INTERNAL_PRIMER_PENALTY"] > self.MAX_PENALTY_FOR_UNINTENDED_PROBE_BINDING:
                    continue
            if effective_pair["AMPLICON_ASSEMBLY_ID"] not in include_covered:
                exclude_not_excluded.append(effective_pair["AMPLICON_ASSEMBLY_ID"])

        if len(exclude_not_excluded) > 0:
            print(f"Specificity check failed for {identifier}")
            return False, effective_pairs
        for include in self.include_assembly_id_list:
            if include not in include_covered:
                print(f"Specificity check failed for {identifier}")
                return False, effective_pairs
        print(f"Specificity check passed for {identifier}")
        return True, effective_pairs

    def makePrimerReport(self, primer_id, status, amplicons, primer):
        output_path = self.primer_report_file_template.format(primer_id)
        with open(output_path, "w") as f:
            f.truncate()
            f.write(f"# Report for {primer_id}\n")
            f.write("## Basic Primer Information\n")
            f.write(f"- status: {'PASSED_SPECIFICITY_CHECK' if status else 'SPECIFICITY_CHECK_FAILED'}\n")
            f.write(f"- probe picked: {self.pick_probe}\n")
            primer_columns = ["type", "penalty", "sequence", "GC%", "TM", "self complementary", "end stability"]
            if self.pick_probe is False:
                primer_table = [
                    primer_columns,
                    [
                        "forward",
                        primer["PRIMER_LEFT_PENALTY"],
                        primer["PRIMER_LEFT_SEQUENCE"],
                        primer["PRIMER_LEFT_GC_PERCENT"],
                        primer["PRIMER_LEFT_TM"],
                        primer["PRIMER_LEFT_SELF_ANY_TH"],
                        primer["PRIMER_LEFT_END_STABILITY"],
                    ],
                    [
                        "reverse",
                        primer["PRIMER_RIGHT_PENALTY"],
                        primer["PRIMER_RIGHT_SEQUENCE"],
                        primer["PRIMER_RIGHT_GC_PERCENT"],
                        primer["PRIMER_RIGHT_TM"],
                        primer["PRIMER_RIGHT_SELF_ANY_TH"],
                        primer["PRIMER_RIGHT_END_STABILITY"],
                    ],
                ]
            else:
                primer_table = [
                    primer_columns,
                    [
                        "forward",
                        primer["PRIMER_LEFT_PENALTY"],
                        primer["PRIMER_LEFT_SEQUENCE"],
                        primer["PRIMER_LEFT_GC_PERCENT"],
                        primer["PRIMER_LEFT_TM"],
                        primer["PRIMER_LEFT_SELF_ANY_TH"],
                        primer["PRIMER_LEFT_END_STABILITY"],
                    ],
                    ["internal", primer["PRIMER_INTERNAL_PENALTY"], primer["PRIMER_INTERNAL_SEQUENCE"], primer["PRIMER_INTERNAL_GC_PERCENT"], primer["PRIMER_INTERNAL_TM"], primer["PRIMER_INTERNAL_SELF_ANY_TH"], "N/A"],
                    [
                        "reverse",
                        primer["PRIMER_RIGHT_PENALTY"],
                        primer["PRIMER_RIGHT_SEQUENCE"],
                        primer["PRIMER_RIGHT_GC_PERCENT"],
                        primer["PRIMER_RIGHT_TM"],
                        primer["PRIMER_RIGHT_SELF_ANY_TH"],
                        primer["PRIMER_RIGHT_END_STABILITY"],
                    ],
                ]
            f.write(tabulate(primer_table, headers="firstrow", tablefmt="pipe"))
            f.write("\n\n")
            f.write("## Amplicons\n")
            for index, amplicon in enumerate(amplicons):
                f.write(f"Amplicon {index}:\n")
                f.write(f"- amplicon size: {amplicon['AMPLICON_SIZE']}\n")
                f.write(f"- organism: {amplicon['AMPLICON_ORGANISM']}\n")
                f.write(f"- assembly accession: {amplicon['AMPLICON_ASSEMBLY_ID']}\n")
                f.write(f"- sequence accession: {amplicon['AMPLICON_SEQUENCE_ACCESSION']}\n")
                amplicon_columns = ["type", "identifier", "start", "end", "matched sequence", "subject sequence", "penalty"]
                if self.pick_probe is False:
                    amplicon_table = [
                        amplicon_columns,
                        [
                            "forward",
                            amplicon["LEFT_PRIMER_IDENTIFIER"],
                            amplicon["LEFT_PRIMER_TEMPLATE_START"],
                            amplicon["LEFT_PRIMER_TEMPLATE_END"],
                            amplicon["LEFT_PRIMER_MATCHED"],
                            amplicon["LEFT_PRIMER_TEMPLATE_MATCHED"],
                            amplicon["LEFT_PRIMER_PENALTY"],
                        ],
                        [
                            "reverse",
                            amplicon["RIGHT_PRIMER_IDENTIFIER"],
                            amplicon["RIGHT_PRIMER_TEMPLATE_START"],
                            amplicon["RIGHT_PRIMER_TEMPLATE_END"],
                            amplicon["RIGHT_PRIMER_MATCHED"],
                            amplicon["LEFT_PRIMER_TEMPLATE_MATCHED"],
                            amplicon["RIGHT_PRIMER_PENALTY"],
                        ],
                    ]
                else:
                    amplicon_table = [
                        amplicon_columns,
                        [
                            "forward",
                            amplicon["LEFT_PRIMER_IDENTIFIER"],
                            amplicon["LEFT_PRIMER_TEMPLATE_START"],
                            amplicon["LEFT_PRIMER_TEMPLATE_END"],
                            amplicon["LEFT_PRIMER_MATCHED"],
                            amplicon["LEFT_PRIMER_TEMPLATE_MATCHED"],
                            amplicon["LEFT_PRIMER_PENALTY"],
                        ],
                        [
                            "internal",
                            amplicon["INTERNAL_PRIMER_IDENTIFIER"],
                            amplicon["INTERNAL_PRIMER_TEMPLATE_START"],
                            amplicon["INTERNAL_PRIMER_TEMPLATE_END"],
                            amplicon["INTERNAL_PRIMER_MATCHED"],
                            amplicon["INTERNAL_PRIMER_TEMPLATE_MATCHED"],
                            amplicon["INTERNAL_PRIMER_PENALTY"],
                        ],
                        [
                            "reverse",
                            amplicon["RIGHT_PRIMER_IDENTIFIER"],
                            amplicon["RIGHT_PRIMER_TEMPLATE_START"],
                            amplicon["RIGHT_PRIMER_TEMPLATE_END"],
                            amplicon["RIGHT_PRIMER_MATCHED"],
                            amplicon["RIGHT_PRIMER_TEMPLATE_MATCHED"],
                            amplicon["RIGHT_PRIMER_PENALTY"],
                        ],
                    ]
                f.write(tabulate(amplicon_table, headers="firstrow", tablefmt="pipe"))
                f.write("\n\n")
            f.write("\n")
            f.write("## Specificity\n")
            f.write(f"This set of primer ({primer_id}) is specific to:\n")
            include_covered = list(set([effective_pair["AMPLICON_ASSEMBLY_ID"] for effective_pair in amplicons if effective_pair["AMPLICON_ASSEMBLY_ID"] in self.include_assembly_id_list]))
            exclude_covered = list(set([effective_pair["AMPLICON_ASSEMBLY_ID"] for effective_pair in amplicons if effective_pair["AMPLICON_ASSEMBLY_ID"] in self.exclude_assembly_id_list]))
            other_covered = list(set([effective_pair["AMPLICON_ASSEMBLY_ID"] for effective_pair in amplicons if effective_pair["AMPLICON_ASSEMBLY_ID"] not in list(self.exclude_assembly_id_list) + list(self.include_assembly_id_list)]))
            cursor = self.db.cursor.execute("SELECT ID,ORGANISM_NAME FROM ASSEMBLY")
            names = {}
            for row in cursor:
                names[row[0]] = row[1]
            for accession in include_covered:
                f.write(f"include: {accession}({names[accession]})\n")
            for accession in other_covered:
                f.write(f"other: {accession}({names[accession]})\n")
            for accession in exclude_covered:
                f.write(f"exclude: {accession}({names[accession]})\n")

            f.write("\n")

    def makeHomologousGroupReport(self, template_id, aligned_group, mask, status) -> None:
        row_length = 240
        print(f"Writing template report for template {template_id}")
        output_path = self.homologous_group_report_file_template.format(template_id)
        if not os.path.exists(os.path.dirname(output_path)):
            os.makedirs(os.path.dirname(output_path))
        primer_found = 0
        for value in status.values():
            if value:
                primer_found += 1

        template_assembly_id = template_id.replace("TEMPLATE", "").split("-")[0]
        template_gene_id = template_id.replace("TEMPLATE", "").replace(template_assembly_id + "_", "")
        template_assembly_info = None
        cursor = self.db.cursor.execute(f'SELECT LEVEL,SUBMISSION_DATE,ORGANISM_NAME FROM ASSEMBLY WHERE ID = "{template_assembly_id}"')
        for row in cursor:
            template_assembly_info = {"level": row[0], "submission_date": row[1], "organism_name": row[2]}
        names = {}
        cursor = self.db.cursor.execute("SELECT ID, ORGANISM_NAME FROM ASSEMBLY")
        for row in cursor:
            names[row[0]] = row[1]
        assert template_assembly_info is not None

        template_gene_info = None
        for record in aligned_group:
            if str(record.id)[0:8] == "TEMPLATE":
                template_gene_info = parseNCBIFastaDescription(record.description)
        assert template_gene_info is not None

        with open(output_path, "w") as f:
            f.truncate()
            f.write(f"# Template report for homologous group of {template_id}\n")
            f.write("-" * row_length + "\n\n")
            f.write("## Basic Template Information\n")
            f.write(f"- sequence accession: {template_gene_info['chromosome']}\n")  # type: ignore
            f.write(f"- start: {template_gene_info['start']}\n")  # type: ignore
            f.write(f"- end: {template_gene_info['end']}\n")  # type: ignore
            f.write(f"- strand: {'plus' if template_gene_info['strand'] == '+' else 'minus'}\n")  # type: ignore
            f.write(f"- gene id: {template_gene_id}\n")  # type: ignore
            f.write(f"- organism: {template_assembly_info['organism_name']}\n")  # type: ignore
            f.write(f"- genome assembly accession: {template_assembly_id}\n")  # type: ignore
            f.write(f"- gene biotype: {template_gene_info['gene_biotype']}\n")  # type: ignore
            f.write(f"- protein product: {template_gene_info['product']}\n")  # type: ignore
            f.write(f"- number of primers found: {primer_found}\n")
            f.write("\n\n")
            f.write("## Alignment & Mask\n")
            f.write('"[assembly accession]": sequence of template assembly\n')
            f.write('"{assembly accession}": sequence of other include assembly\n')
            f.write('" assembly accession ": sequence of other assembly\n')
            f.write('"|assembly accession|": sequence of excluded assembly\n')
            f.write("\n")
            index = 0
            ljust = max(len("consistency mask"), len("Template"))
            include_sequences = {}
            exclude_sequences = {}
            other_sequences = {}
            template_sequence = ""
            include_species = {}
            exclude_species = {}
            template_species = ""
            other_species = {}

            for record in aligned_group:
                description = parseNCBIFastaDescription(record.description)
                assembly_id = description["assembly_id"]
                ljust = max(ljust, len(str(len(assembly_id))) + 2)
                if assembly_id == template_assembly_id:  # type: ignore
                    template_sequence = str(record.seq)
                    template_species = names[assembly_id]
                elif assembly_id in self.include_assembly_id_list:
                    include_sequences[f"{{{assembly_id}}}"] = str(record.seq)
                    include_species[f"{{{assembly_id}}}"] = names[assembly_id]
                elif assembly_id in self.exclude_assembly_id_list:
                    exclude_sequences[f"|{assembly_id}|"] = str(record.seq)
                    exclude_species[f"|{assembly_id}|"] = names[assembly_id]
                else:
                    other_sequences[f" {assembly_id} "] = str(record.seq)
                    other_species[f" {assembly_id} "] = names[assembly_id]

            ljust = ceil((ljust + 2) / 4) * 4
            length = len(str(aligned_group[0].seq))

            f.write(f"[{template_assembly_id}]".ljust(ljust, " "))  # type: ignore

            f.write(template_species)
            f.write("\n")

            for header, sequence in include_sequences.items():
                f.write(header.ljust(ljust, " "))
                f.write(include_species[header])
                f.write("\n")

            for header, sequence in other_sequences.items():
                f.write(header.ljust(ljust, " "))
                f.write(other_species[header])
                f.write("\n")

            for header, sequence in exclude_sequences.items():
                f.write(header.ljust(ljust, " "))
                f.write(exclude_species[header])
                f.write("\n")

            f.write("\n")

            while index < length:
                going_to_add = min(row_length - ljust, length - index)
                f.write(f"[{template_assembly_id}]".ljust(ljust, " "))  # type: ignore
                f.write(template_sequence[index : index + going_to_add])
                f.write("\n")
                for header, sequence in include_sequences.items():
                    f.write(header.ljust(ljust, " "))
                    f.write(sequence[index : index + going_to_add])
                    f.write("\n")
                for header, sequence in other_sequences.items():
                    f.write(header.ljust(ljust, " "))
                    f.write(sequence[index : index + going_to_add])
                    f.write("\n")
                for header, sequence in exclude_sequences.items():
                    f.write(header.ljust(ljust, " "))
                    f.write(sequence[index : index + going_to_add])
                    f.write("\n")
                f.write("Mask".ljust(ljust, " "))
                f.write("".join(mask[index : index + going_to_add]).replace("t", ".").replace("f", " "))
                f.write("\n")
                f.write("\n")
                index += going_to_add
        return


if __name__ == "__main__":
    Database.build(input_dir="/home/hewc/programming/python_project_linux/0-PRIMER_DESIGN_TOOLS/example_raw_dataset/", output_dir="/home/hewc/programming/python_project_linux/0-PRIMER_DESIGN_TOOLS/example_database/")
    db = Database(config_path="/home/hewc/programming/python_project_linux/0-PRIMER_DESIGN_TOOLS/example_database/config.json")
    find(
        db=db,
        include=[
            "Cryptococcus gattii",
        ],
        exclude=["Cryptococcus neoformans"],
        workers=12,
        pick_probe=True,
        reference_id="GCF_000185945.1",
    )

    # status, amplicons = pd.checkPrimers("G&N", left_primer="GTGGGAGCAAGTTCTCTATCCA", right_primer="AGTGTTGGAGGTGAACGTGG")
    # status, amplicons = pd.checkPrimers("G",left_primer="AGCTATCCGCAATCTTCCCG",right_primer="GATTTGTGAGGCCGTCTGGA")
    """status, amplicons = pd.checkPrimers("N", left_primer="ACACTTACAGTTGGCGGCT", right_primer="GGATAGCGTTGTTTGGCACG")
    table_columns = ["Assembly Accession", "Organism name", "Sequence Accession", "Left Primer Start", "Left Primer End", "Left Primer Mismatch", "Right Primer Start", "Right Primer End", "Right Primer Mismatch"]
    table = [table_columns]
    for amplicon in amplicons:
        table.append(
            [
                amplicon["AMPLICON_ASSEMBLY_ID"],
                amplicon["AMPLICON_ORGANISM"],
                amplicon["AMPLICON_SEQUENCE_ACCESSION"],
                amplicon["LEFT_PRIMER_TEMPLATE_START"],
                amplicon["LEFT_PRIMER_TEMPLATE_END"],
                amplicon["LEFT_PRIMER_MISMATCH"],
                amplicon["RIGHT_PRIMER_TEMPLATE_START"],
                amplicon["RIGHT_PRIMER_TEMPLATE_END"],
                amplicon["RIGHT_PRIMER_MISMATCH"],
            ]
        )
    print(tabulate(table, headers="firstrow", tablefmt="pipe"))
        find(
        db=db,
        include=["Cryptococcus gattii"],
        exclude=[
            "Candida albicans",
            "Candida tropicalis",
            "Candida parapsilosis",
            "Candida kruseii",
            "Candida glabrata",
            "Rhizopus microsporus",
            "Rhizopus oryzae",
            "Rhizopus arrhizus",
            "Rhizopus stolonifer",
            "Apophysomyces variabilis",
            "Cunninghamella bertholletiae",
            "Lichtheimia ramosa",
            "Lichtheimia corymbifera",
            "Saksenaea vasiformis",
            "Basidiobolus ranarum",
            "Conidiobolus coronatus",
            "Cryptococcus neoformans",
            "Mucor fuscus",
            "Mucor lanceolatus",
            "Mucor racemosus",
            "Mucor circinelloides",
            "Mucor endophyticus",
            "Aspergillus calidostus",
            "Aspergillus flavus",
            "Aspergillus fumigatus",
            "Aspergillus niger",
            "Aspergillus terreus",
            "Aspergillus nidulans",
        ],
    )

    find(
        db=db,
        include=[
            "Mycobacterium" "Mycobacterium avium",
            "Mycobacterium gordonae",
            "Mycobacterium kansasii",
            "Mycobacterium scrofulaceum",
            "Mycobacterium xenopi",
            "Mycobacterium tuberculosis",
            "Mycobacterium mungi",
            "Mycobacterium suricattae",
        ],
        exclude=[
            "Acinetobacter",
            "Alcaligenes",
            "Bordetella",
            "Borrelia",
            "Brucella",
            "Burkholderia",
            "Chlamydia",
            "Clostridium",
            "Enterobacter",
            "Enterococcus",
            "Escherichia",
            "Haemophilus",
            "Klebsiella",
            "Legionella",
            "Leptospira",
            "Moraxella",
            "Morganella",
            "Mycoplasma",
            "Neisseria",
            "Nocardia",
            "Proteus",
            "Pseudomonas",
            "Serratia",
            "Staphylococcus",
            "Stenotrophomonas",
            "Streptococcus",
            "Treponema",
        ],
        workers=12,
    )

"""
