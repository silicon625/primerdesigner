import sqlite3
import subprocess
import os
import json
import jsonlines
from utils import jsonString, parseGFFAttributes
from urllib.parse import unquote
import platform


class Database:
    """
    Utlities associated with database
    Only works on linux

    non-pythonic dependencies:
        ncbi-blast


    APIs:

        1. __init__(): load pre-built database using path to config file
        2. Database.build(input_dir, output_dir): gather information and build database from a merged ncbi dataset
        3. Database.merge(input_dir,output_dir): find all ncbi datasets from input folder and then merge them to one dataset


    A dataset is a group of files downloaded from ncbi dataset website: https://www.ncbi.nlm.nih.gov/data-hub/genome/
    Users can also create their own dataset using genomic sequence fasta and gff annotation, but they must let the dataset have the same structure with the official one
    A dataset must have the following structure:

    root_dir
    ├── assembly_data_report.jsonl
    ├── dataset_catalog.json
    ├── <first assembly folder named with assembly accession>
    │   ├── files
    │   └── ...
    ├── <second assembly folder named with assembly accession>
    │   ├── files
    │   └── ...
    ├── <other assembly folders named with assembly accession>
    └──...

    dataset_catlog.json contains the mapping of all files in the folder, this file is used for detecting a dataset
    assembly_data_report.jsonl contains the meta data for all assemblies


    If users want to build database using multiple datsets, they should merge the dataset into one folder and modified the 2 files mentioned above
    The static method: Database.merge() is designed for this purpose
    This method will check all children folders in the input folder recursively and try to find all valid datasets, and then merge them into a single dataset


    The static method: Database.build() is used for building database from a valid dataset
    Database used for primer design consists of 2 parts:
        1. BLAST database (used for blast)
        2. sqlite database (store gff information and assembly details)

    User should perform queries using sql query script

    example:

        Database.build("input_dir_path","output_dir_path")
        db = Database("config_path")
        cursor = db.cursor.execute("SELECT ... FROM ...")
        for row in cursor:
            print(row[0])

    sql tables:
        1. ASSEMBLY: metadata for each assembly
            Fields:
                ID: assembly accession with version
                LEVEL: assembly level
                SUBMISSION_DATE: the date of submission
                ORGANISM_NAME: full name for the organism
                SCIENTIFIC_NAME: Genus species
                ORGANISM_ID: taxonomy id for the organism
                ORGANISM_GENUS: genus of the organism
                ORGANISM_SPECIES: species of the organism

                # fields below are values used for assess the quality of the assembly
                CONTIG_L50
                CONTIG_N50
                SCAFFOLD_L50
                SCAFFOLD_N50
                TOTAL_SEQUENCE_LENGTH
                TOTAL_UNGAPPED_LENGTH

        2. GENE: all annotated genes and pseudo_genes
            Fields:
                ID: gene id
                GENE_BIOTYPE: gene or pseudo gene
                PRODUCT: the function of the coded protein
                ASSEMBLY: the parent assembly of the annotated feature
                CHROMOSOME: the accession of the full sequence
                START
                END
                STRAND: + or -

        3. FILE: file mapping for all files in the original database
            Fields:
                PATH: absolute file path
                ASSEMBLY: assembly accession
                TYPE: file type

    """

    # Macros for default storage path
    MERGED_GENOME_FILE = "__temp__merged.fna"
    SQLITE_DB_FILE = "__database.sqlite"
    BLAST_DB_PREFIX = "__blastdb"
    DATABASE_CONFIG_FILE = "config.json"

    @staticmethod
    def checkDependency():
        sys = platform.system()
        if sys != "Linux":
            raise RuntimeError(f"This program is not supported on this platform: {sys}")
        dependency_valid = 1
        try:
            makeblastdb = subprocess.run(["which", "makeblastdb"], capture_output=True, encoding="utf-8")
            assert makeblastdb.stdout != "" and makeblastdb.stderr == ""
            print(f"Found makeblastdb at {makeblastdb.stdout[:-1]}")
        except BaseException:
            print("BLAST is not installed, this program is required for homologous group identification, please check their official website: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
            dependency_valid = 0

        if not dependency_valid:
            raise RuntimeError("Dependency invalid")

    def __init__(self, config_path) -> None:
        """
        read pre-built database from database_config
        """

        self.checkDependency()
        # check if input path is valid
        if not os.path.isfile(config_path):
            raise FileNotFoundError(f"Can not found database config at given path: {config_path}")
        try:
            with open(config_path, "r") as config_file:
                config = json.load(config_file)
        except json.JSONDecodeError:
            raise RuntimeError(f"Can not parse config at given path: {config_path}")

        try:
            assert "db_version" in config.keys()
        except AssertionError:
            raise ValueError(f"The format of the config can't be recognized at: {config_path}")

        # connect database
        db_version = config["db_version"]
        if db_version == "1.2":
            self.version = 1.2
            self.sqlite_db_path = config["sqlite_db_path"]
            self.blast_db_prefix = config["blast_db_prefix"]
            self.merged_fasta = config["merged_fasta"]

            try:
                self.conn = sqlite3.connect(self.sqlite_db_path)
                self.cursor = self.conn.cursor()
                print("Database loaded successfully")
                print("Version: 1.2")
                print(f"BLAST database: {self.blast_db_prefix}")
                print(f"SQLite database: {self.sqlite_db_path}")

                total_assembly_number = self.cursor.execute("SELECT COUNT(*) FROM ASSEMBLY").fetchone()[0]
                total_gene_number = self.cursor.execute("SELECT COUNT(*) FROM GENE").fetchone()[0]
                print(f"Total assembly number: {total_assembly_number}")
                print(f"Total annotated gene: {total_gene_number}")
            except ConnectionError:
                raise ConnectionError(f"Failed to connect database at path: {self.sqlite_db_path}")
        else:
            raise ValueError(f"This version of database is not supported currently: {db_version}, please update")

    # Not in use temporarily
    def assemblyData(self, assembly_id):
        """
        Query assembly meta data by assembly_id (accession)
        """
        if self.version == 1.2:
            return self._queryAssemblyData_V1d2(assembly_id=assembly_id)

    # Not in use temporarily
    def _queryAssemblyData_V1d2(self, assembly_id):
        """
        Internal method for query assembly meta data by assembly_id (accession)
        version 1.2
        """
        result = self.cursor.execute(
            f'SELECT ID,LEVEL,SUBMISSION_DATE,ORGANISM_NAME,SCIENTIFIC_NAME,ORGANISM_ID,ORGANISM_GENUS,ORGANISM_SPECIES,CONTIG_L50,CONTIG_N50,SCAFFOLD_L50,SCAFFOLD_N50,TOTAL_SEQUENCE_LENGTH,TOTAL_UNGAPPED_LENGTH from ASSEMBLY where ID="{assembly_id}"'
        )
        try:
            cnt = 0
            for row in result:
                cnt += 1
            assert cnt == 1
        except AssertionError:
            raise ValueError("Found duplicated accession")
        return dict(
            zip(
                [
                    "id",
                    "level",
                    "submission_date",
                    "organism_name",
                    "organism_id",
                    "organism_genus",
                    "organism_species",
                    "contig_l50",
                    "contig_n50",
                    "scaffold_l50",
                    "scaffold_n50",
                    "total_sequence_length",
                    "total_ungapped_length",
                ],
                result,
            )
        )

    @staticmethod
    def _metaDataEntry(entry, cursor):
        """
        Internal method used for putting the meta data of an assembly to database
        """
        id = f'"{entry["assemblyInfo"]["assemblyAccession"]}"'
        level = f'"{entry["assemblyInfo"]["assemblyLevel"]}"'
        submission_date = f'"{entry["assemblyInfo"]["submissionDate"]}"'
        organism_name = f'"{entry["organismName"]}"'
        organism_id = f'"{entry["taxId"]}"'
        organism_genus = f'"{entry["organismName"].split(" ")[0]}"'
        scientific_name = f'"{entry["organismName"].split(" ")[0]+" "+entry["organismName"].split(" ")[1]}"' if entry["organismName"].split(" ")[1] != "sp." or entry["organismName"].split(" ")[1] != "spp." else "NULL"
        organism_species = f'"{entry["organismName"].split(" ")[1]}"' if entry["organismName"].split(" ")[1] != "sp." or entry["organismName"].split(" ")[1] != "spp." else "NULL"
        contig_l50 = f'"{entry["assemblyStats"]["contigL50"]}"'
        contig_n50 = f'"{entry["assemblyStats"]["contigN50"]}"'
        scaffold_l50 = f'"{entry["assemblyStats"]["scaffoldL50"]}"'
        scaffold_n50 = f'"{entry["assemblyStats"]["scaffoldN50"]}"'
        total_sequence_length = f'"{entry["assemblyStats"]["totalSequenceLength"]}"'
        total_ungapped_length = f'"{entry["assemblyStats"]["totalUngappedLength"]}"'
        command = f"INSERT INTO ASSEMBLY (ID,LEVEL,SUBMISSION_DATE,ORGANISM_NAME,SCIENTIFIC_NAME,ORGANISM_ID,ORGANISM_GENUS,ORGANISM_SPECIES,CONTIG_L50,CONTIG_N50,SCAFFOLD_L50,SCAFFOLD_N50,TOTAL_SEQUENCE_LENGTH,TOTAL_UNGAPPED_LENGTH) \
                VALUES ({id},{level},{submission_date},{organism_name},{scientific_name},{organism_id},{organism_genus},{organism_species},{contig_l50},{contig_n50},{scaffold_l50},{scaffold_n50},{total_sequence_length},{total_ungapped_length})"
        cursor.execute(command)

    @staticmethod
    def _GFFDataEntry(path, assembly_id, cursor):  # , file_mapping):
        """
        Internal method used for putting the meta data of an gff file to database
        """
        products = {}
        # genome = []
        # for file in file_mapping[assembly_id]:
        #    if file["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA":
        #        genome.append(Fasta(file["filePath"]))
        with open(path, "r") as gff:
            for line in gff.readlines():
                if line[0] == "#":
                    continue
                # parsed_line is a dict having following fields
                parsed_line = dict(zip(["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"], line.split("\t")))
                attributes = parseGFFAttributes(parsed_line["attributes"])
                # has parent means this line is not top level
                if "Parent" in attributes.keys():
                    if "product" in attributes.keys():
                        if attributes["Parent"] not in products.keys() or attributes["Parent"] == "" or attributes["Parent"] is None:
                            products[attributes["Parent"]] = unquote(attributes["product"])

        with open(path, "r") as gff:
            for line in gff.readlines():
                if line[0] == "#":
                    continue
                parsed_line = dict(zip(["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"], line.split("\t")))
                attributes = parseGFFAttributes(parsed_line["attributes"])
                if "Parent" not in attributes.keys() and parsed_line["type"] == "gene":
                    id = f'"{attributes["ID"]}"'
                    gene_biotype = f'"{attributes["gene_biotype"]}"'
                    product = f'"{products[attributes["ID"]]}"' if attributes["ID"] in products.keys() else "NULL"
                    assembly = f'"{assembly_id}"'
                    chromosome = f'"{parsed_line["seqid"]}"'
                    start = f'"{parsed_line["start"]}"'
                    end = f'"{parsed_line["end"]}"'
                    strand = f'"{parsed_line["strand"]}"'
                    # for file in genome:
                    #    if parsed_line["seqid"] not in file.keys():
                    #        continue
                    #    sequence = file[parsed_line["seqid"]][int(parsed_line["start"]) - 1 : int(parsed_line["end"])]
                    # if parsed_line["strand"] == "-":
                    #    sequence = sequence.reverse.complement
                    # sequence = f'"{str(sequence)}"'
                    # cursor.execute(
                    #    f"INSERT INTO GENE (ID,GENE_BIOTYPE,PRODUCT,ASSEMBLY,CHROMOSOME,START,END,STRAND,SEQUENCE) \
                    #    VALUES ({id},{gene_biotype},{product},{assembly},{chromosome},{start},{end},{strand},{sequence})"
                    # )
                    cursor.execute(
                        f"INSERT INTO GENE (ID,GENE_BIOTYPE,PRODUCT,ASSEMBLY,CHROMOSOME,START,END,STRAND) \
                        VALUES ({id},{gene_biotype},{product},{assembly},{chromosome},{start},{end},{strand})"
                    )

    @staticmethod
    def _fileDataEntry(file_mapping, cursor):
        """
        Internal method used for putting file mapping into database
        """
        for assembly_id, files in file_mapping.items():
            for file in files:
                path = f'"{file["filePath"]}"'
                file_type = f'"{file["fileType"]}"'
                cursor.execute(f'INSERT INTO FILE (PATH,TYPE,ASSEMBLY) VALUES ({path},{file_type},"{assembly_id}")')

    @staticmethod
    def build(input_dir, output_dir):
        """
        Build database using information from input_path
        """
        Database.checkDependency()
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        except BaseException:
            raise RuntimeError(f"Can not create database at given path: {output_dir}")
        merged_genome_file = os.path.abspath(os.path.join(output_dir, Database.MERGED_GENOME_FILE))
        sqlite_db_file = os.path.abspath(os.path.join(output_dir, Database.SQLITE_DB_FILE))
        blast_db_prefix = os.path.abspath(os.path.join(output_dir, Database.BLAST_DB_PREFIX))
        database_config_file = os.path.abspath(os.path.join(output_dir, Database.DATABASE_CONFIG_FILE))
        # Deal with existing database
        # let the user decide whether to delete the previous database
        # raise RuntimeError when database can't be removed
        if os.path.exists(merged_genome_file) or os.path.exists(sqlite_db_file) or os.path.exists(blast_db_prefix + ".njs") or os.path.exists(database_config_file):
            print(f"Found existing database at: {output_dir}")
            c = input("Delete existing database? (Y/n): ")
            while c != "Y" and c != "n":
                print(f"Can't recognize input character: {c}")
                c = input("Delete existing database? (Y/n): ")
            if c == "Y":
                print("Deleting existing database ...", end=" ")
                subprocess.run(["rm", "-rf", merged_genome_file])
                subprocess.run(["rm", "-rf", sqlite_db_file])
                subprocess.run(["rm", "-rf", database_config_file])
                subprocess.run(["rm", "-rf", blast_db_prefix + "*"])
                print("Done")
            elif c == "n":
                print("Aborted")
                return

        # Initalize database
        try:
            conn = sqlite3.connect(sqlite_db_file)
            cursor = conn.cursor()
            cursor.execute(
                """CREATE TABLE ASSEMBLY
                (
                ID varchar(50),
                LEVEL varchar(20),
                SUBMISSION_DATE DATE,
                ORGANISM_NAME varchar(100),
                SCIENTIFIC_NAME varchar(100),
                ORGANISM_ID INT,
                ORGANISM_GENUS varchar(50),
                ORGANISM_SPECIES varchar(50),
                CONTIG_L50 INT,
                CONTIG_N50 INT,
                SCAFFOLD_L50 INT,
                SCAFFOLD_N50 INT,
                TOTAL_SEQUENCE_LENGTH INT,
                TOTAL_UNGAPPED_LENGTH INT
                )"""
            )
            cursor.execute(
                """CREATE TABLE FILE
                (
                PATH TEXT,
                ASSEMBLY varchar(50),
                TYPE varchar(50)
                )
                """
            )
            cursor.execute(
                """CREATE TABLE GENE
                (
                ID varchar(50),
                GENE_BIOTYPE varchar(50),
                PRODUCT varchar(50),
                ASSEMBLY varchar(50),
                CHROMOSOME varchar(50),
                START INT,
                END INT,
                STRAND varchar(2)
                )
                """
            )
            conn.commit()
        except ConnectionError:
            raise RuntimeError(f"Can not create database at given path: {output_dir}")

        # verify whether input_dir is valid
        # one file named dataset_catalog_path should exists in this folder, containing the naming convention
        try:
            assert os.path.isdir(input_dir)
            dataset_catalog_path = os.path.join(input_dir, "dataset_catalog.json")
            assert os.path.isfile(dataset_catalog_path)
        except AssertionError:
            raise RuntimeError(f"Input directory can not be recognized as parseable database at path: {input_dir}")

        # verify whether dataset_catalog is in right format
        # the file should have the naming convention for data report
        try:
            with open(dataset_catalog_path) as catalog_file:
                catalog = json.load(catalog_file)
                assert "assemblies" in catalog.keys()
                found_data_report = 0
                data_report_path = ""
                assembly_file_mapping = {}
                for assembly in catalog["assemblies"]:
                    assert "files" in assembly.keys()
                    if "accession" in assembly.keys():
                        if assembly["accession"] in assembly_file_mapping.keys():
                            raise KeyError(f"Found duplicated accession: {assembly['accession']} in catalog file")
                        assembly_file_mapping[assembly["accession"]] = assembly["files"]
                        for file in assembly_file_mapping[assembly["accession"]]:
                            file["filePath"] = os.path.join(input_dir, file["filePath"])
                    else:
                        for file in assembly["files"]:
                            assert "filePath" in file.keys()
                            assert "fileType" in file.keys()
                            if file["fileType"] == "DATA_REPORT":
                                found_data_report += 1
                                data_report_path = os.path.join(input_dir, file["filePath"])
                                assert os.path.exists(data_report_path)
                assert found_data_report == 1
        except AssertionError:
            raise RuntimeError(f"Can not parse dataset catalog at path: {dataset_catalog_path}")

        # meta data entry
        print("Prasing meta data ... ", end="")
        with open(data_report_path) as report:
            for item in jsonlines.Reader(report):
                # if the format of assembly_data.jsonl changes in the future, this function should be modified, parser for newer format should be added
                assembly_id = item["assemblyInfo"]["assemblyAccession"]
                Database._metaDataEntry(entry=item, cursor=cursor)
        conn.commit()

        Database._fileDataEntry(file_mapping=assembly_file_mapping, cursor=cursor)
        conn.commit()
        print("Done")

        # gene data entry
        print("Prasing gff files ... ")
        for assembly_id, files in assembly_file_mapping.items():
            for file in files:
                if file["fileType"] == "GFF3":
                    print(f"Dealing with {assembly_id} ... ", end="")
                    Database._GFFDataEntry(path=file["filePath"], assembly_id=assembly_id, cursor=cursor)  # , file_mapping=assembly_file_mapping)
                    print("Done")
        print("Prasing gff files Done")
        conn.commit()
        conn.close()

        print("Merging gemonic fasta ... ", end="")
        with open(merged_genome_file, "w") as f:
            f.truncate()
        with open(merged_genome_file, "a") as o:
            for assembly_id, files in assembly_file_mapping.items():
                for file in files:
                    if file["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA":
                        with open(file["filePath"], "r") as i:
                            for line in i.readlines():
                                if line[0] == ">":
                                    line = ">" + assembly_id + "-" + line[1:]
                                o.write(line)
        print("Done")

        subprocess.run(["makeblastdb", "-dbtype", "nucl", "-input_type", "fasta", "-parse_seqids", "-in", merged_genome_file, "-out", blast_db_prefix])

        # subprocess.run(["rm", "-f", merged_genome_file])

        with open(database_config_file, "w") as f:
            config = {}
            config["db_version"] = "1.2"
            config["sqlite_db_path"] = sqlite_db_file
            config["blast_db_prefix"] = blast_db_prefix
            config["merged_fasta"] = merged_genome_file
            f.write(jsonString(config))

        print("Database building finished successfully")
        print("\n\n")

    @staticmethod
    def merge(input_dir, output_dir):
        """
        find and merge all datasets in the input folder
        """

        # check input
        if isinstance(input_dir, str):
            input_dir = [input_dir]
        elif not isinstance(input_dir, list):
            raise TypeError("Input directory should either be a path or a list of path")

        try:
            for file in input_dir:
                assert os.path.isdir(file)
        except AssertionError:
            raise RuntimeError("Can open input directory")

        # check output
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

        input_dir = [os.path.abspath(file) for file in input_dir]
        output_dir = os.path.abspath(output_dir)

        # create file at output
        output_report_path = os.path.join(output_dir, "assembly_data_report.jsonl")
        output_catalog_path = os.path.join(output_dir, "dataset_catalog.json")
        output_catalog = {}

        with open(output_report_path, "w") as f:
            f.truncate()

        for file in input_dir:
            for root, ds, fs in os.walk(file):
                for f in fs:
                    if f == "dataset_catalog.json":
                        new_database_root_dir = os.path.abspath(root)
                        print(f"Found database at path: {new_database_root_dir}")
                        new_catalog_path = os.path.join(new_database_root_dir, "dataset_catalog.json")
                        new_report_path = None
                        with open(new_catalog_path, "r") as f:
                            new_catalog = json.load(f)
                            if len(output_catalog.keys()) == 0:
                                output_catalog = new_catalog
                            else:
                                new_catalog_assemblies = []
                                for elm in new_catalog["assemblies"]:
                                    if "accession" not in elm.keys():
                                        for file in elm["files"]:
                                            if file["fileType"] == "DATA_REPORT":
                                                new_report_path = os.path.join(new_database_root_dir, file["filePath"])
                                                break
                                        continue
                                    new_catalog_assemblies.append(elm)
                                    for file in elm["files"]:
                                        source = os.path.join(new_database_root_dir, file["filePath"])
                                        destination = os.path.join(output_dir, file["filePath"])
                                        destination_dir = os.path.dirname(destination)
                                        try:
                                            assert (not os.path.exists(destination_dir)) or os.path.isdir(destination_dir)
                                            if not os.path.exists(destination_dir):
                                                os.makedirs(destination_dir)
                                        except AssertionError:
                                            raise RuntimeError(f"Can not copy file to output directory, path: {destination}")
                                        cp_cmd = ["cp", source, destination]
                                        subprocess.run(cp_cmd)
                                output_catalog["assemblies"] += new_catalog_assemblies
                                assert new_report_path is not None
                                with open(output_report_path, "a") as o:
                                    with open(new_report_path, "r") as i:
                                        o.write(i.read())

        with open(output_catalog_path, "w") as f:
            f.truncate()
            f.write(jsonString(output_catalog))


if __name__ == "__main__":
    # Database.merge(input_dir="/home/hewc/Program/datasets", output_dir="/home/hewc/Program/merged")
    # Database.build("/home/hewc/programming/python_project_linux/0-PRIMER_DESIGN_TOOLS/common_bacteria_in_lung/", "/home/hewc/programming/python_project_linux/0-PRIMER_DESIGN_TOOLS/database/common_bacteria_in_lung/")
    # Database(config_path="/home/hewc/programming/python_project_linux/0-PRIMER_DESIGN_TOOLS/database/common_bacteria_in_lung/config.json")
    Database.merge(input_dir="/mnt/e/Download/genome update", output_dir="/mnt/e/Download/genome update/merged")