from typing import Dict, List, Optional, Iterable
import subprocess
import os


class MSA:
    def __init__(self, input_fasta_path: str, output_fasta_path: str, workers: int) -> None:
        """
        non-pythonic dependency:

            mafft

        Encapsulation of multiple sequence aligner: mafft
        """
        self.input_path = input_fasta_path
        if not os.path.exists(self.input_path):
            raise NameError(f"File {input_fasta_path} not found")
        self.output_path = output_fasta_path
        if not os.path.exists(self.output_path):
            try:
                os.makedirs(output_fasta_path)
            except BaseException:
                raise ValueError(f"Path {output_fasta_path} is not valid")
        # use mafft to align the sequence
        os.system(f"mafft --thread {workers} --quiet {self.input_path} > {self.output_path}")

        description_mapping = {}
        with open(input_fasta_path, "r") as f:
            for line in f.readlines():
                if line[0] == ">":
                    header = line[1:].split(" ")[0]
                    description = line[1:].replace(header + " ", "")
                    description_mapping[header] = description

        output_lines = []
        with open(output_fasta_path, "r") as f:
            for line in f.readlines():
                if line[0] == ">":
                    header = line[1:].split(" ")[0]
                    description = description_mapping[header]
                    modified_line = ">" + header + " " + description
                    output_lines.append(modified_line)
                else:
                    output_lines.append(line)
        with open(output_fasta_path, "w") as f:
            f.truncate()
            for line in output_lines:
                f.write(line)


# not in use
class TaxonNode:
    def __init__(self, rank: str, taxonomy_id: int, sciname: str) -> None:
        """
        Internal class used by TaxonTree, store basic information about a TaxonNode

        """
        self.__rank = rank
        self.__taxonomy_id = taxonomy_id
        self.__sciname = sciname
        self.__father = None
        self.__childrens = []

    @property
    def rank(self) -> str:
        """
        The rank of current taxonNode, eg: species/genus/order/...
        """
        return self.__rank

    @property
    def taxonomy_id(self) -> int:
        """
        The taxonomy id of current taxonNode
        """
        return self.__taxonomy_id

    @property
    def sciname(self) -> str:
        """
        The scientific name of current taxonNode
        """
        return self.__sciname

    @property
    def father(self) -> Optional["TaxonNode"]:
        """
        The father node of current TxonNode
        """
        return self.__father

    @property
    def childrens(self) -> List["TaxonNode" or None]:
        """
        The children nodes of current TaxonNode in list
        """
        return self.__childrens

    def setFather(self, father: "TaxonNode") -> None:
        """
        Set the father of current Taxonode to another Taxonode

        Parameters
        ----------
        father: TaxonNode
         father node of current node

        Returns
        -------
        None

        """
        self.__father = father

    def addChildren(self, children: "TaxonNode" or Iterable["TaxonNode"], *args: List["TaxonNode"]) -> None:
        """
        Add child (or children) Taxonode to a Taxonode

        Parameters
        ----------
        TaxonNode or List of TaxonNode or multiple TaxonNode

        Returns
        -------
        None

        Raises
        ------
        ValueError
         when input type are not valid
        """
        try:
            children_list: List
            assert isinstance(children, List) or isinstance(children, TaxonNode)
            if isinstance(children, TaxonNode):
                children_list = [children]
                for extra_children in args:
                    assert isinstance(extra_children, TaxonNode)
                    children_list.append(extra_children)
            elif isinstance(children, Iterable):
                children_list = children
        except AssertionError:
            raise ValueError("Input type invalid")
        for children in children_list:
            self.__childrens += children_list


# not in use
class TaxonTree:
    def __init__(self, taxon_list: List[int] = [], root_node: TaxonNode = TaxonNode(rank="root", taxonomy_id=-1, sciname="root")) -> None:
        """
        Use taxonkit to build taxonomy tree, and than query subtree
        command: taxonkit lineage ...

        Parameters
        ----------
        taxon_list: List
         List of taxonomy id (int)

        root_node: TaxonNode
         the root node of the TaxonTree, default a placeholder

        Raises
        ------
        RuntimeError
         when unable to parse taxonkit output


        """
        self.__root_node = root_node
        self.__node: Dict[int, "TaxonNode"] = {-1: self.__root_node}
        try:
            assert isinstance(taxon_list, List)
            for taxon_id in taxon_list:
                assert isinstance(taxon_id, int)

        except BaseException:
            raise ValueError("Input type not valid")

        if len(taxon_list) > 0:
            progress = subprocess.run(["taxonkit", "lineage", "-R", "-t"], encoding="utf-8", input="\n".join([str(taxon_id) for taxon_id in taxon_list]), capture_output=True)
            if progress.stderr != "":
                raise RuntimeError(progress.stderr)
            for line in progress.stdout.split("\n")[:-1]:
                splited_line = line.split("\t")
                try:
                    assert len(splited_line) == 4
                    id = int(splited_line[0])
                    lineage_name_list = splited_line[1].split(";")
                    lineage_taxon_list = splited_line[2].split(";")
                    lineage_rank_list = splited_line[3].split(";")
                    assert len(lineage_name_list) == len(lineage_rank_list) == len(lineage_taxon_list)
                    previous_taxon_id = self.__root_node.taxonomy_id
                    for name, rank, id in zip(lineage_name_list, lineage_rank_list, lineage_taxon_list):
                        assert previous_taxon_id in self.__node.keys()
                        if int(id) not in self.__node.keys():
                            curnode = TaxonNode(str(rank), int(id), str(name))
                            previous = self.__node[previous_taxon_id]
                            curnode.setFather(previous)
                            previous.addChildren(curnode)
                            self.__node[int(id)] = curnode
                        previous_taxon_id = int(id)

                except AssertionError:
                    raise RuntimeError("Can not parse line {line}")

    def __queryNewickRecursive(self, cur_node: "TaxonNode", taxon_id_list: List):
        """
        Internal method for querying Newick style tree string recursively

        Parameters
        ----------
        cur_node:
         current node in recursion

        taxon_id_list:
         the taxonomy id list for subtree

        Returns
        -------
        partial answer begins with current node

        """
        children_answer = []
        if cur_node.taxonomy_id in taxon_id_list:
            children_answer.append(str(f"({cur_node.taxonomy_id})"))
        for child_node in cur_node.childrens:
            child_answer = self.__queryNewickRecursive(child_node, taxon_id_list)
            if child_node.taxonomy_id in taxon_id_list or child_answer != "":
                children_answer.append(child_answer)
        if len(children_answer) == 1:
            return children_answer[0]
        if len(children_answer) > 1:
            return f'({", ".join(children_answer)})'
        else:
            return ""

    def toNewick(self, sub_tree_taxon_id_list: List):
        """
        Method for querying newick style string of the subtree of a TaxonTree

        Parameters
        ----------
        taxon_id_list:
         the taxonomy id list for subtree

        Returns
        -------
        answer, the string

        """

        for id in sub_tree_taxon_id_list:
            if not isinstance(id, int):
                raise ValueError("Input must be integers")
            if id not in self.__node.keys():
                raise KeyError(f"Node {id} of subtree not found in the whole tree")

        return self.__queryNewickRecursive(self.__root_node, sub_tree_taxon_id_list)


# not in use
class Phast:
    """
    Class for calling phastCons, phyloFit and parsing output, phastCons is a tool for producing conservation score using phylo-HMM model
    """

    def __init__(self, alignment_file_path, output_path: str, model_path: str, newick_tree: str) -> None:
        try:
            assert isinstance(alignment_file_path, str)
            assert os.path.exists(alignment_file_path)
            assert isinstance(output_path, str)
            assert isinstance(model_path, str)
            assert isinstance(newick_tree, str)
            if not os.path.exists(os.path.abspath(os.path.dirname(output_path))):
                os.makedirs(os.path.abspath(os.path.dirname(output_path)))
            if not os.path.exists(os.path.abspath(os.path.dirname(model_path))):
                os.makedirs(os.path.abspath(os.path.dirname(model_path)))
        except AssertionError:
            raise ValueError("Input of Phast invalid, please check again")
        except BaseException:
            raise RuntimeError("Runtime Error when making dirs for Phast input, please check input")
        self.input_path = alignment_file_path
        self.model_path = model_path
        self.output_path = output_path
        self.newick = newick_tree
        PHYLOFIT_COMMAND = ["phyloFit", "--precision", "LOW", "--tree", newick_tree, "--out-root", model_path, alignment_file_path]
        PHASTCON_COMMAND = ["phastCons", "--target-coverage", "0.7", "--expected-length", "25", "--refidx", "0", alignment_file_path, model_path + ".mod"]

        try:
            print("Generating HMM model", end=" ")
            phylofit = subprocess.run(PHYLOFIT_COMMAND, encoding="utf-8", capture_output=True)
            print("Done.")
            print("Calculating conservativity", end=" ")
            with open(output_path, "w") as f:
                phastcon = subprocess.run(PHASTCON_COMMAND, encoding="utf-8", capture_output=True)
                f.write(phastcon.stdout)
            print("Done.")
            assert phylofit.returncode == 0 and phastcon.returncode == 0
        except AssertionError:
            raise RuntimeError("Runtime Error occur when running phastcon, which should not happen")
        self.letter_annotations = []
        with open(output_path) as f:
            for line in f.readlines()[1:]:
                self.letter_annotations.append(float(line.replace("\n", "")))
