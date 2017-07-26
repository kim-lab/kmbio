# """Test some of the modifications since biopython"""
# import os.path as op
# import tempfile
# import urllib.request
#
# import pytest
#
# from kmbio.PDB import FastMMCIFParser, MMCIFParser, PDBParser
#
#
# @pytest.mark.parametrize('filename', ['1A8O.pdb', '2BEG.pdb', '2XHE.pdb'])
# def test_pdb_structure_id(filename):
#     """Test to make sure the PDBParser correctly parses the structure id"""
#     file_path = op.join(op.abspath(op.dirname(__file__)), 'PDB', filename)
#     parser = PDBParser()
#     structure = parser.get_structure(file_path)
#     assert structure.id == op.basename(op.splitext(filename)[0])
#
#
# @pytest.mark.parametrize('filename', [
#     '1A8O.cif', '1LCD.cif', '2OFG.cif', '3JQH.cif', '4CUP.cif', '4ZHL.cif'])
# def test_mmcif_structure_id(filename):
#     """Test to make sure the MMCIFParser correctly parses the structure id"""
#     file_path = op.join(op.abspath(op.dirname(__file__)), 'PDB', filename)
#     parser = MMCIFParser()
#     structure = parser.get_structure(file_path)
#     assert structure.id == op.basename(op.splitext(filename)[0])
#
#
# @pytest.mark.parametrize('filename', [
#     '1A8O.cif', '1LCD.cif', '2OFG.cif', '3JQH.cif', '4CUP.cif', '4ZHL.cif'])
# def test_fastmmcif_structure_id(filename):
#     """Test to make sure the FastMMCIFParser correctly parses the structure id"""
#     file_path = op.join(op.abspath(op.dirname(__file__)), 'PDB', filename)
#     parser = FastMMCIFParser()
#     structure = parser.get_structure(file_path)
#     assert structure.id == op.basename(op.splitext(filename)[0])
#
#
# @pytest.mark.parametrize('pdb_id, chain_ids_auth, chain_ids_label', [
#     ('1A8O', ('A'), ('A', 'B')),
#     ('1LCD', ('B', 'C', 'A'), ('A', 'B', 'C', 'D', 'E', 'F', 'G')),
#     ('2OFG', ('X'), ('A')),
#     ('3JQH', ('A'), ('A', 'B')),
#     ('4CUP', ('A'), ('A', 'B', 'C', 'D', 'E', 'F')),
#     ('4ZHL', ('U', 'P'), ('A', 'B', 'C', 'D')),
# ])
# def test_mmcif_with_label_ids(pdb_id, chain_ids_auth, chain_ids_label):
#     tmpdir = tempfile.mkdtemp()
#
#     pdb_file = op.join(tmpdir, pdb_id + '.pdb')
#     urllib.request.urlretrieve('http://files.rcsb.org/download/' + pdb_id + '.pdb', pdb_file)
#
#     cif_file = op.join(tmpdir, pdb_id + '.cif')
#     urllib.request.urlretrieve('http://files.rcsb.org/download/' + pdb_id + '.cif', cif_file)
#
#     # pdb_parser = PDBParser()
#     # pdb_structure = pdb_parser.get_structure(pdb_file)
#     # assert tuple(chain.id for chain in pdb_structure.get_chains()) == chain_ids_auth
#     #
#     # cif_auth_parser = MMCIFParser(use_label_ids=False)
#     # cif_auth_structure = cif_auth_parser.get_structure(cif_file)
#     # assert tuple(chain.id for chain in cif_auth_structure.get_chains()) == chain_ids_auth
#
#     cif_label_parser = MMCIFParser(use_label_ids=True)
#     cif_label_structure = cif_label_parser.get_structure(cif_file)
#     assert tuple(chain.id for chain in cif_label_structure.get_chains()) == chain_ids_label
#
#     tmpdir.clear()
