"""Test some of the modifications since biopython"""
import os.path as op

import pytest

from kmbio.PDB import FastMMCIFParser, MMCIFParser, PDBParser


@pytest.mark.parametrize('filename', ['1A8O.pdb', '2BEG.pdb', '2XHE.pdb'])
def test_pdb_structure_id(filename):
    file_path = op.join(op.abspath(op.dirname(__file__)), 'PDB', filename)
    parser = PDBParser()
    structure = parser.get_structure(file_path)
    assert structure.id == op.basename(op.splitext(filename)[0])


@pytest.mark.parametrize('filename', [
    '1A8O.cif', '1LCD.cif', '2OFG.cif', '3JQH.cif', '4CUP.cif', '4ZHL.cif'])
def test_mmcif_structure_id(filename):
    file_path = op.join(op.abspath(op.dirname(__file__)), 'PDB', filename)
    parser = MMCIFParser()
    structure = parser.get_structure(file_path)
    assert structure.id == op.basename(op.splitext(filename)[0])


@pytest.mark.parametrize('filename', [
    '1A8O.cif', '1LCD.cif', '2OFG.cif', '3JQH.cif', '4CUP.cif', '4ZHL.cif'])
def test_fastmmcif_structure_id(filename):
    file_path = op.join(op.abspath(op.dirname(__file__)), 'PDB', filename)
    parser = FastMMCIFParser()
    structure = parser.get_structure(file_path)
    assert structure.id == op.basename(op.splitext(filename)[0])
