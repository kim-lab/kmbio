# Copyright 2012 Lenna X. Peterson (arklenna@gmail.com).
# All rights reserved.
#
# Tests adapted from test_PDB.py
#
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the MMCIF portion of the Bio.PDB module."""

import unittest

from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

from kmbio.PDB import PPBuilder, CaPPBuilder
from kmbio.PDB.MMCIFParser import MMCIFParser, FastMMCIFParser


class ParseReal(unittest.TestCase):
    """Testing with real CIF file(s)."""

    def test_parsers(self, ignore_authorId=False):
        """Extract polypeptides from 1A80."""

        parser = MMCIFParser()
        fast_parser = FastMMCIFParser()

        structure = parser.get_structure("PDB/1A8O.cif", "example")
        f_structure = fast_parser.get_structure("PDB/1A8O.cif", "example")


        self.assertEqual(len(structure), 1)
        self.assertEqual(len(f_structure), 1)

        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            # ==========================================================
            # Check that serial_num (model column) is stored properly
            self.assertEqual(structure[0].serial_num, 1)
            self.assertEqual(f_structure[0].serial_num, structure[0].serial_num)

            # First try allowing non-standard amino acids,
            polypeptides = ppbuild.build_peptides(structure[0], False)
            f_polypeptides = ppbuild.build_peptides(f_structure[0], False)

            self.assertEqual(len(polypeptides), 1)
            self.assertEqual(len(f_polypeptides), 1)

            pp = polypeptides[0]
            f_pp = f_polypeptides[0]

            # Check the start and end positions
            self.assertEqual(pp[0].id[1], 151)
            self.assertEqual(pp[-1].id[1], 220)

            self.assertEqual(f_pp[0].id[1], 151)
            self.assertEqual(f_pp[-1].id[1], 220)

            # Check the sequence
            s = pp.get_sequence()
            f_s = f_pp.get_sequence()

            self.assertEqual(s, f_s)  # enough to test this

            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)

            # Here non-standard MSE are shown as M
            self.assertEqual("MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQ"
                             "NANPDCKTILKALGPGATLEEMMTACQG", str(s))

            # ==========================================================
            # Now try strict version with only standard amino acids
            # Should ignore MSE 151 at start, and then break the chain
            # at MSE 185, and MSE 214,215
            polypeptides = ppbuild.build_peptides(structure[0], True)
            self.assertEqual(len(polypeptides), 3)

            # First fragment
            pp = polypeptides[0]
            self.assertEqual(pp[0].id[1], 152)
            self.assertEqual(pp[-1].id[1], 184)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("DIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNW", str(s))

            # Second fragment
            pp = polypeptides[1]
            self.assertEqual(pp[0].id[1], 186)
            self.assertEqual(pp[-1].id[1], 213)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("TETLLVQNANPDCKTILKALGPGATLEE", str(s))

            # Third fragment
            pp = polypeptides[2]
            self.assertEqual(pp[0].id[1], 216)
            self.assertEqual(pp[-1].id[1], 220)
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("TACQG", str(s))

    def testModels(self):
        """Test file with multiple models"""

        parser = MMCIFParser()
        f_parser = FastMMCIFParser()
        structure = parser.get_structure("PDB/1LCD.cif", "example")
        f_structure = f_parser.get_structure("PDB/1LCD.cif", "example")

        self.assertEqual(len(structure), 3)
        self.assertEqual(len(f_structure), 3)

        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            # ==========================================================
            # Check that serial_num (model column) is stored properly
            self.assertEqual(structure[0].serial_num, 1)
            self.assertEqual(structure[1].serial_num, 2)
            self.assertEqual(structure[2].serial_num, 3)
            # First try allowing non-standard amino acids,
            polypeptides = ppbuild.build_peptides(structure[0], False)
            self.assertEqual(len(polypeptides), 1)
            pp = polypeptides[0]
            # Check the start and end positions
            self.assertEqual(pp[0].id[1], 1)
            self.assertEqual(pp[-1].id[1], 51)
            # Check the sequence
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            # Here non-standard MSE are shown as M
            self.assertEqual("MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR",
                             str(s))
            # ==========================================================
            # Now try strict version with only standard amino acids
            polypeptides = ppbuild.build_peptides(structure[0], True)
            self.assertEqual(len(polypeptides), 1)
            pp = polypeptides[0]
            # Check the start and end positions
            self.assertEqual(pp[0].id[1], 1)
            self.assertEqual(pp[-1].id[1], 51)
            # Check the sequence
            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual("MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNR",
                             str(s))

        # This structure contains several models with multiple lengths.
        # The tests were failing.
        structure = parser.get_structure("PDB/2OFG.cif", "example")
        self.assertEqual(len(structure), 3)

    def test_insertions(self):
        """Test file with residue insertion codes"""
        parser = MMCIFParser()
        structure = parser.get_structure("PDB/4ZHL.cif", "example")
        for ppbuild in [PPBuilder(), CaPPBuilder()]:
            # First try allowing non-standard amino acids,
            polypeptides = ppbuild.build_peptides(structure[0], False)
            self.assertEqual(len(polypeptides), 2)
            pp = polypeptides[0]
            # Check the start and end positions (first segment only)
            self.assertEqual(pp[0].id[1], 16)
            self.assertEqual(pp[-1].id[1], 244)
            # Check the sequence
            refseq = "IIGGEFTTIENQPWFAAIYRRHRGGSVTYVCGGSLISPCWVISATHCFIDYPKKEDYIVYLGR" \
                     "SRLNSNTQGEMKFEVENLILHKDYSADTLAYHNDIALLKIRSKEGRCAQPSRTIQTIALPSMY" \
                     "NDPQFGTSCEITGFGKEQSTDYLYPEQLKMTVVKLISHRECQQPHYYGSEVTTKMLCAADPQW" \
                     "KTDSCQGDSGGPLVCSLQGRMTLTGIVSWGRGCALKDKPGVYTRVSHFLPWIRSHTKE"

            s = pp.get_sequence()
            self.assertTrue(isinstance(s, Seq))
            self.assertEqual(s.alphabet, generic_protein)
            self.assertEqual(refseq, str(s))

    def test_filehandle(self):
        """Test if the parser can handle file handle as well as filename"""
        parser = MMCIFParser()
        structure = parser.get_structure("PDB/1A8O.cif", "example")
        self.assertEqual(len(structure), 1)

        structure = parser.get_structure(open("PDB/1A8O.cif"), "example")
        self.assertEqual(len(structure), 1)

    def test_point_mutations_main(self):
        """Test if MMCIFParser parse point mutations correctly."""
        self._run_point_mutation_tests(MMCIFParser())

    def test_point_mutations_fast(self):
        """Test if FastMMCIFParser can parse point mutations correctly."""
        self._run_point_mutation_tests(FastMMCIFParser())

    def _run_point_mutation_tests(self, parser):
        """Common test code for testing point mutations."""
        structure = parser.get_structure("PDB/3JQH.cif", "example")

        # Residue 1 and 15 should be disordered.
        res_1 = structure[0]["A"][1]
        res_15 = structure[0]["A"][15]

        # Cursory check -- this would be true even if the residue just
        # contained some disordered atoms.
        self.assertTrue(res_1.disordered, "Residue 1 is disordered")
        self.assertTrue(res_15.disordered, "Residue 15 is disordered")

        # Check a non-mutated residue just to be sure we didn't break the
        # parser and cause everyhing to be disordered.
        self.assertFalse(
            structure[0]["A"][13].disordered,
            "Residue 13 is not disordered")

        # Check that the residue types were parsed correctly.
        self.assertSetEqual(
            set(res_1.disordered_get_id_list()),
            {"PRO", "SER"},
            "Residue 1 is proline/serine")
        self.assertSetEqual(
            set(res_15.disordered_get_id_list()),
            {"ARG", "GLN", "GLU"},
            "Residue 15 is arginine/glutamine/glutamic acid")

        # Quickly check that we can switch between residues and that the
        # correct set of residues was parsed.
        res_1.disordered_select('PRO')
        self.assertAlmostEqual(
            res_1["CA"].occupancy,
            0.83, 2, "Residue 1 proline occupancy correcy")

        res_1.disordered_select('SER')
        self.assertAlmostEqual(
            res_1["CA"].occupancy,
            0.17, 2, "Residue 1 serine occupancy correcy")




if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
