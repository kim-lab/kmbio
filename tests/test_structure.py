
import ast
from pathlib import Path

import pandas as pd
import pytest

import kmbio.PDB
from kmbio.PDB.core import Structure


@pytest.fixture(
    params=list(Path("PDB").glob("*.pdb")) + list(Path("PDB").glob("*.cif")), ids=lambda p: p.name
)
def structure(request):
    structure = kmbio.PDB.load(request.param)
    # === Populate `structure_dfs` folder ===
    # df = structure.to_dataframe()
    # df.to_csv(
    #     Path("structure_dfs").joinpath(request.param.name).with_suffix(".tsv"),
    #     sep="\t",
    #     index=False,
    # )
    return structure


@pytest.fixture(params=list(Path("structure_dfs").glob("*.tsv")), ids=lambda p: p.name)
def df(request):
    df = pd.read_csv(request.param, sep="\t", na_values=[""], keep_default_na=False)
    df["atom_extra_bonds"] = df["atom_extra_bonds"].apply(ast.literal_eval)
    return df


def test_s_df_s(structure):
    """Test Structure -> DataFrame -> Structure."""
    df = structure.to_dataframe()
    structure_ = Structure.from_dataframe(df)
    assert structure == structure_


def test_df_s_df(df):
    """Test DataFrame -> Structure -> DataFrame."""
    structure = Structure.from_dataframe(df)
    df_ = structure.to_dataframe()
    assert df.equals(df_)
