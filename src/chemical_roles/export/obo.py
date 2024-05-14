# -*- coding: utf-8 -*-

"""Export to OBO."""

from typing import Dict, Iterable, Mapping, Optional, Tuple

from pyobo import Obo, Reference, Term, TypeDef
from pyobo.struct import make_ad_hoc_ontology
from tqdm import tqdm
import click
from chemical_roles.export.utils import get_relations_df
from chemical_roles.constants import directory_option
import os

__all__ = [
    "get_obo",
]


def get_obo() -> Obo:
    """Get Chemical Roles as OBO."""
    return make_ad_hoc_ontology(
        _name="Chemical Roles Graph",
        _ontology="crog",
        terms=iter_terms(),
    )


def iter_terms() -> Iterable[Term]:
    df = get_relations_df()
    it = tqdm(df.dropna().values, total=len(df.index), desc="mapping to OBO", unit_scale=True)
    ref_term = {}
    for (
        source_db,
        source_id,
        source_name,
        modulation,
        target_type,
        target_db,
        target_id,
        target_name,
    ) in it:
        source = Reference(prefix=source_db.upper(), identifier=source_id, name=source_name)
        term = ref_term.get(source)
        if term is None:
            term = ref_term[source] = Term(reference=source)

        typedef = _get_typedef(target_db=target_db, target_type=target_type, modulation=modulation)
        if typedef is not None:
            term.append_relationship(typedef, Reference(prefix=target_db.upper(), identifier=target_id, name=target_name))

    for term in ref_term.values():
        if len(term.relationships) > 0:
            yield term


# TODO fill out rest
_typedefs: Dict[Tuple[str, str, str], TypeDef] = {
    ("go", "biological process", "activator"): TypeDef(
        Reference(prefix="RO", identifier="0002213", name="positively regulates")
    ),
    ("go", "biological process", "inhibitor"): TypeDef(
        Reference(prefix="RO", identifier="0002212", name="negatively regulates")
    ),
    ("go", "biological process", "modulator"): TypeDef(
        Reference(prefix="RO", identifier="0002211", name="regulates")
    ),

}

for p in ["uniprot", "fplx", "hgnc", "pr" , "hgnc.genefamily"]:
    _typedefs.update({
        (p, "protein", "activator"): TypeDef(
            Reference(prefix="RO", identifier="0002450", name="directly positively regulates activity of")
        ),
        (p, "protein", "inhibitor"): TypeDef(
            Reference(prefix="RO", identifier="0002449", name="directly negatively regulates activity of")
        ),
        (p, "protein", "agonist"): TypeDef(
            Reference(prefix="RO", identifier="0018027", name="is agonist of")
        ),
        (p, "protein", "antagonist"): TypeDef(
            Reference(prefix="RO", identifier="0018029", name="is antagonist of")
        ),
        (p, "protein", "antagonist"): TypeDef(
            Reference(prefix="RO", identifier="0018028", name="is inverse agonist of")
        ),
    })



_logged = set()


def _get_typedef(target_db, target_type, modulation) -> Optional[TypeDef]:
    t = (target_db, target_type, modulation)
    rv = _typedefs.get(t)
    if rv is not None:
        return rv
    if t not in _logged:
        _logged.add(t)
        tqdm.write(f"no strategy for: {target_db} {target_type} {modulation}")


@click.command()
@directory_option
def main(directory):
    """Write OBO export."""
    o = get_obo()
    o.write_obo(os.path.join(directory, "crog.obo"))
    o.write_obonet_gz(os.path.join(directory, "crog.obonet.json.gz"))


if __name__ == "__main__":
    main()
