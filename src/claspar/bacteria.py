"""
Module containing functions needed to parse, filter and create bacteria classifications analysis table.
"""

###########
# Imports #

import pandas as pd
from epiweeks import Week
from taxaplease import TaxaPlease

#########
# Setup #

# Kraken filters
READ_THRESHOLD: int = 10
GENUS_RANK_THRESHOLD: int = 3
GENUS_READ_PCT_THRESHOLD: int = 20

# Sylph filters
CONTAINMENT_INDEX_THRESHOLD: float = 0.2
EFFECTIVE_COVERAGE_THRESHOLD: float = 1.0


# Functions


# Kraken:
def _process_collection_date(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Handle dates in metadata dataframe. If the collection_date column is null, it is replaced with the received_date.
    """
    metadata["collection_date"] = pd.to_datetime(metadata["collection_date"])
    metadata["collection_date_merged"] = metadata["collection_date"]
    metadata["received_date"] = pd.to_datetime(metadata["received_date"])
    metadata.loc[
        (metadata["collection_date_merged"].isnull()), "collection_date_merged"
    ] = metadata.loc[(metadata["collection_date_merged"].isnull()), "received_date"]
    metadata["collection_date_month"] = metadata["collection_date_merged"].dt.strftime(
        "%Y-%b"
    )
    metadata["collection_date_epi_week"] = metadata.collection_date_merged.apply(
        Week.fromdate
    ).astype(str)

    return metadata


def _get_parent_taxonomy(
    taxon_id: int, taxaplease_instance: TaxaPlease = None
) -> tuple[bool, dict | None, str | None]:
    """
    Use taxaplease to get 'is bacteria', 'parent record' and 'parent human-readable'.

    :param taxon_id: int, taxon ID
    :param taxaplease_instance: class instance of taxaplease.
    :return: boolean ('is bacteria'), dict ('parent record') and str ('parent human-readable').
    """
    tp = taxaplease_instance if taxaplease_instance else TaxaPlease()

    if taxon_id == 0:
        is_bacteria = False
        parent_record = {}
        parent_human_readable = ""
        return is_bacteria, parent_record, parent_human_readable

    else:
        is_bacteria = tp.isBacteria(taxon_id)
        parent_record = tp.get_parent_record(taxon_id)
        parent_human_readable = parent_record.get("name", None) if parent_record else ""
        return is_bacteria, parent_record, parent_human_readable


def _get_kraken_confidence_rating(
    *, count_descendants, order_in_genus, pct_genus_reads
):
    """
    Apply the thresholds to determine the kraken confidence rating.
    """
    if (
        count_descendants >= READ_THRESHOLD
        and order_in_genus <= GENUS_RANK_THRESHOLD
        and pct_genus_reads >= GENUS_READ_PCT_THRESHOLD
    ):
        return "high"
    else:
        return "low"


def process_kraken(
    classifier_results: pd.DataFrame, taxaplease_instance: TaxaPlease = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Process the kraken classifications. Returns two dataframes:

    species_df: contains all the 'S' bacterial taxa and counts, with info added with respect to the parent genus. A
    confidence level (high or low) is applied based on filters.
    genus_df: contains all of 'G' bacterial taxa and counts, with added info about the species, such as the number that
    pass the filter.

    Also modifies the classifier_results dataframe by adding the parent taxonomy information.

    :params classifier_results: directly takes the kraken report as a pandas dataframe.
    :params taxaplease_instance: instance of taxaplease.
    :returns: species_df dataframe and genus_df dataframe
    """

    tp = taxaplease_instance if taxaplease_instance else TaxaPlease()

    classifier_results = classifier_results.copy()  # Don't edit the original dataframe
    # Add parent taxonomy info
    classifier_results[["is_bacteria", "parent_record", "parent_human_readable"]] = (
        classifier_results.apply(
            lambda row: _get_parent_taxonomy(row["taxon_id"], tp),
            axis=1,
            result_type="expand",
        )
    )

    # Get all genus level rows
    genus_df = classifier_results[
        (classifier_results["raw_rank"] == "G") & (classifier_results["is_bacteria"])
    ].copy()
    # Calculate proportion of reads at the species level compared to genus level
    genus_df["prop_species"] = 1 - (
        genus_df["count_direct"] / genus_df["count_descendants"]
    )

    # get all species level rows for bacterial taxa
    species_df = classifier_results[
        (classifier_results["raw_rank"] == "S") & (classifier_results["is_bacteria"])
    ].copy()
    # add genus ID
    species_df["genus_id"] = species_df["taxon_id"].apply(tp.get_genus_taxid)
    # count number of species
    total_species = (
        species_df.groupby("genus_id")
        .size()
        .reset_index(name="total_species_identified")
        .fillna(0)
    )

    genus_df = pd.merge(
        genus_df, total_species, left_on="taxon_id", right_on="genus_id"
    )

    # count number of species with >= 10 reads (READ_THRESHOLD value)
    filtered_species = (
        species_df[(species_df["count_descendants"] >= READ_THRESHOLD)]
        .groupby("genus_id")
        .size()
        .reset_index(name="filtered_species_identified")
        .fillna(0)
    )
    genus_df = pd.merge(genus_df, filtered_species, on="genus_id", how="left")

    # merge genus info rows with species info
    genus_merge = genus_df[
        [
            "genus_id",
            "count_descendants",
            "total_species_identified",
            "filtered_species_identified",
        ]
    ].rename(columns={"count_descendants": "genus_level_reads"})

    species_df = pd.merge(species_df, genus_merge, on="genus_id")

    # calculate percentage of genus level reads
    species_df["pct_genus_reads"] = (
        species_df["count_descendants"] / species_df["genus_level_reads"] * 100
    )  # Proportion of the species reads from total genus reads
    order_in_genus = (
        species_df.sort_values("pct_genus_reads", ascending=False)
        .groupby("genus_id", sort=False)
        .cumcount()
        .shift()
        .rename("order_in_genus")
    )

    species_df = pd.merge(species_df, order_in_genus, left_index=True, right_index=True)

    species_df["kraken_confidence"] = species_df.apply(
        lambda row: _get_kraken_confidence_rating(
            row["count_descendants"], row["order_in_genus"], row["pct_genus_reads"]
        ),
        axis=1,
    )

    return species_df, genus_df


# Sylph:
def _process_sylph_rank(row: pd.Series) -> tuple[int | None, str | None]:
    """
    Get the species taxon ID and the species name from row, using taxaplease.
    :params row: pd.series, row of a dataframe (used with an apply).
    :return: list of two; taxon ID and human-readable species name, or [None, None] if rank is not species or strain.
    """
    if row["taxon_rank"] == "species":
        return row["taxon_id"], row["human_readable"]
    elif row["taxon_rank"] == "strain":
        tp = TaxaPlease()
        species = tp.get_record(tp.get_species_taxid(int(row["taxon_id"])))
        return species["taxid"], species["name"]
    else:
        print(
            f"Taxon ID returned rank other than species or strain: {row['taxon_id'], row['taxon_rank']}"
        )
        return None, None


def _get_sylph_confidence_rating(
    *, containment_index: float, effective_coverage: float
) -> str:
    """
    Get the confidence rating for sylph results using thresholds.
    """
    if (
        containment_index >= CONTAINMENT_INDEX_THRESHOLD
        and effective_coverage >= EFFECTIVE_COVERAGE_THRESHOLD
    ):
        return "high"
    else:
        return "low"


def process_sylph(
    sylph_out: pd.DataFrame, taxaplease_instance: TaxaPlease = None
) -> pd.DataFrame:
    """
    Process the sylph outputs and apply filters. Normalise to species level (sylph uses reference genomes which could be
    species or strain level), and add a confidence rating using the filters.
    """

    tp = taxaplease_instance if taxaplease_instance else TaxaPlease()
    # Get taxonomic rank of the sylph results
    sylph_out["taxon_rank"] = sylph_out["taxon_id"].apply(
        lambda x: tp.get_record(x)["rank"]
    )

    sylph_out[["species_id", "species_human_readable"]] = sylph_out.apply(
        lambda x: _process_sylph_rank(x), axis=1, result_type="expand"
    )

    sylph_out["genus_id"] = sylph_out["species_id"].apply(
        lambda x: tp.get_genus_taxid(x) if x is not None else None
    )  # Note that genus is not always the parent (unclassified, complexes)!

    sylph_out["cont_ind_eval"] = sylph_out["containment_index"].apply(eval)

    sylph_out["sylph_confidence"] = sylph_out.apply(
        lambda row: _get_sylph_confidence_rating(
            containment_index=row["cont_ind_eval"],
            effective_coverage=row["effective_coverage"],
        ),
        axis=1,
    )
    return sylph_out
