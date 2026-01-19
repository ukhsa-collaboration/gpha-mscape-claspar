"""
Module containing functions needed to parse, filter and create bacteria classifications analysis table.
"""

###########
# Imports #
import logging

import pandas as pd
from taxaplease import TaxaPlease

from claspar import handle_tables


##########
# Kraken #
##########
class KrakenBacteria:
    """
    Class for parsing the bacterial classifications from Kraken.

    The main methods are:
    get_kraken_bacteria_analysis_table - get the analysis table (returns instance of the onyx analysis helper class).
    save_outputs_to_csv - method to save the outputs (saved in instance attributes) to file (returns None)

    Atrributes:
    taxaplease - instance of TaxaPlease.
    classifier_results - the original classifier outputs from Scylla.
    thresholds - dict; the thresholds to filter on.
    sample_id - str; climb-id
    server - str; mscape or synthscape.
    kraken_species_results - pd.DataFrame: All the species kraken identified for the sample, plus the genus id and reads
     at genus level, total species in genus identified (and species that pass the filters), the proportion of total
     genus reads, the rank of that in its genus and the kraken confidence (high or low).
    kraken_genus_results - pd.DataFrame: All the genera kraken identified for the sample, plus some info about the
     species within the genus.
    headline_results - str; the main result, automatically generated to include the final number of taxa that were
     assigned high confidence.
    results - dict; the kraken_species_results dataframe filtered to high confidence species as a dict.
    analysis_table - oa.OnyxAnalysis; instance of the analysis table from the helper, containing all the relevant info.


    :param sample_id: str, climb-id
    :param original_classifier_df: pandas dataframe, the original results from scylla.
    :param kraken_bacteria_thresholds_dict: dict, containing the thresholds to filter.
    :param server: str, mscape or synthscape (database server for Onyx to connect to - will be validated.)
    """

    def __init__(
        self,
        sample_id: str,
        original_classifier_df: pd.DataFrame,
        kraken_bacteria_thresholds_dict: dict[str, int | float],
        server: str = "mscape",
    ):
        """
        On instantiation, populate the thresholds attribute and original dataframe attribute.
        """
        self.taxaplease = TaxaPlease()

        self.classifier_results: pd.DataFrame = original_classifier_df
        self.thresholds: dict = kraken_bacteria_thresholds_dict
        self.sample_id: str = sample_id
        self.server: str = server

        self.headline_result: str
        self.result: dict

        self.kraken_species_results: pd.DataFrame
        self.kraken_genus_results: pd.DataFrame

        self.headline_result, self.results, self.kraken_species_results, self.kraken_genus_results = (
            self._get_kraken_results()
        )

    def _get_parent_taxonomy(self, taxon_id: int) -> tuple[bool, dict | None, str | None]:
        """
        Use taxaplease to get 'is bacteria', 'parent record' and 'parent human-readable'.

        :param taxon_id: int, taxon ID
        :return: boolean ('is bacteria'), dict ('parent record') and str ('parent human-readable').
        """
        tp = self.taxaplease

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

    def _get_kraken_confidence_rating(self, *, count_descendants, order_in_genus, pct_genus_reads):
        """
        Apply the thresholds to determine the kraken confidence rating.
        """
        if (
            count_descendants >= self.thresholds["READ_THRESHOLD"]
            and order_in_genus <= self.thresholds["GENUS_RANK_THRESHOLD"]
            and pct_genus_reads >= self.thresholds["GENUS_READ_PCT_THRESHOLD"]
        ):
            return "high"
        else:
            return "low"

    def _process_kraken(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Process the kraken classifications. Returns two dataframes:

        species_df: contains all the 'S' bacterial taxa and counts, with info added with respect to the parent genus. A
        confidence level (high or low) is applied based on filters.
        genus_df: contains all of 'G' bacterial taxa and counts, with added info about the species, such as the number
        that pass the filter.

        Also modifies the classifier_results dataframe by adding the parent taxonomy information.

        :params classifier_results: directly takes the kraken report as a pandas dataframe.
        :returns: species_df dataframe and genus_df dataframe. If there are no Kraken results (empty dataframe), then
        species and genus df are returned empty.
        """
        tp = self.taxaplease

        classifier_results = self.classifier_results.copy()  # Don't edit the original dataframe

        if classifier_results.empty:
            # If the classifier is empty, just return 2 empty dataframes for species and genus dfs, return early
            return pd.DataFrame(), pd.DataFrame()

        # Add parent taxonomy info
        classifier_results[["is_bacteria", "parent_record", "parent_human_readable"]] = classifier_results.apply(
            lambda row: self._get_parent_taxonomy(row["taxon_id"]),
            axis=1,
            result_type="expand",
        )

        # Get all genus level rows
        genus_df = classifier_results[
            (classifier_results["raw_rank"] == "G") & (classifier_results["is_bacteria"])
        ].copy()
        # Calculate proportion of reads at the species level compared to genus level
        genus_df["prop_species"] = 1 - (genus_df["count_direct"] / genus_df["count_descendants"])

        # get all species level rows for bacterial taxa
        species_df = classifier_results[
            (classifier_results["raw_rank"] == "S") & (classifier_results["is_bacteria"])
        ].copy()
        # add genus ID
        species_df["genus_id"] = species_df["taxon_id"].apply(tp.get_genus_taxid)
        # count number of species
        total_species = species_df.groupby("genus_id").size().reset_index(name="total_species_identified").fillna(0)

        genus_df = pd.merge(genus_df, total_species, left_on="taxon_id", right_on="genus_id")

        # count number of species with >= 10 reads (READ_THRESHOLD value)
        filtered_species = (
            species_df[(species_df["count_descendants"] >= self.thresholds["READ_THRESHOLD"])]
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
            lambda row: self._get_kraken_confidence_rating(
                count_descendants=row["count_descendants"],
                order_in_genus=row["order_in_genus"],
                pct_genus_reads=row["pct_genus_reads"],
            ),
            axis=1,
        )

        return species_df, genus_df

    def _get_kraken_results(self) -> tuple[str, dict, pd.DataFrame, pd.DataFrame]:
        """
        Get the headline result and the results from Kraken for bacteria.
        :param sample_id: string of climb_id.
        :param original_kraken_results: pandas dataframe containing the original kraken results from Scylla.
        :param kraken_thresholds_dict: dictionary containing the filter thresholds for kraken classifications.
        :param taxaplease_instance: instance of taxaplease, default is none and in this case will create a new instance.
        :return: tuple; headline_result (str), result (dict), list of tables to write to csv (all pandas
        dataframes).
        """

        kraken_species, kraken_genus = self._process_kraken()

        if kraken_species.empty:
            high_confidence_species = pd.DataFrame()
            results = {}
        else:
            high_confidence_species = kraken_species.loc[kraken_species["kraken_confidence"] == "high"].reset_index()
            results = high_confidence_species[["human_readable", "taxon_id", "raw_rank"]].to_dict(orient="index")

        headline_result = (
            f"Sample {self.sample_id} has {high_confidence_species.shape[0]} high confidence bacterial "
            f"species classified by Kraken (out of {kraken_species.shape[0]} total assignments)."
        )

        return headline_result, results, kraken_species, kraken_genus

    # Make analysis tables
    def get_kraken_bacteria_analysis_table(self):
        """
        Pull together all the class attributes into the analysis table.
        """
        analysis_table, exitcode = handle_tables.create_analysis_fields(
            domain="bacteria",
            classifier="kraken",
            record_id=self.sample_id,
            thresholds=self.thresholds,
            headline_result=self.headline_result,
            results=self.results,
            server=self.server,
        )

        self.analysis_table = analysis_table
        if exitcode == 1:
            self.exitcode = exitcode

        return analysis_table


##########
# Sylph: #
##########
class SylphBacteria:
    """
    Class for parsing the bacterial classifications from Sylph.

    The main methods are:
    get_sylph_analysis_table - get the analysis table (returns instance of the onyx analysis helper class).
    save_outputs_to_csv - method to save the outputs (saved in instance attributes) to file (returns None)

    Atrributes:
    taxaplease - instance of TaxaPlease.
    classifier_results - the original sylph outputs from Scylla.
    thresholds - dict; the thresholds to filter on.
    sample_id - str; climb-id
    server - str; mscape or synthscape.
    sylph_filtered_results - pd.DataFrame: All the taxa sylph identified for the sample, plus the confidence.
    headline_results - str; the main result, automatically generated to include the final number of taxa that were
     assigned high confidence.
    results - dict; the sylph_filtered_results dataframe filtered to high confidence species as a dict.
    analysis_table - oa.OnyxAnalysis; instance of the analysis table from the helper, containing all the relevant info.


    :param sample_id: str, climb-id
    :param original_classifier_df: pandas dataframe, the original results from scylla.
    :param kraken_bacteria_thresholds_dict: dict, containing the thresholds to filter.
    :param server: str, mscape or synthscape (database server for Onyx to connect to - will be validated.)
    """

    def __init__(
        self,
        sample_id: str,
        original_sylph_df: pd.DataFrame,
        sylph_bacteria_thresholds_dict: dict[str, int | float],
        server: str = "mscape",
    ):
        """
        On instantiation, populate the thresholds attribute and original dataframe attribute.
        """
        self.taxaplease = TaxaPlease()

        self.sylph: pd.DataFrame = original_sylph_df
        self.thresholds: dict = sylph_bacteria_thresholds_dict
        self.sample_id: str = sample_id
        self.server: str = server

        self.headline_result: str
        self.result: dict
        self.sylph_processed_df: pd.DataFrame

        self.headline_result, self.results, self.sylph_processed_df = self._get_sylph_results()

    def _process_sylph_rank(self, row: pd.Series) -> tuple[int | None, str | None]:
        """
        Get the species taxon ID and the species name from row, using taxaplease.

        :params row: pd.series, row of a dataframe (used with an apply).
        :return: list of two; taxon ID and human-readable species name, or [None, None] if rank is not species or strain.
        """
        if row["taxon_id"] is None:
            logging.error(f"Taxon id is not found for {row}")
            return None, None

        if row["taxon_rank"] == "species":
            return row["taxon_id"], row["human_readable"]
        elif row["taxon_rank"] == "strain":
            tp = self.taxaplease
            taxon_id = int(row["taxon_id"])
            species = tp.get_record(tp.get_species_taxid(taxon_id))  # type: ignore
            if species:
                return species["taxid"], species["name"]
            else:
                return None, None
        else:
            print(f"Taxon ID returned rank other than species or strain: {row['taxon_id'], row['taxon_rank']}")
            return None, None

    def _get_sylph_confidence_rating(
        self,
        *,
        containment_index: float,
        effective_coverage: float,
    ) -> str:
        """
        Get the confidence rating for sylph results using thresholds.
        """
        if (
            containment_index >= self.thresholds["CONTAINMENT_INDEX_THRESHOLD"]
            and effective_coverage >= self.thresholds["EFFECTIVE_COVERAGE_THRESHOLD"]
        ):
            return "high"
        else:
            return "low"

    def _process_sylph(
        self,
    ) -> pd.DataFrame:
        """
        Process the sylph outputs and apply filters. Normalise to species level (sylph uses reference genomes which
        could be species or strain level), and add a confidence rating using the filters.
        :return: dataframe, with extra columns for filters. If there are no sylph results, return empty df.
        """

        tp = self.taxaplease

        # Get taxonomic rank of the sylph results
        sylph_df = self.sylph.copy()

        if sylph_df.empty:
            return pd.DataFrame()

        sylph_df["taxon_rank"] = sylph_df["taxon_id"].apply(lambda x: tp.get_record(x)["rank"])  # type: ignore

        sylph_df[["species_id", "species_human_readable"]] = sylph_df.apply(
            lambda x: self._process_sylph_rank(x), axis=1, result_type="expand"
        )

        sylph_df["genus_id"] = sylph_df["species_id"].apply(
            lambda x: tp.get_genus_taxid(x) if x is not None else None
        )  # Note that genus is not always the parent (unclassified, complexes)!

        sylph_df["cont_ind_eval"] = sylph_df["containment_index"].apply(eval)

        sylph_df["sylph_confidence"] = sylph_df.apply(
            lambda row: self._get_sylph_confidence_rating(
                containment_index=row["cont_ind_eval"], effective_coverage=row["effective_coverage"]
            ),
            axis=1,
        )
        return sylph_df

    def _get_sylph_results(self) -> tuple[str, dict, pd.DataFrame]:
        """
        Get the headline result and the results from Sylph.

        :return: tuple; headline_result (str), result (dict), final table from sylph filtering to write to csv. If there
        are no sylph results, return headline result (automatically produced), results as empty dict and
        sylph_processed_df as empty dataframe.
        """
        sylph_processed_df = self._process_sylph()
        if sylph_processed_df.empty:
            high_confidence_species = pd.DataFrame()
            results = {}

        else:
            high_confidence_species = sylph_processed_df.loc[sylph_processed_df["sylph_confidence"] == "high"]
            results = (
                high_confidence_species[["human_readable", "taxon_id", "taxon_rank"]]
                .reset_index(drop=True)
                .to_dict(orient="index")
            )

        headline_result = (
            f"Sample {self.sample_id} has {high_confidence_species.shape[0]} high confidence bacterial "
            f"(and archaeal) species classified by Sylph."
        )

        return headline_result, results, sylph_processed_df

    def get_sylph_analysis_table(self):
        """
        Pull together all the class attributes into the analysis table.
        """
        analysis_table, exitcode = handle_tables.create_analysis_fields(
            domain="bacteria",
            classifier="sylph",
            record_id=self.sample_id,
            thresholds=self.thresholds,
            headline_result=self.headline_result,
            results=self.results,
            server=self.server,
        )

        self.analysis_table = analysis_table
        if exitcode == 1:
            self.exitcode = exitcode

        return analysis_table
