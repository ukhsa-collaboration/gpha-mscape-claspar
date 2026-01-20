"""
Module for handling the viral aligner outputs. Filters are applied to return the taxa of interest.
"""

import logging
import os

import pandas as pd
from onyx_analysis_helper import onyx_analysis_helper_functions as oa

from claspar import handle_tables


class VirusClasPar:
    """
    Class for parsing the viral aligner classifications.

    The main methods are:
    get_virus_analysis_table - get the analysis table (returns instance of the onyx analysis helper class).
    save_outputs_to_csv - method to save the outputs (saved in instance attributes) to file (returns None)

    Attributes:
    exitcode - int; updated after init and get_virus_analysis_table.
    data_input - pd.DataFrame; original viral aligner results from scylla.
    thresholds - dict; the thresholds for filtering
    sample_id - str; the climb-id
    server - str; mscape or synthscape
    filtered_data - pd.DataFrame; the viral aligner results after filtering
    headline_results - str; the main result, automatically generated to include the final number of taxa that remained
    after filtering
    results - dict; the filtered dataframe as a dict.
    analysis_table - oa.OnyxAnalysis; instance of the analysis table from the helper, containing all the relevant info.

    :param sample_id: str, climb-id
    :param original_viral_aligner_df: pandas dataframe, the original results from scylla.
    :param virus_thresholds_dict: dict, containing the thresholds to filter.
    :param server: str, mscape or synthscape (database server for Onyx to connect to - will be validated.)
    """

    def __init__(
        self,
        sample_id: str,
        original_viral_aligner_df: pd.DataFrame,
        virus_thresholds_dict: dict[str, int | float],
        server: str = "mscape",
    ):
        """
        On instantiation, populate the thresholds attribute and original dataframe attribute.
        """
        self.data_input: pd.DataFrame = original_viral_aligner_df
        self.thresholds: dict = virus_thresholds_dict
        self.sample_id: str = sample_id
        self.server: str = server

        self.filtered_data: pd.DataFrame
        self.exitcode: int
        self.exitcode, self.filtered_data = self._filter_viral_aligner()

        self.headline_results: str
        self.results: dict
        self.headline_results, self.results = self._get_viral_aligner_results()

    def _filter_viral_aligner(self) -> tuple[int, pd.DataFrame]:
        """
        Process the viral aligner data by applying the filters. If the input is empty, return empty dataframe.

        :param viral_aligner_df: dataframe of the viral aligner outputs straight from scylla.
        :param thresholds_dict: dictionary of the filters to apply to the dataframe.
        :return: tuple of exitcode, filtered df.
        Filtered df is pandas dataframe with the taxa that pass the filter. (Empty df returned if input is also empty.)
        """
        # check if empty:
        if self.data_input.empty:
            return 0, pd.DataFrame()
        # apply filter
        else:
            try:
                filtered_viral_aligner_df = self.data_input.loc[
                    (self.data_input["evenness_value"] >= self.thresholds["EVENNESS_VALUE"])
                    & (self.data_input["coverage_1x"] >= self.thresholds["COVERAGE_1X"])
                    & (self.data_input["uniquely_mapped_reads"] >= self.thresholds["UNIQUELY_MAPPED_READS"])
                    & (self.data_input["mean_read_identity"] >= self.thresholds["MEAN_READ_IDENTITY"])
                    & (self.data_input["mean_alignment_length"] >= self.thresholds["MEAN_ALIGNMENT_LENGTH"])
                ]
                return 0, filtered_viral_aligner_df
            except KeyError as k:
                logging.error("The viral alignment data from Scylla is missing expected column:%s", k)
                return 1, pd.DataFrame()

    def _get_viral_aligner_results(self) -> tuple[str, dict]:
        """
        Get the headline result and the results from Sylph.

        :return: tuple; headline_result (str), result (dict). If no alignment data, headline_result is string and
        results is empty dict.
        """
        headline_result = (
            f"Sample {self.sample_id} has {self.filtered_data.shape[0]} viral taxa classified by the Viral Aligner "
            f"that passed the filters (out of a total of {self.data_input.shape[0]})."
        )
        if self.filtered_data.empty:
            return headline_result, {}
        else:
            results = self.filtered_data[["human_readable", "taxon_id"]].reset_index(drop=True).to_dict(orient="index")

        return headline_result, results

    def get_virus_analysis_table(self) -> oa.OnyxAnalysis:
        """
        Pull together all the class attributes into the analysis table.
        """

        analysis_table, exitcode = handle_tables.create_analysis_fields(
            domain="virus",
            classifier="viral aligner",
            record_id=self.sample_id,
            thresholds=self.thresholds,
            headline_result=self.headline_results,
            results=self.results,
            server=self.server,
        )
        # Check the exitcode attribute - it might have broken elsewhere...
        if exitcode == 1:
            logging.error("There was an error creating the analysis table. Check the Onyx Analysis Helper.")
            self.exitcode = 1

        # Populate the instance attribute:
        self.analysis_table: oa.OnyxAnalysis = analysis_table

        return analysis_table

    def save_outputs_to_csv(self, results_dir: str | os.PathLike) -> None:
        """
        Save the final results to csv.
        :param filename: str, name of file to save to.
        :param results_dir: str or path to directory to save to.
        """
        unique_filename = f"{self.sample_id}_filtered_viral_aligner_results"
        handle_tables.write_df_to_csv(df=self.filtered_data, filename=unique_filename, results_dir=results_dir)
