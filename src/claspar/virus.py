"""
Module for handling the viral aligner outputs. Filters are applied to return the taxa of interest.
"""

import pandas as pd
from onyx_analysis_helper import onyx_analysis_helper_functions as oa

from claspar import handle_tables


class VirusClasPar:
    """ """

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

        self.filtered_data: pd.DataFrame = self._filter_viral_aligner()

        self.headline_results: str = ""
        self.results: dict = {}
        self.headline_results, self.results = self._get_viral_aligner_results()

    def _filter_viral_aligner(self) -> pd.DataFrame:
        """
        Process the viral aligner data by applying the filters.

        :param viral_aligner_df: dataframe of the viral aligner outputs straight from scylla.
        :param thresholds_dict: dictionary of the filters to apply to the dataframe.
        :return: pandas dataframe with the taxa that pass the filter.
        """
        # apply filter
        filtered_viral_aligner_df = self.data_input.loc[
            (self.data_input["evenness_value"] >= self.thresholds["EVENNESS_VALUE"])
            & (self.data_input["coverage_1x"] >= self.thresholds["COVERAGE_1X"])
            & (self.data_input["uniquely_mapped_reads"] >= self.thresholds["UNIQUELY_MAPPED_READS"])
            & (self.data_input["mean_read_identity"] >= self.thresholds["MEAN_READ_IDENTITY"])
            & (self.data_input["mean_alignment_length"] >= self.thresholds["MEAN_ALIGNMENT_LENGTH"])
        ]
        return filtered_viral_aligner_df

    def _get_viral_aligner_results(self) -> tuple[str, dict]:
        """
        Get the headline result and the results from Sylph.

        :return: tuple; headline_result (str), result (dict), list of tables to write to csv (all pandas
        dataframes - in this case just one).
        """
        headline_result = (
            f"Sample {self.sample_id} has {self.filtered_data.shape[0]} viral taxa classified by the Viral Aligner "
            f"that passed the filters."
        )
        results = self.filtered_data[["human_readable", "taxon_id"]].to_dict(orient="index")

        return headline_result, results

    def get_virus_analysis_table(self) -> tuple[oa.OnyxAnalysis, int]:
        """
        Pull together all the class attributes into the analysis table.
        """
        exitcode = 0

        analysis_table, exitcode = handle_tables.create_analysis_fields(
            domain="virus",
            classifier="viral_aligner",
            record_id=self.sample_id,
            thresholds=self.thresholds,
            headline_result=self.headline_results,
            results=self.results,
            server=self.server,
        )
        return analysis_table, exitcode
