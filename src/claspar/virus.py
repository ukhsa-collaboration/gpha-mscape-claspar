"""
Module for handling the viral aligner outputs. Filters are applied to return the taxa of interest.
"""

import pandas as pd


def process_viral_aligner(
    viral_aligner_df: pd.DataFrame, thresholds_dict: dict
) -> pd.DataFrame:
    """
    Process the viral aligner data by applying the filters.

    :param viral_aligner_df: dataframe of the viral aligner outputs straight from scylla.
    :param thresholds_dict: dictionary of the filters to apply to the dataframe.
    :return: pandas dataframe with the taxa that pass the filter.
    """
    # apply filter
    filtered_viral_aligner_df = viral_aligner_df.loc[
        (viral_aligner_df["evenness_value"] >= thresholds_dict["EVENNESS_VALUE"])
        & (viral_aligner_df["coverage_1x"] >= thresholds_dict["COVERAGE_1X"])
        & (
            viral_aligner_df["uniquely_mapped_reads"]
            >= thresholds_dict["UNIQUELY_MAPPED_READS"]
        )
        & (
            viral_aligner_df["mean_read_identity"]
            >= thresholds_dict["MEAN_READ_IDENTITY"]
        )
        & (
            viral_aligner_df["mean_alignment_length"]
            >= thresholds_dict["MEAN_ALIGNMENT_LENGTH"]
        )
    ]
    return filtered_viral_aligner_df


def get_viral_aligner_results(
    *,
    sample_id: str,
    original_viral_aligner_results: pd.DataFrame,
    viral_aligner_thresholds_dict: dict,
) -> tuple[str, dict, list[pd.DataFrame]]:
    """
    Get the headline result and the results from Sylph.
    :param sample_id: string of climb_id.
    :param original_viral_aligner_results: pandas dataframe containing the original viral aligner results from Scylla.
    :param viral_aligner_thresholds_dict: dictionary containing the filter thresholds for viral aligner classifications.
    :return: tuple; headline_result (str), result (dict), list of tables to write to csv (all pandas
    dataframes - in this case just one).
    """
    viral_aligner_filtered_df = process_viral_aligner(
        original_viral_aligner_results, viral_aligner_thresholds_dict
    )
    headline_result = (
        f"Sample {sample_id} has {viral_aligner_filtered_df.shape[0]} viral taxa classified by the Viral Aligner "
        f"that passed the filters."
    )
    results = viral_aligner_filtered_df[["human_readable", "taxon_id"]].to_dict(
        orient="index"
    )

    return headline_result, results, [viral_aligner_filtered_df]
