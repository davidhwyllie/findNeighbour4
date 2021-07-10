#!/usr/bin/env python
""" 
A component of a findNeighbour4 server which provides relatedness information for bacterial genomes.
Provides a persistence layer for the output of PCA, including modelling of category frequencies over time.

Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without tcen the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
"""

import pandas as pd
import warnings
import statsmodels.api as sm
from patsy import dmatrices


class ModelCounts:
    """fits negative binomial and Poission models to count data, where each pc_cat in a pc is represented"""

    def __init__(
        self,
        earliest_date,
        latest_date,
        counts,
        denominators,
        pcas_int_id,
        pc_cat,
        show_bar=True,
    ):
        """
        Stores data to fit linear Poisson models to count data

        Parameters:
        -----------
        earliest_date: the earliest date to analyse
        latest_date: the latest date to analyse
        counts:  a pandas data frame containing counts
        pc_cat: the pc_cat analysed.  Is not required, but pcas_int_id is
        pcas_int_id: the primary key to the pca_summary table row in which the population/pc_cat combination is specified
        denominators: a pandas data frame containing denominators

        Suitable input can be generated by the  PCADatabaseManager.pcas_count_table(format = 2 ) method

        Returns:
        --------
        None

        """

        # check that the count data frame has the expected columns
        if not (
                    set(counts.keys()),
                    set(["sample_date","pc_cat", "n"]),
                ):
                raise KeyError("Passed data frame has {0} columns when {1} are expected.".format(set(counts.keys()),set(["sample_date","pc_cat", "n"])))

        self.show_bar = show_bar
        self.earliest_date = earliest_date
        self.latest_date = latest_date
        self.date_range = pd.DataFrame(
            {
                "sample_date": [
                    x.date()
                    for x in pd.date_range(earliest_date, latest_date, freq="D")
                ]
            }
        )
        self.interval_analysed_days = len(self.date_range)

        # chop off any data after latest date
        self.counts = counts[counts["sample_date"] <= latest_date]

        denominator_subset = denominators[denominators["sample_date"] <= latest_date]
        self.denominators = denominator_subset.rename(columns={"n": "n_total"})
        self.pc_cat = pc_cat
        self.pcas_int_id = pcas_int_id

        self.date_fields = ["sample_date",'pc_cat']

        to_model = self.date_range.merge(
            self.denominators, how="left", on=['sample_date']
        )
        if 'pc_cat' not in self.counts.columns.to_list():
            self.counts['pc_cat'] = pc_cat
            
        pc_cats = pd.DataFrame(
            {
                "pc_cat":self.counts['pc_cat'].unique()
            }
        
        )
        
        to_model = to_model.merge(pc_cats, how='cross')
        to_model['n_total'] = to_model['n_total'].fillna(0)

        self.to_model = to_model.merge(self.counts, how="left", on=self.date_fields)
        self.to_model['n'] = self.to_model['n'].fillna(0)
        

    def data_to_model(self, drop_days_with_no_samples=True):
        """
        returns a table consisting of
        date,
        pc_cat,
        count,
        denominator count

        for each date between self.start_date and self.end_date, including dates with zeros
        (if there is any data)
        or None (if there is no data)

        Parameters:
        drop_days_with_no_samples: if there are no denominator samples on a day, drop this day (required for Poisson model with offset)

        Returns: a pandas DataFrame
        """

        fit_df = self.to_model.copy()

        fit_df["t"] = [x.days for x in fit_df["sample_date"] - self.latest_date]
        fit_df["day_of_week"] = [x.weekday() for x in fit_df["sample_date"]]
        fit_df['sample_dow'] = ['dow'+str(x) for x in fit_df['day_of_week']]
        if drop_days_with_no_samples is True:
            fit_df = fit_df[fit_df["n_total"] > 0]

        return fit_df
    
    def fit_nb(self, raise_error_for_unittest=False):
        """fits negative binomial model with offset"""

        cnt_df = self.data_to_model(drop_days_with_no_samples = True)
        analysis_status = "started"
        errors_returned = None
        coeff_df = None
        # R-style formula
        formula = """n ~ t*C(pc_cat, Treatment(reference='{0}'))""".format(self.pc_cat)
        warnings.simplefilter("error", category=RuntimeWarning)

        try:       # debug disabled
            if raise_error_for_unittest:
                raise ZeroDivisionError("Test error")
            
            # generate model inputs
            
            y, X = dmatrices(formula, cnt_df, return_type="dataframe")

            # fit model.  exposure = X just means offset = log(x).

            # model raises some divide by zero and failure to converge as warnings.
            # for the purpose of tracking modelling , we regard all of this as failure to fit
            # poisson_fit = sm.GLM(
            #    y, X, exposure=cnt_df["n_total"], family=sm.families.Poisson()
            #).fit()

            poisson_fit = sm.GLM(
                y, X, exposure=cnt_df["n_total"], family=sm.families.NegativeBinomial()
            ).fit()

            cnt_df["pred"] = poisson_fit.predict(X)
            coeff_df = poisson_fit.summary2().tables[1]
            
            # select the intercept and the time effect for the relevant pc_cat, which is set to the reference ccategory
            coeff_df = coeff_df.loc[["Intercept", "t"]]
        
            # rename the data to match the StatisticalModelFit table

            coeff_df["param"] = coeff_df.index
            coeff_df = coeff_df.rename(
                columns={
                    "Coef.": "estimate",
                    "Std.Err.": "std_err",
                    "[0.025": "estimate_ci_low",
                    "0.975]": "estimate_ci_high",
                    "P>|z|": "p_value",
                }
            )

            coeff_df["estimate2naturalspace"] = "exp"
            coeff_df["is_reference"] = False
            coeff_df["has_ci"] = False
            coeff_df["comments"] = "Poisson regression on a subset of data"
            coeff_df["param_desc"] = "-"
            coeff_df.at[
                "t", "param_desc"
            ] = "Estimated lineage change over the analysis period"
            coeff_df.at[
                "Intercept", "param_desc"
            ] = "Modelled rate at end of time period, per sample tested"
            analysis_status = "completed"

        except Exception as e:     # debug disabled
            errors_returned = str(e)
            analysis_status = "failed"
            cnt_df["pred"] = None
        warnings.resetwarnings()

        # reference to the PCASummary table
        retVal = dict(
            statistical_model=dict(
                pcas_int_id=self.pcas_int_id,
                date_start=self.earliest_date,
                date_end=self.latest_date,
                formula=formula,
                analysis_readable_title="Negative binomial regression",
                analysis_type="GLM:NegativeBinomial",
                readable_description="NegativeBinomial regression, estimating rate of events relative to total number of samples.  Samples are analysed for one pc, allowing different rates for different pc_cats",
                analysis_software="python:statmodels.GLM with NegBinomial family",
                interval_analysed_days=self.interval_analysed_days,
                analysis_status=analysis_status,
                errors_returned=errors_returned,
                run_parameters="default",
            ),
            modelled_data=cnt_df,
            coefficients=coeff_df,
        )

   
        return retVal

    def fit_poisson(self, raise_error_for_unittest=False):
        """fits poisson model with offset"""

        cnt_df = self.data_to_model(drop_days_with_no_samples = True)
        cnt_df = cnt_df[cnt_df['pc_cat']==self.pc_cat]

        analysis_status = "started"
        errors_returned = None
        coeff_df = None
        # R-style formula
        formula = """n ~ t """
        warnings.simplefilter("error", category=RuntimeWarning)

        try:

            if raise_error_for_unittest:
                raise ZeroDivisionError("Test error")
            
            # generate model inputs
            
            y, X = dmatrices(formula, cnt_df, return_type="dataframe")

            # fit model.  exposure = X just means offset = log(x).

            # model raises some divide by zero and failure to converge as warnings.
            # for the purpose of tracking modelling , we regard all of this as failure to fit
            poisson_fit = sm.GLM(
                y, X, exposure=cnt_df["n_total"], family=sm.families.Poisson()
            ).fit()

            cnt_df["pred"] = poisson_fit.predict(X)

            # coeffs to data frame
            coeff_df = poisson_fit.summary2().tables[1]

            # select the intercept and the time effect
            coeff_df = coeff_df.loc[["Intercept", "t"]]

            # rename the data to match the StatisticalModelFit table

            coeff_df["param"] = coeff_df.index
            coeff_df = coeff_df.rename(
                columns={
                    "Coef.": "estimate",
                    "Std.Err.": "std_err",
                    "[0.025": "estimate_ci_low",
                    "0.975]": "estimate_ci_high",
                    "P>|z|": "p_value",
                }
            )

            coeff_df["estimate2naturalspace"] = "exp"
            coeff_df["is_reference"] = False
            coeff_df["has_ci"] = False
            coeff_df["comments"] = "Poisson regression on a subset of data"
            coeff_df["param_desc"] = "-"
            coeff_df.at[
                "t", "param_desc"
            ] = "Estimated lineage change over the analysis period"
            coeff_df.at[
                "Intercept", "param_desc"
            ] = "Modelled rate at end of time period, per sample tested"
            analysis_status = "completed"

        except Exception as e:
            errors_returned = str(e)
            analysis_status = "failed"
            cnt_df["pred"] = None
            cnt_df['pc_cat'] = self.pc_cat
        warnings.resetwarnings()
        # reference to the PCASummary table
        retVal = dict(
            statistical_model=dict(
                pcas_int_id=self.pcas_int_id,
                date_start=self.earliest_date,
                date_end=self.latest_date,
                formula=formula,
                analysis_readable_title="Poisson regression, pc_cat wise",
                analysis_type="GLM:Poisson regression",
                readable_description="Poisson regression, estimating rate of events relative to total number of samples.  Samples are analysed in strata, pc_cat by pc_cat",
                analysis_software="python:statmodels.GLM with Poisson family",
                interval_analysed_days=self.interval_analysed_days,
                analysis_status=analysis_status,
                errors_returned=errors_returned,
                run_parameters="default",
            ),
            modelled_data=cnt_df,
            coefficients=coeff_df,
        )

        return retVal

