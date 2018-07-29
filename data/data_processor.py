import contextlib
import datetime as dt
import os
import random
import re
import sys
import time
from itertools import product
from math import exp

import numpy as np
import pandas as pd
from scipy import stats


class DataProcessor:
    """
    Keyword arguments:
    file_pahts: A list of filepaths to the files containing the time series
    skip_headers: Whether to ignore the headers
    delta: Length of the interval in seconds between two adjacent observations
    regime: Start time and end time of the trading regime
    """

    def __init__(self):
        """
        Initialises attributes and starts the counter.
        """
        self.start_time = time.clock()

    def __del__(self):
        print(f"\n\nFinished in {time.clock() - self.start_time} seconds.\n\n")

    @staticmethod
    def _peek_line(file):
        """
        Reads a single csv line of a file without moving the file's pointer.
        Returns the line string.
        """
        pos = file.tell()
        line = file.readline()
        file.seek(pos)
        return line

    @staticmethod
    def _read_line(file):
        """
        Reads a single csv line of a file and returns a list of values.
        """
        line = file.readline()

        if not line:
            return False

        line = line[:-1]
        return line.split(",")

    @staticmethod
    def _parse_line(line):
        """
        Parses a single csv line, returns a (timestamp, value) tuple.
        """
        pattern = r"^([\d/]+,[\d:.]+),([\d.]+)"
        match = re.search(pattern, line).groups()
        timestamp = dt.datetime.strptime(match[0], "%m/%d/%Y,%H:%M:%S.%f")
        value = float(match[1])
        return (timestamp, value)

    @staticmethod
    def _round_datetime_up(t, delta):
        """
        Rounds a timestamp t down to the nearest delta-interval and
        returns the rounded timestamp.
        """
        seconds = (t - t.min).seconds
        rounding = (seconds + delta) // delta * delta
        return t + dt.timedelta(0, rounding - seconds, -t.microsecond)

    @staticmethod
    def _write_line(x, file):
        """
        Concats list x into a csv line and writes it into the file.
        """
        line = ",".join(map(str, x)) + "\n"
        file.write(line)

    @staticmethod
    def _update_progress(files, file_sizes):
        """
        Updates the progress in the console.
        """
        print("", end="\r")
        for i, f in enumerate(files):
            perc = round(f.tell() / file_sizes[i] * 100, 2)
            print("File {}: {}%".format(i + 1, perc), end=" ")
            sys.stdout.flush()

    @staticmethod
    def skip_lines(files, n):
        """
        Skips n lines in each of the files.
        """
        for f in files:
            for _ in range(n):
                f.readline()

    def previous_tick(self, file_paths, skip_headers, delta, regime):
        """
        Implementation of Previous Tick synchronisation scheme as used by
        Zhang, L. in Estimating covariation: Epps effect, microstructure
        noise (2011).

        Synchronises multiple time series from files object using the
        Previous Tick scheme. The last known observation before each timestamp
        is used. Observations are equally spaced with delta seconds inbetween
        them.

        Keyword arguments:

        """
        files = [stack.enter_context(open(fp)) for fp in file_paths]
        file_sizes = [os.stat(f.fileno()).st_size for f in files]

        if skip_headers:
            self.skip_lines(files, 1)

        regime_start = dt.datetime.strptime(regime[0], "%H:%M").time()
        regime_end = dt.datetime.strptime(regime[1], "%H:%M").time()

        # get the shared initial timestamp t1
        obs = [self._peek_line(f) for f in files]
        obs = [self._parse_line(l) for l in obs]
        t1 = max([o[0] for o in obs])
        t1 = self._round_datetime_up(t1, delta)
        t1 = max(t1, dt.datetime.combine(t1, regime_start))

        print("\n\n# Previous Tick Synchronisation #\n")

        for f in files:
            eof = False
            new_timestamp = t1
            last_value = None
            nobs = 0

            # generate the output filename
            f_name, f_ext = os.path.splitext(f.name)
            out_path = f_name + '_pt' + str(delta) + 's' + f_ext

            with open(out_path, "w") as out_file:
                while not eof:
                    line = f.readline()
                    if line:
                        (timestamp, new_value) = self._parse_line(line)
                        # timestamp must be within the trading regime
                        if (timestamp.time() < regime_start or
                                timestamp.time() > regime_end or
                                timestamp < t1 - dt.timedelta(0, delta)):
                            continue
                        value = new_value
                        # imputing the missing observations and writing results
                        while new_timestamp <= timestamp:
                            obs = [new_timestamp, last_value, nobs]
                            self._write_line(obs, out_file)
                            # update timestamp
                            new_timestamp += dt.timedelta(0, delta)
                            # if over the regime -> switch to the next day
                            if new_timestamp.time() > regime_end:
                                new_timestamp += dt.timedelta(days=1)
                                new_timestamp = dt.datetime.combine(
                                    new_timestamp, regime_start) + dt.timedelta(0, delta)
                            # record progress
                            if nobs > 0:
                                self._update_progress(files, file_sizes)
                            nobs = 0
                        # the last observation in the period not reached yet
                        last_value = value
                        nobs += 1
                    else:
                        eof = True
                        obs = [new_timestamp, last_value, nobs]
                        self._write_line(obs, out_file)
                        self._update_progress(files, file_sizes)

    def nelson_siegel(self, series, lmbda, out_path):
        """
        Goes trhough the synchronised time series, computes yields for each series and
        then computes the Nelson-Siegel factor parameters.

        Keyword arguments:
        series: data about the time series
        lmbda: lambda parameter from the nelson-siegel model
        out_path: file path to the merged output file
        """

        file_paths = [s['path'] for s in series]
        conv_factors = [s['conv_factor'] for s in series]
        contr_mults = [s['contr_mult'] for s in series]
        maturities = [s['maturity'] for s in series]
        fvs = [s['fv'] for s in series]

        # X matrix for the Nelson-Siegel factor regression
        X_ns = np.array(
            [[1, (1 - exp(-lmbda * m)) / (lmbda * m), (1 - exp(-lmbda * m)) / (lmbda * m) - exp(-lmbda * m)]
             for m in maturities]
        )

        # create output file header
        tickers = [s['ticker'] for s in series]
        variables = ["close", "nobs", "yield"]
        out_header = "timestamp," + ",".join(["_".join(x)
                                              for x in product(tickers, variables)])
        out_header += ",ns_beta0,ns_beta1,ns_beta2"

        print("\n\n# Nelson-Siegel Merging #\n")

        # write the merged csv
        with open(out_path, "w") as out_file:
            files = [stack.enter_context(open(fp)) for fp in file_paths]
            file_sizes = [os.stat(f.fileno()).st_size for f in files]

            out_file.write(out_header + "\n")
            eof = False

            while not eof:
                lines = []

                # read a line in each file
                for f in files:
                    line = f.readline()
                    if line:
                        line = line.replace("\n", "").split(",")
                        lines.append(line)
                    else:
                        eof = True

                if not eof:
                    # all timestamps mnust be equal
                    timestamps = [l[0] for l in lines]
                    if len(set(timestamps)) == 1:
                        timestamp = timestamps[0]
                    else:
                        raise ValueError(f"Unequal timestamps: {timestamps}.")
                    w_line = timestamp + ","

                    # get a close, nobs and yield for each series
                    close = [float(l[1]) for l in lines]
                    nobs = [int(l[2]) for l in lines]
                    yields = [(fv / (c * cf * cm))**(1 / mat) - 1
                              for c, cf, cm, mat, fv in zip(close, conv_factors, contr_mults, maturities, fvs)]
                    w_line += ",".join([",".join([str(y) for y in x])
                                        for x in zip(close, nobs, yields)]) + ","

                    # compute Nelson-Siegel Factors
                    y_ns = np.array(yields)
                    coef_ns = np.linalg.lstsq(X_ns, y_ns, rcond=None)[0]
                    w_line += ",".join([str(x) for x in coef_ns])

                    # write the line into the output file
                    out_file.write(w_line + "\n")
                    self._update_progress(files, file_sizes)

    def realised_volatility(self, file_path, cols_rvs, cols_lobs, cols_sums, out_path):
        f = stack.enter_context(open(file_path))
        f_size = os.stat(f.fileno()).st_size
        eof = False
        timestamp = None

        headers = self._read_line(f)
        idx_rvs = [i for i, x in enumerate(headers) if x in cols_rvs]
        idx_lobs = [i for i, x in enumerate(headers) if x in cols_lobs]
        idx_sums = [i for i, x in enumerate(headers) if x in cols_sums]

        rvs = [0] * len(idx_rvs)
        rvs_last_obs = False
        lobs = [0] * len(idx_lobs)
        sums = [0] * len(idx_sums)

        print("\n\n# Realised Volatility Aggregation #\n")

        with open(out_path, "w") as out_file:
            out_header = "timestamp," + \
                ",".join([x + '_rv' for x in cols_rvs]) + "," + \
                ",".join([x + '_lobs' for x in cols_lobs]) + "," + \
                ",".join([x + '_sum' for x in cols_sums]) + "\n"
            out_file.write(out_header)

            while not eof:
                line = self._read_line(f)

                if not line:
                    eof = True
                else:
                    new_timestamp = line[0][0:10]

                # first pass only
                if timestamp is None:
                    timestamp = new_timestamp

                # when there is new timestamp, write an observation and reset variables
                if new_timestamp != timestamp or not line:
                    # write to out_file, only when there were any observations
                    if not all(x == 0 for x in sums):
                        w_line = [timestamp] + rvs + lobs + sums
                        self._write_line(w_line, out_file)
                        self._update_progress([f], [f_size])
                    # reset variables
                    rvs = [0] * len(idx_rvs)
                    rvs_last_obs = False
                    lobs = [0] * len(idx_lobs)
                    sums = [0] * len(idx_sums)
                    timestamp = new_timestamp

                if not eof:
                    # get realised variances
                    rvs_curr_obs = [line[i] for i in idx_rvs]
                    if rvs_last_obs:
                        rvs_add = [(float(x) - float(y)) ** 2
                                   for x, y in zip(rvs_curr_obs, rvs_last_obs)]
                        rvs = [sum(x) for x in zip(rvs, rvs_add)]
                    rvs_last_obs = rvs_curr_obs
                    # get last observations
                    lobs = [line[i] for i in idx_lobs]
                    # get sums
                    for i, idx in enumerate(idx_sums):
                        sums[i] += int(line[idx])

    def subset_data(self, file_path, out_path):
        """
        Removes public holidays, weekends and incomplete observations
        from the dataset.

        Keyword arguments:
        file_path: path to the original dataset
        out_path: path where to save the subset of the original dataset
        """

        print("\n\n# Reading the Original Dataset #\n")
        df = pd.read_csv(file_path, index_col="timestamp")

        # remove US federal holidays from the dataset
        print("# Removing Federal Holidays #\n")
        holidays = pd.read_csv("usholidays.csv", usecols=["Date"],
                               squeeze=True).values
        df = df.loc[df.index.difference(holidays)]

        # remove weekends from the dataset
        print("# Removing Weekends #\n")

        def is_weekday(x):
            dow = int(dt.datetime.strptime(x, '%Y-%m-%d').strftime('%w'))
            return dow not in {0, 6}

        df = df.loc[df.index.map(is_weekday)]

        print("# Saving the Result #\n")
        df.to_csv(out_path)

    def lambda_optim(self, file_path):
        # maturities to test
        t_maturities = [2, 5, 10, 25]

        # get the lambda values to test
        t_lambdas = pd.read_csv("testlambdas.csv", index_col=["maturity"])

        # read the final dataset
        df = pd.read_csv(file_path, index_col="timestamp")
        df = df[["TU_yield_lobs", "FV_yield_lobs",
                 "TY_yield_lobs", "US_yield_lobs"]]

        # computes square error of a regression
        def compute_sq_error(x, X_ns):
            y_ns = x.values
            coef_ns = np.linalg.lstsq(X_ns, y_ns, rcond=None)[0]
            resid = y_ns - np.dot(X_ns, coef_ns)
            ssr = np.sum(resid**2)
            return ssr

        with open("lambda_mse.csv", "w") as out_file:
            out_file.write("maturity,lambda,mse\n")
            # go through each maturity-lambda pair and compute RMSE
            for m, lmbda in t_lambdas.itertuples():
                # X matrix for the Nelson-Siegel factor regression
                X_ns = np.array(
                    [[1, (1 - exp(-lmbda * m)) / (lmbda * m), (1 - exp(-lmbda * m)) / (lmbda * m) - exp(-lmbda * m)]
                     for m in t_maturities]
                )
                sq_errors = df.apply(func=compute_sq_error,
                                     axis=1, args=(X_ns,))
                mse = sq_errors.sum() / (df.shape[0] * 4)
                out_file.write(str(m) + "," + str(lmbda) +
                               "," + str(mse) + "\n")
                print(f"maturity: {m}, lambda: {lmbda}, mse: {mse}")


if __name__ == "__main__":
    with contextlib.ExitStack() as stack:
        dp = DataProcessor()

        dp.previous_tick(
            file_paths=[
                "tickTU.csv",
                "tickFV.csv",
                "tickUS.csv",
                "tickTY.csv"
            ],
            skip_headers=True,
            delta=300,
            regime=("07:20", "14:00")
        )

        dp.nelson_siegel(
            series=[
                {
                    "path": "tickTU_pt300s.csv",
                    "ticker": "TU",
                    "maturity": 2,
                    "fv": 200000,
                    "conv_factor": 0.8885,
                    "contr_mult": 2000
                },
                {
                    "path": "tickFV_pt300s.csv",
                    "ticker": "FV",
                    "maturity": 5,
                    "fv": 100000,
                    "conv_factor": 0.7441,
                    "contr_mult": 1000
                },
                {
                    "path": "tickTY_pt300s.csv",
                    "ticker": "TY",
                    "maturity": 10,
                    "fv": 100000,
                    "conv_factor": 0.5537,
                    "contr_mult": 1000
                },
                {
                    "path": "tickUS_pt300s.csv",
                    "ticker": "US",
                    "maturity": 25,
                    "fv": 100000,
                    "conv_factor": 0.2281,
                    "contr_mult": 1000
                }
            ],
            # lmbda=0.0609,  # used by Diebold and Li 
            lmbda=0.6329,  # maximises curvature at 2.83 years which minimises mse of resulting fits
            out_path="nelson_siegel-merged.csv"
        )

        dp.realised_volatility(
            file_path="nelson_siegel-merged.csv",
            cols_rvs=["TU_yield", "FV_yield", "TY_yield", "US_yield",
                      "ns_beta0", "ns_beta1", "ns_beta2"],
            cols_lobs=["TU_close", "TU_yield",
                       "FV_close", "FV_yield",
                       "TY_close", "TY_yield",
                       "US_close", "US_yield",
                       "ns_beta0", "ns_beta1", "ns_beta2"],
            cols_sums=["TU_nobs", "FV_nobs", "TY_nobs", "US_nobs"],
            out_path="realised_volatility.csv"
        )

        dp.subset_data(
            file_path="realised_volatility.csv",
            out_path="realised_volatility-subset.csv"
        )

        # dp.lambda_optim(
        #     file_path="realised_volatility-subset.csv"
        # )

        del dp
