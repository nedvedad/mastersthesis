import os
import sys
import time
import contextlib
import re
import random
import datetime as dt

# input file paths
FILES = [
    "/Volumes/flashstorage/thesis_data/tickTU.csv",
    "/Volumes/flashstorage/thesis_data/tickFV.csv",
    "/Volumes/flashstorage/thesis_data/tickTY.csv",
    "/Volumes/flashstorage/thesis_data/tickUS.csv"
]
# output file path
OUT_FILE = "refresh-time.csv"
# a function returning a statistic of the subset of the data being synchronised
STATISTIC = random.choice
# number of seconds between equally spaced observations
DELTA = 60


class TimeSeriesSyncr:
    """
    Keyword arguments:
    files: A list of File objects, containing the time series
           to be synchronised
    """

    def __init__(self, files):
        """
        Initialises attributes and starts the counter.
        """
        self.files = files
        self.files_size = [os.stat(f.fileno()).st_size for f in self.files]
        self.start_time = time.clock()

    def __del__(self):
        print("\nFinished in {} seconds.".format(time.clock()-self.start_time))

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
    def _parse_line(line):
        """
        Parses a single csv line, returns a (timestamp, value) tuple.
        """
        pattern = "^([\d/]+,[\d:.]+),([\d.]+)"
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
        # rounding = seconds // delta * delta
        return t + dt.timedelta(0, rounding-seconds, -t.microsecond)

    @staticmethod
    def _write_line(x, file):
        """
        Concats list x into a csv line and writes it into the file.
        """
        line = ",".join(map(str, x))+"\n"
        file.write(line)

    def _update_progress(self):
        """
        Updates the progress in the console.
        """
        print("", end="\r")
        for i, f in enumerate(self.files):
            perc = round(f.tell()/self.files_size[i]*100, 2)
            print("File #{}: {}%".format(i+1, perc), end=" ")
            sys.stdout.flush()

    def skip_lines(self, n):
        """
        Skips n lines in each of the files.
        """
        for f in self.files:
            for _ in range(n):
                f.readline()

    def previous_tick(self, delta, out_path):
        """
        Implementation of Previous Tick synchronisation scheme as used by
        Zhang, L. in Estimating covariation: Epps effect, microstructure
        noise (2011).

        Synchronises multiple time series from files object using the
        Previous Tick scheme. The last known observation before each timestamp
        is used. Observations are equally spaced with delta seconds inbetween
        them.

        Keyword arguments:
        delta: The length of a period between two observations in seconds,
        out_path: A path to the output file
        """
        # get the shared initial timestamp t1
        obs = [self._peek_line(f) for f in self.files]
        obs = [self._parse_line(l) for l in obs]
        t1 = max([o[0] for o in obs])
        t1 = self._round_datetime_up(t1, delta)

        for i, f in enumerate(self.files):
            eof = False
            new_timestamp = t1
            nobs = 0
            path = os.path.splitext(out_path)
            path = path[0] + str(i).zfill(2) + path[1]
            with open(path, "w") as out_file:
                while not eof:
                    line = f.readline()
                    if line:
                        (timestamp, value) = self._parse_line(line)
                        # imputing the missing observations and writing results
                        while new_timestamp <= timestamp:
                            obs = [new_timestamp, last_value, nobs]
                            self._write_line(obs, out_file)
                            # update timestamp and progress
                            new_timestamp += dt.timedelta(0, delta)
                            if nobs > 0:
                                self._update_progress()
                            nobs = 0
                        # the last observation in the period not reached yet
                        last_value = value
                        nobs += 1
                    else:
                        eof = True
                        obs = [new_timestamp, last_value, nobs]
                        self._write_line(obs, out_file)
                        self._update_progress()

    def gst(self, fun, out_path):
        """
        Implementation of Generalised Sampling Time proposed
        by Aït-Sahalia, Y., Fan, J., and Xiu, D. in High-Frequency Covariance
        Estimates With Noisy and Asynchronous Financial Data (2010).

        Synchronises multiple time series from files object using the
        Generalised Sampling Time with fun as the statistic of the
        data. The output is written to a file specified in out_path.

        Keyword arguments:
        fun: A function returning a desired statistic of the subset of the time
             series being synchronised.
        out_path: A path to the output file.
        """
        eof = False  # end of file indicator

        with open(out_path, "w") as out_file:
            while not eof:
                obs = []
                # new timestamp is the maximum of the next observations
                lines = [f.readline() for f in self.files]
                new_obs = [self._parse_line(l) for l in lines]
                new_timestamp = max([o[0] for o in new_obs])

                for i, f in enumerate(self.files):
                    # save the first observation
                    values = [new_obs[i][1]]

                    # get all observations before the new timestamp
                    while True:
                        line = f.readline()
                        if line:
                            (timestamp, value) = self._parse_line(line)
                            if timestamp <= new_timestamp:
                                values.append(value)
                            else:
                                # return the file's pointer to the prev. line
                                f.seek(f.tell()-len(line)-1)
                                break
                        else:
                            eof = True
                            break

                    # append the new obs computed using the desired statistic
                    if not eof:
                        obs.append(fun(values))

                # write output
                if not eof:
                    self._write_line([new_timestamp] + obs, out_file)

                self._update_progress()


if __name__ == "__main__":
    with contextlib.ExitStack() as stack:
        files = [stack.enter_context(open(file_path)) for file_path in FILES]
        ts = TimeSeriesSyncr(files)
        ts.skip_lines(1)  # skip the headers
        ts.gst(lambda x: x[-1], OUT_FILE)
        # ts.previous_tick(delta=DELTA, out_path=OUT_FILE)
        del ts
