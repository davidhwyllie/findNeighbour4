""" produces a human readable, interactive HTML/js depiction of findneighbour4 server status

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

# import libraries
import os
import logging
import pandas as pd
import numpy as np
import dateutil.parser
import datetime
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, CDSView, IndexFilter, DataTable, TableColumn
from bokeh.models.widgets import Div, Panel, Tabs
from bokeh.plotting import figure


class MakeHumanReadable:
    """converts server monitoring codes into human readable formats"""

    def __init__(self):
        """creats a new converter"""
        self.tag2readable = {
            "content|activity|guid": "sequence identifier",
            "server|mstat|": "fn4 memory",
            "server|catwalk|": "catwalk comparison engine ",
            "server|pcstat|n": "in-memory sequences (preComparer) ",
            "server|scstat|n": "in-memory sequences (py_seqComparer) ",
            "interval": "Interval between events (seconds) ",
            "dstats|guid2neighbour|": "db, guid2neighbour collection (stores similar sequence pairs), ",
            "dstats|guid2meta|": "db, guid2metadata collection (stores meta data about sequences), ",
            "dstats|server_monitoring|": "db, server monitoring collection (stores server state information over time), ",
            "avgObjSize": "Average record size (bytes)",
            "count": "Number of records",
            "storageSize": "Storage size used (bytes)",
            "totalIndexSize": "Total index size",
            "available": "Available RAM (bytes)",
            "free": "Free RAM (bytes)",
            "slab": "slab RAM (bytes)",
            "percent": "Percent RAM Used",
            "total": "Total system (bytes)",
            "used": "Used RAM (bytes)",
        }

    def convert(self, column_header):
        """generates human readable columns for plotting"""
        result = column_header
        for item in self.tag2readable.keys():
            if item in column_header:
                result = result.replace(item, self.tag2readable[item])
        return result


class DepictServerStatus:
    """depicts server status using data in the format cached by the server in the server_monitoring collection
    in mongodb"""

    def _set_server_info(self):
        """sets server summary  information"""
        data = [
            ["Server URL", self.server_url],
            ["Server Port", self.server_port],
            ["Server Description", self.server_description],
            ["Time report generated", datetime.datetime.now().isoformat()],
        ]
        if "recent_server" in self.monitoring_data.keys():
            if len(self.monitoring_data["recent_server"].index) > 0:
                for item in self.monitoring_data["recent_server"].index:
                    for column_name in self.monitoring_data[
                        "recent_server"
                    ].columns.values.tolist():
                        if column_name.startswith("server|"):
                            data.append(
                                [
                                    self.mhr.convert(column_name),
                                    self.monitoring_data["recent_server"].loc[
                                        item, column_name
                                    ],
                                ]
                            )
                        if column_name == "t":
                            data.append(
                                [
                                    "Note: this data is based on the most recent server status report, which was at",
                                    self.monitoring_data["recent_server"]
                                    .loc[item, column_name]
                                    .isoformat(),
                                ]
                            )
                    break
        self.monitoring_data["server_info"] = pd.DataFrame(
            columns=["Description", "Value"], data=data
        )
        self._create_cds("server_info")

    def __init__(
        self,
        logfiles=None,
        server_url="Not specified",
        server_port="Not specified",
        server_description="Not specified",
        max_interval=600,
    ):
        """initialises DepictServerStatus object

        input:
        logfile:    the server log file.  This will be read, but not written to.
        server_url: the server url.  only used to display info
        server_port:server port.  only used to display info
        server_description: server description.  only used as display info
        max_interval: don't plot intervals more than this (seconds)
        returns:
        nothing
        """

        self.logfiles = logfiles
        self.server_url = server_url
        self.server_port = server_port
        self.server_description = server_description
        self.mhr = MakeHumanReadable()
        self.source = {}
        self.source_columns = {}
        self.max_interval = max_interval
        self.monitoring_data = {}

        self._set_server_info()

    def logfile_tail(self, process, nlines=100):
        """returns the last nlines of a file logfile, if it exists

        input:
        nlines:    the number of lines to return.
        process:    one of the keys of self.logfiles: currently one of server, dbmonitor, ..

        returns:
        the tail of the file, in html format (line breaks converted to <br>)
        """

        if self.logfiles is None:
            return "No Log file was specified. "
        if not os.path.exists(self.logfiles[process]):
            return "No log file exists.  The specified {0} log file was {1}".format(
                process, self.logfiles[process]
            )
        # https://stackoverflow.com/questions/2301789/read-a-file-in-reverse-order-using-python/26747854
        try:
            with open(self.logfiles[process], "rt"):
                with open(self.logfiles[process]) as qfile:
                    qfile.seek(0, os.SEEK_END)
                    endf = position = qfile.tell()
                    linecnt = 0
                    while position >= 0:
                        qfile.seek(position)
                        next_char = qfile.read(1)
                        if next_char == "\n" and position != endf - 1:
                            linecnt += 1

                        if linecnt == nlines:
                            break
                        position -= 1

                    if position < 0:
                        qfile.seek(0)

                    return qfile.read().replace("\n", "<br />")

        except PermissionError:
            return "Log file specified cannot be accessed; Permission denied; another process may be using it."

    def read_json(self, input_json, data_tag="all"):
        """loads json data source, optionally filtering it

        input_json: json object as generated by calls to /api/v2/server_memory_usage or to
            PERSIST.recent_server_monitoring(max_reported= max_reported) where PERSIST is an
            fn3persistence object.  See below for examples.
        data_tag: a string identifying this data set

        """

        self.monitoring_data[data_tag] = pd.DataFrame.from_records(input_json)
        nrows_recovered = len(self.monitoring_data[data_tag].index)
        if nrows_recovered == 0:
            self._create_cds(data_tag)
            return

        # generate x axis and labels for tooltips
        self.monitoring_data[data_tag]["order"] = -self.monitoring_data[data_tag].index
        self.monitoring_data[data_tag]["t"] = [
            dateutil.parser.parse(x)
            for x in self.monitoring_data[data_tag]["context|time|time_now"].tolist()
        ]
        self.monitoring_data[data_tag]["_id"] = list(
            range(len(self.monitoring_data[data_tag].index))
        )
        # sort by t

        self.monitoring_data[data_tag].sort_values(
            "t", 0, ascending=False, inplace=True
        )
        self.monitoring_data[data_tag]["label_t"] = self.monitoring_data[data_tag][
            "context|time|time_now"
        ].tolist()
        self.monitoring_data[data_tag]["label_event"] = self.monitoring_data[data_tag][
            "context|info|message"
        ].tolist()
        try:
            self.monitoring_data[data_tag]["guid"] = self.monitoring_data[data_tag][
                "content|activity|guid"
            ].tolist()
        except KeyError:
            self.monitoring_data[data_tag]["guid"] = None
        # create a column data source
        self._create_cds(data_tag)

    def _create_cds(self, data_tag):
        """creates interval data & columnar data source for the data set with data_tag"""
        if data_tag not in self.monitoring_data.keys():
            # no data
            self.source[data_tag] = ColumnDataSource()
            return

        if len(self.monitoring_data[data_tag].index) > 0:

            # reindex
            self.monitoring_data[data_tag].reset_index(drop=True, inplace=True)

            # compute timings between events ('interval')
            if "t" in self.monitoring_data[data_tag].columns.tolist():
                self.monitoring_data[data_tag].loc[:, "interval"] = None
                n = 0
                previous = None
                for ix in self.monitoring_data[data_tag].index:
                    n += 1
                    if n == 1:
                        interval = None
                    else:
                        # compute interval

                        try:
                            interval = (
                                self.monitoring_data[data_tag].loc[previous, "t"]
                                - self.monitoring_data[data_tag].loc[ix, "t"]
                            ).total_seconds()
                        except TypeError:
                            interval = None
                    if interval is not None:
                        if interval > self.max_interval:
                            interval = None  # don't show very high intervals.  These normally result from the server waiting for new data
                    self.monitoring_data[data_tag].loc[ix, "interval"] = interval
                    previous = ix
        columns = []
        for item in self.monitoring_data[data_tag].columns.tolist():
            columns.append(TableColumn(field=item, title=item))
        self.source_columns[data_tag] = columns
        self.source[data_tag] = ColumnDataSource(data=self.monitoring_data[data_tag])

    def subset_data(
        self,
        from_data_tag,
        to_data_tag,
        column_name="context|info|message",
        cell_values=["About to insert"],
    ):
        """subsets data, selecting only rows containing cell_value in column_name

        input:
        from_data_tag: a data_tag for an existing data set, as loaded by read_json
        to_data_tag:   a data tag for the new dataset
        column_name:   a column name to select against
        cell_values:   a list of valid cell values.  Any cells matching these will be selected"""

        try:
            self.monitoring_data[to_data_tag] = self.monitoring_data[from_data_tag][
                self.monitoring_data[from_data_tag][column_name].isin(cell_values)
            ]
            self.monitoring_data[to_data_tag].reset_index(drop=True, inplace=True)
        except KeyError:  # fails if empty, or invalid column proposed
            pass

        # create intervals & a column data source
        self._create_cds(to_data_tag)

    def most_recent_guids(self, data_tag, n=10):
        """returns the most recent nGuids guids from the dataset defined by data_tag

        input:
        data_tag: a data_tag for an existing data set, as loaded by read_json
        n:   the number of most recent guids to select"""

        guids = set()

        for ix in self.monitoring_data[data_tag].index:
            try:
                guids.add(
                    self.monitoring_data[data_tag].loc[ix, "content|activity|guid"]
                )
                if len(guids) >= n:
                    return guids
            except KeyError:
                # fails if empty data base with no reports
                pass
        return guids

    def show_underlying_data(self, data_tag="all", tab_title="No title"):
        """returns a bokeh data table containing the source data described in self.monitoring_data[data_tag]"""

        data_table = DataTable(
            source=self.source[data_tag],
            columns=self.source_columns[data_tag],
            width=1200,
            height=800,
            editable=False,
        )
        return Panel(child=data_table, title=tab_title)

    def server_info_tab(self, tab_title="Synopsis"):
        """returns a bokeh data table containing the source data described in self.monitoring_data[data_tag]"""
        data_table = DataTable(
            source=self.source["server_info"],
            columns=self.source_columns["server_info"],
            width=1200,
            height=800,
            editable=False,
        )
        return Panel(child=data_table, title=tab_title)

    def depict(
        self,
        data_tag="all",
        tab_title="No title",
        metrics=None,
        colorpoints=None,
        x_axis_label="Order status recorded (Oldest --> Most recent)",
    ):
        """produces depiction of server activity using a json representation of the server's status.
        Depicts one data point at each point that the server inserts.

        input:

        metrics:    what to output
                    if None, reports all columns in the data supplied

                    if a string, reports columns starting with that string
                    if a list, reports those columns exactly matching the list elements
        data_tag:   an identifier for a dataset, as loaded using read_json
        tab_title:  the title for a Bokeh tab
        x_axis_label: the x axis label in the plots
        colourpoints: a categorical variable in the data frame which is
        returns:
        output:     a Bokeh tab containing graphics of server function
        """

        if data_tag not in self.monitoring_data.keys():
            # no data: a special case
            text = "No data for {0}".format(data_tag)
            div = Div(text=text, render_as_text=False, width=1000, height=800)
            tab1 = Panel(child=div, title=text)
            return tab1

        self.server_info_tab()
        # specify the tools used to plot and a container for plots
        tools = [
            "box_select",
            "hover",
            "wheel_zoom",
            "box_zoom",
            "pan",
            "save",
            "reset",
        ]
        single_plots = {"order": {}}
        multi_plots = []
        TOOLTIPS = [
            ("t", "@label_t"),
            ("event", "@label_event"),
            ("sampleId", "@guid"),
            ("interval/secs", "@interval"),
        ]

        ## select columns to output
        # if we are not supplied any metrics to output, we output them all
        if metrics is None:
            report_columns = self.monitoring_data[data_tag].columns.tolist()
        elif isinstance(metrics, str):
            # then we only report columns which start with metrics
            report_columns = []
            for item in self.monitoring_data[data_tag].columns.tolist():
                if item.startswith(metrics):
                    report_columns.append(item)
        elif isinstance(metrics, list):
            report_columns = metrics
        report_columns.insert(0, "interval")
        # plot time vs. reporting order
        item = "t"
        nnan = []
        for ix in self.monitoring_data[data_tag].index:
            try:
                if not isinstance(
                    self.monitoring_data[data_tag].loc[ix, item], np.datetime64
                ):
                    nnan.append(ix)
            except TypeError:
                pass

        view = CDSView(source=self.source[data_tag], filters=[IndexFilter(nnan)])

        single_plots["order"][item] = figure(
            plot_height=150,
            title="Time event recorded",
            plot_width=1000,
            y_axis_type="datetime",
            tools=tools,
            tooltips=TOOLTIPS,
        )
        single_plots["order"][item].xaxis.axis_label = ""
        single_plots["order"][item].yaxis.axis_label = "Time"
        single_plots["order"][item].circle(
            x="order",
            y=item,
            size=2,
            hover_color="red",
            source=self.source[data_tag],
            view=view,
        )
        multi_plots.append([single_plots["order"][item]])

        # plot everything else
        # for each column, identify numerical rows
        last_plot = None
        for item in report_columns:  # for the metrics we're supposed to start with
            if (
                item in self.monitoring_data[data_tag].columns.tolist()
            ):  # it's a valid column
                nnan = []
                for ix in self.monitoring_data[data_tag].index:
                    try:
                        if not np.isnan(self.monitoring_data[data_tag].loc[ix, item]):
                            nnan.append(ix)
                    except TypeError:
                        pass
                view = CDSView(
                    source=self.source[data_tag], filters=[IndexFilter(nnan)]
                )

                if len(nnan) > 0:
                    # plot by event
                    single_plots["order"][item] = figure(
                        plot_height=150,
                        title=self.mhr.convert(item),
                        plot_width=1000,
                        tools=tools,
                        tooltips=TOOLTIPS,
                    )
                    single_plots["order"][item].xaxis.axis_label = ""
                    single_plots["order"][item].yaxis.axis_label = ""
                    single_plots["order"][item].circle(
                        x="order",
                        y=item,
                        size=2,
                        hover_color="red",
                        source=self.source[data_tag],
                        view=view,
                    )

                    # line is not supported, see https://github.com/bokeh/bokeh/issues/9388
                    # single_plots['order'][item].line(x='order', y=item, hover_color="red", source=self.source[data_tag], view=view)

                    multi_plots.append([single_plots["order"][item]])

                    last_plot = item
                else:
                    logging.info("No data for {0}".format(item))

        if last_plot is not None:  # add x axis label
            single_plots["order"][last_plot].xaxis.axis_label = x_axis_label

        if len(multi_plots) > 0:
            tab1 = Panel(child=gridplot(multi_plots), title=tab_title)
            return tab1

        else:
            ## be clear something went wrong
            text = "? asked for metrics which are not present.<br>No data for metrics {0} -> columns {1}".format(
                metrics, report_columns
            )
            div = Div(text=text, render_as_text=False, width=1000, height=800)
        tab1 = Panel(child=div, title=tab_title)
        return tab1

    def make_report(self, all_result, recent_result, db_result):
        """makes a Bokeh report on the server's health and activity

        input:
        all_result:     result of calling PERSIST.server_monitoring
        recent_result:  as above - subsets can be different.
        db_result:      as above - selecting output of  dbManager
        returns: bokeh document
        """

        doc = {}

        # try to read json.
        self.read_json(all_result, data_tag="pre_insert")
        self.read_json(recent_result, data_tag="recent_all")
        self.read_json(db_result, data_tag="recent_db")

        if (
            len(self.monitoring_data["pre_insert"].index) == 0
            or len(self.monitoring_data["recent_all"].index) == 0
        ):
            doc["Report"] = Div(text="No data is available to plot")
            return doc

        tab_server = self.depict(
            data_tag="pre_insert",
            tab_title="In RAM sequence",
            metrics="server|pcstat",
            x_axis_label="Order sequences added (Oldest --> Most recent)",
        )
        tab_memory = self.depict(
            data_tag="pre_insert",
            tab_title="RAM usage",
            metrics="server|mstat",
            x_axis_label="Order sequences added (Oldest --> Most recent)",
        )

        tab_g2n = self.depict(
            data_tag="pre_insert",
            tab_title="Neighbours (on insert)",
            metrics="dstats|guid2neighbour",
            x_axis_label="Order sequences added (Oldest --> Most recent)",
        )
        tab_g2m = self.depict(
            data_tag="pre_insert",
            tab_title="Meta (on insert)",
            metrics="dstats|guid2meta",
            x_axis_label="Order sequences added (Oldest --> Most recent)",
        )
        tab_sm = self.depict(
            data_tag="pre_insert",
            tab_title="ServerMon",
            metrics="dstats|server_monitoring",
            x_axis_label="time (results noted every 25 repackings) (Oldest --> Most recent)",
        )

        # get details of the most recent guids
        n = 10
        recent_guids = self.most_recent_guids("recent_all", n=n)
        self.subset_data(
            from_data_tag="recent_all",
            to_data_tag="recent_server",
            column_name="content|activity|guid",
            cell_values=[x for x in recent_guids],
        )

        tab_rserver = self.depict(
            data_tag="recent_server",
            tab_title="Last {0}: Seqs".format(n),
            metrics="server|pcstat",
        )
        tab_rmemory = self.depict(
            data_tag="recent_server",
            tab_title="Last {0}: RAM".format(n),
            metrics="server|mstat",
        )

        # get details of database monitor
        tab_dbmon = self.depict(
            data_tag="recent_db",
            tab_title="Neighbours (all)",
            metrics="dstats|guid2neighbour",
        )

        self._set_server_info()
        s1 = self.server_info_tab()

        # logging
        n_latest_lines = 100
        tab_list = [
            s1,
            tab_rserver,
            tab_rmemory,
            tab_server,
            tab_memory,
            tab_dbmon,
            tab_g2n,
            tab_g2m,
            tab_sm,
        ]
        for process in self.logfiles.keys():
            res = self.logfile_tail(process, n_latest_lines)
            div = Div(
                text="[last {0} lines of log file are shown]<br/>".format(
                    n_latest_lines
                )
                + res.replace("\n", "<br />"),
                render_as_text=False,
                width=1000,
                height=800,
            )
            tab_list.append(Panel(child=div, title="Log:{0}".format(process)))

        doc["Report"] = Tabs(tabs=tab_list)
        return doc
