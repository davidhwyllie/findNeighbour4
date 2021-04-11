""" unit tests for depictStatus """
import os
import json
import unittest
from bokeh.io import output_file
from bokeh.models.widgets import Div, Panel, Tabs
from findn.depictStatus import DepictServerStatus

class test_init(unittest.TestCase):
    """ tests init method of DepictServerStatus """
    def runTest(self):
        """ tests init """
        
        DepictServerStatus()

class test_logfile(unittest.TestCase):
    def runTest(self):
        """ tests tailing of a file """
        
        # no parameters except SNV threshold
        dss = DepictServerStatus()
        self.assertEqual(dss.logfile_tail(None), 'No Log file was specified. ')

        dss = DepictServerStatus(logfiles={'nil':'nonexisting_file.log'})
        self.assertEqual(dss.logfile_tail('nil'), 'No log file exists.  The specified nil log file was nonexisting_file.log')

        dss = DepictServerStatus(logfiles={'test':os.path.join( "testdata","monitoring","logfile_small.log")})        
        res = dss.logfile_tail(process ='test')        
        self.assertTrue("2018-10-30 12:02:09,964" in res)    

        dss = DepictServerStatus(logfiles={'test':os.path.join( "testdata","monitoring","logfile_big.log")})
        res = dss.logfile_tail(process = 'test')
        self.assertTrue("2018-11-02 10:33:41,831" in res)    
 
class test_depict_1(unittest.TestCase):
    def runTest(self):
        """ tests depiction """
        inputfile = os.path.join( "testdata","monitoring","m50.json")
        with open(inputfile,'rt') as f:
            res = json.load(f)
            
        logfiles = {'big':os.path.join( "testdata","monitoring","logfile_big.log"),
                    'small':os.path.join( "testdata","monitoring","logfile_small.log")}
        
        dss = DepictServerStatus(logfiles= logfiles)
        dss.read_json(res, data_tag='all')
        tab_g2n_all = dss.depict(data_tag='all', tab_title="Guid2neighbour", metrics='dstats|guid2neighbour',  x_axis_label = 'Order sequences added (Oldest --> Most recent)')

        dss.subset_data(from_data_tag='all', to_data_tag='pre_insert', column_name="context|info|message", cell_values=["About to insert"])
        tab_server = dss.depict(data_tag='pre_insert', tab_title="In RAM sequence", metrics='server|pcstat',  x_axis_label = 'Order sequences added (Oldest --> Most recent)')
        tab_memory = dss.depict(data_tag='pre_insert', tab_title="RAM usage", metrics='server|mstat', x_axis_label = 'Order sequences added (Oldest --> Most recent)')
        tab_g2n = dss.depict(data_tag='pre_insert', tab_title="Db: guid->neighbours", metrics='dstats|guid2neighbour',  x_axis_label = 'Order sequences added (Oldest --> Most recent)')

        tab_g2m = dss.depict(data_tag='pre_insert', tab_title="Db: guid->metadata", metrics='dstats|guid2meta',  x_axis_label = 'Order sequences added (Oldest --> Most recent)')
        tab_sm = dss.depict(data_tag='pre_insert', tab_title="Db: server monitor", metrics='dstats|server_monitoring',  x_axis_label = 'Order sequences added (Oldest --> Most recent)')
           
        # get details of the most recent 5 guids
        n=5
        recent_guids = dss.most_recent_guids('all',n=n)
        dss.subset_data(from_data_tag='all', to_data_tag='recent_guids', column_name='content|activity|guid', cell_values=[x for x in recent_guids])
        tab_rserver = dss.depict(data_tag='recent_guids', tab_title="Last {0} RAM sequences".format(n), metrics='server|pcstat')
        tab_rmemory = dss.depict(data_tag='recent_guids', tab_title="Last {0} RAM usage".format(n), metrics='server|mstat')
        
        # get tail of logfile
        n_latest_lines = 200

        # render
        tab_list = [tab_server, tab_memory, tab_g2n, tab_g2n_all, tab_g2m, tab_sm, tab_rserver, tab_rmemory]
        for process in logfiles.keys():
            res = dss.logfile_tail(process,n_latest_lines)
            div = Div(text= "[last {0} lines of log file are shown]<br/>".format(n_latest_lines) + res.replace('\n','<br />'), render_as_text=False, width=1000, height=800)
            tab_list.append( Panel(child=div, title='{0} log'.format(process))) 


        doc = Tabs(tabs=tab_list)
        self.assertTrue(doc is not None)
        #show(doc)

class test_depict_2(unittest.TestCase):
    def runTest(self):
        """ tests depiction of guid2neighbour size"""
        inputfile = os.path.join( "testdata","monitoring","m50.json")
        with open(inputfile,'rt') as f:
            res = json.load(f)
            
        logfiles = {'big':os.path.join( "testdata","monitoring","logfile_big.log"),
                    'small':os.path.join( "testdata","monitoring","logfile_small.log")}
        
        dss = DepictServerStatus(logfiles= logfiles)
        dss.read_json(res, data_tag='all')
        tab_g2n_all = dss.depict(data_tag='all', tab_title="Db: guid->neighbours", metrics='dstats|guid2neighbour',  x_axis_label = 'Order sequences added (Oldest --> Most recent)')
 
          
        
        # render
        tab_list = [tab_g2n_all]
        
        doc = Tabs(tabs=tab_list)
        self.assertTrue(doc is not None)
        #show(doc)

class test_tail(unittest.TestCase):
    def runTest(self):
        """ tests tailing of a file """
        
        # no parameters except SNV threshold
        dss = DepictServerStatus(logfiles = {'test':os.path.join("testdata","monitoring","logfile_big.log")})
        res = dss.logfile_tail('test',100)
        self.assertTrue("2018-11-02 10:33:41,831" in res)    
 
        output_file("div.html")

        div = Div(text= res.replace('\n','<br />'), render_as_text=False, width=1000, height=800)
        tab1 = Panel(child=div, title='Log tail')        
        doc = Tabs(tabs=[ tab1 ])

        self.assertTrue(doc is not None)
        #show(doc)