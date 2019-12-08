import os, io, random
import string
import numpy as np
import pandas as pd
import time
import faker
import collections

from Bio.Seq import Seq
from Bio import AlignIO, SeqIO

from bokeh.plotting import figure, save, output_file
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, Span, Panel, Tabs
from bokeh.models.glyphs import Text, Rect
from bokeh.models.widgets import Div, Button
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn
from bokeh.models.callbacks import CustomJS
from bokeh.events import ButtonClick

from bokeh.layouts import grid
from bokeh.embed import file_html
from bokeh.plotting import show
from bokeh.resources import CDN

import unittest

class SimulateSequenceData():
    """ makes simulated sequence data for testing """
    def make_seq(self, alignment_width=40):    
        return ''.join([random.choice(['A','C','T','G']) for i in range(alignment_width)])

    def mutate_seq(self, seq, iupac=False):
        """
        Mutate a sequence randomly.  
        If iupac is true, includes all iupac codes.  
        Otherwise, just includes A,C,G,T,M,N.
        """
        if not iupac:
           options = ['A','C','T','G','M','N']
        else:
           options = ['A','C','T','G','R','Y','S','M','K','B','D','H','V','N']

        seq = list(seq)
        pos = np.random.randint(1,len(seq),6)    
        for i in pos:
            seq[i] = random.choice(options)
        return ''.join(seq)

    def make_msa(self, nSeqs= 10, alignment_width = 40):

        f = faker.Faker('en_UK')
        seq = self.make_seq(alignment_width)
        msa_dict = {}
        for seq_id in range(nSeqs):
            seq_name = "seq_{0}".format(seq_id)
            msa_dict[seq_name]= {'aligned_seq':self.mutate_seq(seq, iupac=True)}
            msa_dict[seq_name]['aligned_mseq'] = msa_dict[seq_name]['aligned_seq']
            msa_dict[seq_name]['alignN'] = sum([bool(x=='N') for x in msa_dict[seq_name]['aligned_seq']])
            msa_dict[seq_name]['alignM'] = sum([bool(x in ['R','Y','S','M','K']) for x in msa_dict[seq_name]['aligned_seq']])
            msa_dict[seq_name]['alignN_or_M'] = msa_dict[seq_name]['alignN']+msa_dict[seq_name]['alignM']
            msa_dict[seq_name]['allN'] = np.random.binomial(4.4e6, 1e-3)
            msa_dict[seq_name]['allM'] = np.random.binomial(4.4e6, 1e-4)
            msa_dict[seq_name]['allN_or_M'] = msa_dict[seq_name]['allN']+msa_dict[seq_name]['allM']
            msa_dict[seq_name]['Surname'] = f.last_name()
            msa_dict[seq_name]['Forename'] = f.first_name()
            msa_dict[seq_name]['Patient_id'] = random.randint(1,1e9)
            msa_dict[seq_name]['mlp_value1'] = random.randint(0,8)
            msa_dict[seq_name]['mlp_value2'] = random.randint(0,8)
            msa_dict[seq_name]['mlp_value3'] = random.randint(0,8)
            msa_dict[seq_name]['mlp_value4'] = random.randint(0,8)
            msa_dict[seq_name]['p_value1'] = 10^-msa_dict[seq_name]['mlp_value1']
            msa_dict[seq_name]['p_value2'] = 10^-msa_dict[seq_name]['mlp_value3']
            msa_dict[seq_name]['p_value3'] = 10^-msa_dict[seq_name]['mlp_value3']
            msa_dict[seq_name]['p_value4'] = 10^-msa_dict[seq_name]['mlp_value4']

        return pd.DataFrame.from_dict(msa_dict, orient= 'index')
            


class DepictMSA():
       """ depicts a multi sequence alignment """
       def __init__(self, msar_df, iupac=None, positions_analysed=None, identify_sequence_by=None):
           """ loads an MSA result object's dataframe.

               parameters:
               msar_df: a pandas DataFrame containing the multisequence alignment.
               The output of SimulateSequenceData.make_msa() has the correct format, as does the 
               .df property of a MultiSequencingAlignmentResult object.
               iupac:  what kind of characters reflect mixtures
                       if True, assumes IUPAC codes encode mixtures https://www.bioinformatics.org/sms/iupac.html
                       if False, expects mixtures to be represented by an M
                       if None, autodetects.
               positions_analysed: a list, the same length as the multisequence alignment, which is used as labels 
                       in the MSA.  Typically, the MSA consists of positions varying between sequences - so position 1
                       of the MSA may correspond to position 10002 of the underlying reference.  positions_analysed
                       allows data about which positions, genes, codons etc are contributing to the variation.
                       If none, the labels are an integer series increasing from 1.
               identify_sequence_by: a list of columns in msar_df, or None.  If not None, the data is presented ordere
                        by identify_sequence_by and the values in identify_sequence_by are included in the sequence identifier.
               returns:
                        None
                        
               errors:
                        ValueError, if sequences contain non-iupac characters"""
           self.df = msar_df
           self.iupac = iupac
           self.positions_analysed = positions_analysed
	   
           # todo: qc on data passed

           # note: by default, the index will become the y-labels.  To display in a different order, reindex.
           self.df['msa_id']=self.df.index.tolist()
           if identify_sequence_by is not None:
               # check if any identify_sequence_by fields which are not in the column headers.
               missing_columns = set(identify_sequence_by)-set(self.df.columns.tolist())
               if len(missing_columns)>0:
                   raise ValueError("asked to identify sequences by missing fields {0}".format(missing_columns))
               # construct a new index
               new_identifiers = []

               for ix in self.df.index:
                   new_id = []
                   for col in identify_sequence_by:
                       new_id.append(self.df.at[ix, col])
                   new_id.append('[{0}]'.format(ix))
                   new_identifiers.append(" ".join(new_id))

               self.df['msa_id']=new_identifiers
               self.df.sort_values(by='msa_id', inplace=True, ascending=False)
               self.df.set_index('msa_id', drop=True, inplace=True)
           # sort

           # ensure sequences are upper case
           self.df['aligned_mseq'] = [x.upper() for x in self.df['aligned_mseq'].tolist()]        # ensure uppercase

           # store number of sequences and alignment width
           self.nSeqs = len(self.df.index)
           if self.nSeqs==0:     # no data
               self.align_width=0
           else:
               first_seq_name = self.df.index.tolist()[0] 
               self.align_width = len(self.df.loc[first_seq_name, 'aligned_mseq'])
        
           # compute composition of sequences
           self.composition = collections.Counter()
           for seq in self.df['aligned_mseq'].tolist():
              for nucl in list(seq):
                   self.composition.update(nucl)

           # compute composition of sequences
           self.compositions = {}
           for ix in self.df.index:
               seq = self.df.at[ix,'aligned_seq']
               self.compositions[ix]=collections.Counter()
               for nucl in list(seq):
                   self.compositions[ix].update(nucl)
           self.compositions = pd.DataFrame.from_dict(self.compositions, orient='index')
           self.compositions['msa_id'] = self.compositions.index
           self.compositions.fillna(0, inplace=True)
           # check iupac character composition.  The keys of the valid_characters dictionary
           # reflect whether all iupac characters are expected.
           valid_characters = {False:set(['A','C','G','T','M','N','-']),
                                True:set(['A','C','G','T','M','N','-','R','Y','S','W','K','B','D','H','V'])
                              }
           if self.iupac is None:
               if len(self.composition.keys() - valid_characters[False]) > 0:            # there are iupac, or other disallowed characters in the sequences
                   self.iupac = True
               else:
                   self.iupac = False

           # check there are no invalid characters
           if self.iupac:
               illegal_characters = self.composition.keys() - valid_characters[True]
               if len(illegal_characters) > 0:               # illegal characters
                   raise ValueError("Illegal characters present in sequences {0}".format(illegal_characters))

           # define what the iupac code mean https://www.bioinformatics.org/sms/iupac.html
           base2base = {
                     'A':['A'],
                     'C':['C'],
                     'T':['T'],
                     'G':['G'],
                     '-':['-'], 
                     'N':['N'],
                     'R':['A','G'],
                     'Y':['C','T'],
                     'S':['G','C'],
                     'W':['A','T'],
                     'K':['G','T'],
                     'M':['A','C'],
                     'K':['G','T'],
                     'B':['C','G','T'],
                     'D':['A','G','T'],
                     'H':['A','C','T'],
                     'V':['A','C','G']
                     }

           if self.iupac is False:
               # M reflects any kind of mixture
               base2base['M'] = 'M'

           clrs =  {'A':'green','T':'red','G':'black','C':'lightblue','-':'white','N':'white', 'M':'yellow'}
       
           # build a data frame describing the rectangles which will be drawn for each type of sequence
           ix = 0
           rectangle_map = {}
           for nucl in valid_characters[self.iupac]:
               n_mapped_nucleotides = len(base2base[nucl])      # how many nucleotides comprise the mixture, if any
               for i,mapped_to in enumerate(base2base[nucl]):
                   ix +=1
                   rectangle_map[ix] = {'base_y':i/n_mapped_nucleotides,
                                        'height':1/n_mapped_nucleotides,
                                        'nucl':nucl, 
                                        'component':mapped_to, 
                                        'full_fill_colour':clrs[mapped_to],
                                        'limited_fill_colour':clrs[mapped_to],                                        
                                        'is_mixed':len(base2base[nucl])>1
                                        }
                   if len(base2base[nucl])>1:
                          rectangle_map[ix]['limited_fill_colour']='yellow'     # all mixed positions yellow
                   rectangle_map[ix]['fill_colour']=rectangle_map[ix]['limited_fill_colour']     # all mixed positions yellow
 
           self.rectangles = pd.DataFrame.from_dict(rectangle_map, orient='index')

           # compute fasta
           fasta_str = []
           for ix in self.df.index:
               fasta_str.append(">"+ix)
               fasta_str.append(self.df.at[ix,'aligned_mseq'])
           self.fasta = "<br>".join(fasta_str)

       def render_msa(self, cluster_name=None):
           """ depict msa """
           # prepare view
           if cluster_name is None:
               cluster_name = ""
               
           info = Div(text="Cluster {0} with {1} variant bases in alignment".format(cluster_name, self.align_width))
           
           # fasta
           div = Div(text=self.fasta)
           p1 = Panel(child=div, title = "Fasta")


           # table of compositions
           comp_source = ColumnDataSource(self.compositions)
           target_columns = self.compositions.columns.to_list()
           columns = []
           for item in target_columns:
               columns.append(TableColumn(field=item, title=item))
           data_table = DataTable(source=comp_source, columns=columns, width=1000, height =500)
           pC = Panel(child=data_table, title = 'Composition')

    
           # table of base counts and and -logp values
           grid_source = ColumnDataSource(self.df)
           target_columns = ['Surname','Forename','Patient_id','msa_id','mlp_value1','mlp_value2','mlp_value3','mlp_value4','alignN','allN','alignM','allM','alignN_or_M','allN_or_M']
           columns = []
           for item in target_columns:
               if item in self.df.columns:
                    columns.append(TableColumn(field=item, title=item))
           data_table = DataTable(source=grid_source, columns=columns, width=1000, height =500)
           p0 = Panel(child=data_table, title = 'Summary')
           
           # grid
           
           # https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner
 
           # create a one-base per row data frame laying out the positions of each character in the msa
           x = np.arange(0,self.align_width)
           y = np.arange(0,self.nSeqs)      #self.df.index.tolist()
           xx, yy = np.meshgrid(x,y)
           gx = xx.ravel()
           gy = yy.flatten()
           h = 1/self.nSeqs
           nucl = [i for s in self.df['aligned_mseq'].tolist() for i in s]

           # create a dataframe and link to colours           
           perBase = pd.DataFrame.from_dict({'x':gx, 'y':gy, 'nucl':nucl})
           to_depict = perBase.merge(self.rectangles, how='left', left_on='nucl', right_on='nucl')
           no_color_assigned = to_depict[to_depict['full_fill_colour'].isnull()]

           # sanity check: everything should be mapped to a colour.
           if len(no_color_assigned.index)>0:
               raise KeyError("Not all cells were mapped to a colour; {0}".format(no_color_assigned))

           # apply offsets as the centre, not the LH corner, of the rectangle is required
           to_depict['base_y'] = to_depict['y']+to_depict['base_y']
           to_depict['centre_y'] = to_depict['base_y']  + (to_depict['height']/2)
           to_depict['centre_x'] = to_depict['x']       +0.5

           # construct column data source for plotting for plotting
           # we plot on a linear axis but apply text tick labels
           source = ColumnDataSource(to_depict)
           plot_height = (self.nSeqs*15)+  50
           x_range = Range1d(0,10+self.align_width, bounds='auto')
           y_ticks = [x+0.5 for x in range(self.nSeqs)]
           y_tick_relabelling = dict(zip(y_ticks, self.df.index.tolist()))
           x_ticks = [x+0.5 for x in range(self.align_width)]
           default_x_labels = [str(x) for x in range(self.align_width)]
           if self.positions_analysed is not None:
               x_tick_relabelling = dict(zip(x_ticks, [str(x) for x in sorted(list(self.positions_analysed))]))
           else:
               x_tick_relabelling = dict(zip(x_ticks, default_x_labels))

           # configure view length
           if self.align_width > 100:
               viewlen =0
           else:
               viewlen = self.align_width

           view_range = (0,viewlen)
           tools = "xpan, xwheel_zoom, reset, save"
           
           # construct figure, add rectangles
           p = figure(title=None, plot_width= 800, plot_height=600, 
                      tools=tools, min_border =0, toolbar_location='below',
                      x_range= x_range)
           rects = Rect(x='centre_x',y='centre_y',line_width=0, 
                        fill_color= 'fill_colour', width=1, height= 'height', 
                        fill_alpha=1)
           p.add_glyph(source, rects)
           p.grid.visible = False

           # configure y-axis labels
           p.yaxis.ticker = y_ticks
           p.yaxis.major_label_overrides = y_tick_relabelling
           
           # configure x-axis labels
           if self.positions_analysed is not None:
               p.xaxis.ticker = x_ticks
               p.xaxis.major_label_overrides = x_tick_relabelling
           p.xaxis.major_label_orientation = "vertical"

           # lines to divide up sequences
           for division in y_ticks:
               divide_sequence = Span(location=division-0.5, dimension='width', line_color='white', line_width=0.5)
               p.add_layout(divide_sequence)

           p.add_glyph(source, rects)

           callback1 = CustomJS(args=dict(source=source),
                                code=""" var data = source.data;
           console.log('Clicked callback 1')
           x = data['full_fill_colour'];
           y = data['fill_colour'];
           for (i = 0; i < x.length; i++) {
              y[i] = x[i];
           }
           source.change.emit();
           """)
           callback2 = CustomJS(args=dict(source=source),
                                code=""" var data = source.data;
           console.log('Clicked callback 1')
           x = data['limited_fill_colour'];
           y = data['fill_colour'];
           for (i = 0; i < x.length; i++) {
              y[i] = x[i];
           }
           source.change.emit();
           """)
           button1 = Button(label="Show mixed bases' components", button_type="success")
           button2 = Button(label="Show mixed bases in yellow", button_type="success")
                     
           button1.js_on_event(ButtonClick, callback1)
           button2.js_on_event(ButtonClick, callback2)
 
           p2 = Panel(child = grid([[button1,button2],p]), title='Image')
           # add text stating name of cluster and alignment width
           g = grid([info,Tabs(tabs=[p0, pC, p1, p2 ]) ])
           #doc = Tabs(tabs=[p0, p1, p2 ])
           html = file_html(g, CDN, "Multisequence alignment")

           return html

class Test_SimulateSequenceData_1(unittest.TestCase):
    """ tests make_seq """
    def runTest(self):
        ssd = SimulateSequenceData()
        
        seq = ssd.make_seq(30)
        self.assertIsInstance(seq , str)
        self.assertEqual(30, len(seq))  
class Test_SimulateSequenceData_2(unittest.TestCase):
    """ tests mutate_seq """
    def runTest(self):
        ssd = SimulateSequenceData()
        
        seq = ssd.make_seq(100)
        mseq = ssd.mutate_seq(seq)

        self.assertIsInstance(mseq , str)
        self.assertEqual(100, len(mseq))

        seq = ssd.make_seq(100)
        mseq = ssd.mutate_seq(seq, iupac=True)

        self.assertIsInstance(mseq , str)
        self.assertEqual(100, len(mseq))
class Test_SimulateSequenceData_3(unittest.TestCase):
    """ tests msa generation """
    def runTest(self):
        ssd = SimulateSequenceData()
        msa = ssd.make_msa(10,40)
        self.assertIsInstance(msa, pd.DataFrame)
        self.assertEqual(len(msa.index), 10)

class Test_DepictMSA_1(unittest.TestCase):
    """ tests msa object """
    def runTest(self):
        ssd = SimulateSequenceData()
        nSeqs = 20
        alignLen = 200

        x_labels = ["pos{0}".format(x) for x in range(alignLen)]

        msa = ssd.make_msa(nSeqs, alignLen)
        dep_msa = DepictMSA(msa)
        self.assertEqual(dep_msa.align_width,alignLen)
        self.assertEqual(dep_msa.nSeqs, nSeqs)
        self.assertIsInstance(dep_msa.composition, collections.Counter)
        self.assertIsInstance(dep_msa.rectangles, pd.DataFrame)
        retVal = dep_msa.render_msa()

        self.assertIsInstance(retVal, str)

        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)

        msa = ssd.make_msa(nSeqs, alignLen)
        dep_msa = DepictMSA(msa)

        self.assertEqual(dep_msa.align_width,alignLen)
        self.assertEqual(dep_msa.nSeqs, nSeqs)
        self.assertIsInstance(dep_msa.composition, collections.Counter)
        self.assertIsInstance(dep_msa.rectangles, pd.DataFrame)
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)

        dep_msa = DepictMSA(msa, positions_analysed  = x_labels)
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)

        with self.assertRaises(ValueError):
            dep_msa = DepictMSA(msa, identify_sequence_by = ['no field'])
            retVal = dep_msa.render_msa()
            self.assertIsInstance(retVal, str)

        dep_msa = DepictMSA(msa, identify_sequence_by = ['Surname','Forename'])
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)
        self.assertTrue('</html>' in retVal)
        

