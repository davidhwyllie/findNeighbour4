""" provides depictions of multisequence alignments using interactive static html/javascript pages 

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
"""

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
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d, Span, Panel, Tabs, RangeTool
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
       def __init__(self, msar_df, iupac=None, positions_analysed=None, identify_sequence_by=None, max_elements_in_plot=300000):
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
               identify_sequence_by: a list of columns in msar_df, or None.  If not None, the data is presented ordered
                        by identify_sequence_by and the values in identify_sequence_by are included in the sequence identifier.
               max_elements_in_plot: the maximum size of heatmap to produce.  More than about 250k render very slowly.
               returns:
                        None
                        
               errors:
                        ValueError, if sequences contain non-iupac characters"""
           self.df = msar_df
           self.iupac = iupac
           self.positions_analysed = positions_analysed
           self.max_elements_in_plot = max_elements_in_plot
	   
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
           self.compositions.fillna(0, inplace=True)
           for compulsory_column in ['A','C','G','T']:
                if not compulsory_column in self.compositions.columns.tolist():
                    self.compositions[compulsory_column]=0
           self.df = self.df.merge(self.compositions, left_index=True,right_index=True)
           self.compositions['msa_id'] = self.compositions.index
 
           self.df.sort_values(by=['p_value3','A','C','G','T'], inplace=True, ascending=False)  # sort by mixture detection statistic
           self.df.set_index('msa_id', drop=True, inplace=True)
           self.df['y_order']=list(range(len(self.df.index)))
    
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
           line_colour_map = {True:'purple', False:'white'}
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
                                        'is_mixed':len(base2base[nucl])>1,
                                        'line_colour':line_colour_map[len(base2base[nucl])>1]
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

           # if positions_analysed is supplied, produce data for what is effectively a qqplot
           # depicting the rank of the positions analysed vs. the proportion of the genome (estimated from the highest position analysed)
           qq_dict = {}
           if self.positions_analysed is not None:
               try:
                   pa = sorted([int(x) for x in self.positions_analysed])
                   pa_is_integer = True
                   if len(pa) == 0:
                      pa_is_integer = False     # nothing to display
               except ValueError:       # not an integer; some other code used
                   pa_is_integer = False
               if pa_is_integer:  # numeric positions analysed are supplied and there's >=1 of them
                   max_pa = max(pa)
                   for i,pos in enumerate(pa):
                       qq_dict[i]={'order':i,'order_prop':i/self.align_width, 'pos':pos, 'pos_prop':pos/max_pa}
                   self.qq = pd.DataFrame.from_dict(qq_dict, orient='index')
               else:
                  self.qq = None 
           else:
               self.qq = None
       def render_msa(self, cluster_name=None):
           """ depict msa """
           # prepare view
           if cluster_name is None:
               cluster_name = ""
               
           info = Div(text="Cluster {0} with {1} sequences and {2} variant bases in alignment. ".format(cluster_name, len(self.df.index), self.align_width))
           
           # fasta
           d = Div(text=self.fasta)
           p1 = Panel(child=d, title = "Fasta")


           # table of compositions
           comp_source = ColumnDataSource(self.compositions)
           target_columns = self.compositions.columns.to_list()
           columns = []
           for item in target_columns:
               columns.append(TableColumn(field=item, title=item))
           data_table = DataTable(source=comp_source, columns=columns, width=1200, height =500)
           pC = Panel(child=data_table, title = 'Composition of variant bases')

    
           # table of base counts and and -logp values
           self.df['sample_id'] = self.df.index.tolist()
           grid_source = ColumnDataSource(self.df)
           target_columns = ['Surname','Forename','Patient_id','sample_id','p_value1','p_value2','p_value3','p_value4','alignN','allN','alignM','allM','alignN_or_M','allN_or_M']
           columns = []
           for item in target_columns:
               if item in self.df.columns:
                    columns.append(TableColumn(field=item, title=item))
           data_table = DataTable(source=grid_source, columns=columns, width=1200, height =500)
           p0 = Panel(child=data_table, title = 'Summary')
           
           # grid - provided it is not too large
           if self.nSeqs * self.align_width > self.max_elements_in_plot:
               d = Div(text = "The alignment is too large to show with {0} elements.  A maximum of {1} elements are allowed.  The cutoff, which is in place for performance reasons, can be altered by setting max_elements_in_plot when constructing a MSAViewer object.".format(self.nSeqs * self.align_width , self.max_elements_in_plot))
               p2 = Panel(child = d, title='No interactive MSA')
           else:          
               # create graphical view of msa.  cf. https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner
     
               # create a one-base per row data frame laying out the positions of each character in the msa
               x = np.arange(0,self.align_width)
               y = self.df['y_order']  #np.arange(0,self.nSeqs)      #self.df.index.tolist()
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
               plot_height = (self.nSeqs*15)+  200
               plot_width =  1200

               # configure view length
               if self.align_width > 100:
                   viewlen =100
               else:
                   viewlen = self.align_width

               x_range = (0,viewlen)  
               y_ticks = [x+0.5 for x in self.df['y_order'].tolist()]
               y_tick_relabelling = dict(zip(y_ticks, self.df.index.tolist()))
               x_ticks = [x+0.5 for x in range(self.align_width)]
               default_x_labels = [x for x in range(self.align_width)]
               if self.positions_analysed is not None:
                   x_tick_relabelling = dict(zip(x_ticks, [str(x) for x in sorted(list(self.positions_analysed))]))
               else:
                   x_tick_relabelling = dict(zip(x_ticks, default_x_labels))

               tools = "xpan, xwheel_zoom, reset, save"
               
               # construct figure, add rectangles
               p = figure(title=None, plot_width= plot_width, plot_height=plot_height, 
                          tools=tools, min_border =0, toolbar_location='below',
                          x_range= x_range, output_backend='webgl')
               rects = Rect(x='centre_x',y='centre_y',line_width=1, 
                            fill_color= 'fill_colour', line_color= 'line_colour', width=1, height= 'height', 
                            fill_alpha=1)
               p.add_glyph(source, rects)
               p.grid.visible = False

               # configure y-axis labels
               p.yaxis.ticker = y_ticks
               p.yaxis.major_label_overrides = y_tick_relabelling
               p.yaxis.major_label_text_font = 'Courier'

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

               select = figure(title="Drag the middle and edges of the selection box to see multisequence alignment",
                    plot_height=self.nSeqs+100, plot_width=1200, y_range=p.y_range,
                    y_axis_type=None,x_axis_type=None,
                    tools="", toolbar_location=None, output_backend='webgl')

               range_tool = RangeTool(x_range=p.x_range)
               range_tool.overlay.fill_color = "navy"
               range_tool.overlay.fill_alpha = 0.2

               select.circle('centre_x', 'centre_y', size=1, fill_color='fill_colour',line_color=None, source=source)
               select.ygrid.grid_line_color = None
               select.add_tools(range_tool)
               select.toolbar.active_multi = range_tool

               if self.qq is not None:
                   if len(self.qq.index)>0:
                       qq_source = ColumnDataSource(self.qq)
                       qq_plot = figure(width=1200, height=200,title="QQ plot: non linearity is evidence for non-random distribution of variants. Yellow=expected, Blue= observed.")
                       qq_plot.line(x=[0,1],y=[0,1], color='yellow')
                       qq_plot.line(x='order_prop',y='pos_prop', source=qq_source)
                       qq_plot.xaxis.axis_label = 'Prop. variants'
                       qq_plot.yaxis.axis_label = 'Prop. genome'
                   else:
                       qq_plot = Div(text = "The alignment is of zero length.")
                       
               else:
                   qq_plot = Div(text = "Positions of variant bases were not supplied.")

               p2 = Panel(child = grid([qq_plot, select, [button1,button2],p]), title='Interactive MSA')
               # add text stating name of cluster and alignment width
           g = grid([info,Tabs(tabs=[p2, p0, pC, p1]) ])
           html = file_html(g, CDN, "Multisequence alignment")
           # show(g)
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

        x_labels = ["seq{0}".format(x) for x in range(alignLen)]
        x_labels_int = [x for x in range(alignLen)]

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

        dep_msa = DepictMSA(msa, positions_analysed  = x_labels_int)
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
        
        dep_msa = DepictMSA(msa, identify_sequence_by = ['Surname','Forename'], max_elements_in_plot=20)
        retVal = dep_msa.render_msa()
        self.assertIsInstance(retVal, str)
        self.assertTrue('</html>' in retVal)

