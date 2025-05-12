import warnings
warnings.simplefilter("ignore")

import os
os.environ["KIVY_NO_CONSOLELOG"] = "1"

from numpy.core.numeric import ones
from pandas._config.config import set_option
from kivy.app import App
from kivy.uix.floatlayout import FloatLayout
from kivy.factory import Factory
from kivy.properties import ObjectProperty, StringProperty, AliasProperty
from kivy.uix.popup import Popup
from kivy.uix.checkbox import CheckBox 
from kivy.uix.label import Label 
from kivy.uix.widget import Widget
from kivy.uix.scrollview import ScrollView
from kivy.uix.textinput import TextInput
from kivy.uix.colorpicker import ColorPicker
from kivy.clock import Clock
from kivy.event import EventDispatcher

# Loading this in is slow here, so might be nice to put elsewhere...
# Can I bring up over python programs to laod this in?
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as pl
from matplotlib import rcParams
from matplotlib import rc
import pandas
import scipy
from time import time
from matplotlib.colors import LinearSegmentedColormap
from sklearn.utils import resample
from matplotlib.lines import Line2D

from aims_immune import aims_loader as aimsLoad
from aims_immune import aims_analysis as aims
from aims_immune import aims_classification as classy

from kivy.uix.button import Button
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput

from kivy.uix.screenmanager import Screen
from os.path import dirname, join
from kivy.lang import Builder
from kivy.properties import NumericProperty, StringProperty, BooleanProperty,\
    ListProperty

import re

# DEFINE A CUSTOM COLORMAP FOR OUR NICE encode_mats:
import matplotlib as mpl
upper = mpl.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
encode_cmap = np.vstack(( lower, upper ))
encode_cmap = mpl.colors.ListedColormap(encode_cmap, name='myColorMap', N=encode_cmap.shape[0])

# Before we do some special stuff, save the directory the user is calling the
# command from... need them to be able to easily find their files if they want
startDir = os.getcwd()

# More special stuff for the pip version of the script....
import aims_immune
# The -11 is because the filepath includes '__init__.py'
# So we need to remove that to get our data path.
datPath = aims_immune.__file__[:-11]
# Need to change the directory to here to see all of our nice 
# screens that are preloaded in...
os.chdir(datPath)

# Evidently these lines will remove a weird multi-touch emulator from kivy:
from kivy.config import Config
Config.set('input', 'mouse', 'mouse,multitouch_on_demand')
# Jesus, this thing really was built for phone apps...

# Hopefully defining this outside everything lets it be called whenever I want:
# Lets try to create our own color gradient:
# code from: https://medium.com/@BrendanArtley/matplotlib-color-gradients-21374910584b
def get_color_gradient(c1, c2,c3, n):
    assert n > 1
    
    c1_rgb = np.array(c1)
    c2_rgb = np.array(c2)
    c3_rgb = np.array(c3)
    mix_pcts = [x/(n-1) for x in range(n)]
    rgb_colors1 = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
    rgb_colors2 = [((1-mix)*c2_rgb + (mix*c3_rgb)) for mix in mix_pcts]
    rgb_colors=rgb_colors1+rgb_colors2
    return (rgb_colors)

class Root(Screen):
    # REALLY Sucks to use global variables, but I'm unsure how else to
    # pass all of these across functions
    global N; N = 4
    global LFile; LFile = ['']
    global LFile_cut; LFile_cut = ['']
    fullscreen = BooleanProperty(False)

    def on_pre_enter(self):
        global loadLabel
        global loadButton
        if 'loadLabel' not in globals():
            global LFile; LFile = ['']
            global LFile_cut; LFile_cut = ['']
            N = 4
            a=0
            # Need to re-read through this shit to see what's going on...
            for j in np.arange(int(N)):
                
                # No idea what this lambda function does, but it lets us bind show load to the button
                if molecule == 'ig' or molecule == 'pep':
                    xxname = 'File '
                    xname = ''
                else:
                    xxname = 'FASTA '
                    xname = 'FASTA '
                
                if LFile == ['']:
                    if j > 0:
                        button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                        on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                        background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0),
                        disabled = True)
                    else:
                        button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                        on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                        background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
                else:
                    if j > len(LFile):
                        button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                        on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                        background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0),
                        disabled = True)
                    else:
                        button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                        on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                        background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
                
                # What an absolute fucking nightmare solution the line above... works though
                # Need to make sure we don't overwrite the labels every time we load shit in
                if j >= len(LFile):
                    label = Label(text=xname +'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')
                else:
                    if LFile[int(j)] != '':
                        # Shouldnt worry about entering IF on LFile then showing LFile_cut,
                        # because they are defined together
                        label = Label(text=LFile_cut[int(j)], size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')
                    else:
                        label = Label(text=xname+'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')

                if a == 0:
                    loadButton = [button]
                    loadLabel = [label]
                else:
                    loadButton = loadButton + [button]
                    loadLabel = loadLabel + [label]
                a = a + 1

            for i in loadButton:
                FloatLayout.add_widget(self, i)
            for k in loadLabel:
                FloatLayout.add_widget(self, k)

    def get_path(self):
        return(startDir)
        #return(os.getcwd())

    def dismiss_popup(self):
        self._popup.dismiss()

    # This is getting more and more confusing. Apologies for whoever has to go through this.
    # this is basically just a dummy function to keep track of which FASTA file is where...
    def do_thing(self,win2 = 1):
        global FASTA_L
        # What a wild way for all of this to work... BUT IT DOES
        if str(type(win2)) == "<class 'kivy.uix.button.Button'>":
            FASTA_L = int(win2.text[-1])-1
        else:
            FASTA_L = win2

    # Basically just recreate the entire screen every time we want to 
    # add more fastas. Kind of a pain, but it works.
    def check_one(self):
        global onlyONE
        global LFile
        if self.tbutt1.state == 'down':
            for i in loadButton[1:]:
                i.disabled = True
            onlyONE = True
            if LFile[0] == '':
                return
            else:
                if len(LFile) > 1:
                    LFile = [LFile[0]]
                    LFile_cut = [LFile_cut[0]]
                self.next1_1.disabled = False
        elif self.tbutt2.state == 'down':
            onlyONE = False
            if LFile != ['']:
                loadButton[1].disabled = False
            if len(LFile) >= 2:
                return
            else:
                self.next1_1.disabled = True
            return

    # win is how we will try to keep track of using the right button...
    def show_load(self, win = 2):
        
        content = LoadDialog(load=self.load, cancel=self.dismiss_popup, fas1 = self.do_thing(win2 = win))
        self._popup = Popup(title="Load file", content=content,
                            size_hint=(0.75, 0.9))
        self._popup.open()

    def load(self, path, filename):
        global LFile
        # LFile_cut will be a more visually appealing path
        # Users will hopefully be able to see the full LFile by clicking a button
        global LFile_cut
        path1 = os.path.join(path, filename[0])
        # So FASTA_L should tell you WHERE
        # The loadfile is coming from
        while FASTA_L+1 > len(LFile):
            LFile = LFile + ['']
            LFile_cut = LFile_cut + ['']
        LFile[FASTA_L] = path1
        # Break up path1 a little bit
        # rfind is "find the last one"
        path1_cut = path1[path1.rfind("/")+1:]
        LFile_cut[FASTA_L] = path1_cut
        # Need to have two separate options because we move from Kivy defined buttons to
        # python defined buttons. I handle those slightly differently.
        loadLabel[FASTA_L].text = path1_cut
        if len(loadButton) > FASTA_L+1: 
            loadButton[FASTA_L+1].disabled = False
        self.check_one()

        if len(LFile) >= 2:
            self.next1_1.disabled = False 

        self.dismiss_popup()

    def make_path(self):
        global dir_name
        self.text1 = self.v_input1.text
        dir_name = self.text1

    def reset_loadScrn(self):
        # need to do something to let people know they need to re-enter everything...
        global N
        global loadButton
        global loadLabel
        global LFile
        self.lessF.disabled = True
        self.next1_1.disabled = True

        for i in loadButton:
            FloatLayout.remove_widget(self,i)
        for j in loadLabel:
            FloatLayout.remove_widget(self,j)

        # Remove any extra entries, just in case...
        while len(LFile) > N:
            LFile = LFile[:-1]
            LFile_cut = LFile_cut[:-1]


    def more_fastas(self):
        global N
        global loadButton
        global loadLabel
        # Alright so this works...
        self.lessF.disabled = False
        N = N + 1
        a = 0
        for i in loadButton:
            FloatLayout.remove_widget(self,i)
        for j in loadLabel:
            FloatLayout.remove_widget(self,j)

        for j in np.arange(int(N)):

            # No idea what this lambda function does, but it lets us bind show load to the button
            if molecule == 'ig' or molecule == 'pep':
                xxname = 'File '
                xname = ''
            else:
                xxname = 'FASTA '
                xname = 'FASTA '

            if LFile == ['']:
                if j > 0:
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0),
                    disabled = True)
                else:
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
            else:
                if j > len(LFile):
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0),
                    disabled = True)
                else:
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
            # What an absolute fucking nightmare solution the line above... works though
            # Need to make sure we don't overwrite the labels every time we load shit in
            if j >= len(LFile):
                label = Label(text=xname +'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')
            else:
                if LFile[int(j)] != '':
                    label = Label(text=LFile_cut[int(j)], size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')
                else:
                    label = Label(text=xname+'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')

            if a == 0:
                loadButton = [button]
                loadLabel = [label]
            else:
                loadButton = loadButton + [button]
                loadLabel = loadLabel + [label]

            a = a + 1
        for i in loadButton:
            FloatLayout.add_widget(self, i)
        for k in loadLabel:
            FloatLayout.add_widget(self, k)

    def less_fastas(self):
        global N
        global loadButton
        global loadLabel
        # Alright so this works...
        N = N - 1
        a = 0
        # WHAT A MESS THIS IS, BUT IT WORKS!

        if N == 4:
            self.lessF.disabled = True

        for i in loadButton:
            FloatLayout.remove_widget(self,i)
        for j in loadLabel:
            FloatLayout.remove_widget(self,j)

        for j in np.arange(int(N)):

            # No idea what this lambda function does, but it lets us bind show load to the button
            if molecule == 'ig' or molecule == 'pep':
                xxname = 'File '
                xname = ''
            else:
                xxname = 'FASTA '
                xname = 'FASTA '

            if LFile == ['']:
                if j > 0:
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0),
                    disabled = True)
                else:
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
            else:
                if j > len(LFile):
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0),
                    disabled = True)
                else:
                    button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.25, 'center_y':.75-int(a)*0.6/N},
                    on_release=lambda x = int(j):self.show_load(win = x),background_down='app_data/butt_down.png', 
                    background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
            
            # What an absolute fucking nightmare solution the line above... works though
            # Need to make sure we don't overwrite the labels every time we load shit in
            if j >= len(LFile):
                label = Label(text=xname +'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')
            else:
                if LFile[int(j)] != '':
                    label = Label(text=LFile_cut[int(j)], size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')
                else:
                    label = Label(text=xname+'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.62, 'center_y':.75-int(a)*0.6/N},font_name='app_data/Poppins-Light.ttf')
            if a == 0:
                loadButton = [button]
                loadLabel = [label]
            else:
                loadButton = loadButton + [button]
                loadLabel = loadLabel + [label]

            a = a + 1
        for i in loadButton:
            FloatLayout.add_widget(self, i)
        for k in loadLabel:
            FloatLayout.add_widget(self, k)

    def add_widget(self, *args):
        if 'content' in self.ids:
            return self.ids.content.add_widget(*args)
        return super(Root, self).add_widget(*args)

# Alright for our alignment entries, let's try breaking all of 
# this code up into another class...
class aligner(Screen):
    # Need to redefine a different type of loading from the wild one above.
    def dismiss_popup(self):
        self._popup.dismiss()

    def show_load(self, *args):
        content = LoadDialog(load=self.load, cancel=self.dismiss_popup)
        self._popup = Popup(title="Load file", content=content,size_hint=(0.9, 0.9))
        self._popup.open()

    def drop_buts(self):
        x = self.align_button1; FloatLayout.remove_widget(self, x)
        x = self.align_button2; FloatLayout.remove_widget(self, x)
        global align_label
        global alignEnt
        global thing_loaded
        align_label = Label(text = '')
        alignEnt = ''
        thing_loaded = False
    # This is our special load function just for this page, so 
    # do everything inside this one function
    #... probably if I did this more often, I wouldn't need to assign so many
    # global variables... whoops.
    def load(self, path, filename):
        global align_label
        global alignEnt
        global thing_loaded
        pos_path = os.path.join(path, filename[0])
        matPos = pandas.read_csv(pos_path,header=0).values
        # Check that they have at least two entries
        if len(matPos) < 2:
            if align_label.text == 'ERROR: Matrix only has 1 row. Must have at least 2':
                FloatLayout.remove_widget(self,align_label) 
            align_label = Label(text='ERROR: Matrix only has 1 row. Must have at least 2', size_hint=(0.4, 0.2),
                    pos_hint={'center_x':.5, 'center_y':.5}, color = (1, 0, 0, 1))
            FloatLayout.add_widget(self, align_label)
            button_moved = Button(text='Load Matrix', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.4, 'center_y':0.1}, on_release=self.show_load,
            background_down='app_data/butt_down.png',background_normal='app_data/butt_up.png',
            color=(0, 0.033, 0.329, 1), border=(0, 0, 0, 0))
            button2_moved = Button(text='Resize Entries', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.6, 'center_y':0.1}, on_release=self.load_aligners,
            background_down='app_data/butt_down.png',background_normal='app_data/butt_up.png',
            color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
            FloatLayout.add_widget(self, button2_moved)
            FloatLayout.add_widget(self, button_moved)
        # Check that they have the right number of columns, or else the program will crash
        # LEAVE THIS BE FOR NOW, BUT ACTUALLY IMPLEMENT ASAP
        #elif len(matPos[0]) != 6:
        #    if align_label.text == 'ERROR: Matrix only has 1 row. Must have at least 2':
        #        FloatLayout.remove_widget(self,align_label) 
        #    align_label = Label(text='ERROR: Matrix only has 1 row. Must have at least 2', size_hint=(0.4, 0.2),
        #            pos_hint={'center_x':.5, 'center_y':.5}, color = (1, 0, 0, 1))
        else:
            if align_label.text == 'ERROR: Matrix only has 1 row. Must have at least 2':
                FloatLayout.remove_widget(self,align_label)
                align_label.text = ''
            # REALIZING NOW I DONT HAVE ANY FUNCTIONS TO REMOVE THESE BUTTONS,
            # PROBABLY DOESN'T MATTER, BUT THEORETICALLY COULD CAUSE A CRASH?
            button_moved = Button(text='Load Matrix', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.4, 'center_y':0.1}, on_release=self.show_load,
            background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
            button2_moved = Button(text='Resize Entries', size_hint=(None, None),background_down='app_data/butt_down.png',
            background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0),
            size =('150dp','48dp'),pos_hint={'center_x':.6, 'center_y':0.1}, on_release=self.load_aligners)
            FloatLayout.add_widget(self, button2_moved)
            FloatLayout.add_widget(self, button_moved)

            # If we're running this script again, we've almost certainly gotta remove it...
            if alignEnt != '':
                for entry in alignEnt:
                    for j in entry:
                        FloatLayout.remove_widget(self, j)

            for j in np.arange(len(matPos)):
                if matPos[j,0] != '':
                    name = matPos[j,0]
                else:
                    name = 'FASTA ' + str(j+1)
                textinput0 = TextInput(text=name, multiline=False, size_hint = (None, None),write_tab =False,
                pos_hint = {'center_x': 0.1, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='125dp')
                textinput1 = TextInput(text=str(matPos[j,1]), multiline=False, size_hint = (None, None),write_tab = False,
                pos_hint = {'center_x': 0.25, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
                textinput2 = TextInput(text=str(matPos[j,2]), multiline=False, size_hint = (None, None),write_tab = False,
                pos_hint = {'center_x': 0.4, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
                textinput3 = TextInput(text=str(matPos[j,3]), multiline=False, size_hint = (None, None),write_tab = False,
                pos_hint = {'center_x': 0.55, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
                textinput4 = TextInput(text=str(matPos[j,4]), multiline=False, size_hint = (None, None),write_tab = False,
                pos_hint = {'center_x': 0.7, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
                textinput5 = TextInput(text=str(matPos[j,5]), multiline=False, size_hint = (None, None),write_tab = False,
                pos_hint = {'center_x': 0.85, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
                
                if int(j) == 0:
                    alignEnt = [[textinput0,textinput1,textinput2,textinput3,textinput4,textinput5]]
                else:
                    alignEnt = alignEnt + [[textinput0,textinput1,textinput2,textinput3,textinput4,textinput5]]
            # Before adding in all of the new ones
            for entry in alignEnt:
                for j in entry:
                    FloatLayout.add_widget(self, j)
            thing_loaded = True

        self.dismiss_popup()

    def load_aligners(self, *args):
        global align_label
        global alignEnt
        global thing_loaded
        self.next1_2.disabled = False
        # Literally just here so we can have these messages or not
        if align_label.text == 'ERROR: MUST LOAD AT LEAST TWO FILES':
            FloatLayout.remove_widget(self,align_label)
            align_label.text = ''
        button_moved = Button(text='Load Matrix', size_hint=(None, None),
        size =('150dp','48dp'),pos_hint={'center_x':.4, 'center_y':0.1}, on_release=self.show_load,
        background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
        button2_moved = Button(text='Resize Entries', size_hint=(None, None),
        size =('150dp','48dp'),pos_hint={'center_x':.6, 'center_y':0.1}, on_release=self.load_aligners,
        background_normal='app_data/butt_up.png',color=(0, 0.033, 0.329, 1),border=(0, 0, 0, 0))
        FloatLayout.add_widget(self, button2_moved)
        FloatLayout.add_widget(self, button_moved)

        # If we're running this script again, we've almost certainly gotta remove it...
        if alignEnt != '':
            for entry in alignEnt:
                for j in entry:
                    FloatLayout.remove_widget(self, j)

        for j in np.arange(len(LFile)):
            if LFile[j] != '':
                startname = [y.end() for y in re.finditer('/',LFile[j])]
                endname = [y.start() for y in re.finditer('.fasta',LFile[j])]
                if len(endname) == 0:
                    endname = [int(len(LFile[j])-1)]
                name = LFile[j][int(startname[len(startname)-1]):int(endname[0])]
            else:
                name = 'File ' + str(j+1)
            textinput0 = TextInput(text=name, multiline=False, size_hint = (None, None),write_tab =False,
            pos_hint = {'center_x': 0.1, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='125dp')
            textinput1 = TextInput(text='', multiline=False, size_hint = (None, None),write_tab = False,
            pos_hint = {'center_x': 0.25, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
            textinput2 = TextInput(text='', multiline=False, size_hint = (None, None),write_tab = False,
            pos_hint = {'center_x': 0.4, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
            textinput3 = TextInput(text='', multiline=False, size_hint = (None, None),write_tab = False,
            pos_hint = {'center_x': 0.55, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
            textinput4 = TextInput(text='', multiline=False, size_hint = (None, None),write_tab = False,
            pos_hint = {'center_x': 0.7, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
            textinput5 = TextInput(text='', multiline=False, size_hint = (None, None),write_tab = False,
            pos_hint = {'center_x': 0.85, 'center_y': 0.60-int(j)*0.5/N}, height = '32dp',width='42dp')
            
            if int(j) == 0:
                alignEnt = [[textinput0,textinput1,textinput2,textinput3,textinput4,textinput5]]
            else:
                alignEnt = alignEnt + [[textinput0,textinput1,textinput2,textinput3,textinput4,textinput5]]
        # Before adding in all of the new ones
        for entry in alignEnt:
            for j in entry:
                FloatLayout.add_widget(self, j)
        thing_loaded = True

    # So now we're done with our loading stuff, let's wrap things up and get ready to move on...
    def make_path(self):
        global dir_name
        self.text1 = self.input1.text
        dir_name = self.text1
    
    def check_aligns(self):
        global mat_coords
        global labels
        global check_run
        # Real quick check if users want to drop degenerate sequences or not.
        global exp_drop
        exp_drop = self.degen_drop.active
        # Really idiot-proof this thing, in case someone is just buzzing through the screens
        if 'thing_loaded' in globals():
            if thing_loaded:
                # Convert our weird little kivy objects into a numpy array
                x,y = np.shape(alignEnt)
                for row in np.arange(x):
                    for column in np.arange(y):
                        if column == 0:
                            mat_pre = alignEnt[row][column].text
                        else:
                            mat_pre = np.hstack((mat_pre,alignEnt[row][column].text))
                    if row == 0:
                        mat_coords = mat_pre
                    else:
                        mat_coords = np.vstack((mat_coords,mat_pre))
        if molecule == 'mhc' and onlyONE:
            labels = mat_coords[0]
        else:
            labels = mat_coords[:,0]
            # All of this here is to fix issues in naming
            # this fix will also break if users for some reason end filenames
            # with "___". Better hope they don't do that
            dum_len = len(labels)
            for i in np.arange(dum_len):
                for j in np.arange(dum_len):
                    if i == j:
                        continue
                    if labels[i].find(labels[j]) == 0:
                        labels[j] = labels[j]+'___'
        # When you go back this far, make sure we recreate the checkbox screen
        if 'check_run' in globals():
            # So intelligently, if you define something as a global variable
            # but don't "fill" it, it won't register as being in globals.
            del check_run
    
class LoadDialog(FloatLayout):
    load = ObjectProperty(size_hint=(0.2,0.2))
    cancel = ObjectProperty(size_hint=(0.2,0.2))
    fas1 = ObjectProperty(None)

    def get_path(self):
        return(startDir)
        #return(os.getcwd())

# let's try to have all of my nice python functions in 
# the "analysis" screen instance:
class Analysis(Screen):
    def get_matrix(self):
        # Ok so here is where I will actually call my python code
        # For now, print all the shit
        #########################################
        # DEFINE PLOT PARAMETERS:
        font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 16}
        COLOR = 'black'
        rcParams['text.color'] = 'black'
        rcParams['axes.labelcolor'] = COLOR
        rcParams['xtick.color'] = COLOR
        rcParams['ytick.color'] = COLOR

        rc('font', **font)

        # Custom colormap code from: https://stackoverflow.com/questions/49367144/modify-matplotlib-colormap
        import matplotlib as mpl
        upper = mpl.cm.jet(np.arange(256))
        lower = np.ones((int(256/4),4))
        for i in range(3):
            lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])

        global cmap
        cmap = np.vstack(( lower, upper ))
        cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])
        ##########################################
        global seq_MIf
        global seqNameF
        global seq_final
        global labels
        #this_dir = os.getcwd()
        this_dir = startDir
        paths = LFile
        # mat_coords is predefined from above function

        ######### PRE LOAD DATA SO TESTING CAN GO FASTER!!!! #############
        # x1 = np.array(['cd1','HLA-A','UFA','UAA','UDA']); x2 = np.array(['124','170','22','2','2'])
        # x3 = np.array(['167','210','66','49','49']); x4 = np.array(['209','260','105','93','93'])
        # x5 = np.array(['262','306','158','152','152']);  x6 = np.array(['303','348','199','193','193'])
        # mat_coords = np.transpose(np.vstack((x1,x2,x3,x4,x5,x6)))
        # paths = ['/Users/boughter/Desktop/AIMS/gui/AIMS/mhc_testData/cd1_seqs.fasta',
        # '/Users/boughter/Desktop/AIMS/gui/AIMS/mhc_testData/hlaA_seqs.fasta',
        # '/Users/boughter/Desktop/AIMS/gui/AIMS/mhc_testData/cd1_ufa_genes.fasta',
        # '/Users/boughter/Desktop/AIMS/gui/AIMS/mhc_testData/UAA_seqs.fasta',
        # '/Users/boughter/Desktop/AIMS/gui/AIMS/mhc_testData/UDA_seqs.fasta']
        #################### DEBUG VERSION ONLY #########################

        AA_num_key = aims.get_props()[1]
        AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        # Labels will always be predefined for Abs
        # Predefined, but poorly formatted... reformat here
        global align
        global labels
        global LOOPnum
        global exp_drop
        global mat_size
        global ID
        if molecule == 'mhc':
            if onlyONE:
                data = mat_coords[1:]
            else:
                data = mat_coords[:,1:]

        for i in np.arange(len(paths)):
            if molecule == 'mhc':
                #if any(data[:,i] == ['','','','','']):
                #    self.img1.source = this_dir + '/app_data/error_pos.png'
                #    return()
                # turn data into an integer.
                if onlyONE:
                    int_dat = [int(x) for x in data]
                    seq,seq_key = aimsLoad.mhc_loader(paths[i],int_dat,labels[i],drop_dups = exp_drop)
                else:
                    int_dat = [int(x) for x in data[i]]
                    seq,seq_key = aimsLoad.mhc_loader(paths[i],int_dat,labels[i],drop_dups = exp_drop)
            elif molecule == 'ig':
                # My first ever error handling! Works pretty well (for now)
                # Probably need more "exceptions" here for formatting errors
                try:
                    seq = aimsLoad.Ig_loader(paths[i],labels[i],loops=LOOPnum,drop_degens=exp_drop)
                except pandas.errors.ParserError:
                    # Alrighty we want a popup here on this screen
                    popup = Popup(title='ERROR (Click Anywhere to Dismiss)',
                    content=Label(text='Wrong # loops, go back and redefine'),
                    size_hint=(None, None), size=(600, 600))
                    popup.open()
                    return
            elif molecule == 'msa':
                # ADD IN MSA LOADING
                seq = aimsLoad.msa_loader(paths[i],label=labels[i],drop_dups = exp_drop)
            elif molecule == 'pep':
                # Add in pep loading
                # Should probably throw in an option for a length cutoff...
                # DONT have a cutoff for now...
                seq = aimsLoad.pep_loader(paths[i],label=labels[i],
                                                         drop_degens = exp_drop,len_cutoff=cutoff)
            ID_pre = np.ones(np.shape(seq)[1])
            if i == 0:
                seq_final = seq
                seq_size = np.shape(seq)[1]
                seqNameF = labels[i]
                ID = i*ID_pre
                if molecule == 'mhc':
                    seq_keyF = seq_key
                mat_size = aims.get_sequence_dimension(np.array(seq))[0]
            else:
                seq_final = pandas.concat([seq_final,seq],axis = 1)
                seqNameF = np.vstack((seqNameF,labels[i]))
                seq_size = np.vstack((seq_size,np.shape(seq)[1]))
                ID = np.hstack((ID, i*ID_pre))
                if molecule == 'mhc':
                    seq_keyF = np.hstack((seq_keyF,seq_key))
                mat_size2 = aims.get_sequence_dimension(np.array(seq))[0]
                if type(mat_size) != int:
                    max_lenp=np.zeros(len(mat_size))
                    for i in np.arange(len(mat_size)):
                        max_lenp[i]=int(max(mat_size[i],mat_size2[i]))
                else:
                    max_lenp = int(max(mat_size,mat_size2))
                mat_size = max_lenp
        # Alright if we have gotten this far, bring the "next" button back online
        self.next1_3.disabled = False
        # Obviously need to define this somewhere in the software
        if molecule == 'ig' or molecule == 'pep':
            if self.alnc.active:
                align = 'center'
            elif self.alnb.active:
                align = 'bulge'
            elif self.alnl.active:
                align = 'left'
            elif self.alnr.active:
                align = 'right'
        else:
            align = 'center'

        if molecule == 'msa':
            AA_key_dash = AA_key + ['-']
            AA_num_key_dash = np.hstack((AA_num_key,[0]))
            seq_MI = aims.gen_MSA_matrix(np.array(seq_final),AA_key_dash = AA_key_dash, key = AA_num_key_dash, giveSize = mat_size)
        if molecule == 'pep':
            seq_MI = aims.gen_tcr_matrix(np.array(seq_final),key = AA_num_key, giveSize = mat_size, alignment = align, bulge_pad = bulge_pad)
        else:
            seq_MI = aims.gen_tcr_matrix(np.array(seq_final),key = AA_num_key, giveSize = mat_size, alignment = align)
        # Convert our MI matrix to a pandas dataframe
        seq_MIf = pandas.DataFrame(np.transpose(seq_MI),columns = seq_final.columns)
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
        ax[0,0].imshow(seq_MI, interpolation='nearest', aspect='auto',cmap=encode_cmap)
        pl.close()

        # Need to try to draw lines separating the distinct groups...
        seq_locs = 0
        if onlyONE:
            pass
        else:
            for i in np.arange(len(seq_size)-1):
                seq_locs = seq_locs + seq_size[i]
                ax[0,0].plot(np.arange(len(np.transpose(seq_MI))),np.ones(len(np.transpose(seq_MI)))*seq_locs,'black',linewidth = 3)

        # Alright now we want to change xlabel to actually talk about the features...
        ax[0,0].set_ylabel('Sequence Number')
        global xtick_loc
        if type(mat_size) != int:
            for i in np.arange(len(mat_size)):
                if i == 0:
                    xtick_loc = [mat_size[i]/2]
                else:
                    pre_loc = sum(mat_size[:i])
                    xtick_loc = xtick_loc + [mat_size[i]/2 + pre_loc]
            ax[0,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
        else:
            xtick_loc = [mat_size/2]
            ax[0,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
        if molecule == 'mhc':
            ax[0,0].set_xticklabels(['Strand 1','Helix 1','Strand 2','Helix 2'])
        if molecule == 'ig':
            if LOOPnum == 6:
                ax[0,0].set_xticklabels(['CDR1L','CDR2L','CDR3L','CDR1H','CDR2H','CDR3H'])
            elif LOOPnum == 2:
                ax[0,0].set_xticklabels(['CDR3H','CDR3L'])
            elif LOOPnum == 3:
                ax[0,0].set_xticklabels(['CDR1','CDR2','CDR3'])
            elif LOOPnum == 1:
                ax[0,0].set_xticklabels(['CDR Loop'])
            
            Numclones = np.shape(seq_MIf)[1]
            if type(mat_size) != int:
                for i in np.arange(len(mat_size)-1):
                    ax[0,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(Numclones),np.arange(Numclones),'black',linewidth = 3)

        # NOTE need to make sure that I do something so that files are not overwritten...
        if os.path.exists(this_dir +'/' + dir_name):
            None
        else:
            os.mkdir(this_dir +'/' + dir_name)
        fig.savefig(this_dir + '/' + dir_name + '/matrix.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/matrix.png',format='png',dpi=500)
        # save raw data at all steps now...
        np.savetxt(this_dir + '/' + dir_name + '/raw_matrix.dat',seq_MIf,fmt='%.3f')
        # Received an explicit request to copy the encoding BACK into sequences...
        # I think I have a script for this, just need to save it into a file
        seqT = np.transpose(seq_MIf)
        fin_convert = []
        for i in np.arange(len(seqT)):
            tt = aims.decode_mat(seqT.iloc[i],num_key_AA=AA_num_key,key_AA=AA_key)
            fin_convert = fin_convert + [*tt]
        fff = np.array(fin_convert).reshape(len(seqT),len(tt))
        pandas.DataFrame(fff).to_csv(this_dir + '/' + dir_name + '/AIMS_encodedSeqs.csv')
        if molecule == 'mhc':
            np.savetxt(this_dir + '/' + dir_name + '/sequence_key.txt',seq_keyF,fmt='%s')

        self.img1.source = this_dir + '/' + dir_name + '/matrix.png'

    def reduce_dim(self):
        self.clust_butt.disabled = False
        global chosen_dset
        full_big.index = seq_MIf.columns

        if self.dat1.active:
            dChoice = 'net'
        elif self.dat2.active:
            dChoice = 'parsed'

        if dChoice == 'full':
            chosen_dset = full_big
        elif dChoice == 'parsed':
            chosen_dset = parsed_mat
        elif dChoice == 'net':
            first = True
            for i in seqNameF:
                pre_label = i[0]
                index = [column for column in seq_MIf.columns if pre_label in column]
                pre_array = np.array(full_big.loc[index])
                dat_size = len(index)
                seq_bigF = pre_array.reshape(dat_size,61,len(seq_MIf))
                if first:
                    big_reshape = seq_bigF
                    first = False
                else:
                    big_reshape = np.vstack((big_reshape,seq_bigF))
            chosen_dset = pandas.DataFrame(np.average(big_reshape,axis=2))
            chosen_dset.columns = prop_names

        chosen_dset.index = seq_MIf.columns
        final_chosen = np.transpose(chosen_dset)
        #### For datasets with > 2 entries need to allow those to pick and choose:
        global plot_label
        if onlyONE:
            plot_label = seqNameF
            ID_new = ID
            pass
        elif len(labels) > 2:
            aaa = 0; bbb = 0
            for i in range(len(checked)):
                # Don't love this solution but it works... so...
                pre_label = seqNameF[i][0]
                index = [column for column in seq_MIf.columns if pre_label in column]
                # Alright so our new matrix, checked is "i" rows by three
                # 0 is the "group 1", 1 is the "group 2"
                if checked[i,0]:
                    if aaa == 0:
                        for_clust = final_chosen[index]
                        ID_new = aaa * np.ones(np.shape(for_clust)[1])
                        plot_label = pre_label
                        aaa = 1
                    else:
                        for_clust = pandas.concat([for_clust,final_chosen[index]],axis=1)
                        ID_new = np.hstack((ID_new, aaa * np.ones(np.shape(final_chosen[index])[1])))
                        plot_label = np.hstack((plot_label,pre_label))
                        aaa = aaa + 1
                elif checked[i,1]:
                    if bbb == 0:
                        no_clust = final_chosen[index]
                        bbb = 1
                    else:
                        no_clust = pandas.concat([no_clust,final_chosen[index]],axis=1)
                        bbb = bbb + 1
            # Then just redefine "for_clust" as "chosen_dset" for consistency with labels <=2
            chosen_dset = np.transpose(for_clust)
        else:
            plot_label = np.hstack((seqNameF[0][0],seqNameF[1][0]))
            ID_new = ID

        global dim
        if self.pca2.active:
            reduce = 'pca'; dim = '2d'
        elif self.pca3.active:
            reduce = 'pca'; dim = '3d'
        elif self.umap2.active:
            reduce = 'umap'; dim = '2d'
        elif self.umap3.active:
            reduce = 'umap'; dim = '3d'
        else:
            return
        if reduce == 'pca':
            from sklearn.decomposition import PCA
            pca = PCA(n_components=3, svd_solver='full')
            final=pca.fit_transform(chosen_dset)
            transform = pandas.DataFrame(np.transpose(final),columns = chosen_dset.index)
            # Alright if using PCA, output the explained variance and top10 components
            explain_var = pca.explained_variance_ratio_
            comp1_df = np.transpose(pandas.DataFrame([np.abs(pca.components_[0]),chosen_dset.columns]))
            comp2_df = np.transpose(pandas.DataFrame([np.abs(pca.components_[1]),chosen_dset.columns]))
            comp3_df = np.transpose(pandas.DataFrame([np.abs(pca.components_[2]),chosen_dset.columns]))
            comp1_top10 = comp1_df.sort_values(0,ascending=False).values[0:10]
            comp2_top10 = comp2_df.sort_values(0,ascending=False).values[0:10]
            comp3_top10 = comp3_df.sort_values(0,ascending=False).values[0:10]
            # Save all the good stuff
            #this_dir = os.getcwd()
            this_dir = startDir
            np.savetxt(this_dir + '/' + dir_name + '/pca_explainVar.dat',explain_var,fmt='%.3f')
            np.savetxt(this_dir + '/' + dir_name + '/pca_comp1_top10.dat',comp1_top10,fmt='%s')
            np.savetxt(this_dir + '/' + dir_name + '/pca_comp2_top10.dat',comp2_top10,fmt='%s')
            np.savetxt(this_dir + '/' + dir_name + '/pca_comp3_top10.dat',comp3_top10,fmt='%s')
        elif reduce == 'umap':
            import umap
            reducer = umap.UMAP(n_components=3, n_neighbors = 25, n_jobs=1, random_state = 47)
            final = reducer.fit_transform(chosen_dset)
            transform = pandas.DataFrame(np.transpose(final),columns = chosen_dset.index)
            #this_dir = os.getcwd()
            this_dir = startDir
        elif reduce == 'tsne':
            from sklearn.manifold import TSNE
            tsne = TSNE(n_components = 3, random_state = 47, n_jobs = 1)
            final=tsne.fit_transform(chosen_dset)
            transform = pandas.DataFrame(np.transpose(final),columns = chosen_dset.index)
            #this_dir = os.getcwd()
            this_dir = startDir
        
        global clust_input
        clust_input = np.array(np.transpose(transform))

        # A bit of an odd solution to how we add in a legend to the PCA/UMAP
        from matplotlib.lines import Line2D
        cmap2 = pl.get_cmap('rainbow')
        global cmap_discrete
        cmap_discrete = cmap2(np.linspace(0, 1, len(seqNameF)))
        for i in np.arange(len(seqNameF)):
            if i == 0:
                custom_lines = [Line2D([0], [0], color='w', marker='o',markerfacecolor=cmap_discrete[i], markersize= 10)]
            else:
                custom_lines = custom_lines + [Line2D([0], [0], color='w', marker='o',markerfacecolor = cmap_discrete[i], markersize= 10)]

        if dim == '2d':
            fig = pl.figure(figsize = (14, 10))
            pl.scatter(clust_input[:,0],clust_input[:,1],c = ID_new, cmap = 'rainbow')
            pl.legend(custom_lines, plot_label)
            pl.xlabel('AX1'); pl.ylabel('AX2')
        elif dim == '3d':
            from mpl_toolkits import mplot3d
            fig = pl.figure(figsize = (14, 10))
            ax3d = fig.add_subplot(111, projection='3d')
            ax3d.scatter(clust_input[:,0],clust_input[:,1],clust_input[:,2],c = ID_new, cmap ='rainbow')
            pl.legend(custom_lines, plot_label)
            ax3d.set_xlabel('AX1',labelpad=20)
            ax3d.set_ylabel('AX2',labelpad=20)
            ax3d.set_zlabel('AX3',labelpad=10)
        fig.savefig(this_dir + '/' + dir_name + '/'+reduce+'_'+dim+'_fig.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/'+reduce+'_'+dim+'_fig.png',format='png',dpi=500)
        np.savetxt(this_dir + '/' + dir_name + '/'+reduce+'_'+dim+'.dat',clust_input,fmt='%.3f')
        #NOTE DELETE THIS BELOW LINE ASAP
        self.img6.source = this_dir + '/' + dir_name + '/'+reduce+'_'+dim+'_fig.png'
        pl.close()

    def cluster_seqs(self):
        global cluster_dset
        import sklearn.cluster as cluster
        global clust
        global NClust
        self.next1_5.disabled = False
        if self.kmean.active:
            clust = 'kmean'
            NClust = int(self.kclust.text)
            clusts = cluster.KMeans(n_clusters=NClust).fit_predict(clust_input)
        elif self.optics.active:
            clust = 'optics'
            OPsamp = int(self.opticnum.text)
            clusts = cluster.OPTICS(min_samples=OPsamp).fit_predict(clust_input)
        elif self.dbscan.active:
            clust = 'dbscan'
            DBrad = float(self.dbrad.text)
            clusts = cluster.DBSCAN(eps=DBrad).fit_predict(clust_input)
        else:
            return

        cluster_dset = pandas.DataFrame(clusts,columns=['cluster'])
        # Want the min cluster to be white for non kmean cluster, as these are unclustered
        if dim == '2d':
            fig = pl.figure(figsize = (14, 10))
            if clust == 'kmean':
                pl.scatter(clust_input[:,0],clust_input[:,1],c = clusts, cmap = 'rainbow')
            else:
                pl.scatter(clust_input[:,0],clust_input[:,1],c = clusts, cmap = cmap)
            pl.xlabel('AX1'); pl.ylabel('AX2')
        elif dim == '3d':
            from mpl_toolkits import mplot3d
            fig = pl.figure(figsize = (14, 10))
            ax3d = fig.add_subplot(111, projection='3d')
            if clust == 'kmean':
                ax3d.scatter(clust_input[:,0],clust_input[:,1],clust_input[:,2],c = clusts, cmap='rainbow')
            else:
                ax3d.scatter(clust_input[:,0],clust_input[:,1],clust_input[:,2],c = clusts, cmap=cmap)
            ax3d.set_xlabel('AX1',labelpad=20)
            ax3d.set_ylabel('AX2',labelpad=20)
            ax3d.set_zlabel('AX3',labelpad=10)

        #this_dir = os.getcwd()
        this_dir = startDir
        fig.savefig(this_dir + '/' + dir_name + '/'+clust+'_'+dim+'_fig.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/'+clust+'_'+dim+'_fig.png',format='png',dpi=500)
        np.savetxt(this_dir + '/' + dir_name + '/'+clust+'_'+dim+'.dat',clusts,fmt='%.3f')
        pl.close()
        self.img6.source = this_dir + '/' + dir_name + '/'+clust+'_'+dim+'_fig.png'

    def analyze_clusts(self):
        self.next1_6.disabled = False
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,10))
        global fin_clustL
        fin_clustL = len(cluster_dset['cluster'].drop_duplicates())
        cluster_bins = np.zeros((fin_clustL,len(seqNameF)))
        kmean = True
        for i in np.sort(cluster_dset['cluster'].drop_duplicates()):
            pre_clust = seq_MIf[chosen_dset.index[cluster_dset[cluster_dset['cluster'] == i].index]]
            clustID = np.transpose(pandas.DataFrame(i*np.ones(np.shape(pre_clust)[1])))
            clustID.columns = pre_clust.columns
            pre_clustF = pandas.concat([pre_clust,clustID],axis=0)
            for j in np.arange(len(seqNameF)):
                index = [column for column in clustID.columns if seqNameF[j][0] in column]
                cluster_bins[i,j] = np.shape(clustID[index])[1]
            # Do not plot the unclustered data
            if i == -1:
                kmean = False
                continue
            if i == 0:
                clustered = pre_clustF
            else:
                clustered = pandas.concat([clustered, pre_clustF],axis = 1)
            ax[0,0].plot(np.arange(len(seq_MIf)),np.ones(len(seq_MIf))*(np.shape(clustered)[1]),'black',linewidth = 3)

        ax[0,0].imshow(np.transpose(np.array(clustered))[:,:-1], interpolation='nearest', aspect='auto',cmap=encode_cmap)
        ax[0,0].set_ylabel('Sequence Number')
        #this_dir = os.getcwd()
        this_dir = startDir
        fig.savefig(this_dir + '/' + dir_name + '/'+clust+'_mat_fig.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/'+clust+'_mat_fig.png',format='png',dpi=500)
        self.img7.source = this_dir + '/' + dir_name + '/'+clust+'_mat_fig.png'
        pl.close()
        
        # Skip all of this stuff if only looking at ONE file... can't get cluster membership then
        if onlyONE:
            return
        else:
            fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,10))
            # Same deal, need to ignore the unclustered (last entry)
            clustmax = 0
            if kmean:
                for i in np.arange(len(labels)):
                    newmax = max(cluster_bins[:,i])
                    if newmax > clustmax:
                        clustmax = newmax
                    pad = i*0.8/len(labels)
                    pl.bar(np.arange(fin_clustL)+pad,cluster_bins[:,i],width=0.8/len(labels),alpha = 0.5, color = cmap_discrete[i])
            else:
                for i in np.arange(len(labels)):
                    newmax = max(cluster_bins[:-1,i])
                    if newmax > clustmax:
                        clustmax = newmax
                    pad = i*0.8/len(labels)
                    pl.bar(np.arange(fin_clustL-1)+pad,cluster_bins[:-1,i],width=0.8/len(labels),alpha = 0.5, color = cmap_discrete[i])
            pl.xlabel('Cluster #')
            pl.ylabel('Count')
            pl.legend(labels)
            if clust == 'kmean':
                subt_clust = 1
            else:
                subt_clust = 2
            for i in np.arange(fin_clustL-subt_clust):
                place = int(i) + 0.9 - 0.4/len(labels)
                pl.plot(place*np.ones(int(clustmax)),np.arange(int(clustmax)),'black')
            unclust_str = ('Number Unclustered: '+str(sum(cluster_bins[-1,:])))
            fig.savefig(this_dir + '/' + dir_name + '/'+clust+'_bars_fig.pdf',format='pdf',dpi=500)
            fig.savefig(this_dir + '/' + dir_name + '/'+clust+'_bars_fig.png',format='png',dpi=500)
            self.img7_5.source = this_dir + '/' + dir_name + '/'+clust+'_bars_fig.png'
            pl.close()

    def prep7(self):
        #this_dir = os.getcwd()
        this_dir = startDir
        # Select out our clusters of interest for downstream analysis
        global sub1_MI; global sub2_MI
        global sub1_seqs; global sub2_seqs
        global sel1; global sel2
        global skip7
        # Define these so you can reset user-defined cmaps
        global user_cmap1; global user_cmap2
        user_cmap1 = ['']; user_cmap2 = ['']
        global grad_cmap1; global grad_cmap2
        grad_cmap1 = ['']; grad_cmap2 = ['']
        # This is defined to be sure that 
        if clust == 'kmean':
            subt_clust = 1
        else:
            subt_clust = 2
        if self.ori_binary.active:
            if onlyONE:
                self.ori_binary.active = False
                sel1 = int(self.clust1.text)
                sel2 = int(self.clust2.text)
                if sel1 > fin_clustL-subt_clust or sel2 > fin_clustL-subt_clust:
                    sel1 = 0; sel2 = 1
                sub1_MI = seq_MIf[seq_MIf.columns[cluster_dset[cluster_dset['cluster'] == sel1].index]]
                sub2_MI = seq_MIf[seq_MIf.columns[cluster_dset[cluster_dset['cluster'] == sel2].index]]

                # Alright now get the sequences
                sub1_seqs = np.transpose(seq_final[sub1_MI.columns])
                sub2_seqs = np.transpose(seq_final[sub2_MI.columns])

                sub1_seqs.to_csv(this_dir + '/' + dir_name +'/clust2seq_'+str(sel1)+'.txt',header=None,index=None)
                sub2_seqs.to_csv(this_dir + '/' + dir_name +'/clust2seq_'+str(sel2)+'.txt',header=None,index=None)
                # AVOID the binary clusters screen
                skip7 = True
                popup = Popup(title='ERROR (Click Anywhere to Dismiss)',
                    content=Label(text='Why would you click that button? Revert to cluster#'),
                    size_hint=(None, None), size=(800, 800))
                popup.open()
                return

            # Bring back our "binary clusters" screen
            skip7 = False
            return
        else:
            sel1 = int(self.clust1.text)
            sel2 = int(self.clust2.text)
            if sel1 > fin_clustL-subt_clust or sel2 > fin_clustL-subt_clust:
                sel1 = 0; sel2 = 1
                popup = Popup(title='ERROR (Click Anywhere to Dismiss)',
                content=Label(text='Selected cluster out of range, default to clust0, clust1'),
                size_hint=(None, None), size=(800, 800))
                popup.open()
            # AVOID the binary clusters screen
            skip7 = True
        
        sub1_MI = seq_MIf[seq_MIf.columns[cluster_dset[cluster_dset['cluster'] == sel1].index]]
        sub2_MI = seq_MIf[seq_MIf.columns[cluster_dset[cluster_dset['cluster'] == sel2].index]]

        # Alright now get the sequences
        sub1_seqs = np.transpose(seq_final[sub1_MI.columns])
        sub2_seqs = np.transpose(seq_final[sub2_MI.columns])

        sub1_seqs.to_csv(this_dir + '/' + dir_name +'/clust2seq_'+str(sel1)+'.txt',header=None,index=None)
        sub2_seqs.to_csv(this_dir + '/' + dir_name +'/clust2seq_'+str(sel2)+'.txt',header=None,index=None)

    def get_pos_props(self):
        self.next1_8.disabled = False
        property1 = str(self.prop1_sel.text)
        property2 = str(self.prop2_sel.text)
        global sel1; global sel2

        #this_dir = os.getcwd()
        this_dir = startDir
        full_big.index = seq_MIf.columns
        global sels
        global labels_new
        # This will be the color scheme used for figures downstream of clustering
        global cmap_discrete_fin
        global cmap_discrete
        # Skip7 is our sign of "use the original labels"
        if skip7:
            # define the IDs by the datasets already subsected
            sel1_pre = seq_MIf.columns[cluster_dset[cluster_dset['cluster'] == sel1].index]
            sel2_pre = seq_MIf.columns[cluster_dset[cluster_dset['cluster'] == sel2].index]
            sel_pre = np.hstack((sel1_pre,sel2_pre))
            # So I believe this should just make a big vector with IDs that matches the vector length
            ID_vect = [sel1] * len(sel1_pre) + [sel2] * len(sel2_pre)
            sels = np.transpose(pandas.DataFrame((sel_pre,ID_vect)))
            sels.columns = ['selection','ID']

            labels_new =['Cluster '+str(sel1), 'Cluster '+str(sel2)]

            # If we're skipping 7, we need to redefine our colormap to have more entries...
            # Redefine as cmap_discrete_fin so it doesn't mess with the previously
            # defined cmap_discrete (based on # datasets)
            #if clust == 'kmean':
            
            # I think that no matter what, we want this...
            cmap2 = pl.get_cmap('rainbow')
            cmap_discrete_fin = cmap2(np.linspace(0, 1, fin_clustL))
            #else:
                #cmap_discrete_fin = cmap(np.linspace(0, 1, fin_clustL))
        else:
            cmap2 = pl.get_cmap('rainbow')
            cmap_discrete_fin = cmap2(np.linspace(0, 1, len(seqNameF)))
            # try to sort of cheat in creating the labels here...
            labels_new = [''] * len(group_a_id)
            first = True; second = False
            for i in np.arange(len(group_a_id)):
                # lda_checked is now a "dump this data" variable
                # This madness with first and second has to do with 
                # the color selection options... will need to clean up eventually...
                if len(group_a_id) > 2:
                    if lda_checked[i]:
                        continue
                if first:
                    sel_pre = [column for column in seq_MIf.columns if seqNameF[i][0] in column]
                    ID_vect = [int(group_a_id[i][0])] * len(sel_pre)
                    labels_new[int(group_a_id[i][0])] = seqNameF[i][0]
                    first = False
                    second = True
                    continue
                if second:
                    temp_sel = [column for column in seq_MIf.columns if seqNameF[i][0] in column]
                    sel_pre = np.hstack((sel_pre,temp_sel))
                    ID_vect = ID_vect + [int(group_a_id[i][0])] * len(temp_sel)
                    if labels_new[int(group_a_id[i][0])] == '':
                        labels_new[int(group_a_id[i][0])] = seqNameF[i][0]
                    else:
                        labels_new[int(group_a_id[i][0])] = labels_new[int(group_a_id[i][0])] + ' + ' + seqNameF[i][0]
                    second = False
                else:
                    temp_sel = [column for column in seq_MIf.columns if seqNameF[i][0] in column]
                    sel_pre = np.hstack((sel_pre,temp_sel))
                    ID_vect = ID_vect + [int(group_a_id[i][0])] * len(temp_sel)
                    if labels_new[int(group_a_id[i][0])] == '':
                        labels_new[int(group_a_id[i][0])] = seqNameF[i][0]
                    else:
                        labels_new[int(group_a_id[i][0])] = labels_new[int(group_a_id[i][0])] + ' + ' + seqNameF[i][0]
            sels = np.transpose(pandas.DataFrame((sel_pre,ID_vect)))
            sels.columns = ['selection','ID']
            labels_new = [a for a in labels_new if a != '']

        # Let users define their own colormaps. This should also allow for an easy reset (I hope)
        if len(user_cmap1) != 1:
            cmap_discrete_fin[sel1] = user_cmap1
        if len(user_cmap2) != 1:
            cmap_discrete_fin[sel2] = user_cmap2

        # Also let users define the properties
        prop_dict = {'Charge':1,'Hydropathy':2, 'Flexibility':4, 'Bulk':3}
        if property1 == 'Property1':
            prop1 = 1
            property1 = 'Charge'
        else: 
            prop1 = prop_dict[property1]
        if property2 == 'Property2':
            prop2 = 2
            property2 = 'Hydropathy'
        else:
            prop2 = prop_dict[property2]

        # Alright we need to change the way we plot these too:
        fig, ax = pl.subplots(2, 1,squeeze=False,figsize=(16,8))
        for j in sels['ID'].drop_duplicates():
            findex = sels[sels['ID'] == j]['selection']
            pre_array = np.array(full_big.loc[findex])
            dat_size = len(findex)
            seq_bigF = pre_array.reshape(dat_size,61,len(seq_MIf))
            # If we need the sequences, we can call them this same way but
            # with the variable seq_final rather than full_big
            if self.boot_pos_sen.active and self.vis_indi.active:
                popup = Popup(title='ERROR (Click Anywhere to Dismiss)',
                    content=Label(text='Can only click one of the boxes on this page'),
                    size_hint=(None, None), size=(600, 600))
                popup.open()
                return()
            elif self.boot_pos_sen.active:
                # For now just make everyone do 1000 bootstrap replicates...
                # Should eventually make this tunable...
                boots = 1000
                prop_avg1 = []; prop_avg2 = []
                for i in np.arange(boots):
                    re_big = resample(seq_bigF)
                    prop_avg1.append(np.average(re_big[:,prop1,:],axis=0))
                    prop_avg2.append(np.average(re_big[:,prop2,:],axis=0))
                fin_avg1 = np.average(prop_avg1,axis=0); fin_std1 = np.std(prop_avg1,axis=0)
                fin_avg2 = np.average(prop_avg2,axis=0); fin_std2 = np.std(prop_avg2,axis=0)
                ax[0,0].plot(fin_avg1,marker='o',linewidth=2.5,color=cmap_discrete_fin[j])
                ax[0,0].fill_between(np.arange(len(fin_avg1)),fin_avg1+fin_std1,fin_avg1-fin_std1,alpha=0.3,color=cmap_discrete_fin[j])
                ax[1,0].plot(fin_avg2,marker='o',linewidth=2.5,color=cmap_discrete_fin[j])
                ax[1,0].fill_between(np.arange(len(fin_avg2)),fin_avg2+fin_std2,fin_avg2-fin_std2,alpha=0.3,color=cmap_discrete_fin[j])
            elif self.vis_indi.active:
                for oneLine in np.arange(len(seq_bigF[:,prop1,:])):
                    ax[0,0].plot(seq_bigF[oneLine,prop1,:],linewidth=2.5,color=cmap_discrete_fin[j],alpha=0.3)
                    ax[1,0].plot(seq_bigF[oneLine,prop2,:],linewidth=2.5,color=cmap_discrete_fin[j],alpha=0.3)
                plotProp1 = np.average(seq_bigF[:,prop1,:],axis = 0)
                ax[0,0].plot(plotProp1,marker='o',linewidth=2.5,color=cmap_discrete_fin[j])
                plotProp2 = np.average(seq_bigF[:,prop2,:],axis = 0)
                ax[1,0].plot(plotProp2,marker='o',linewidth=2.5,color=cmap_discrete_fin[j])
            else:
                plotProp1 = np.average(seq_bigF[:,prop1,:],axis = 0)
                ax[0,0].plot(plotProp1,marker='o',linewidth=2.5,color=cmap_discrete_fin[j])
                plotProp2 = np.average(seq_bigF[:,prop2,:],axis = 0)
                ax[1,0].plot(plotProp2,marker='o',linewidth=2.5,color=cmap_discrete_fin[j])
            np.savetxt(this_dir + '/' + dir_name + '/position_sensitive_mat_clust'+str(j)+'.dat',pre_array,fmt='%.3f')

        ax[0,0].set_ylabel(property1)
        ax[1,0].set_ylabel(property2)
        ax[1,0].set_xlabel('Sequence Position')
        if self.vis_indi.active or self.boot_pos_sen.active:
            legend_elements=[]
            for j in np.arange(len(labels_new)):
                element = [Line2D([0], [0], marker='o', color='w', label=labels_new[j],markerfacecolor=cmap_discrete_fin[j], markersize=10)]
                legend_elements = legend_elements+element
            pl.legend(handles=legend_elements,ncol=len(labels_new))
        else:
            pl.legend(labels_new)

        if molecule == 'mhc':
            ax[0,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
            ax[0,0].set_xticklabels(['Strand 1','Helix 1','Strand 2','Helix 2'])
            ax[1,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
            ax[1,0].set_xticklabels(['Strand 1','Helix 1','Strand 2','Helix 2'])
        elif molecule == 'ig':
            if LOOPnum == 6:
                ax[0,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[0,0].set_xticklabels(['CDR1L','CDR2L','CDR3L','CDR1H','CDR2H','CDR3H'])
                ax[1,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[1,0].set_xticklabels(['CDR1L','CDR2L','CDR3L','CDR1H','CDR2H','CDR3H'])
            elif LOOPnum == 3:
                ax[0,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[0,0].set_xticklabels(['CDR1','CDR2','CDR3'])
                ax[1,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[1,0].set_xticklabels(['CDR1','CDR2','CDR3'])
            elif LOOPnum == 2:
                ax[0,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[0,0].set_xticklabels(['CDR3H','CDR3L'])
                ax[1,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[1,0].set_xticklabels(['CDR3H','CDR3L'])
            elif LOOPnum == 1:
                ax[0,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[0,0].set_xticklabels(['CDR Loop'])
                ax[1,0].set_xticks(np.array(xtick_loc).reshape(len(xtick_loc)))
                ax[1,0].set_xticklabels(['CDR Loop'])
        # NOTE NOTE: For now lets not annotate any of the peptide/MSA axes
        # Maybe at some point in the future I want to change that, but not now.

        fig.savefig(this_dir + '/' + dir_name + '/pos_prop'+str(prop1)+str(prop2)+'.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/pos_prop'+str(prop1)+str(prop2)+'.png',format='png',dpi=500)

        self.img8.source = this_dir + '/' + dir_name + '/pos_prop'+str(prop1)+str(prop2)+'.png'
        pl.close()

    def get_clone_pos_props(self):
        self.next1_9.disabled = False
        property1 = str(self.prop1_sel1.text)
        #this_dir = os.getcwd()
        this_dir = startDir
        # Also let users define the properties
        prop_dict = {'Charge':1,'Hydropathy':2, 'Flexibility':4, 'Bulk':3}
        if property1 == 'Property':
            prop1 = 1
            property1 = 'Charge'
        else: 
            prop1 = prop_dict[property1]
        # Generate the position sensitive charge across all clones in the dataset
        num_figs = int(np.ceil(len(sels['ID'].drop_duplicates())/2))
        fig, axs = pl.subplots(num_figs, 2,squeeze=False,figsize=(20,8))
        # Lets get some more control over plot colors
        global grad_cmap1
        global grad_cmap2
        global gradMap
        if grad_cmap1== [''] or grad_cmap2 == ['']:
            gradMap = cm.PiYG
        else:
            # Need to take off the transparency
            # pretty sure this is all 1 on the wheel anyway
            c1=grad_cmap1[:3]
            c2=[1,1,1]
            c3=grad_cmap2[:3]
            grad_list = get_color_gradient(c1,c2,c3,n=1000)
            gradMap = LinearSegmentedColormap.from_list('myGrad', grad_list, N=1000)

        fig_track = 0; track2 = 0
        aa = 0
        for j in sels['ID'].drop_duplicates():
            findex = sels[sels['ID'] == j]['selection']
            pre_array = np.array(full_big.loc[findex])
            dat_size = len(findex)
            seq_bigF = pre_array.reshape(dat_size,61,len(seq_MIf))
            # Some tricks to make sure 0 is at the center of the colorbar
            min_temp = np.min(seq_bigF[:,prop1,:])
            max_temp = np.max(seq_bigF[:,prop1,:])
            if abs(min_temp) > abs(max_temp):
                propMin = min_temp
                propMax = -min_temp
            elif abs(max_temp) > abs(min_temp):
                propMin = -max_temp
                propMax = max_temp
            else:
                propMin = min_temp
                propMax = max_temp
            # If we need the sequences, we can call them this same way but
            # with the variable seq_final rather than full_big
            x = axs[fig_track,track2].imshow(seq_bigF[:,prop1,:],interpolation='nearest', aspect='auto',
            cmap=gradMap,vmin=propMin,vmax=propMax)
            # NEED TO CHANGE THIS SO IT ISNT DEFINED FOR EVERY FIGURE
            axs[fig_track,track2].set_title(labels_new[aa] + ' - '+property1)
            aa += 1
            fig.colorbar(x, ax=axs[fig_track,track2])
            if fig_track == num_figs:
                axs[fig_track,track2].set_xlabel('Sequence Position')
            axs[fig_track,track2].set_ylabel('Sequence Number')

            if track2 == 0:
                track2 = 1
            else:
                fig_track += 1
                track2 = 0

        fig.savefig(this_dir + '/' + dir_name + '/clone_pos_'+property1+'.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/clone_pos_'+property1+'.png',format='png',dpi=500)
        self.img9.source = this_dir + '/' + dir_name + '/clone_pos_'+property1+'.png'
        pl.close()

    def get_props(self):
        do_stats = True
        self.next1_10.disabled = False
        #this_dir = os.getcwd()
        this_dir = startDir
        # Generate the position sensitive charge across all clones in the dataset
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
        x_axis = np.array([-0.2,0.9,2,3.1])
        # Need to have some kind of color wheel to replace this...
        # We want to exclude prop0 (the simple 1-21 AA representation entries)
        aa = 0
        for j in sels['ID'].drop_duplicates():
            findex = sels[sels['ID'] == j]['selection']
            pre_array = np.array(full_big.loc[findex])
            dat_size = len(findex)
            seq_bigF = pre_array.reshape(dat_size,61,len(seq_MIf))
            # If we need the sequences, we can call them this same way but
            # with the variable seq_final rather than full_big
            if self.boot_prop.active:
                # For now just make everyone do 1000 bootstrap replicates...
                # Should eventually make this tunable...
                boots = 1000
                prop_avg1 = []; prop_avg2 = []
                prop_avg3 = []; prop_avg4 = []
                for i in np.arange(boots):
                    re_big = resample(seq_bigF)
                    prop_avg1.append(np.average(re_big[:,1,:]))
                    prop_avg2.append(np.average(re_big[:,2,:]))
                    prop_avg3.append(np.average(re_big[:,3,:]))
                    prop_avg4.append(np.average(re_big[:,4,:]))
                plotProp1 = np.average(prop_avg1,axis=0); stdProp1 = np.std(prop_avg1,axis=0)
                plotProp2 = np.average(prop_avg2,axis=0); stdProp2 = np.std(prop_avg2,axis=0)
                plotProp3 = np.average(prop_avg3,axis=0); stdProp3 = np.std(prop_avg3,axis=0)
                plotProp4 = np.average(prop_avg4,axis=0); stdProp4 = np.std(prop_avg4,axis=0)
            else:
                plotProp1 = np.average(np.average(seq_bigF[:,1,:],axis = 1))
                plotProp2 = np.average(np.average(seq_bigF[:,2,:],axis = 1))
                plotProp3 = np.average(np.average(seq_bigF[:,3,:],axis = 1))
                plotProp4 = np.average(np.average(seq_bigF[:,4,:],axis = 1))
                stdProp1 = np.std(np.average(seq_bigF[:,1,:],axis = 1))
                stdProp2 = np.std(np.average(seq_bigF[:,2,:],axis = 1))
                stdProp3 = np.std(np.average(seq_bigF[:,3,:],axis = 1))
                stdProp4 = np.std(np.average(seq_bigF[:,4,:],axis = 1))
            compile_avg = [plotProp1,plotProp2,plotProp3,plotProp4]
            compile_std = [stdProp1,stdProp2,stdProp3,stdProp4]
            if aa == 0:
                compile_names = ['Charge','Hydrophobicity','Bulkiness','Flexibility']
                save_these = [compile_names,compile_avg,compile_std]
            else:
                save_these = save_these + [compile_avg,compile_std]
            plotIT = np.hstack((plotProp1, plotProp2,plotProp3,plotProp4))
            stdIT = np.hstack((stdProp1, stdProp2,stdProp3,stdProp4))
            ax[0,0].bar(x_axis+aa*1/len(labels_new), plotIT,
                        yerr = stdIT,alpha = 0.5, width = 1/len(labels_new),color=cmap_discrete_fin[j])
            aa += 1

        if do_stats:
            first = True
            # Calculate Welch's t-test for the statistics in the plot
            # NOTE: This is essentially Student's T test, but when
            # we cannot assume equal variance or sample size
            # The null hypothesis here is that the means are equal
            temp_labels = sels['ID'].drop_duplicates().values
            for a in np.arange(len(sels['ID'].drop_duplicates())):
                for b in np.arange(len(sels['ID'].drop_duplicates())):
                    # Only have a trailing b
                    if a >= b:
                        continue
                    findex1 = sels[sels['ID'] == temp_labels[a]]['selection']
                    findex2 = sels[sels['ID'] == temp_labels[b]]['selection']
                    dat_size1 = len(findex1); dat_size2 = len(findex2)
                    for cc in np.arange(4):
                        # temp_prop will cycle through charge, phob, bulk, flex
                        # in that particular order
                        temp_prop = np.transpose(save_these)[cc]
                        # So it seems like the data is saved linearly. avg, std, avg, std, etc
                        mean1 = float(temp_prop[a*2+1]); mean2 = float(temp_prop[b*2+1])
                        std1 = float(temp_prop[a*2+2]); std2 = float(temp_prop[b*2+2])
                        t = (mean1-mean2)/np.sqrt(std1**2/dat_size1 + std2**2/dat_size2)
                        
                        # Two sided t-test here, null: "the means are not equal"
                        # This scipy is just the t-statistic lookup table
                        # Here using dof = Nsample-1 (Nsample=2)
                        fin_p = scipy.stats.t.sf(abs(t), df=1)*2
                        
                        if first:
                            fin_ptest = pandas.DataFrame([compile_names[cc],labels_new[a],labels_new[b],fin_p])
                            first = False
                        else:
                            fin_ptest = pandas.concat([fin_ptest,pandas.DataFrame([compile_names[cc],labels_new[a],labels_new[b],fin_p])],axis=1)
            final_ptest = np.transpose(fin_ptest)
            final_ptest.columns = [['Property','Dset1','Dset2','p-value']]
            final_ptest.to_csv(this_dir + '/' + dir_name + '/avg_props_stats.csv',index=False)
        ax[0,0].legend(labels_new)
        ax[0,0].set_xticks([0.2,1.3,2.4,3.5])
        ax[0,0].set_xticklabels(['Charge','Hydrophobicity','Bulkiness','Flexibility'])
        ax[0,0].set_xlabel('Biophysical Property')
        ax[0,0].set_ylabel('Normalized Property Value')

        fig.savefig(this_dir + '/' + dir_name + '/avg_props.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/avg_props.png',format='png',dpi=500)
        self.img10.source = this_dir + '/' + dir_name + '/avg_props.png'
        np.savetxt(this_dir + '/' + dir_name + '/avg_props.dat',save_these,fmt='%s')
        pl.close()

    def do_lda(self):
        #this_dir = os.getcwd()
        this_dir = startDir
        numVects = int(self.inputLDA.text)

        findex1 = sels[sels['ID'] == sels['ID'].drop_duplicates().values[0]]['selection']
        findex2 = sels[sels['ID'] == sels['ID'].drop_duplicates().values[1]]['selection']
        pre_array1 = np.array(full_big.loc[findex1])
        pre_array2 = np.array(full_big.loc[findex2])
        sub1_seqs = np.transpose(seq_final[findex1])
        sub2_seqs = np.transpose(seq_final[findex2])

        pg1 = np.transpose(sub1_seqs.values); num1 = np.shape(pg1)[1]
        pg2 = np.transpose(sub2_seqs.values); num2 = np.shape(pg2)[1]
        total_mat = np.vstack((pre_array1,pre_array2))
        new_big, weights, acc_all, mda_all, parsed_mat, top_names = aims.compile_MP(total_mat, pg1, pg2, final_size = numVects, cat = False)

        # Seaborn plots look nicer for these LDA figures
        import seaborn as sns
        fig = pl.figure(figsize = (12, 12))
        dset = ["Linear Discriminant Analysis" for x in range(num1+num2)]
        reacts = [labels_new[0] for x in range(num1)] + [labels_new[1] for x in range(num2)]

        d1 = {'Dataset': dset, 'Linear Discriminant 1': mda_all.reshape(len(mda_all)),
            'Class' : reacts}
        df1 = pandas.DataFrame(data=d1)
        sns.set(style="white", color_codes=True,font_scale=1.5)
        sns.swarmplot(x="Dataset", y="Linear Discriminant 1", data=df1, hue = 'Class', palette = "Dark2")

        fig.savefig(this_dir + '/' + dir_name + '/lda.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/lda.png',format='png',dpi=500)
        pl.close()
        # And save the raw data
        np.savetxt(this_dir + '/' + dir_name + '/lda_data.dat',mda_all,fmt='%.3f')

        # Need to put acc_all somewhere on the screen
        indices_topProps = np.argsort(abs(weights))
        compile_Props = np.vstack((top_names,weights))
        final_pc1 = compile_Props[:,indices_topProps[0][-10:]]
        print(acc_all)
        np.savetxt(this_dir + '/' + dir_name + '/lda_weights.dat',compile_Props,fmt='%s')

        # For running multiple times with vectors less than 10, need to reset screen...
        self.aname10.text = ''; self.aname9.text = ''; self.aname8.text = ''; self.aname7.text = ''
        self.aname6.text = ''; self.aname5.text = ''; self.aname4.text = ''; self.aname3.text = ''
        self.aname2.text = ''; self.aname1.text = ''
        self.aweight10.text = ''; self.aweight9.text = ''; self.aweight8.text = ''; self.aweight7.text = ''
        self.aweight6.text = ''; self.aweight5.text = ''; self.aweight4.text = ''; self.aweight3.text = ''
        self.aweight2.text = ''; self.aweight1.text = ''
        # Update the weights and labels onto the screen:
        self.aname10.text =str(final_pc1[0,0]); self.aweight10.text= str(round(final_pc1[1,0],4))
        if numVects > 1:
            self.aname9.text = str(final_pc1[0,1]); self.aweight9.text = str(round(final_pc1[1,1],4))
        if numVects > 2:
            self.aname8.text = str(final_pc1[0,2]); self.aweight8.text = str(round(final_pc1[1,2],4))
        if numVects > 3:
            self.aname7.text = str(final_pc1[0,3]); self.aweight7.text = str(round(final_pc1[1,3],4))
        if numVects > 4:
            self.aname6.text = str(final_pc1[0,4]); self.aweight6.text = str(round(final_pc1[1,4],4))
        if numVects > 5:
            self.aname5.text = str(final_pc1[0,5]); self.aweight5.text = str(round(final_pc1[1,5],4))
        if numVects > 6:
            self.aname4.text = str(final_pc1[0,6]); self.aweight4.text = str(round(final_pc1[1,6],4))
        if numVects > 7:
            self.aname3.text = str(final_pc1[0,7]); self.aweight3.text = str(round(final_pc1[1,7],4))
        if numVects > 8:
            self.aname2.text = str(final_pc1[0,8]); self.aweight2.text = str(round(final_pc1[1,8],4))
        if numVects > 9:
            self.aname1.text = str(final_pc1[0,9]); self.aweight1.text = str(round(final_pc1[1,9],4))

        self.img14.source = this_dir + '/' + dir_name + '/lda.png'

    def get_freq(self):
        self.next1_13.disabled = False
        #this_dir = os.getcwd()
        this_dir = startDir
        global freq1,freq2
        global mat_size
        AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        freq1 = freq_hold[0]
        freq2 = freq_hold[1]
        # Lets get some more control over plot colors
        global grad_cmap1
        global grad_cmap2
        global gradMap
        if grad_cmap1== [''] or grad_cmap2 == ['']:
            gradMap = cm.PiYG
        else:
            # Need to take off the transparency
            # pretty sure this is all 1 on the wheel anyway
            c1=grad_cmap1[:3]
            c2=[1,1,1]
            c3=grad_cmap2[:3]
            grad_list = get_color_gradient(c1,c2,c3,n=1000)
            gradMap = LinearSegmentedColormap.from_list('myGrad', grad_list, N=1000)

        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
        freqMax = np.max(freq1[:,1:]-freq2[:,1:]); freqMin = np.min(freq1[:,1:]-freq2[:,1:])
        freqBound = max(abs(freqMax),abs(freqMin))
        x=ax[0,0].pcolormesh(freq1[:,1:]-freq2[:,1:],cmap=gradMap,vmin = -freqBound, vmax = freqBound)
        pl.colorbar(x)
        pl.ylabel('Sequence Position')
        xax=pl.setp(ax,xticks=np.arange(20)+0.5,xticklabels=AA_key)
        pl.title(labels_new[0]+ ' AA Freq - ' + labels_new[1] + ' AA Freq')

        place=0
        if type(mat_size) != int:
            for i in mat_size:
                place += i
                pl.plot(np.arange(21),place*np.ones(21),'black')
        
        # Since there's only two now, just easier to hard code these...
        fig.savefig(this_dir + '/' + dir_name + '/frequency.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/frequency.png',format='png',dpi=500)
        # And save the raw data
        np.savetxt(this_dir + '/' + dir_name + '/frequency_mat1.dat',freq1,fmt='%.3f')
        np.savetxt(this_dir + '/' + dir_name + '/frequency_mat2.dat',freq2,fmt='%.3f')

        self.img13.source = this_dir + '/' + dir_name + '/frequency.png'
        pl.close()

    def get_shannon(self):
        self.next1_11.disabled = False
        #this_dir = os.getcwd()
        this_dir = startDir
        global freq_hold
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
        shannon_hold = np.zeros((len(labels_new),len(seq_MIf)))
        freq_hold = np.zeros((len(labels_new),len(seq_MIf),21))
        aa = 0
        for j in sels['ID'].drop_duplicates():
            findex = sels[sels['ID'] == j]['selection']
            sub_MI = seq_MIf[findex]
            shannon_hold[aa],freq_hold[aa],coverage = aims.calculate_shannon(np.transpose(np.array(sub_MI)))
            pl.plot(shannon_hold[aa], marker='o',linewidth=2.5,color=cmap_discrete_fin[j])
            aa += 1

        pl.legend(labels_new); pl.xlabel('Position'); pl.ylabel('Shannon Entropy (Bits)')
        # Entropy is a rare case where we know the exact bounds of the values...
        # Comment this out if your max entropy is substantially lower
        pl.ylim([0,4.12])
        # Guide the eyes with these lines
        if type(mat_size) != int:
            for i in np.arange(len(mat_size)-1):
                ax[0,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(0,4.2,100),'black',linewidth = 3)
        
        global LOOPnum
        global xtick_loc
        ax[0,0].legend(labels_new)
        # Since there's only two now, just easier to hard code these...
        pl.ylabel('Shannon Entropy (Bits)')
        fig.savefig(this_dir + '/' + dir_name + '/shannon.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/shannon.png',format='png',dpi=500)
        # And save the raw data

        self.img11.source = this_dir + '/' + dir_name + '/shannon.png'
        pl.close()
        # Don't let them go on if not a binary problem
        if len(labels_new) > 2:
            self.next1_11.disabled = True

    def get_MI(self):
        ylim_min = float(self.min_mi.text)
        ylim_max = float(self.max_mi.text)
        if ylim_min > ylim_max:
            popup = Popup(title='ERROR (Click Anywhere to Dismiss)',
                    content=Label(text='MI min > MI max'),
                    size_hint=(None, None), size=(600, 600))
            popup.open()
            return()
        self.next1_12.disabled = False
        #this_dir = os.getcwd()
        this_dir = startDir
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(10,8))
        poses = len(seq_MIf)
        # Lets get some more control over plot colors
        global grad_cmap1
        global grad_cmap2
        global gradMap
        if grad_cmap1== [''] or grad_cmap2 == ['']:
            gradMap = cm.PiYG
        else:
            # Need to take off the transparency
            # pretty sure this is all 1 on the wheel anyway
            c1=grad_cmap1[:3]
            c2=[1,1,1]
            c3=grad_cmap2[:3]
            grad_list = get_color_gradient(c1,c2,c3,n=1000)
            gradMap = LinearSegmentedColormap.from_list('myGrad', grad_list, N=1000)
        # So to get this far we are REQUIRING Binary entries
        findex1 = sels[sels['ID'] == sels['ID'].drop_duplicates().values[0]]['selection']
        sub1_MI = seq_MIf[findex1]
        findex2 = sels[sels['ID'] == sels['ID'].drop_duplicates().values[1]]['selection']
        sub2_MI = seq_MIf[findex2]
            
        MI1,entropy_cond1,counted1=aims.calculate_MI(np.transpose(np.array(sub1_MI)))
        MI2,entropy_cond2,counted2=aims.calculate_MI(np.transpose(np.array(sub2_MI)))
        x = pl.imshow(MI1 - MI2, cmap=gradMap, vmin = ylim_min, vmax = ylim_max)
        pl.colorbar(x); pl.title(labels_new[0]+ ' MI - ' + labels_new[1] + ' MI')

        # Help Guide the eyes a bit
        if type(mat_size) != int:
            for i in np.arange(len(mat_size)-1):
                ax[0,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(0,poses,100),'black',linewidth = 3)
                ax[0,0].plot( np.linspace(0,poses,100), (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100) ,'black',linewidth = 3)
        
        pl.xlabel('Sequence Position'); pl.ylabel('Sequence Position')
        fig.savefig(this_dir + '/' + dir_name + '/MI.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/MI.png',format='png',dpi=500)
        # And save the raw data
        np.savetxt(this_dir + '/' + dir_name + '/MI_mat1.dat',MI1,fmt='%.3f')
        np.savetxt(this_dir + '/' + dir_name + '/MI_mat2.dat',MI2,fmt='%.3f')
        self.img12.source = this_dir + '/' + dir_name + '/MI.png'
        pl.close()

    def show_colorWheel1(self):
        def on_color(instance, value):
            global user_cmap1
            user_cmap1 = value
        content = ColorPicker()
        self._popup = Popup(title="Pick Color1, Click Anywhere to Exit", content=content,
                            size_hint=(0.8, 0.8))
        content.bind(color=on_color)
        self._popup.open()

    def show_colorWheel2(self):
        def on_color(instance, value):
            global user_cmap2
            user_cmap2 = value
        content = ColorPicker()
        self._popup = Popup(title="Pick Color1, Click Anywhere to Exit", content=content,
                            size_hint=(0.8, 0.8))
        content.bind(color=on_color)
        self._popup.open()

    def reset_userColor(self):
        global user_cmap1
        global user_cmap2
        user_cmap1 = ['']
        user_cmap2 = ['']

    # Define all the gradient functions separately so we can
    # separately define and reset the color options
    def show_gradWheel1(self):
        def on_color(instance, value):
            global grad_cmap1
            grad_cmap1 = value
        content = ColorPicker()
        self._popup = Popup(title="Pick Color1, Click Anywhere to Exit", content=content,
                            size_hint=(0.8, 0.8))
        content.bind(color=on_color)
        self._popup.open()

    def show_gradWheel2(self):
        def on_color(instance, value):
            global grad_cmap2
            grad_cmap2 = value
        content = ColorPicker()
        self._popup = Popup(title="Pick Color1, Click Anywhere to Exit", content=content,
                            size_hint=(0.8, 0.8))
        content.bind(color=on_color)
        self._popup.open()

    def reset_gradColor(self):
        global gradMap
        global grad_cmap1
        global grad_cmap2
        grad_cmap1 = ['']
        grad_cmap2 = ['']
        gradMap=cm.PiYG

class checker(Screen):
    # Need to redefine a different type of loading from the wild one above.

    def on_pre_enter(self, *args):
        global check_status
        #### NOTE, BOX 3 has been removed... no longer can exclude points from the plot.
        if 'check_run' not in globals():
            check_L1 = Label(text='Include in Clustering',size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.15, 'center_y': 0.6},font_name='app_data/Poppins-Medium.ttf')
            FloatLayout.add_widget(self, check_L1)
            # only add in these you have MORE than 2 options
            if onlyONE:
                pass
            elif len(labels) > 2:
                check_L2 = Label(text='Exclude from Clustering',size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.15, 'center_y': 0.45},font_name='app_data/Poppins-Medium.ttf')
                FloatLayout.add_widget(self, check_L2)
                #check_L3 = Label(text='Exclude from Plot',size_hint = (None, None), height = '48dp',
                #    pos_hint = {'center_x': 0.15, 'center_y': 0.3},font_name='app_data/Poppins-Medium.ttf')
                #check_L4 = Label(text='(Optional)',size_hint = (None, None), height = '48dp',
                #    pos_hint = {'center_x': 0.15, 'center_y': 0.25},font_name='app_data/Poppins-Medium.ttf')
                #FloatLayout.add_widget(self, check_L3)
                #FloatLayout.add_widget(self, check_L4)
            if onlyONE:
                names = Label(text=labels, size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.3+int(0)*0.6/N, 'center_y': 0.70},font_name='app_data/Poppins-Medium.ttf')
                box1 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(0),state = 'down',
                pos_hint = {'center_x': 0.3+int(0)*0.6/N, 'center_y': 0.6}, allow_no_selection = False)
                check_status = [[names,box1]]
                # Before adding in all of the new ones
            else:
                for j in np.arange(len(labels)):
                    names = Label(text=labels[j], size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.70},font_name='app_data/Poppins-Medium.ttf')
                    if len(labels) > 2:
                        box1 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j),state = 'down',
                        pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.6},allow_no_selection = False)
                        box2 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j),
                        pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.45},allow_no_selection = False)
                        #box3 = CheckBox(size_hint = (None,None), height = '48dp',
                        #pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.3},allow_no_selection = False)
                    else:
                        box1 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j),state = 'down',
                        pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.6},allow_no_selection = False)
                    if int(j) == 0:
                        if len(labels) > 2:
                            check_status = [[names,box1,box2]]#,box3]]
                        else:
                            check_status = [[names,box1]]
                    else:
                        if len(labels) > 2:
                            check_status = check_status + [[names,box1,box2]]#,box3]]
                        else:
                            check_status = check_status + [[names,box1]]
                    # Before adding in all of the new ones
            for entry in check_status:
                for j in entry:
                    FloatLayout.add_widget(self, j)
            # Make sure these boxes are deleted and recreated...
            global check_run
            check_run = True

    def check_checks(self):
        global checked
        global check_status
        # Convert our weird little kivy objects into a numpy array
        x,y = np.shape(check_status)
        for row in np.arange(x):
            for column in np.arange(y):
                # Exclude the name column
                if column == 0:
                    continue
                if column == 1:
                    check_pre = check_status[row][column].active
                else:
                    check_pre = np.hstack((check_pre,check_status[row][column].active))
            if row == 0:
                checked = check_pre
            else:
                checked = np.vstack((checked,check_pre))

    def get_bigMat(self):
        self.next1_4.disabled = False
        dsetF = seq_final.values
        global parsed_mat
        global full_big
        global prop_names
        bigass = classy.get_bigass_matrix(dsetF, giveSize = mat_size, alignment = align )
        # This definition will be important if I get parallel processing into the app... see notebook for EX
        total_mat = bigass
        prop_list_old = ['Phobic1','Charge','Phobic2','Bulk','Flex','Kid1','Kid2','Kid3','Kid4','Kid5','Kid6','Kid7','Kid8','Kid9','Kid10']
        prop_list_new = ['Hot'+str(b+1) for b in range(46)]

        prop_names = prop_list_old + prop_list_new
        num_locs = int(np.shape(total_mat)[1]/61)
        Bigass_names = []
        for i in prop_names:
            for j in np.arange(num_locs):
                Bigass_names = Bigass_names + [ i + '-' + str(j) ]

        # AND THEN GO ON WITH THE REST OF THE dropping of correlated vectors
        full_big = pandas.DataFrame(total_mat,columns = Bigass_names)
        drop_zeros = [column for column in full_big.columns if all(full_big[column] == 0 )]
        y = full_big.drop(full_big[drop_zeros], axis=1)
        z_pre = np.abs(np.corrcoef(np.transpose(y)))
        z = pandas.DataFrame(z_pre,columns=y.columns,index=y.columns)
        # Select upper triangle of correlation matrix
        upper = z.where(np.triu(np.ones(z.shape), k=1).astype(bool))
        to_drop = [column for column in upper.columns if ( any(upper[column] > 0.75) ) ]

        parsed_mat = y.drop(y[to_drop], axis=1)

class aligner_ab(Screen):
    
    def make_path(self):
        global dir_name
        self.text1 = self.input1.text
        dir_name = self.text1
    
    def on_pre_enter(self, *args):
        global ab_check
        global alignAb
        global labels
        if 'ab_check' not in globals():
            for j in np.arange(len(LFile)):
                name = 'File '+str(j)
                y_val = float(0.50 - 0.1*np.floor(int(j)/5))
                x_val = float(0.1 + int(j)%5*0.9/(5))

                textinputAb = TextInput(text=name, multiline=False, size_hint = (None, None),write_tab =False,
                    pos_hint = {'center_x': x_val, 'center_y': y_val}, height = '32dp',width='125dp')
                
                if int(j) == 0:    
                    alignAb = [textinputAb]
                else:
                    alignAb = alignAb + [textinputAb]

            # Before adding in all of the new ones
            for entry in alignAb:
                    FloatLayout.add_widget(self, entry)
        else:
            for entry in alignAb:
                    FloatLayout.remove_widget(self, entry)
            # HERE, why are we going back to "File" when clearly this page
            # has already been run? change this to read in "labels" or something
            for j in np.arange(len(LFile)):
                if 'labels' in globals():
                    if len(labels) > j:
                        name = labels[j]
                    else:
                        name = 'File '+str(j)
                else:
                    name = 'File '+str(j)
                y_val = float(0.50 - 0.1*np.floor(int(j)/5))
                x_val = float(0.1 + int(j)%5*0.9/(5))

                textinputAb = TextInput(text=name, multiline=False, size_hint = (None, None),write_tab =False,
                    pos_hint = {'center_x': x_val, 'center_y': y_val }, height = '32dp',width='125dp')

                if int(j) == 0:    
                    alignAb = [textinputAb]
                else:
                    alignAb = alignAb + [textinputAb]

            # Before adding in all of the new ones
            for entry in alignAb:
                    FloatLayout.add_widget(self, entry)
            # Make sure these boxes are deleted and recreated...
        global ab_check
        ab_check = True

    def check_aligns(self):
        global labels
        global ab_check
        x = len(alignAb)
        for row in np.arange(x):
            if row == 0:
                labels = alignAb[row].text
            else:
                labels = np.vstack((labels,alignAb[row].text))
        # Get labels out of the GD weird format.
        if onlyONE:
            pass
        else:
            labels = [i[0] for i in labels]
            # All of this here is to fix issues in naming
            # this fix will also break if users for some reason end filenames
            # with "___". Better hope they don't do that
            dum_len = len(labels)
            for i in np.arange(dum_len):
                for j in np.arange(dum_len):
                    if i == j:
                        continue
                    if labels[i].find(labels[j]) == 0:
                        labels[j] = labels[j]+'___'
        # Figure out what number of loops we are doing
        global LOOPnum
        if self.cdr1.active == True:
            LOOPnum = 1
        elif self.cdr2.active == True:
            LOOPnum = 2
        elif self.cdr3.active == True:
            LOOPnum = 3
        elif self.cdr6.active == True:
            LOOPnum = 6
        global exp_drop
        exp_drop = self.degen_drop.active

# Try to make a general aligner screen for MSA/Peptide
# Since we won't need to do all that much extra for these...
class aligner_gen(Screen):
    
    def make_path(self):
        global dir_name
        self.text1 = self.input1.text
        dir_name = self.text1
    
    def on_pre_enter(self, *args):
        global gen_check
        global alignGen
        global labels
        if 'gen_check' not in globals():
            for j in np.arange(len(LFile)):
                name = 'File '+str(j)
                y_val = float(0.50 - 0.1*np.floor(int(j)/5))
                x_val = float(0.1 + int(j)%5*0.9/(5))

                textinputGen = TextInput(text=name, multiline=False, size_hint = (None, None),write_tab =False,
                    pos_hint = {'center_x': x_val, 'center_y': y_val}, height = '32dp',width='125dp')
                
                if int(j) == 0:    
                    alignGen = [textinputGen]
                else:
                    alignGen = alignGen + [textinputGen]

            # Before adding in all of the new ones
            for entry in alignGen:
                    FloatLayout.add_widget(self, entry)
        else:
            for entry in alignGen:
                    FloatLayout.remove_widget(self, entry)
            # HERE, why are we going back to "File" when clearly this page
            # has already been run? change this to read in "labels" or something
            for j in np.arange(len(LFile)):
                if 'labels' in globals():
                    if len(labels) > j:
                        name = labels[j]
                    else:
                        name = 'File '+str(j)
                else:
                    name = 'File '+str(j)
                y_val = float(0.50 - 0.1*np.floor(int(j)/5))
                x_val = float(0.1 + int(j)%5*0.9/(5))

                textinputGen = TextInput(text=name, multiline=False, size_hint = (None, None),write_tab =False,
                    pos_hint = {'center_x': x_val, 'center_y': y_val }, height = '32dp',width='125dp')

                if int(j) == 0:    
                    alignGen = [textinputGen]
                else:
                    alignGen = alignGen + [textinputGen]

            # Before adding in all of the new ones
            for entry in alignGen:
                    FloatLayout.add_widget(self, entry)
            # Make sure these boxes are deleted and recreated...
        global gen_check
        gen_check = True

    def check_aligns(self):
        global labels
        global gen_check
        x = len(alignGen)
        for row in np.arange(x):
            if row == 0:
                labels = alignGen[row].text
            else:
                labels = np.vstack((labels,alignGen[row].text))
        # Get labels out of the GD weird format.
        if onlyONE:
            pass
        else:
            labels = [i[0] for i in labels]
            # All of this here is to fix issues in naming
            # this fix will also break if users for some reason end filenames
            # with "___". Better hope they don't do that
            dum_len = len(labels)
            for i in np.arange(dum_len):
                for j in np.arange(dum_len):
                    if i == j:
                        continue
                    if labels[i].find(labels[j]) == 0:
                        labels[j] = labels[j]+'___'
        # Figure out what number of loops we are doing
        global cutoff
        global bulge_pad
        if molecule == 'pep':
            cutoff = int(self.lenCut.text)
            bulge_pad = int(self.bulger.text)
        global exp_drop
        exp_drop = self.degen_drop.active

################## NOW WE MOVE INTO THE LDA/BINARY COMPARISON SECTION ###############
class lda_binary(Screen):
    # Need to redefine a different type of loading from the wild one above.

    def on_pre_enter(self, *args):
        global lda_status
        if 'check_lda' not in globals():
            lda_L1 = Label(text='Group  ID',size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.15, 'center_y': 0.6},font_name='app_data/Poppins-Medium.ttf')
            #lda_L2 = Label(text='Binary Class 2',size_hint = (None, None), height = '48dp',
            #    pos_hint = {'center_x': 0.15, 'center_y': 0.45},font_name='app_data/Poppins-Medium.ttf')
            FloatLayout.add_widget(self, lda_L1)
            #FloatLayout.add_widget(self, lda_L2)

            if len(labels) > 2:
                lda_L3 = Label(text='Exclude from',size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.15, 'center_y': 0.45},font_name='app_data/Poppins-Medium.ttf')
                lda_L4 = Label(text='Analysis',size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.15, 'center_y': 0.40},font_name='app_data/Poppins-Medium.ttf')
                FloatLayout.add_widget(self, lda_L3)
                FloatLayout.add_widget(self, lda_L4)
            
            for j in np.arange(len(labels)):
                names = Label(text=labels[j], size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.3+int(j)*0.6/len(labels), 'center_y': 0.70},font_name='app_data/Poppins-Medium.ttf')
                if j == 0:
                    box1 = TextInput(text=str(j), multiline=False, size_hint = (None, None),write_tab =False,
                    pos_hint = {'center_x': 0.3+int(j)*0.6/len(labels), 'center_y': 0.6}, height = '32dp',width='32dp')

                else:
                    box1 = TextInput(text=str(j), multiline=False, size_hint = (None, None),write_tab =False,
                    pos_hint = {'center_x': 0.30+int(j)*0.6/len(labels), 'center_y': 0.6}, height = '32dp',width='32dp')

                if len(labels) > 2:
                    box3 = CheckBox(size_hint = (None,None), height = '48dp',
                    pos_hint = {'center_x': 0.3+int(j)*0.6/len(labels), 'center_y': 0.45},allow_no_selection = True)
                if int(j) == 0:
                    if len(labels) > 2:
                        lda_status = [[names,box1,box3]]
                    else:
                        lda_status = [[names,box1]]
                else:
                    if len(labels) > 2:
                        lda_status = lda_status + [[names,box1,box3]]
                    else:
                        lda_status = lda_status + [[names,box1]]
                # Before adding in all of the new ones
            for entry in lda_status:
                for j in entry:
                    FloatLayout.add_widget(self, j)
            # Make sure these boxes are deleted and recreated...
            global check_lda
            check_lda = True

    def check_checks(self):
        global lda_checked
        global group_a_id
        global sel1; global sel2
        # Convert our weird little kivy objects into a numpy array
        x,y = np.shape(lda_status)
        for row in np.arange(x):
            for column in np.arange(y):
                # Exclude the name column
                if column == 0:
                    continue
                # Read in the text column
                if column == 1:
                    text_pre_l = lda_status[row][column].text
                else:
                    check_pre_l = lda_status[row][column].active
            if row == 0:
                sel1 = int(text_pre_l)
                group_a_id = text_pre_l
                if column > 1:
                    lda_checked = check_pre_l
            elif row == 1:
                sel2 = int(text_pre_l)
                group_a_id = np.vstack((group_a_id,text_pre_l))
                if column > 1:
                    lda_checked = np.vstack((lda_checked,check_pre_l))
            else:
                group_a_id = np.vstack((group_a_id,text_pre_l))
                if column > 1:
                    lda_checked = np.vstack((lda_checked,check_pre_l))

class intro(Screen):
    global skip7
    skip7 = False
    def junk(self):
        x=1
#################### THEN THE APP, CLOSE OUT ANY ANALYSIS #######################
class AIMSApp(App):
    # ALRIGHT FINALLY MAYBE GETTING THIS
    #v_label2 = StringProperty('')
    #v_label3 = StringProperty('')
    #v_label4 = StringProperty('')

    index = NumericProperty(-1)
    current_title = StringProperty()
    time = NumericProperty(0)
    screen_names = ListProperty([])
    hierarchy = ListProperty([])

    def build(self):
        self.title = 'AIMS - Automated Immune Molecule Separator'
        # Alright since we 
        Clock.schedule_interval(self._update_clock, 1 / 60.)
        self.screens = {}

        # Honestly we can probably just hardcode this... That way we can control order for now.
        self.available_screens = ['app_data/screens/main_screen.kv','app_data/screens/mhc_1.kv',
        'app_data/screens/mhc_2.kv','app_data/screens/mhc_3.kv','app_data/screens/mhc_4.kv',
        'app_data/screens/mhc_5.kv','app_data/screens/mhc_6.kv','app_data/screens/mhc_7.kv',
        'app_data/screens/mhc_8.kv','app_data/screens/mhc_9.kv','app_data/screens/mhc_10.kv',
        'app_data/screens/mhc_11.kv','app_data/screens/mhc_12.kv','app_data/screens/mhc_13.kv',
        'app_data/screens/mhc_14.kv',
        'app_data/screens/ab_1.kv','app_data/screens/ab_2.kv','app_data/screens/ab_3.kv',
        'app_data/screens/ab_4.kv','app_data/screens/ab_5.kv','app_data/screens/ab_6.kv',
        'app_data/screens/ab_7.kv','app_data/screens/ab_8.kv','app_data/screens/ab_9.kv',
        'app_data/screens/ab_10.kv','app_data/screens/ab_11.kv','app_data/screens/ab_12.kv',
        'app_data/screens/ab_13.kv','app_data/screens/ab_14.kv',
        'app_data/screens/msa_1.kv','app_data/screens/msa_2.kv','app_data/screens/msa_3.kv',
        'app_data/screens/msa_4.kv','app_data/screens/msa_5.kv','app_data/screens/msa_6.kv',
        'app_data/screens/msa_7.kv','app_data/screens/msa_8.kv','app_data/screens/msa_9.kv',
        'app_data/screens/msa_10.kv','app_data/screens/msa_11.kv','app_data/screens/msa_12.kv',
        'app_data/screens/msa_13.kv','app_data/screens/msa_14.kv',
        'app_data/screens/pep_1.kv','app_data/screens/pep_2.kv','app_data/screens/pep_3.kv',
        'app_data/screens/pep_4.kv','app_data/screens/pep_5.kv','app_data/screens/pep_6.kv',
        'app_data/screens/pep_7.kv','app_data/screens/pep_8.kv','app_data/screens/pep_9.kv',
        'app_data/screens/pep_10.kv','app_data/screens/pep_11.kv','app_data/screens/pep_12.kv',
        'app_data/screens/pep_13.kv','app_data/screens/pep_14.kv']
        global molecule
        molecule = ''
        self.go_next_screen(num = 0)
    
    def go_main_screen(self,num = 0):
        self.index = num
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='right')
        self.current_title = screen.name

    def go_next_screen(self,num = 0):
        global molecule
        if skip7 and num == 7:
            num = 8        
        if molecule == 'ig':
            # somewhat janky way to do this, but should work.
            # Basically "skip" the mhc screens
            num = num + 14
        if molecule == 'msa':
            num = num + 28
        if molecule == 'pep':
            num = num + 42
        self.index = num
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='left')
        self.current_title = screen.name

    def go_mhc(self):
        # Need to make and delete lots of shit when
        # we get back to this main screen. HARD restart
        global N; N = 4
        global LFile; LFile = ['']
        global LFile_cut; LFile_cut = ['']
        # Delete some global variables so we have a HARD reset.
        global check_run
        if 'check_run' in globals():
            # So intelligently, if you define something as a global variable
            # but don't "fill" it, it won't register as being in globals.
            del check_run
        global check_lda
        if 'check_lda' in globals():
            del check_lda
        global loadLabel
        if 'loadLabel' in globals():
            del loadLabel
        
        # Needed to delete and define everything BEFORE loading screen
        global molecule
        molecule = 'mhc'
        self.index = 1
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='left')
        self.current_title = screen.name

    def go_msa(self):
        # Need to make and delete lots of shit when
        # we get back to this main screen. HARD restart
        global N; N = 4
        global LFile; LFile = ['']
        global LFile_cut; LFile_cut = ['']
        # Delete some global variables so we have a HARD reset.
        global check_run
        if 'check_run' in globals():
            # So intelligently, if you define something as a global variable
            # but don't "fill" it, it won't register as being in globals.
            del check_run
        global check_lda
        if 'check_lda' in globals():
            del check_lda
        global loadLabel
        if 'loadLabel' in globals():
            del loadLabel
        
        # Needed to delete and define everything BEFORE loading screen
        global molecule
        molecule = 'msa'
        self.index = 29
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='left')
        self.current_title = screen.name

    def go_Ig(self):
        # Need to make and delete lots of shit when
        # we get back to this main screen. HARD restart
        global N; N = 4
        global LFile; LFile = ['']
        global LFile_cut; LFile_cut = ['']
        # Delete some global variables so we have a HARD reset.
        global check_run
        if 'check_run' in globals():
            # So intelligently, if you define something as a global variable
            # but don't "fill" it, it won't register as being in globals.
            del check_run
        global check_lda
        if 'check_lda' in globals():
            del check_lda
        global loadLabel
        if 'loadLabel' in globals():
            del loadLabel

         # Needed to delete and define everything BEFORE loading screen
        global molecule
        molecule = 'ig'
        self.index = 15
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='left')
        self.current_title = screen.name

    def go_pep(self):
        # Need to make and delete lots of shit when
        # we get back to this main screen. HARD restart
        global N; N = 4
        global LFile; LFile = ['']
        global LFile_cut; LFile_cut = ['']
        # Delete some global variables so we have a HARD reset.
        global check_run
        if 'check_run' in globals():
            # So intelligently, if you define something as a global variable
            # but don't "fill" it, it won't register as being in globals.
            del check_run
        global check_lda
        if 'check_lda' in globals():
            del check_lda
        global loadLabel
        if 'loadLabel' in globals():
            del loadLabel

         # Needed to delete and define everything BEFORE loading screen
        global molecule
        molecule = 'pep'
        self.index = 43
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='left')
        self.current_title = screen.name

    def go_prev_screen(self,num = 0):
        global molecule
        if skip7 and num == 7:
            num = 6
        if molecule == 'ig':
            # somewhat janky way to do this, but should work.
            # Basically "skip" the mhc screens
            num = num + 14
        if molecule == 'msa':
            num = num + 28
        if molecule == 'pep':
            num = num + 42
        self.index = num
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='right')
        self.current_title = screen.name

    def load_screen(self, index):
        if index in self.screens:
            return self.screens[index]
        screen = Builder.load_file(self.available_screens[index])
        self.screens[index] = screen
        return screen

    def _update_clock(self, dt):
          self.time = time()

Factory.register('Root', cls=Root)
Factory.register('intro', cls=intro)
Factory.register('aligner', cls=aligner)
Factory.register('aligner_ab', cls=aligner_ab)
Factory.register('aligner_gen', cls=aligner_gen)
Factory.register('Analysis', cls=Analysis)
Factory.register('checker', cls=checker)

main = launch_new_instance = AIMSApp().run

if __name__ == '__main__':
    main()
