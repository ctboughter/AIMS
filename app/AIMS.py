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
from time import time

import aims_loader as aimsLoad
import aims_analysis as aims
import aims_classification as classy

from kivy.uix.button import Button
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput

from kivy.uix.screenmanager import Screen
from os.path import dirname, join
from kivy.lang import Builder
from kivy.properties import NumericProperty, StringProperty, BooleanProperty,\
    ListProperty

import os
import re

# Evidently these lines will remove a weird multi-touch emulator from kivy:
from kivy.config import Config
Config.set('input', 'mouse', 'mouse,multitouch_on_demand')
# Jesus, this thing really was built for phone apps...

class Root(Screen):
    # REALLY Sucks to use global variables, but I'm unsure how else to
    # pass all of these across functions
    global N; N = 4
    global LFile; LFile = ['']
    fullscreen = BooleanProperty(False)

    def on_pre_enter(self):
        global loadLabel
        global loadButton
        if 'loadLabel' not in globals():
            global LFile; LFile = ['']
            N = 4
            a=0
            # Need to re-read through this shit to see what's going on...
            for j in np.arange(int(N)):
                
                # No idea what this lambda function does, but it lets us bind show load to the button
                if molecule == 'ig':
                    xxname = 'File '
                    xname = ''
                else:
                    xxname = 'FASTA '
                    xname = 'FASTA '
                button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
                pos_hint={'center_x':.15, 'center_y':.75-int(a)*0.6/N},
                on_release=lambda x = int(j):self.show_load(win = x))
                # What an absolute fucking nightmare solution the line above... works though
                # Need to make sure we don't overwrite the labels every time we load shit in
                if j >= len(LFile):
                    label = Label(text=xname +'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})
                else:
                    if LFile[int(j)] != '':
                        label = Label(text=LFile[int(j)], size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})
                    else:
                        label = Label(text=xname+'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                        pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})

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
        return(os.getcwd())

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

    # win is how we will try to keep track of using the right button...
    def show_load(self, win = 2):
        content = LoadDialog(load=self.load, cancel=self.dismiss_popup, fas1 = self.do_thing(win2 = win))
        self._popup = Popup(title="Load file", content=content,
                            size_hint=(0.9, 0.9))
        self._popup.open()

    def load(self, path, filename):
        global LFile
        path1 = os.path.join(path, filename[0])
        # So FASTA_L should tell you WHERE
        # The loadfile is coming from
        while FASTA_L+1 > len(LFile):
            LFile = LFile + ['']
        LFile[FASTA_L] = path1
        # Need to have two separate options because we move from Kivy defined buttons to
        # python defined buttons. I handle those slightly differently.
        loadLabel[FASTA_L].text = path1

        if len(LFile) >= 2:
            self.next1_1.disabled = False 

        self.dismiss_popup()

    def make_path(self):
        global dir_name
        self.text1 = self.v_input1.text
        dir_name = self.text1

    # Basically just recreate the entire screen every time we want to 
    # add more fastas. Kind of a pain, but it works.

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
            if molecule == 'ig':
                xxname = 'File '
                xname = ''
            else:
                xxname = 'FASTA '
                xname = 'FASTA '
            button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
            pos_hint={'center_x':.15, 'center_y':.75-int(a)*0.6/N},
            on_release=lambda x = int(j):self.show_load(win = x))
            # What an absolute fucking nightmare solution the line above... works though
            # Need to make sure we don't overwrite the labels every time we load shit in
            if j >= len(LFile):
                label = Label(text=xname +'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})
            else:
                if LFile[int(j)] != '':
                    label = Label(text=LFile[int(j)], size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})
                else:
                    label = Label(text=xname+'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})

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
            if molecule == 'ig':
                xxname = 'File '
                xname = ''
            else:
                xxname = 'FASTA '
                xname = 'FASTA '
            button = Button(text='Load '+ xxname + str(j+1), size_hint=(0.2, 0.075),
            pos_hint={'center_x':.15, 'center_y':.75-int(a)*0.6/N},
            on_release=lambda x = int(j):self.show_load(win = x))
            # What an absolute fucking nightmare solution the line above... works though
            # Need to make sure we don't overwrite the labels every time we load shit in
            if j >= len(LFile):
                label = Label(text=xname +'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})
            else:
                if LFile[int(j)] != '':
                    label = Label(text=LFile[int(j)], size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})
                else:
                    label = Label(text=xname+'File ' + str(j+1) + ' Path', size_hint=(0.2, 0.075),
                    pos_hint={'center_x':.65, 'center_y':.75-int(a)*0.6/N})

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
            size =('150dp','48dp'),pos_hint={'center_x':.4, 'center_y':0.1}, on_release=self.show_load)
            button2_moved = Button(text='Resize Entries', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.6, 'center_y':0.1}, on_release=self.load_aligners)
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
            size =('150dp','48dp'),pos_hint={'center_x':.4, 'center_y':0.1}, on_release=self.show_load)
            button2_moved = Button(text='Resize Entries', size_hint=(None, None),
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
        # Literally just here so we can have these messages or not
        if len(LFile) < 2:
            if align_label.text == 'ERROR: MUST LOAD AT LEAST TWO FILES':
                FloatLayout.remove_widget(self,align_label) 
            align_label = Label(text='ERROR: MUST LOAD AT LEAST TWO FILES', size_hint=(0.4, 0.2),
                    pos_hint={'center_x':.5, 'center_y':.5}, color = (1, 0, 0, 1))
            FloatLayout.add_widget(self, align_label)
            button_moved = Button(text='Load Matrix', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.4, 'center_y':0.1}, on_release=self.show_load)
            button2_moved = Button(text='Resize Entries', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.6, 'center_y':0.1}, on_release=self.load_aligners)
            FloatLayout.add_widget(self, button2_moved)
            FloatLayout.add_widget(self, button_moved)
        else:
            if align_label.text == 'ERROR: MUST LOAD AT LEAST TWO FILES':
                FloatLayout.remove_widget(self,align_label)
                align_label.text = ''
            button_moved = Button(text='Load Matrix', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.4, 'center_y':0.1}, on_release=self.show_load)
            button2_moved = Button(text='Resize Entries', size_hint=(None, None),
            size =('150dp','48dp'),pos_hint={'center_x':.6, 'center_y':0.1}, on_release=self.load_aligners)
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

        labels = mat_coords[:,0]
        # When you go back this far, make sure we recreate the checkbox screen
        if 'check_run' in globals():
            # So intelligently, if you define something as a global variable
            # but don't "fill" it, it won't register as being in globals.
            del check_run
    
class LoadDialog(FloatLayout):
    load = ObjectProperty(None)
    cancel = ObjectProperty(None)
    fas1 = ObjectProperty(None)

    def get_path(self):
        return(os.getcwd())

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
        ##########################################
        global seq_MIf
        global seqNameF
        global seq_final
        global labels
        this_dir = os.getcwd()
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

        if len(paths) < 2:
            self.img1.source = this_dir + '/app_data/error_load.png'
            return()
        AA_num_key = aims.get_props()[1]
        # Labels will always be predefined for Abs
        # Predefined, but poorly formatted... reformat here
        global labels
        global LOOPnum
        if molecule == 'mhc':
            data = mat_coords[:,1:]

        if len(labels) != len(paths):
            # need to make a new error message for this.
            self.img1.source = this_dir + '/app_data/error_load.png'
            return()

        for i in np.arange(len(paths)):
            if molecule == 'mhc':
                #if any(data[:,i] == ['','','','','']):
                #    self.img1.source = this_dir + '/app_data/error_pos.png'
                #    return()
                # turn data into an integer.
                int_dat = [int(x) for x in data[i]]
                seq,seq_key = aimsLoad.mhc_loader(paths[i],int_dat,labels[i])
            else:
                seq = aimsLoad.Ig_loader(paths[i],labels[i],loops=LOOPnum)
            if i == 0:
                seq_final = seq
                seq_size = np.shape(seq)[1]
                seqNameF = labels[i]
                if molecule == 'mhc':
                    seq_keyF = seq_key
                mat_size = aims.get_sequence_dimension(np.array(seq))[0]
            else:
                seq_final = pandas.concat([seq_final,seq],axis = 1)
                seqNameF = np.vstack((seqNameF,labels[i]))
                seq_size = np.vstack((seq_size,np.shape(seq)[1]))
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
        seq_MI = aims.gen_tcr_matrix(np.array(seq_final),key = AA_num_key, giveSize = mat_size)
        # Convert our MI matrix to a pandas dataframe
        seq_MIf = pandas.DataFrame(np.transpose(seq_MI),columns = seq_final.columns)
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
        ax[0,0].imshow(seq_MI, interpolation='nearest', aspect='auto',cmap=cm.jet)

        # Need to try to draw lines separating the distinct groups...
        seq_locs = 0
        for i in np.arange(len(seq_size)-1):
            seq_locs = seq_locs + seq_size[i]
            ax[0,0].plot(np.arange(len(np.transpose(seq_MI))),np.ones(len(np.transpose(seq_MI)))*seq_locs,'white',linewidth = 3)

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
            ax[0,0].set_xticks(xtick_loc)
        else:
            xtick_loc = [mat_size/2]
            ax[0,0].set_xticks(xtick_loc)
        if molecule == 'mhc':
            ax[0,0].set_xticklabels(['Strand 1','Helix 1','Strand 2','Helix 2'])
        else:
            if LOOPnum == 6:
                ax[0,0].set_xticklabels(['CDR1L','CDR2L','CDR3L','CDR1H','CDR2H','CDR3H'])
            elif LOOPnum == 3:
                ax[0,0].set_xticklabels(['CDR1','CDR2','CDR3'])
            elif LOOPnum == 1:
                ax[0,0].set_xticklabels(['CDR Loop'])

        # NOTE need to make sure that I do something so that files are not overwritten...
        if os.path.exists(this_dir +'/' + dir_name):
            None
        else:
            os.mkdir(this_dir +'/' + dir_name)
        fig.savefig(this_dir + '/' + dir_name + '/matrix.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/matrix.png',format='png',dpi=500)
        # save raw data at all steps now...
        np.savetxt(this_dir + '/' + dir_name + '/raw_matrix.dat',seq_MIf,fmt='%.3f')
        if molecule == 'mhc':
            np.savetxt(this_dir + '/' + dir_name + '/sequence_key.txt',seq_keyF,fmt='%s')

        self.img1.source = this_dir + '/' + dir_name + '/matrix.png'

    def get_props(self):
        global pcaF
        pca_props = aims.gen_clone_props(np.array(np.transpose(seq_MIf)))
        pcaF = pandas.DataFrame(pca_props,columns = seq_MIf.columns)
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
        x_axis = x_axis = np.array([-0.2,0.9,2,3.1])
        # Need to have some kind of color wheel to replace this...
        #colors = ['purple','green','black','orange']
        for i in np.arange(len(seqNameF)):
            index = [column for column in pcaF.columns if seqNameF[i][0] in column]
            plotThis = np.array(pcaF[index])
            # Specifically plot the first 5 properties here... Make sure that I 
            # have some way soon to actually select these...
            ax[0,0].bar(x_axis+i*1/len(seqNameF), np.average(plotThis[1:5,:],axis = 1),
                        yerr = np.std(plotThis[1:5,:],axis = 1),alpha = 0.5, width = 1/len(seqNameF))
        ax[0,0].legend(labels)
        ax[0,0].set_xticks([0.2,1.3,2.4,3.5])
        ax[0,0].set_xticklabels(['Charge','Hydrophobicity','Flexibility','Bulkiness'])
        ax[0,0].set_xlabel('Biophysical Property')
        ax[0,0].set_ylabel('Normalized Property Value')
        this_dir = os.getcwd()
        fig.savefig(this_dir + '/' + dir_name + '/physical_props.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/physical_props.png',format='png',dpi=500)
        np.savetxt(this_dir + '/' + dir_name + '/property_masks.dat',pcaF,fmt='%.3f')
        #NOTE DELETE THIS BELOW LINE ASAP
        self.img4.source = this_dir + '/' + dir_name + '/physical_props.png'
    
    def do_pca(self):
        from sklearn.decomposition import PCA
        global checked
        this_dir = os.getcwd()
        # The "things" are here to track number of entries in PCA groups
        thing1 = True
        for i in np.arange(len(seqNameF)):
            index = [column for column in pcaF.columns if seqNameF[i][0] in column]
            # Alright so our new matrix, checked is "i" rows by three
            if checked[i,0]:
                if thing1:
                    pg1 = np.array(pcaF[index])
                    thing1 = False
                else:
                    pg1 = np.hstack((pg1,pcaF[index]))

        pca_Doit = np.transpose(pg1)
        pca = PCA(n_components=3, svd_solver='full')
        pca.fit(pca_Doit)

        explain_var = pca.explained_variance_ratio_
        comps1=pca.components_[0]
        comps2=pca.components_[1]
        comps3=pca.components_[2]
        prop_names = pandas.read_csv('app_data/property_names.txt',header=None).values
        
        all_comps = np.vstack((prop_names.reshape(len(prop_names)),comps1,comps2,comps3))
        # Alrightn now we need to sort on the pca components
        indices_topProps = np.argsort(abs(comps1))
        compile_Props = np.vstack((prop_names.reshape(len(prop_names)),comps1))
        final_pc1 = compile_Props[:,indices_topProps[-10:]]

        self.exp1.text = str(round(explain_var[0],3))
        self.exp2.text = str(round(explain_var[1],3))
        self.exp3.text = str(round(explain_var[2],3))

        # We square the weights so that they then sum to 1.
        self.name10.text = str(final_pc1[0,0]); self.weight10.text = str(round(final_pc1[1,0]**2,4))
        self.name9.text = str(final_pc1[0,1]); self.weight9.text = str(round(final_pc1[1,1]**2,4))
        self.name8.text = str(final_pc1[0,2]); self.weight8.text = str(round(final_pc1[1,2]**2,4))
        self.name7.text = str(final_pc1[0,3]); self.weight7.text = str(round(final_pc1[1,3]**2,4))
        self.name6.text = str(final_pc1[0,4]); self.weight6.text = str(round(final_pc1[1,4]**2,4))
        self.name5.text = str(final_pc1[0,5]); self.weight5.text = str(round(final_pc1[1,5]**2,4))
        self.name4.text = str(final_pc1[0,6]); self.weight4.text = str(round(final_pc1[1,6]**2,4))
        self.name3.text = str(final_pc1[0,7]); self.weight3.text = str(round(final_pc1[1,7]**2,4))
        self.name2.text = str(final_pc1[0,8]); self.weight2.text = str(round(final_pc1[1,8]**2,4))
        self.name1.text = str(final_pc1[0,9]); self.weight1.text = str(round(final_pc1[1,9]**2,4))
        # Alright so let me try to have labels encoded in to the MHC_6 kivy file
        # and then just populate them with things every time we run this function

        # THIS SECTION IS WHERE YOU WOULD PUT PLOTS THAT YOU WANNA SEE
        from mpl_toolkits import mplot3d
        fig3d = pl.figure(figsize = (10, 10))
        ax3d = fig3d.add_subplot(111, projection='3d')
        legName = []
        for i in np.arange(len(seqNameF)):
            # What a dumb way to do this but I'm unsure how else to reverse this...
            if len(seqNameF) > 2:
                if checked[i,2] == False:
                    index = [column for column in pcaF.columns if seqNameF[i][0] in column]
                    group_transform = pca.transform(np.transpose(pcaF[index]))
                    legName = legName + [labels[i]]
                    ax3d.scatter(group_transform[:,0],group_transform[:,1],group_transform[:,2],marker='o',alpha=0.5,s=50)
            else:
                index = [column for column in pcaF.columns if seqNameF[i][0] in column]
                group_transform = pca.transform(np.transpose(pcaF[index]))
                legName = legName + [labels[i]]
                ax3d.scatter(group_transform[:,0],group_transform[:,1],group_transform[:,2],marker='o',alpha=0.5,s=50)

        #pl.legend(['seq1','seq2','seq3'])
        ax3d.set_xlabel('PC1',labelpad=20)
        ax3d.set_ylabel('PC2',labelpad=20)
        ax3d.set_zlabel('PC3',labelpad=10)
        ax3d.legend(legName)
        fig3d.savefig(this_dir + '/' + dir_name + '/pca_3D.pdf',format='pdf',dpi=500)
        fig3d.savefig(this_dir + '/' + dir_name + '/pca_3D.png',format='png',dpi=500)

        fig, ax2 = pl.subplots(1, 1,squeeze=False,figsize=(10,10))
        for i in np.arange(len(seqNameF)):
            if len(seqNameF) > 2:
                if checked[i,2] == False:
                    index = [column for column in pcaF.columns if seqNameF[i][0] in column]
                    group_transform = pca.transform(np.transpose(pcaF[index]))
                    pl.scatter(group_transform[:,0],group_transform[:,1])
            else:
                index = [column for column in pcaF.columns if seqNameF[i][0] in column]
                group_transform = pca.transform(np.transpose(pcaF[index]))
                pl.scatter(group_transform[:,0],group_transform[:,1])
        ax2[0,0].set_xlabel('PC1')
        ax2[0,0].set_ylabel('PC2')
        fig.legend(legName)
        fig.savefig(this_dir + '/' + dir_name + '/pca_2D.pdf',format='pdf',dpi=500)
        # And save the raw data
        np.savetxt(this_dir + '/' + dir_name + '/pca_transformed_data.dat',group_transform,fmt='%.3f')
        np.savetxt(this_dir + '/' + dir_name + '/pca_component_weights.dat',all_comps,fmt='%s')
        np.savetxt(this_dir + '/' + dir_name + '/pca_explained_variance.dat',explain_var,fmt='%.3f')

        self.img6.source = this_dir + '/' + dir_name + '/pca_3D.png'

    def get_pos_props(self):
        this_dir = os.getcwd()
        global lda_checked
        # The "things" are here to track number of entries in PCA groups
        thing1 = True
        thing2 = True
        for i in np.arange(len(seqNameF)):
            index = [column for column in pcaF.columns if seqNameF[i][0] in column]
            # Alright so our new matrix, checked is "i" rows by three
            # 0 is the "group 1", 1 is the "group 2"
            if lda_checked[i,0]:
                if thing1:
                    pg1 = np.array(seq_MIf[index])
                    thing1 = False
                else:
                    pg1 = np.hstack((pg1,seq_MIf[index]))
            elif lda_checked[i,1]:
                if thing2:
                    pg2 = np.array(seq_MIf[index])
                    thing2 = False
                else:
                    pg2 = np.hstack((pg2,seq_MIf[index]))

        pos_sens1=aims.gen_dset_props(np.array(np.transpose(pg1)),stdev=False)
        pos_sens2=aims.gen_dset_props(np.array(np.transpose(pg2)),stdev=False)

        fig, ax = pl.subplots(2, 1,squeeze=False,figsize=(14,10))
        global LOOPnum
        global xtick_loc

        for prop in np.arange(2):
            if prop == 0:
                x = 0; y = 0
            elif prop == 1:
                x = 1; y = 0
            plotThis1 = pos_sens1[prop+2]
            plotThis2 = pos_sens2[prop+2]
            # Specifically plot the first 5 properties here... Make sure that I 
            # have some way soon to actually select these...
            ax[x,y].plot(np.arange(len(plotThis1)), plotThis1)
            ax[x,y].plot(np.arange(len(plotThis2)), plotThis2)
            ax[x,y].legend(['Group1', 'Group2'])

            if molecule == 'mhc':
                ax[x,y].set_xticks(xtick_loc)
                ax[0,0].set_xticklabels(['Strand 1','Helix 1','Strand 2','Helix 2'])
            else:
                if LOOPnum == 6:
                    ax[x,y].set_xticks(xtick_loc)
                    ax[x,y].set_xticklabels(['CDR1L','CDR2L','CDR3L','CDR1H','CDR2H','CDR3H'])
                if LOOPnum == 3:
                    ax[x,y].set_xticks(xtick_loc)
                    ax[x,y].set_xticklabels(['CDR1','CDR2','CDR3'])
                if LOOPnum == 1:
                    ax[x,y].set_xticks(xtick_loc)
                    ax[x,y].set_xticklabels(['CDR Loop'])
        # Since there's only two now, just easier to hard code these...
        ax[0,0].set_ylabel('Normalized Charge')
        ax[1,0].set_ylabel('Normalized Hydrophobicity')

        fig.savefig(this_dir + '/' + dir_name + '/pos_prop.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/pos_prop.png',format='png',dpi=500)
        # And save the raw data
        np.savetxt(this_dir + '/' + dir_name + '/position_sensitive_mat1.dat',pos_sens1,fmt='%.3f')
        np.savetxt(this_dir + '/' + dir_name + '/position_sensitive_mat2.dat',pos_sens2,fmt='%.3f')

        self.img8.source = this_dir + '/' + dir_name + '/pos_prop.png'

    def do_lda(self):
        this_dir = os.getcwd()
        # The "things" are here to track number of entries in PCA groups
        thing1 = True
        thing2 = True
        for i in np.arange(len(seqNameF)):
            index = [column for column in pcaF.columns if seqNameF[i][0] in column]
            # Alright so our new matrix, checked is "i" rows by three
            # 0 is the "group 1", 1 is the "group 2"
            if lda_checked[i,0]:
                if thing1:
                    pg1 = np.array(seq_final[index])
                    thing1 = False
                else:
                    pg1 = np.hstack((pg1,seq_final[index]))
            elif lda_checked[i,1]:
                if thing2:
                    pg2 = np.array(seq_final[index])
                    thing2 = False
                else:
                    pg2 = np.hstack((pg2,seq_final[index]))

        # Should probably make NumVects an option you can change...
        numVects = int(self.inputLDA.text)
        # Here we also need to get a pre-defined matrix size...

        num1 = np.shape(pg1)[1]
        num2 = np.shape(pg2)[1]
        #x,y,MatrixSize = aims.gen_tcr_matrix(pg1,pre_mono=pg2,binary=True,return_Size=True)
        acc_all,weights,cols,indices,mda_all,dropped = classy.do_linear_split(pg1, pg2, matSize = numVects)
        # Seaborn plots look nicer for these LDA figures
        import seaborn as sns
        fig = pl.figure(figsize = (12, 12))
        dset = ["Linear Discriminant Analysis" for x in range(num1+num2)]
        reacts = ["Group 1" for x in range(num1)] + ["Group 2" for x in range(num2)]

        d1 = {'Dataset': dset, 'Linear Discriminant 1': mda_all.reshape(len(mda_all)),
            'Reactivity' : reacts}
        df1 = pandas.DataFrame(data=d1)
        sns.set(style="white", color_codes=True,font_scale=1.5)
        sns.swarmplot(x="Dataset", y="Linear Discriminant 1", data=df1, hue = 'Reactivity', palette = "Dark2")

        fig.savefig(this_dir + '/' + dir_name + '/lda.pdf',format='pdf',dpi=500)
        fig.savefig(this_dir + '/' + dir_name + '/lda.png',format='png',dpi=500)
        # And save the raw data
        np.savetxt(this_dir + '/' + dir_name + '/lda_data.dat',mda_all,fmt='%.3f')

        # Need to put acc_all somewhere on the screen
        indices_topProps = np.argsort(abs(weights))
        compile_Props = np.vstack((cols[indices],weights))
        final_pc1 = compile_Props[:,indices_topProps[0][-10:]]
        np.savetxt(this_dir + '/' + dir_name + '/lda_weights.dat',compile_Props,fmt='%s')

        # For running multiple times with vectors less than 10, need to reset screen...
        self.aname10.text = ''; self.aname9.text = ''; self.aname8.text = ''; self.aname7.text = ''
        self.aname6.text = ''; self.aname5.text = ''; self.aname4.text = ''; self.aname3.text = ''
        self.aname2.text = ''; self.aname1.text = '';
        self.aweight10.text = ''; self.aweight9.text = ''; self.aweight8.text = ''; self.aweight7.text = ''
        self.aweight6.text = ''; self.aweight5.text = ''; self.aweight4.text = ''; self.aweight3.text = ''
        self.aweight2.text = ''; self.aweight1.text = '';
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

        self.img9.source = this_dir + '/' + dir_name + '/lda.png'

class checker(Screen):
    # Need to redefine a different type of loading from the wild one above.

    def on_pre_enter(self, *args):
        global check_status
        if 'check_run' not in globals():
            check_L1 = Label(text='Include in PCA',size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.15, 'center_y': 0.6})
            FloatLayout.add_widget(self, check_L1)
            # only add in these you have MORE than 2 options
            if len(labels) > 2:
                check_L2 = Label(text='Exclude from PCA',size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.15, 'center_y': 0.45})
                FloatLayout.add_widget(self, check_L2)
                check_L3 = Label(text='Exclude from Plot',size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.15, 'center_y': 0.3})
                check_L4 = Label(text='(Optional)',size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.15, 'center_y': 0.25})
                FloatLayout.add_widget(self, check_L3)
                FloatLayout.add_widget(self, check_L4)
            
            for j in np.arange(len(labels)):
                names = Label(text=labels[j], size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.70})
                if len(labels) > 2:
                    box1 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j),state = 'down',
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.6})
                    box2 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j),
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.45})
                    box3 = CheckBox(size_hint = (None,None), height = '48dp',
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.3})
                else:
                    box1 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j),state = 'down',
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.6})
                if int(j) == 0:
                    if len(labels) > 2:
                        check_status = [[names,box1,box2,box3]]
                    else:
                        check_status = [[names,box1]]
                else:
                    if len(labels) > 2:
                        check_status = check_status + [[names,box1,box2,box3]]
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
        labels = [i[0] for i in labels]
        
        # Figure out what number of loops we are doing
        global LOOPnum
        if self.cdr1.active == True:
            LOOPnum = 1
        elif self.cdr3.active == True:
            LOOPnum = 3
        elif self.cdr6.active == True:
            LOOPnum = 6

################## NOW WE MOVE INTO THE LDA/BINARY COMPARISON SECTION ###############
class lda_binary(Screen):
    # Need to redefine a different type of loading from the wild one above.

    def on_pre_enter(self, *args):
        global lda_status
        if 'check_lda' not in globals():
            lda_L1 = Label(text='Binary Class 1',size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.15, 'center_y': 0.6})
            lda_L2 = Label(text='Binary Class 2',size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.15, 'center_y': 0.45})
            FloatLayout.add_widget(self, lda_L1)
            FloatLayout.add_widget(self, lda_L2)

            if len(labels) > 2:
                lda_L3 = Label(text='Exclude from',size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.15, 'center_y': 0.3})
                lda_L4 = Label(text='Binary Analysis',size_hint = (None, None), height = '48dp',
                    pos_hint = {'center_x': 0.15, 'center_y': 0.25})
                FloatLayout.add_widget(self, lda_L3)
                FloatLayout.add_widget(self, lda_L4)
            
            for j in np.arange(len(labels)):
                names = Label(text=labels[j], size_hint = (None, None), height = '48dp',
                pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.70})
                if j == 0:
                    box1 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j+50),state = 'down',
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.6})
                    box2 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j+50),
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.45})
                else:
                    box1 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j+50),
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.6})
                    box2 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j+50),state = 'down',
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.45})
                if len(labels) > 2:
                    box3 = CheckBox(size_hint = (None,None), height = '48dp', group = 'g'+str(j+50),
                    pos_hint = {'center_x': 0.3+int(j)*0.6/N, 'center_y': 0.3})
                if int(j) == 0:
                    if len(labels) > 2:
                        lda_status = [[names,box1,box2,box3]]
                    else:
                        lda_status = [[names,box1,box2]]
                else:
                    if len(labels) > 2:
                        lda_status = lda_status + [[names,box1,box2,box3]]
                    else:
                        lda_status = lda_status + [[names,box1,box2]]
                # Before adding in all of the new ones
            for entry in lda_status:
                for j in entry:
                    FloatLayout.add_widget(self, j)
            # Make sure these boxes are deleted and recreated...
            global check_lda
            check_lda = True

    def check_checks(self):
        global lda_checked
        # Convert our weird little kivy objects into a numpy array
        x,y = np.shape(lda_status)
        for row in np.arange(x):
            for column in np.arange(y):
                # Exclude the name column
                if column == 0:
                    continue
                if column == 1:
                    check_pre_l = lda_status[row][column].active
                else:
                    check_pre_l = np.hstack((check_pre_l,lda_status[row][column].active))
            if row == 0:
                lda_checked = check_pre_l
            else:
                lda_checked = np.vstack((lda_checked,check_pre_l))

class intro(Screen):
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
        'app_data/screens/mhc_8.kv','app_data/screens/mhc_9.kv',
        'app_data/screens/ab_1.kv','app_data/screens/ab_2.kv','app_data/screens/ab_3.kv',
        'app_data/screens/ab_4.kv','app_data/screens/ab_5.kv','app_data/screens/ab_6.kv',
        'app_data/screens/ab_7.kv','app_data/screens/ab_8.kv','app_data/screens/ab_9.kv']
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
        if molecule == 'ig':
            # somewhat janky way to do this, but should work.
            # Basically "skip" the mhc screens
            num = num + 9
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

    def go_Ig(self):
        # Need to make and delete lots of shit when
        # we get back to this main screen. HARD restart
        global N; N = 4
        global LFile; LFile = ['']
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
        self.index = 10
        screen = self.load_screen(self.index)
        sm = self.root.ids.sm
        sm.switch_to(screen, direction='left')
        self.current_title = screen.name

    def go_prev_screen(self,num = 0):
        global molecule
        if molecule == 'ig':
            # somewhat janky way to do this, but should work.
            # Basically "skip" the mhc screens
            num = num + 9
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
Factory.register('Analysis', cls=Analysis)
Factory.register('checker', cls=checker)

if __name__ == '__main__':
    AIMSApp().run()
