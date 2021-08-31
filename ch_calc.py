#!/usr/bin/env python3

""" Copyright © 2020 Borys Olifirov
Channels crosstalk calculator

"""

import sys
import os
import random
import logging
import yaml

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


def wavelen2rgb(Wavelength, MaxIntensity=100):
    """Calculate RGB values given the wavelength of visible light.

    Arguments:
    * Wavelength:  Wavelength in nm.  Scalar floating.
    * MaxIntensity:  The RGB value for maximum intensity.  Scalar 
      integer.

    Returns:
    * 3-element list of RGB values for the input wavelength.  The
      values are scaled from 0 to MaxIntensity, where 0 is the
      lowest intensity and MaxIntensity is the highest.  Integer
      list.

    Visible light is in the range of 380-780 nm.  Outside of this
    range the returned RGB triple is [0,0,0].

    Based on code by Earl F. Glynn II at:
       http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm
    See also:
       http://www.physics.sfasu.edu/astro/color/spectra.html
    whose code is what Glynn's code is based on.

    Example:
    >>> from wavelen2rgb import wavelen2rgb
    >>> waves = [300.0, 400.0, 600.0]
    >>> rgb = [wavelen2rgb(waves[i], MaxIntensity=255) for i in range(3)]
    >>> print rgb
    [[0, 0, 0], [131, 0, 181], [255, 190, 0]]
    """

    def Adjust_and_Scale(Color, Factor, Highest=100):
        """Gamma adjustment.

        Arguments:
        * Color:  Value of R, G, or B, on a scale from 0 to 1, inclusive,
          with 0 being lowest intensity and 1 being highest.  Floating
          point value.
        * Factor:  Factor obtained to have intensity fall off at limits 
          of human vision.  Floating point value.
        * Highest:  Maximum intensity of output, scaled value.  The 
          lowest intensity is 0.  Scalar integer.

        Returns an adjusted and scaled value of R, G, or B, on a scale 
        from 0 to Highest, inclusive, as an integer, with 0 as the lowest 
        and Highest as highest intensity.

        Since this is a helper function I keep its existence hidden.
        See http://www.efg2.com/Lab/ScienceAndEngineering/Spectra.htm and
        http://www.physics.sfasu.edu/astro/color/spectra.html for details.
        """
        Gamma = 0.80

        if Color == 0.0:
            result = 0
        else:
            result = int( round(pow(Color * Factor, Gamma) * round(Highest)) )
            if result < 0:        result = 0
            if result > Highest:  result = Highest

        return result


    if (Wavelength >= 380.0) and (Wavelength < 440.0):
        Red   = -(Wavelength - 440.) / (440. - 380.)
        Green = 0.0
        Blue  = 1.0

    elif (Wavelength >= 440.0) and (Wavelength < 490.0):
        Red   = 0.0
        Green = (Wavelength - 440.) / (490. - 440.)
        Blue  = 1.0

    elif (Wavelength >= 490.0) and (Wavelength < 510.0):
        Red   = 0.0
        Green = 1.0
        Blue  = -(Wavelength - 510.) / (510. - 490.)

    elif (Wavelength >= 510.0) and (Wavelength < 580.0):
        Red   = (Wavelength - 510.) / (580. - 510.)
        Green = 1.0
        Blue  = 0.0

    elif (Wavelength >= 580.0) and (Wavelength < 645.0):
        Red   = 1.0
        Green = -(Wavelength - 645.) / (645. - 580.)
        Blue  = 0.0

    elif (Wavelength >= 645.0) and (Wavelength <= 780.0):
        Red   = 1.0
        Green = 0.0
        Blue  = 0.0

    else:
        Red   = 0.0
        Green = 0.0
        Blue  = 0.0


    #- Let the intensity fall off near the vision limits:

    if (Wavelength >= 380.0) and (Wavelength < 420.0):
        Factor = 0.3 + 0.7*(Wavelength - 380.) / (420. - 380.)
    elif (Wavelength >= 420.0) and (Wavelength < 701.0):
        Factor = 1.0
    elif (Wavelength >= 701.0) and (Wavelength <= 780.0):
        Factor = 0.3 + 0.7*(780. - Wavelength) / (780. - 700.)
    else:
        Factor = 0.0


    #- Adjust and scale RGB values to 0 to MaxIntensity integer range:

    R = Adjust_and_Scale(Red,   Factor, MaxIntensity)
    G = Adjust_and_Scale(Green, Factor, MaxIntensity)
    B = Adjust_and_Scale(Blue,  Factor, MaxIntensity)


    #- Return 3-element list value:

    return [R, G, B]

class fluorophore:
    """ Fluorophore class, creating for each spectra
    
    """
    def __init__(self, fluo_name, fluo_df, laser_list, ch_list):
        logging.info('==> {} spectra uploaded! <=='.format(fluo_name))
        self.fluo_name = fluo_name
        self.ex_spec = fluo_df[['w', 'ex']]
        self.em_spec = fluo_df[['w', 'em']]
        self.laser_list = laser_list
        self.ch_list = ch_list

        self.ex_dict = {}
        for laser in self.laser_list:
            try:
                self.ex_lvl = int(self.ex_spec.loc[self.ex_spec['w'] == laser, 'ex'])
            except ValueError:
                logging.fatal('{} DOESN`T excite at {} nm!'.format(self.fluo_name, laser))
                continue
            self.ex_dict.update({laser: self.ex_lvl})
            logging.info('{}|{} nm = {}'.format(self.fluo_name, laser, self.ex_lvl))

        self.ch_int = {}
        for self.ch_num in self.ch_list:
            self.ch_band = self.ch_list[self.ch_num]
            self.band = self.em_spec.loc[(self.em_spec['w'] >= self.ch_band[0]) & (self.em_spec['w'] <= self.ch_band[1]), 'em']
            self.band_int = int(self.band.sum())
            self.ch_int.update({self.ch_num: self.band_int})

            logging.info('{} integral int.|{}nm = {}'.format(self.fluo_name, self.ch_band, self.band_int))

class calcRes:
    """ Results of integral intensity and crosstalk calculation

    """
    def __init__(self, fluo_name):
        pass

FORMAT = "%(asctime)s| %(levelname)s [%(filename)s: - %(funcName)20s]  %(message)s"
logging.basicConfig(level=logging.INFO,
                    format=FORMAT)

plt.style.use('dark_background')
plt.rcParams['figure.facecolor'] = '#272b30'


# read settings file
with open('settings.yml') as f:
    settings_dict = yaml.safe_load(f)
    logging.info('Settings YAML file uploaded!')

mode = 'ch'

fluo_list = settings_dict['fluo_list']
ex_list = settings_dict['ex_list']
ch_dict = {}
try:
    ch_dict = settings_dict['ch_reg']
except KeyError:
    mode = 'view'
    logging.info('NO band pass selected, view mod!')


# read CSV spectra
mol_dict = {}
for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if file.endswith('.csv') and file.split('.')[0] in fluo_list:
            file_path = os.path.join(root, file)

            raw_csv = pd.read_csv(file_path)
            raw_csv.columns = ['w', 'ex', 'em']

            if raw_csv['em'].max() <= 1:
              raw_csv['ex'] = raw_csv['ex'] *100
              raw_csv['em'] = raw_csv['em'] *100

            mol_dict.update({file.split('.')[0]: fluorophore(file.split('.')[0], raw_csv, ex_list, ch_dict)})


# channels ratio calc
if len(fluo_list) > 1 and len(ch_dict.keys()) > 0:
    logging.info('Double fluorophore mode')
    mol_1 = mol_dict[fluo_list[0]] 
    mol_2 = mol_dict[fluo_list[1]]

    for ch in ch_dict:
        ch_ratio = round(mol_1.ch_int[ch] / mol_2.ch_int[ch], 3)
        logging.info('Ch. {} {} em./{} em. ratio = {}'.format(ch_dict[ch], mol_1.fluo_name, mol_2.fluo_name, ch_ratio))

        for ex_laser in ex_list:
            try:
                A = round(mol_1.ex_dict[ex_laser] / mol_2.ex_dict[ex_laser], 3)  # excitation ratio for two first fluoropheres
                # logging.info('  {}nm A factor ({} ex./{} ex.) = {}'.format(ex_laser, mol_1.fluo_name, mol_2.fluo_name, A))
                int_ratio = round(ch_ratio * A, 3)

                logging.info('  Ch. {} at {}nm corrected ratio = {} (A={})'.format(ch_dict[ch], ex_laser, int_ratio, A))

            except ZeroDivisionError:
                logging.fatal('  {}nm A factor ({} ex./{} ex.) DOESN`t exist, {} ex. is zero!'.format(ex_laser, mol_1.fluo_name, mol_2.fluo_name, mol_2.fluo_name))
                continue


# build plots
for fluo in mol_dict:
    fluo_spectra = mol_dict[fluo]
    plt.plot(fluo_spectra.em_spec.w, fluo_spectra.em_spec.em,  # emission
             label='{}, em'.format(fluo),
             color=wavelen2rgb(fluo_spectra.em_spec.loc[fluo_spectra.em_spec['em'].idxmax(), 'w'], 1))
    plt.plot(fluo_spectra.ex_spec.w, fluo_spectra.ex_spec.ex,  # excitation
             label='{}, ex'.format(fluo),
             linestyle='--',
             color=wavelen2rgb(fluo_spectra.em_spec.loc[fluo_spectra.em_spec['em'].idxmax(), 'w'], 1))

if mode != 'view':
    for ch in ch_dict:
        ch_band = ch_dict[ch]
        plt.fill_between(x=range(ch_band[0], ch_band[1], 1),
                     y1=0,
                     y2=110,
                     alpha=0.35,
                     label=ch)

for laser in ex_list:
    plt.plot([laser, laser], [0, 110],
             label='{} nm'.format(laser),
             linewidth=2,
             color=wavelen2rgb(laser, 1))

plt.xlabel('Wavelength, nm')
plt.ylabel('Intensity, %')
plt.xlim(300,800)  # 300, 800
plt.xticks(range(300,850,25))
plt.yticks(range(0,110,10))
plt.legend()
plt.tight_layout()
plt.show()
