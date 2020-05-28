#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This file contains supportive/backend functions related to modeling
'''

import pandas as pd
import numpy as np
import params
from pyteomics import mass
from IsoSpecPy import IsoSpecPy

class Cycler:
    '''
    Helper class allowing cycle time calculation
    '''
    def __init__(self, parallel):
        '''
        Constructor. Intializes three device queues:
            IS for ion source,
            OT for orbitrap,
            IT for ion trap
        
        Parameters
        ----------
        parallel : bool
            Use parallelization.

        Returns
        -------
        None.

        '''
        self.IS = [0] #ion source queue
        self.OT = [0] #OT queue
        self.IT = [0] #IT queue
        self.parallel = parallel #parallelization

    def whenFree(self, queue):
        '''
        Find timepoint when specific queue is free

        Parameters
        ----------
        queue : string
            Name of the queue, should be one of OT/IT/IS.

        Returns
        -------
        numerical
            time when the queue will be available.

        '''
        try:
            return self.__getattribute__(queue)[-1]
        except:
            raise ValueError("Unknown queue name: {}".format(queue))
    
    def whenAllFree(self):
        '''
        Find timepoint when all queues are free

        Returns
        -------
        numeric
            time when all queues are free.

        '''
        return max(self.IS[-1], self.IT[-1], self.OT[-1])
    
    def pushToQueue(self, queue, start, duration):
        '''
        Add element to queue

        Parameters
        ----------
        queue : string
            name of the queue, should be one of IT/OT/IS.
        start : numeric
            start time of an element.
        duration : numeric
            duration of an element.

        Returns
        -------
        None.

        '''
        try:
            workingQueue = self.__getattribute__(queue)
        except:
            raise ValueError("Unknown queue name: {}".format(queue))
        
        if len(workingQueue) == 1: # first element
            workingQueue.pop()
        
        workingQueue.extend([start, start + duration])
        
    def pushTask(self, task):
        '''
        Add acquisition task to cycle

        Parameters
        ----------
        task : tuple
            task should consist of three elements:
            ion collection time (numeric),
            name of device queue (string), should be IT/OT,
            dwell time (numeric).

        Returns
        -------
        None.

        '''            
        isTime, device, devTime = task
        
        if not device in ['OT', 'IT']:
            raise ValueError('Device {} is not allowed'.format(device))
        
        if self.parallel:
            starttime = max(self.whenFree('IS'), self.whenFree(device) - isTime)
        else:
            starttime = self.whenAllFree()
            
        self.pushToQueue('IS', starttime, isTime)
        self.pushToQueue(device, starttime + isTime, devTime)
        
    def getCycle(self):
        '''
        Get cycle time and queues
        
        Queues are represented as [N, 2] numpy.ndarrays,
        first dimnension - blocks, second dimension - from - to:
        [[start1, end1]
             ...
         [startN, endN]]
        
        Depending on parallelization, cycle time can be shorter, than longest
        device queue

        Returns
        -------
        cycletime : numeric
            length of cycle time.
        dict
            contains content of all three device queues,
            queue name is used as a key.

        '''
        if self.parallel: #parallalelize first ion collection with last dwell time
            firstInjection = self.IS[1]
            devFree = max(self.whenFree('OT'), self.whenFree('IT'))
            cycletime = max(self.whenFree('IS'), devFree - firstInjection)
        else:
            cycletime = self.whenAllFree()
            
        #checking for empty queues
        for queue in ['IS', 'IT', 'OT']:
            if len(self.__getattribute__(queue)) == 1:
                self.__getattribute__(queue).pop()
            
        return cycletime, {'IS': np.array(self.IS).reshape(-1, 2),
                           'IT': np.array(self.IT).reshape(-1, 2),
                           'OT': np.array(self.OT).reshape(-1, 2)}

      
def get_ions(pep_mass, charge):
    '''
    Convert neutral mass to m/z value
    '''
    return  (pep_mass + charge * 1.00727) / charge  

def get_profile_peak(mz, intensity, mz_grid, sigma):
    '''
    Generate gausian peak from position (mz), height (intensity),
    width (sigma), and plot grid (mz_grid)
    '''
    sigma /= 2 * np.sqrt(2 * np.log(2))
    return intensity*np.exp(-(mz_grid-mz)**2/(2*(sigma**2)))

def get_profile_spectrum(mz_intensity_list, r, points=41):        
    '''
    Parameters
    ----------
    mz_intensity_list : 2D numpy array shaped like
                            [[mz0, intensity0],
                             [mz1, intensity1],
                             ....
                             [mzN, intensityN]]
    r : resolution @ 200
    
    Return
        tuple of `numpy.ndarray` (m/z values, intensities)
    '''
    
    if mz_intensity_list.shape[0] == 0: #early termination for empty spectrum
       return np.array([params.low_mass, params.high_mass]), np.array([0, 0])
    
    full_spectrum = np.array([])
    full_grid = np.array([])
    real_r = r * np.sqrt(200 / mz_intensity_list[0][0])
    sigma = mz_intensity_list[0][0] / real_r
    sigma_coef = 2.5
    grid = np.linspace(mz_intensity_list[0][0] - sigma_coef * sigma, 
                       mz_intensity_list[0][0] + sigma_coef * sigma,
                       points)
    full_grid = np.append(full_grid, grid)
    y0 = get_profile_peak(mz_intensity_list[0][0], mz_intensity_list[0][1], full_grid, sigma)
    full_spectrum = np.append(full_spectrum, y0)

    if len(mz_intensity_list) > 1:
        for ind, mz_int in enumerate(mz_intensity_list[1:], 1):
            real_r = r * np.sqrt(200 / mz_int[0])
            sigma = mz_int[0] / real_r
            new_grid = np.linspace(mz_int[0] - sigma_coef * sigma, 
                                   mz_int[0] + sigma_coef * sigma,
                                   points) 
            
            if new_grid[0] < full_grid[-1]:
                start = np.argmin(np.abs(full_grid - new_grid[0]))
                end = np.argwhere( (new_grid - full_grid[-1]) < 0 )[-1][-1]
                full_spectrum[start:] += get_profile_peak(mz_int[0], mz_int[1], full_grid[start:], sigma)
                full_spectrum = np.append(full_spectrum, get_profile_peak(mz_int[0], mz_int[1], new_grid[end+1:], sigma))
                full_grid = np.append(full_grid, new_grid[end+1:])
            else:
                full_grid = np.append(full_grid, new_grid)
                full_spectrum = np.append(full_spectrum, get_profile_peak(mz_int[0], mz_int[1],new_grid, sigma))
                
    return full_grid, full_spectrum

def get_peptides(peptide_collection_size):
    '''
    Randomly selects peptides from generated tryptic digest.
    Parameters
    ----------
    peptide_collection_size : int, number of peptides to be selected.
    
    Return peptides list
    '''
    peptides = pd.read_csv('./assets/peptides.csv')
    peptide_slice = np.random.choice(peptides.index, peptide_collection_size, replace=False)
    return peptides.loc[peptide_slice, :].reset_index(drop=True)

def get_peptide_abundance(distribution, peptide_collection_size):
    '''
    Randomly generates abundances for peptides according to distribution.
    Parameters
    ----------
    distribution : str, one of 'equal', 'lognormal', 'lognormal-major'
    peptide_collection_size : sample size
    
    Return list of abundancies.
    '''
    if distribution == 'equal':
        return np.repeat(1, peptide_collection_size)
    elif distribution == 'lognormal':
        result = np.random.lognormal(9, 2, peptide_collection_size)
        #trimming
        overflow = result > np.exp(14)
        result[overflow] = result[overflow] / np.exp(9)
        underflow = result < np.exp(4)
        result[underflow] = result[underflow] * np.exp(9)
        return result
    elif distribution == 'lognormal-major': #90% of abundance for N Major peaks
        nmajor = int(peptide_collection_size * 0.05)
        result = get_peptide_abundance('lognormal', peptide_collection_size)
        order = np.argsort(result)
        result[order[:-nmajor]] = 0.1 * result[order[-nmajor:]].sum() * result[order[:-nmajor]] /\
                                    result[order[:-nmajor]].sum()
        
        return result
    else:
        raise ValueError('''Distribution should be one of equal|lognormal|lognormal-major;
                         '{}' provided'''.format(distribution))

def get_charge_state_probabilities(peptide_collection_size):
    '''
    Randomly generates intensities for two charge states.
    Parameters
    ----------
    peptide_collection_size : sample size
    
    Return `numpy.array`(sample size, 2).
    '''
    charge_state_2 = np.random.rand(peptide_collection_size,1)
    return [charge_state_2 ,1 - charge_state_2]

def expand_isotopes(peptide, charge_states=[2,3]):
    '''
    Convert peptide to DataFrame of isotopic peaks
    Input
        Series, should contain 'sequence', 'z+' columns, and model columns
    Return
        DataFrame with one row for each isotopic peak
        columns are:
            mz - m/z of ion
            ic_XX - ion abundance acording to XX model
            z  - charge
            sequence - peptide sequence
    '''
    formula=''.join(['{}{}'.format(x, y) for x, y in mass.Composition(peptide['sequence']).items()])
    cluster = IsoSpecPy.IsoSpec.IsoFromFormula(formula, cutoff=0.005, method='threshold_absolute').getConfs()
    mz0 = np.array(cluster[0])
    int0 = np.exp(cluster[1])    
    mz = np.concatenate([get_ions(mz0, z) for z in charge_states])
    ic = np.concatenate([int0 * peptide['{}+'.format(z)] for z in charge_states])
    charge = np.concatenate([np.repeat(z, mz0.shape[0]) for z in charge_states])
    result = pd.DataFrame({'mz': mz, 'ic': ic, 'z': charge})
    result['sequence'] = peptide['sequence']
    for model in params.ion_models:
        result['ic_{}'.format(model)] = result['ic'] * peptide[model]
        
    return result

def get_ion_data(nPeptides):
    '''
    Generate pandas.DataFrame with all ion data
    '''
    peptide_data = get_peptides(nPeptides)
    prob2, prob3 = get_charge_state_probabilities(nPeptides)
    peptide_data['2+'] = prob2
    peptide_data['3+'] = prob3
    for model in params.ion_models:
        peptide_data[model] = get_peptide_abundance(model, nPeptides)
    
    ion_data = pd.concat(peptide_data.apply(expand_isotopes, axis=1).tolist(), ignore_index=True)
    ion_data.sort_values(by='mz', inplace=True)
    
    return ion_data

def normalize_ion_currents(ion_data, low, high):
    '''
    Restrict m/z to (low, high) mass range
    '''
    in_mass = np.logical_and(ion_data['mz'] >= low, ion_data['mz'] <= high)
    ion_data.drop(ion_data.index[~in_mass], axis='index', inplace=True)

def scale_ion_currents(ion_data, tic):
    '''
    Scale ion intensities to get desired TIC
    '''
    for model in params.ion_models:
        ion_data['ic_{}'.format(model)] *= tic / ion_data['ic_{}'.format(model)].sum()

def get_boxes(low, high, nBoxes, nScans, overlap):
    '''
    Generate BoxCar boxes
    low - float - lowest mass
    high - float - highest mass
    nBoxes - int - number of boxes per scan
    nScans - int - number of scans
    overlap - float - overlap between boxes
    
    Return
        list of numpy.arrays
        the number of element of the list corresponds to number of scans
        the elements of the list is 2D arrays of window edges
        [[min1, max1],
         [min2, max2],
         .....
         [minX, maxX]]
    '''
    edges = np.linspace(low, high, nBoxes * nScans + 1)
    edges =  np.repeat(edges,2)[1:-1].reshape((-1,2)) + [[overlap/-2,overlap/2]]
    
    return [edges[s::nScans] for s in range(nScans)]

def add_boxes(ion_data, boxes):
    '''
    Assign ions to boxes
    '''
    for scan in range(len(boxes)):
        for box in range(boxes[scan].shape[0]):
            selector = np.logical_and(ion_data['mz'] > boxes[scan][box][0],
                                      ion_data['mz'] <= boxes[scan][box][1])
            ion_data.loc[selector, 'box'] = box
            ion_data.loc[selector, 'scan'] = scan

def sample_ions(ion_data, distribution, agc_target, max_it):
    '''
    Perform ion sampling using parameters below.
    
    Parameters
    ----------
    ion_data : DataFrame, contains ion currents for ions to be sampled
    distribution : str, one of 'equal', 'lognormal', 'lognormal-major'
    agc_target: float, number of ions to sample
    max_it: float, maximal injection time in milliseconds
    
    Return
        tuple of three elements
        1. 1D array of number of ions sampled
        2. required scan time seconds
        3. acquired number of ions
    '''
    tic = ion_data['ic_{}'.format(distribution)].sum()
    agc = min(agc_target, tic * max_it * 1e-3)
    scan_time = agc / tic
    probabilities = np.cumsum(ion_data['ic_{}'.format(distribution)].values)
    probabilities = np.append([0], probabilities / probabilities[-1])
    intensities = np.histogram(np.random.random(int(agc)), 
                               bins=probabilities)[0].astype(float)
    
    return intensities, scan_time, agc


def get_full_spectrum(ion_data, distribution, agc_target, max_it):
    '''
    Create centroids for full spectrum using parameters below.
    
    Parameters
    ----------
    ion_data : DataFrame, contains ion currents for all ions
    distribution : str, one of 'equal', 'lognormal', 'lognormal-major'
    agc_target: float, number of ions to sample
    max_it: float, maximal injection time in milliseconds
    
    Return
        tuple of six elements
        1. 2D array of mz and intensities
                [[mz0, intensity0],
                 [mz1, intensity1],
                 ....
                 [mzN, intensityN]]        
        2. required scan time in milliseconds
        3. acquired number of ions
        4. set of observed peptide sequences
        5. maximum observed ion intensity under the distrubution
        6. minimum observed ion intensity under the distribution
    '''
    intensities, scan_time, agc = sample_ions(ion_data, distribution, agc_target, max_it)
    dyn_range_filter = intensities > max(intensities.max() * 1e-4, 10)
    mzdata = np.stack((ion_data['mz'].values, intensities), axis=-1)
    scan_ion_data = ion_data[dyn_range_filter]
    if scan_ion_data.shape[0] > 0: #non-empty
        peptides = set(scan_ion_data['sequence'])
        max_int = scan_ion_data['ic_' + distribution].max()
        min_int = scan_ion_data['ic_' + distribution].min()
    else:
        peptides, max_int, min_int = set(), -1, -1
        
    return mzdata[dyn_range_filter, :], scan_time*1000, agc, peptides, max_int, min_int

def get_boxcar_spectra(ion_data, distribution, agc_target, max_it, nBoxes, nScans):
    '''
    Create centroids for boxcar spectra using parameters below.
    
    Parameters
    ----------
    ion_data : DataFrame, contains ion currents for all ions
    distribution : str, one of 'equal', 'lognormal', 'lognormal-major'
    agc_target: float, number of ions to sample per scan
    max_it: float, maximal injection time in milliseconds per scan
    nBoxes: int, number of boxes per scan
    nScans: int, number of scans
    
    Return
        one tuple of six elements per each boxcar scan
        1. 2D array of mz and intensities
                [[mz0, intensity0],
                 [mz1, intensity1],
                 ....
                 [mzN, intensityN]]        
        2. required scan time in milliseconds
        3. acquired number of ions
        4. set of observed peptide sequences
        5. maximum observed ion intensity under the distrubution
        6. minimum observed ion intensity under the distribution
    '''
    BCscans = []
    for scan in range(nScans):
        scan_mz = []
        scan_counts = []
        scan_time = 0
        agc = 0
        for box in range(nBoxes):
            selector = np.logical_and(ion_data['box'] == box, ion_data['scan'] == scan)
            if selector.sum() > 0:
                intensities, box_time, box_agc = sample_ions(ion_data[selector],
                                                             distribution,
                                                             agc_target/nBoxes,
                                                             max_it/nBoxes)
                scan_time += box_time
                agc += box_agc
                scan_mz.append(ion_data.loc[selector, 'mz'].values)
                scan_counts.append(intensities)
            else:
               scan_time += 1e-3 * max_it / nBoxes

        scan_mz = np.concatenate(scan_mz)
        scan_counts = np.concatenate(scan_counts)
        dyn_range_filter = scan_counts > max(scan_counts.max() * 1e-4, 10)
        mzdata = np.stack((scan_mz, scan_counts), axis=-1)
        scan_ion_data = ion_data[ion_data['scan'] == scan][dyn_range_filter]
        
        if scan_ion_data.shape[0] > 0: #non-empty
            peptides = set(scan_ion_data['sequence'])
            max_int = scan_ion_data['ic_' + distribution].max()
            min_int = scan_ion_data['ic_' + distribution].min()
        else:
            peptides, max_int, min_int = set(), -1, -1
            
        BCscans.append((mzdata[dyn_range_filter, :], scan_time*1000, agc,
                        peptides, max_int, min_int))
    
    return np.array(BCscans)

def get_MS_counts(scan_method, acc_time, resolution, topN, ms2resolution,
                  ms2IT, time, parallel=False):
    '''
    Calculate number of MS1 and MS2 scans using parameters below.
    Parameters:
        scan_method, str, one of 'full'|'boxcar'
        acc_time, float (full scan) or iterable of floats (boxcar), ion accumulation times required for the scan
            in case of full scan only one value is provided
            in case of boxcar scan, the iterable has to be, MS1 accumulation time, and acc_times for all boxcar scans
        resolution, int, used resolution (used to calculate transient time)
        topN, float, average TopN
        ms2resolution, int or string, resolution for MS/MS scans
        ms2IT, float, injection time for MS/MS scans
        time, float, the length of the gradient
        parallel, bool, parallelization mode
    Return:
        tuple, (cycle time, number of MS1 scans, number of MS2 scans, scan cycle queues)
    '''
    ms2device = 'IT' if ms2resolution == 'IT' else 'OT' #select ms2 device
    
    cycler = Cycler(parallel)
    
    if scan_method == 'full':
        cycler.pushTask((acc_time, 'OT', params.transients[resolution]))
   
    elif scan_method == 'boxcar':
        for at in acc_time:
            cycler.pushTask((at, 'OT', params.transients[resolution]))
    else:
        raise ValueError('scan_method has to be one of "full"|"boxcar"')
    
    for _ in range(topN):
        cycler.pushTask((ms2IT, ms2device, params.transients[ms2resolution]))
        
    cycletime, queues = cycler.getCycle()
    
    nMS1 = int(60000 * time / cycletime)
    nMS2 = int(topN * nMS1)
    
    return cycletime, nMS1, nMS2, queues
